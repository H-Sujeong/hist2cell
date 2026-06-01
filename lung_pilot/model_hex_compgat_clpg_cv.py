# model_hex_context_pg.py
# HEX(ALF+FDS) baseline with optional:
#   1) concat-context GAT spatial encoder
#   2) residual marker graph correction (PG as output refinement)
#   3) CL projection head support
# Input features are pre-extracted h-optimus-0 vectors: (N, 1024).

import math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.ndimage import gaussian_filter1d
from scipy.signal.windows import triang


# ============================================================
# 1. Graph construction
# ============================================================

def build_knn_graph(coords, k=4, max_distance=None):
    """Build directed kNN graph from coordinates."""
    device = coords.device
    N = coords.shape[0]
    if N <= 1:
        return torch.empty((2, 0), dtype=torch.long, device=device)

    k = min(int(k), N - 1)
    coords_cpu = coords.detach().cpu().float()

    try:
        from sklearn.neighbors import NearestNeighbors
        nbrs = NearestNeighbors(n_neighbors=k + 1, algorithm="auto").fit(coords_cpu.numpy())
        dists, indices = nbrs.kneighbors(coords_cpu.numpy())
        dists, indices = dists[:, 1:], indices[:, 1:]
        dst = torch.arange(N).view(-1, 1).repeat(1, k).reshape(-1)
        src = torch.tensor(indices, dtype=torch.long).reshape(-1)
        dist_flat = torch.tensor(dists, dtype=torch.float32).reshape(-1)
        edge_index = torch.stack([src, dst], dim=0)
        if max_distance is not None and max_distance > 0:
            edge_index = edge_index[:, dist_flat <= float(max_distance)]
        return edge_index.to(device)
    except Exception:
        dist = torch.cdist(coords.float(), coords.float())
        knn_dist, knn = dist.topk(k=k + 1, largest=False)
        knn, knn_dist = knn[:, 1:], knn_dist[:, 1:]
        dst = torch.arange(N, device=device).view(-1, 1).repeat(1, k).reshape(-1)
        src = knn.reshape(-1)
        edge_index = torch.stack([src, dst], dim=0)
        if max_distance is not None and max_distance > 0:
            edge_index = edge_index[:, knn_dist.reshape(-1) <= float(max_distance)]
        return edge_index


# ============================================================
# 2. Lightweight GAT layer
# ============================================================

class EdgeGATLayer(nn.Module):
    def __init__(self, in_dim, out_dim, heads=4, dropout=0.1, concat=True, negative_slope=0.2):
        super().__init__()
        self.out_dim = int(out_dim)
        self.heads = int(heads)
        self.concat = bool(concat)
        self.dropout = float(dropout)
        self.negative_slope = float(negative_slope)

        self.lin = nn.Linear(in_dim, self.heads * self.out_dim, bias=False)
        self.att_src = nn.Parameter(torch.empty(self.heads, self.out_dim))
        self.att_dst = nn.Parameter(torch.empty(self.heads, self.out_dim))
        out_dim_total = self.heads * self.out_dim if self.concat else self.out_dim
        self.out_proj = nn.Linear(out_dim_total, out_dim_total)
        self.reset_parameters()

    def reset_parameters(self):
        nn.init.xavier_uniform_(self.lin.weight)
        nn.init.xavier_uniform_(self.att_src)
        nn.init.xavier_uniform_(self.att_dst)
        nn.init.xavier_uniform_(self.out_proj.weight)
        nn.init.zeros_(self.out_proj.bias)

    def edge_softmax(self, scores, dst, num_nodes):
        if scores.numel() == 0:
            return scores
        _, H = scores.shape
        dst_expand = dst.view(-1, 1).expand(-1, H)
        max_per_dst = torch.full((num_nodes, H), -1e15, device=scores.device, dtype=scores.dtype)
        max_per_dst.scatter_reduce_(0, dst_expand, scores, reduce="amax", include_self=True)
        exp_scores = torch.exp(scores - max_per_dst[dst])
        denom = torch.zeros((num_nodes, H), device=scores.device, dtype=scores.dtype)
        denom.index_add_(0, dst, exp_scores)
        return exp_scores / (denom[dst] + 1e-8)

    def forward(self, x, edge_index):
        N = x.shape[0]
        if edge_index is None or edge_index.numel() == 0:
            h = self.lin(x).view(N, self.heads, self.out_dim)
            out = h.reshape(N, self.heads * self.out_dim) if self.concat else h.mean(dim=1)
            return self.out_proj(out), torch.empty((0,), device=x.device, dtype=x.dtype)

        src, dst = edge_index[0], edge_index[1]
        h = self.lin(x).view(N, self.heads, self.out_dim)
        h_src, h_dst = h[src], h[dst]
        scores = F.leaky_relu(
            (h_src * self.att_src).sum(-1) + (h_dst * self.att_dst).sum(-1),
            negative_slope=self.negative_slope,
        )
        alpha = self.edge_softmax(scores, dst, N)
        alpha = F.dropout(alpha, p=self.dropout, training=self.training)
        msg = h_src * alpha.unsqueeze(-1)
        out = torch.zeros((N, self.heads, self.out_dim), device=x.device, dtype=x.dtype)
        out.index_add_(0, dst, msg)
        out = out.reshape(N, self.heads * self.out_dim) if self.concat else out.mean(dim=1)
        return self.out_proj(out), alpha.mean(dim=1)


class ConcatGATContextEncoder(nn.Module):
    """
    Safer GAT context module.

    Instead of replacing or residual-adding to h, it computes a context vector and
    lets an MLP decide how much of [h, gamma * context] to use:

        context = GAT(h)
        h_out = MLP([h, gamma * context])

    This preserves the HEX feature path more strongly than residual GAT.
    """
    def __init__(self, feature_dim=128, num_layers=1, heads=4, dropout=0.1, gamma_init=0.01):
        super().__init__()
        assert feature_dim % heads == 0, "feature_dim must be divisible by heads"
        self.layers = nn.ModuleList([
            EdgeGATLayer(feature_dim, feature_dim // heads, heads=heads, dropout=dropout, concat=True)
            for _ in range(num_layers)
        ])
        self.norms = nn.ModuleList([nn.LayerNorm(feature_dim) for _ in range(num_layers)])
        self.dropout = float(dropout)
        self.gamma = nn.Parameter(torch.tensor(float(gamma_init)))
        self.fusion = nn.Sequential(
            nn.Linear(feature_dim * 2, feature_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.LayerNorm(feature_dim),
        )

    def forward(self, h, edge_index):
        ctx = h
        alpha_last = None
        for layer, norm in zip(self.layers, self.norms):
            msg, alpha = layer(ctx, edge_index)
            ctx = norm(msg)
            ctx = F.relu(ctx)
            ctx = F.dropout(ctx, p=self.dropout, training=self.training)
            alpha_last = alpha
        h_out = self.fusion(torch.cat([h, self.gamma * ctx], dim=1))
        return h_out, alpha_last


class RelationGATContextEncoder(nn.Module):
    """
    Relation-specific context encoder for compartment-aware GAT.

    edge_dict maps relation names to edge_index tensors. Each relation has its
    own GAT stack and a small learned gamma initialized near zero. The original
    HEX feature path is preserved by fusing [h, gamma_r * context_r ...].
    """
    def __init__(self, feature_dim=128, relation_names=None, num_layers=1, heads=4, dropout=0.1, gamma_init=0.01):
        super().__init__()
        assert feature_dim % heads == 0, "feature_dim must be divisible by heads"
        if relation_names is None:
            relation_names = ["same", "inter"]
        self.relation_names = list(relation_names)
        self.feature_dim = int(feature_dim)
        self.dropout = float(dropout)

        self.rel_layers = nn.ModuleDict()
        self.rel_norms = nn.ModuleDict()
        for rel in self.relation_names:
            self.rel_layers[rel] = nn.ModuleList([
                EdgeGATLayer(feature_dim, feature_dim // heads, heads=heads, dropout=dropout, concat=True)
                for _ in range(num_layers)
            ])
            self.rel_norms[rel] = nn.ModuleList([nn.LayerNorm(feature_dim) for _ in range(num_layers)])

        self.gamma = nn.ParameterDict({
            rel: nn.Parameter(torch.tensor(float(gamma_init))) for rel in self.relation_names
        })
        self.fusion = nn.Sequential(
            nn.Linear(feature_dim * (1 + len(self.relation_names)), feature_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.LayerNorm(feature_dim),
        )

    def forward(self, h, edge_dict):
        contexts = []
        alpha_dict = {}
        empty = torch.empty((2, 0), dtype=torch.long, device=h.device)
        for rel in self.relation_names:
            edge_index = None if edge_dict is None else edge_dict.get(rel, None)
            if edge_index is None:
                edge_index = empty
            ctx = h
            alpha_last = torch.empty((0,), dtype=h.dtype, device=h.device)
            for layer, norm in zip(self.rel_layers[rel], self.rel_norms[rel]):
                msg, alpha = layer(ctx, edge_index)
                ctx = norm(msg)
                ctx = F.relu(ctx)
                ctx = F.dropout(ctx, p=self.dropout, training=self.training)
                alpha_last = alpha
            contexts.append(self.gamma[rel] * ctx)
            alpha_dict[rel] = alpha_last
        h_out = self.fusion(torch.cat([h] + contexts, dim=1))
        return h_out, alpha_dict


def build_compartment_edge_dict(coords, comp, k=8, max_distance=None, mode="binary_comp"):
    """Build relation-specific directed kNN edges from hard compartment labels.

    Compartment convention:
      0 other, 1 tumor, 2 T cell, 3 B cell, 4 myeloid, 5 stroma, 6 vascular
    """
    edge_index = build_knn_graph(coords, k=k, max_distance=max_distance)
    device = coords.device
    if edge_index.numel() == 0:
        names = ["same", "inter"] if mode == "binary_comp" else [
            "same", "tumor_tcell", "tumor_bcell", "tumor_myeloid", "tumor_stroma",
            "tumor_vascular", "immune_myeloid", "immune_vascular", "other_inter"
        ]
        return {name: edge_index for name in names}

    comp = comp.to(device=device, dtype=torch.long)
    src, dst = edge_index[0], edge_index[1]
    cs, cd = comp[src], comp[dst]
    same = cs == cd

    if mode == "binary_comp":
        return {
            "same": edge_index[:, same],
            "inter": edge_index[:, ~same],
        }

    def pair(a, b):
        return ((cs == a) & (cd == b)) | ((cs == b) & (cd == a))

    tumor_tcell = pair(1, 2)
    tumor_bcell = pair(1, 3)
    tumor_myeloid = pair(1, 4)
    tumor_stroma = pair(1, 5)
    tumor_vascular = pair(1, 6)
    immune_myeloid = pair(2, 4) | pair(3, 4)
    immune_vascular = pair(2, 6) | pair(3, 6)
    used = same | tumor_tcell | tumor_bcell | tumor_myeloid | tumor_stroma | tumor_vascular | immune_myeloid | immune_vascular
    other_inter = (~same) & (~used)
    return {
        "same": edge_index[:, same],
        "tumor_tcell": edge_index[:, tumor_tcell],
        "tumor_bcell": edge_index[:, tumor_bcell],
        "tumor_myeloid": edge_index[:, tumor_myeloid],
        "tumor_stroma": edge_index[:, tumor_stroma],
        "tumor_vascular": edge_index[:, tumor_vascular],
        "immune_myeloid": edge_index[:, immune_myeloid],
        "immune_vascular": edge_index[:, immune_vascular],
        "other_inter": edge_index[:, other_inter],
    }


# ============================================================
# 3. Residual marker graph correction
# ============================================================

class ResidualMarkerGraphHead(nn.Module):
    """
    Protein graph as a residual marker correction, not a full decoder.

        raw = Linear(h)                         # HEX-like prediction
        A   = softmax(B + beta * A_prior)       # marker dependency
        pred = raw + eta * raw @ A^T

    This keeps the direct HEX prediction path intact and only allows a small
    marker-dependency correction.
    """
    def __init__(self, in_dim=128, num_markers=19, prior_adj=None, use_prior=True, eta_init=0.01, beta_init=0.1, dropout=0.0):
        super().__init__()
        self.num_markers = int(num_markers)
        self.use_prior = bool(use_prior)
        self.raw_head = nn.Linear(in_dim, num_markers)
        self.adj_logits = nn.Parameter(torch.zeros(num_markers, num_markers))
        self.eta = nn.Parameter(torch.tensor(float(eta_init)))
        self.beta = nn.Parameter(torch.tensor(float(beta_init)))
        self.dropout = float(dropout)
        if prior_adj is None:
            prior_adj = torch.zeros(num_markers, num_markers)
        self.register_buffer("prior_adj", prior_adj.float())
        self.register_buffer("eye", torch.eye(num_markers))

    def get_adj(self):
        scores = self.adj_logits
        if self.use_prior:
            scores = scores + self.beta * self.prior_adj
        # avoid self-correction dominating; raw path already keeps self information
        scores = scores.masked_fill(self.eye.bool(), -1e4)
        A = torch.softmax(scores, dim=-1)
        A = F.dropout(A, p=self.dropout, training=self.training)
        return A

    def forward(self, h):
        raw = self.raw_head(h)
        A = self.get_adj()
        correction = torch.matmul(raw, A.t())
        pred = raw + self.eta * correction
        return pred


# ============================================================
# 4. FDS copied/adapted from HEX-style implementation
# ============================================================

def calibrate_mean_var(matrix, m1, v1, m2, v2, clip_min=0.1, clip_max=10.0):
    if torch.sum(v1) < 1e-10:
        return matrix
    v1p, v2p = torch.clamp(v1, min=0.0), torch.clamp(v2, min=0.0)
    if (v1p == 0.).any():
        valid = (v1p != 0.)
        if valid.any():
            factor = torch.clamp(v2p[valid] / v1p[valid], clip_min, clip_max)
            matrix[:, valid] = (matrix[:, valid] - m1[valid]) * torch.sqrt(factor) + m2[valid]
        return matrix
    factor = torch.clamp(v2p / v1p, clip_min, clip_max)
    return (matrix - m1) * torch.sqrt(factor) + m2


class FDS(nn.Module):
    def __init__(self, feature_dim, bucket_num=50, bucket_start=0, start_update=0, start_smooth=10, kernel='gaussian', ks=9, sigma=2, momentum=0.9):
        super().__init__()
        self.feature_dim, self.bucket_num, self.bucket_start = int(feature_dim), int(bucket_num), int(bucket_start)
        self.half_ks, self.momentum = (int(ks) - 1) // 2, float(momentum)
        self.start_update, self.start_smooth = int(start_update), int(start_smooth)
        B, D = self.bucket_num - self.bucket_start, self.feature_dim
        self.register_buffer('epoch', torch.zeros(1, dtype=torch.long).fill_(self.start_update))
        self.register_buffer('running_mean', torch.zeros(B, D))
        self.register_buffer('running_var', torch.ones(B, D))
        self.register_buffer('running_mean_last_epoch', torch.zeros(B, D))
        self.register_buffer('running_var_last_epoch', torch.ones(B, D))
        self.register_buffer('smoothed_mean_last_epoch', torch.zeros(B, D))
        self.register_buffer('smoothed_var_last_epoch', torch.ones(B, D))
        self.register_buffer('num_samples_tracked', torch.zeros(B))
        self.register_buffer('kernel_window', self._get_kernel_window(kernel, int(ks), sigma))

    @staticmethod
    def _get_kernel_window(kernel, ks, sigma):
        assert kernel in ['gaussian', 'triang', 'laplace']
        half_ks = (ks - 1) // 2
        if kernel == 'gaussian':
            base = np.array([0.] * half_ks + [1.] + [0.] * half_ks, dtype=np.float32)
            win = gaussian_filter1d(base, sigma=sigma)
            win = win / np.sum(win)
        elif kernel == 'triang':
            win = triang(ks)
            win = win / np.sum(win)
        else:
            laplace = lambda x: np.exp(-abs(x) / sigma) / (2. * sigma)
            vals = np.array(list(map(laplace, np.arange(-half_ks, half_ks + 1))), dtype=np.float32)
            win = vals / np.sum(vals)
        return torch.tensor(win, dtype=torch.float32)

    def _get_bucket_idx_vec(self, labels):
        labels = np.clip(np.asarray(labels, dtype=np.float32), 0.0, 1.0)
        buckets = np.clip((labels * self.bucket_num).astype(np.int64), 0, self.bucket_num - 1)
        if self.bucket_start > 0:
            buckets = np.maximum(buckets, self.bucket_start)
        return buckets

    def _update_last_epoch_stats(self):
        self.running_mean_last_epoch.copy_(self.running_mean)
        self.running_var_last_epoch.copy_(self.running_var)
        mean_in = F.pad(self.running_mean_last_epoch.unsqueeze(1).permute(2, 1, 0), (self.half_ks, self.half_ks), mode='reflect')
        var_in = F.pad(self.running_var_last_epoch.unsqueeze(1).permute(2, 1, 0), (self.half_ks, self.half_ks), mode='reflect')
        w = self.kernel_window.view(1, 1, -1)
        self.smoothed_mean_last_epoch.copy_(F.conv1d(mean_in, w).permute(2, 1, 0).squeeze(1))
        self.smoothed_var_last_epoch.copy_(F.conv1d(var_in, w).permute(2, 1, 0).squeeze(1))

    def update_last_epoch_stats(self, epoch):
        if epoch == int(self.epoch.item()) + 1:
            self.epoch += 1
            self._update_last_epoch_stats()

    @torch.no_grad()
    def update_running_stats_from_moments(self, count, sum_feat, sumsq_feat, epoch):
        if epoch < int(self.epoch.item()):
            return
        b0, b1 = int(self.bucket_start), int(self.bucket_num)
        count, sum_feat, sumsq_feat = count.long(), sum_feat.float(), sumsq_feat.float()
        for bucket in range(b0, b1):
            n = int(count[bucket].item())
            if n <= 0:
                continue
            mean = sum_feat[bucket] / float(n)
            if n > 1:
                var = (sumsq_feat[bucket] - (sum_feat[bucket] * sum_feat[bucket]) / float(n)) / float(n - 1)
            else:
                var = torch.zeros_like(mean)
            var = torch.clamp(var, min=0.0)
            factor = 0.0 if epoch == self.start_update else float(self.momentum)
            idx = bucket - b0
            self.num_samples_tracked[idx] += float(n)
            self.running_mean[idx] = (1 - factor) * mean.to(self.running_mean.dtype) + factor * self.running_mean[idx]
            self.running_var[idx] = (1 - factor) * var.to(self.running_var.dtype) + factor * self.running_var[idx]

        present = set(int(b + b0) for b in torch.nonzero(count[b0:b1] > 0, as_tuple=False).view(-1).cpu().tolist())
        for bucket in range(b0, b1):
            idx = bucket - b0
            if bucket in present or float(self.num_samples_tracked[idx].item()) > 0:
                continue
            if bucket == b0:
                self.running_mean[idx] = self.running_mean[min(idx + 1, b1 - b0 - 1)]
                self.running_var[idx] = self.running_var[min(idx + 1, b1 - b0 - 1)]
            elif bucket == b1 - 1:
                self.running_mean[idx] = self.running_mean[idx - 1]
                self.running_var[idx] = self.running_var[idx - 1]
            else:
                self.running_mean[idx] = (self.running_mean[idx - 1] + self.running_mean[idx + 1]) / 2.0
                self.running_var[idx] = (self.running_var[idx - 1] + self.running_var[idx + 1]) / 2.0

    def smooth(self, features, labels, epoch):
        if epoch < self.start_smooth:
            return features
        feat = features.float().clone()
        labels_np = labels.detach().cpu().numpy() if torch.is_tensor(labels) else labels
        buckets = self._get_bucket_idx_vec(labels_np)
        for bucket in np.unique(buckets):
            mask = torch.as_tensor((buckets == bucket).astype(bool), device=feat.device)
            idx = bucket - self.bucket_start
            feat[mask] = calibrate_mean_var(
                feat[mask],
                self.running_mean_last_epoch[idx],
                self.running_var_last_epoch[idx],
                self.smoothed_mean_last_epoch[idx],
                self.smoothed_var_last_epoch[idx],
            )
        return feat.to(features.dtype)


# ============================================================
# 5. Unified model
# ============================================================

class HEXContextModel(nn.Module):
    """
    Baseline:
        h-optimus feature -> HEX MLP 1024->256->128 -> linear output

    Optional modules:
        --use_spatial_encoder: concat-context GAT after HEX MLP
        --use_protein_graph: residual marker graph correction head
        --use_fds: HEX-style FDS during training
        --use_cl: handled in main.py using z_cl returned here
    """
    def __init__(
        self,
        in_dim=1024,
        num_markers=19,
        hidden_dim1=256,
        hidden_dim2=128,
        dropout=0.5,
        use_spatial_encoder=False,
        gat_layers=1,
        gat_heads=4,
        gat_dropout=0.1,
        gat_gamma_init=0.01,
        gat_mode="knn",
        comp_source="oracle",
        num_compartments=7,
        use_comp_head=False,
        use_protein_graph=False,
        marker_graph_eta_init=0.01,
        marker_graph_beta_init=0.1,
        marker_graph_dropout=0.0,
        prior_adj=None,
        use_prior=True,
        use_fds=True,
        fds_active_markers=None,
        fds_bucket_num=50,
        fds_bucket_start=0,
        fds_label_min=-3.0,
        fds_label_max=3.0,
        fds_momentum=0.9,
        fds_start_smooth=10,
        fds_start_update=0,
        fds_kernel='gaussian',
        fds_ks=9,
        fds_sigma=2.0,
        cl_proj_dim=64,
    ):
        super().__init__()
        self.num_markers = int(num_markers)
        self.hidden_dim2 = int(hidden_dim2)
        self.use_spatial_encoder = bool(use_spatial_encoder)
        self.gat_mode = str(gat_mode)
        self.comp_source = str(comp_source)
        self.num_compartments = int(num_compartments)
        self.use_comp_head = bool(use_comp_head)
        self.use_protein_graph = bool(use_protein_graph)
        self.use_fds = bool(use_fds)
        self.training_status = True
        self.fds_label_min = float(fds_label_min)
        self.fds_label_max = float(fds_label_max)
        self.fds_active_markers = list(range(num_markers)) if fds_active_markers is None else [int(x) for x in fds_active_markers]

        self.regression_head = nn.Sequential(
            nn.Linear(in_dim, hidden_dim1),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim1, hidden_dim2),
            nn.ReLU(),
            nn.Dropout(dropout),
        )

        if self.use_spatial_encoder and self.gat_mode in ["binary_comp", "relation_comp"]:
            relation_names = ["same", "inter"] if self.gat_mode == "binary_comp" else [
                "same", "tumor_tcell", "tumor_bcell", "tumor_myeloid", "tumor_stroma",
                "tumor_vascular", "immune_myeloid", "immune_vascular", "other_inter"
            ]
            self.spatial_encoder = RelationGATContextEncoder(
                hidden_dim2, relation_names=relation_names, num_layers=gat_layers,
                heads=gat_heads, dropout=gat_dropout, gamma_init=gat_gamma_init
            )
        elif self.use_spatial_encoder:
            self.spatial_encoder = ConcatGATContextEncoder(hidden_dim2, gat_layers, gat_heads, gat_dropout, gat_gamma_init)
        else:
            self.spatial_encoder = None

        self.comp_head = nn.Sequential(
            nn.Linear(hidden_dim2, hidden_dim2),
            nn.ReLU(),
            nn.Dropout(gat_dropout),
            nn.Linear(hidden_dim2, self.num_compartments),
        ) if self.use_comp_head else None

        if self.use_protein_graph:
            self.output_head = ResidualMarkerGraphHead(
                in_dim=hidden_dim2,
                num_markers=num_markers,
                prior_adj=prior_adj,
                use_prior=use_prior,
                eta_init=marker_graph_eta_init,
                beta_init=marker_graph_beta_init,
                dropout=marker_graph_dropout,
            )
            self.linear_output = None
        else:
            self.linear_output = nn.Linear(hidden_dim2, num_markers)
            self.output_head = None

        # CL is applied on z_cl, not directly on regression h.
        self.cl_projection = nn.Sequential(
            nn.Linear(hidden_dim2, hidden_dim2),
            nn.ReLU(),
            nn.Linear(hidden_dim2, cl_proj_dim),
        )

        self.FDS = nn.ModuleList([
            FDS(hidden_dim2, fds_bucket_num, fds_bucket_start, fds_start_update, fds_start_smooth, fds_kernel, fds_ks, fds_sigma, fds_momentum)
            for _ in range(num_markers)
        ]) if self.use_fds else nn.ModuleList([])

    def _labels_to_fds_range(self, labels):
        y = labels.detach().float()
        y = (y - self.fds_label_min) / (self.fds_label_max - self.fds_label_min + 1e-6)
        return torch.clamp(y, 0.0, 1.0)

    def _apply_output_head(self, h):
        return self.output_head(h) if self.use_protein_graph else self.linear_output(h)

    def forward(
        self,
        features,
        coords=None,
        edge_index=None,
        edge_dict=None,
        comp_target=None,
        knn_k=4,
        max_edge_distance=None,
        labels=None,
        epoch=0,
    ):
        h = self.regression_head(features)

        comp_logits = self.comp_head(h) if self.comp_head is not None else None

        if self.use_spatial_encoder:
            if coords is None:
                raise ValueError("coords are required when use_spatial_encoder=True")
            if self.gat_mode in ["binary_comp", "relation_comp"]:
                if edge_dict is None:
                    if self.comp_source == "predicted":
                        if comp_logits is None:
                            raise ValueError("comp_source='predicted' requires use_comp_head=True")
                        comp_for_graph = comp_logits.detach().argmax(dim=1)
                    else:
                        if comp_target is None:
                            raise ValueError("comp_source='oracle' requires comp_target")
                        comp_for_graph = comp_target
                    edge_dict = build_compartment_edge_dict(
                        coords, comp_for_graph, k=knn_k, max_distance=max_edge_distance, mode=self.gat_mode
                    )
                h, edge_alpha = self.spatial_encoder(h, edge_dict)
                edge_index = torch.empty((2, 0), dtype=torch.long, device=features.device)
            else:
                if edge_index is None:
                    edge_index = build_knn_graph(coords, k=knn_k, max_distance=max_edge_distance)
                h, edge_alpha = self.spatial_encoder(h, edge_index)
        else:
            edge_index = torch.empty((2, 0), dtype=torch.long, device=features.device)
            edge_alpha = torch.empty((0,), dtype=features.dtype, device=features.device)

        # Projection branch for CL. Regression path still uses h directly.
        z_cl = self.cl_projection(h)

        if self.use_fds and self.training and self.training_status and labels is not None and len(self.fds_active_markers) > 0 and epoch >= self.FDS[0].start_smooth:
            labels_fds = self._labels_to_fds_range(labels)

            if self.use_protein_graph:
                # Conservative FDS for residual PG: average active marker-specific smoothed features,
                # then apply residual marker correction head once.
                h_accum = torch.zeros_like(h)
                for j in self.fds_active_markers:
                    h_accum += self.FDS[int(j)].smooth(h.clone(), labels_fds[:, int(j)], epoch)
                h_for_pred = h_accum / max(len(self.fds_active_markers), 1)
                pred = self.output_head(h_for_pred)
            else:
                pred = self.linear_output(h)
                weight, bias = self.linear_output.weight, self.linear_output.bias
                if len(self.fds_active_markers) == 1:
                    j0 = int(self.fds_active_markers[0])
                    pred = self.linear_output(self.FDS[j0].smooth(h.clone(), labels_fds[:, j0], epoch))
                else:
                    for j in self.fds_active_markers:
                        j = int(j)
                        h_smooth = self.FDS[j].smooth(h.clone(), labels_fds[:, j], epoch)
                        pred[:, j] = F.linear(h_smooth, weight[j:j+1], None if bias is None else bias[j:j+1]).squeeze(1)
        else:
            pred = self._apply_output_head(h)

        return {
            "pred": pred,
            "h": h,
            "z_cl": z_cl,
            "edge_index": edge_index,
            "edge_alpha": edge_alpha,
            "comp_logits": comp_logits,
        }

    def create_fds_moment_buffers(self, device):
        if not self.use_fds or len(self.fds_active_markers) == 0:
            return None
        bucket_num = int(self.FDS[0].bucket_num)
        feature_dim = int(self.FDS[0].feature_dim)
        k = len(self.fds_active_markers)
        return {
            "count": torch.zeros(k, bucket_num, device=device, dtype=torch.long),
            "sum": torch.zeros(k, bucket_num, feature_dim, device=device),
            "sumsq": torch.zeros(k, bucket_num, feature_dim, device=device),
        }

    @torch.no_grad()
    def accumulate_fds_moments(self, buffers, h, labels):
        if buffers is None or not self.use_fds:
            return
        h = h.detach().float()
        labels_fds = self._labels_to_fds_range(labels).detach().float()
        bucket_num = int(self.FDS[0].bucket_num)
        for k, j in enumerate(self.fds_active_markers):
            y = labels_fds[:, int(j)]
            idx = torch.clamp((y * bucket_num).long(), 0, bucket_num - 1)
            buffers["count"][k].index_add_(0, idx, torch.ones_like(idx, dtype=torch.long))
            buffers["sum"][k].index_add_(0, idx, h)
            buffers["sumsq"][k].index_add_(0, idx, h * h)

    @torch.no_grad()
    def finalize_fds_epoch(self, buffers, epoch):
        if buffers is None or not self.use_fds:
            return
        for k, j in enumerate(self.fds_active_markers):
            j = int(j)
            self.FDS[j].update_running_stats_from_moments(
                buffers["count"][k],
                buffers["sum"][k],
                buffers["sumsq"][k],
                epoch,
            )
            self.FDS[j].update_last_epoch_stats(epoch + 1)
