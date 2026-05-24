"""Multi-GPU inference for Hist2Cell on a prepared PyG Data file.

Splits the spots into N shards (N = number of GPUs) via mp.spawn. Each
worker binds to one GPU, runs NeighborLoader over its shard with
input_nodes=shard, and writes its partial predictions. Main process
merges shards into predictions.csv + predictions.npy.

Usage:
    python inference/infer.py \
        --data    inference/slide1_085_12/slide1_085_12.pt \
        --weights model_weights/humanlung_cell2location_leave_A50_out.pth \
        --output  inference/slide1_085_12
"""

import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.multiprocessing as mp
import torch.nn as nn
import torchvision.models as tv_models
from torch.nn import Linear
from torch_geometric.loader import NeighborLoader
from torch_geometric.nn import GATv2Conv, LayerNorm
import torch_geometric

torch_geometric.typing.WITH_PYG_LIB = False

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from model.ViT import Mlp, VisionTransformer  # noqa: E402


# ---------------- model ----------------

class Hist2Cell(nn.Module):
    """Re-declaration of the Hist2Cell architecture from
    tutorial_training/training_tutorial.ipynb cell 5. Architecture must
    stay identical for state_dict loading."""

    def __init__(self, cell_dim=80, vit_depth=3):
        super().__init__()
        backbone = tv_models.resnet18(weights=None)
        self.resnet18 = nn.Sequential(*list(backbone.children())[:-1])
        self.embed_dim = 32 * 8
        self.head = 8
        self.dropout = 0.3
        self.conv1 = GATv2Conv(
            in_channels=512,
            out_channels=int(self.embed_dim / self.head),
            heads=self.head,
        )
        self.norm1 = LayerNorm(in_channels=self.embed_dim)
        self.cell_transformer = VisionTransformer(
            num_classes=cell_dim,
            embed_dim=self.embed_dim,
            depth=vit_depth,
            mlp_head=True,
            drop_rate=self.dropout,
            attn_drop_rate=self.dropout,
        )
        self.spot_fc = Linear(in_features=512, out_features=256)
        self.spot_head = Mlp(in_features=256, hidden_features=512 * 2, out_features=cell_dim)
        self.local_head = Mlp(in_features=256, hidden_features=512 * 2, out_features=cell_dim)
        self.fused_head = Mlp(in_features=256, hidden_features=512 * 2, out_features=cell_dim)

    def forward(self, x, edge_index, return_features: bool = False):
        x_spot = self.resnet18(x).squeeze(-1).squeeze(-1)   # [N, 512]
        x_local = self.conv1(x=x_spot, edge_index=edge_index)
        x_local = self.norm1(x_local)
        x_local_seq = x_local.unsqueeze(0)                  # [1, N, 256]
        x_cell = x_local_seq

        x_spot_e = self.spot_fc(x_spot)                     # [N, 256]
        cp_spot = self.spot_head(x_spot_e)
        x_local_flat = x_local_seq.squeeze(0)
        cp_local = self.local_head(x_local_flat)

        cp_global, x_global = self.cell_transformer(x_cell)
        cp_global = cp_global.squeeze(0) if cp_global.dim() == 3 else cp_global
        x_global = x_global.squeeze(0) if x_global.dim() == 3 else x_global

        fused_repr = (x_spot_e + x_local_flat + x_global) / 3.0
        cp_fused = self.fused_head(fused_repr)
        cp = torch.relu((cp_spot + cp_local + cp_global + cp_fused) / 4.0)
        if return_features:
            return cp, {"resnet": x_spot, "fused": fused_repr}
        return cp


# ---------------- worker ----------------

def worker(rank: int, world_size: int, data_path: str, weight_path: str,
           batch_size: int, hop: int, out_dir: str):
    torch.cuda.set_device(rank)
    device = torch.device(f"cuda:{rank}")

    model = Hist2Cell(cell_dim=80, vit_depth=3).to(device).eval()
    sd = torch.load(weight_path, map_location=device)
    if isinstance(sd, dict) and "state_dict" in sd:
        sd = sd["state_dict"]
    missing, unexpected = model.load_state_dict(sd, strict=False)
    if rank == 0:
        if missing:
            print(f"[gpu0] WARN missing keys: {len(missing)} (first 3: {missing[:3]})", flush=True)
        if unexpected:
            print(f"[gpu0] WARN unexpected keys: {len(unexpected)} (first 3: {unexpected[:3]})", flush=True)

    data = torch.load(data_path, map_location="cpu", weights_only=False)
    if hasattr(data, "spot_id"):
        del data.spot_id  # python list w/ len==num_nodes confuses PyG's collate
    n = data.num_nodes
    shard = torch.arange(rank, n, world_size)
    if rank == 0:
        print(f"[gpu0] total nodes = {n}, per-gpu shard = {len(shard)}", flush=True)

    loader = NeighborLoader(
        data,
        num_neighbors=[-1] * hop,
        batch_size=batch_size,
        directed=False,
        shuffle=False,
        input_nodes=shard,
        num_workers=0,
    )

    n_my = len(shard)
    preds = np.zeros((n_my, 80), dtype=np.float32)
    feats_resnet = np.zeros((n_my, 512), dtype=np.float32)
    feats_fused = np.zeros((n_my, 256), dtype=np.float32)
    indices = np.zeros(n_my, dtype=np.int64)
    pos = 0
    t0 = time.time()
    with torch.no_grad():
        for sub in loader:
            sub = sub.to(device, non_blocking=True)
            cp, feats = model(sub.x, sub.edge_index, return_features=True)
            bs = int(sub.batch_size)
            preds[pos:pos + bs] = cp[:bs].cpu().numpy()
            feats_resnet[pos:pos + bs] = feats["resnet"][:bs].cpu().numpy()
            feats_fused[pos:pos + bs] = feats["fused"][:bs].cpu().numpy()
            # PyG 2.7 returns input_id as a *local* index into input_nodes,
            # not the global graph index — map back through the shard.
            indices[pos:pos + bs] = shard[sub.input_id.cpu()].numpy()
            pos += bs
            if rank == 0 and (pos // batch_size) % 50 == 0:
                rate = pos / (time.time() - t0 + 1e-6)
                print(f"[gpu0] {pos}/{n_my}  ({rate:.1f} spots/s)", flush=True)

    out_path = Path(out_dir) / f"shard_{rank}.npz"
    np.savez(out_path,
             indices=indices[:pos],
             preds=preds[:pos],
             feats_resnet=feats_resnet[:pos],
             feats_fused=feats_fused[:pos])
    if rank == 0:
        print(f"[gpu0] done {pos}/{n_my} in {time.time()-t0:.1f}s", flush=True)


# ---------------- main ----------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="path to prepared .pt (PyG Data)")
    ap.add_argument("--weights", required=True, help="path to Hist2Cell weights .pth")
    ap.add_argument("--output", required=True, help="output directory")
    ap.add_argument("--cell-types",
                    default="example_data/humanlung_cell2location/cell_types.pkl",
                    help="pickle with the list of 80 cell type names")
    ap.add_argument("--batch-size", type=int, default=16, help="NeighborLoader batch size per GPU")
    ap.add_argument("--hop", type=int, default=2, help="hop size for NeighborLoader")
    ap.add_argument("--world-size", type=int, default=torch.cuda.device_count(),
                    help="number of GPUs to use; default = all available")
    args = ap.parse_args()

    if args.world_size < 1:
        raise RuntimeError("No CUDA devices visible.")

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Launching multi-GPU inference (world_size={args.world_size})")
    print(f"  data:    {args.data}")
    print(f"  weights: {args.weights}")
    print(f"  output:  {out_dir}")

    mp.spawn(
        worker,
        args=(args.world_size, args.data, args.weights, args.batch_size, args.hop, str(out_dir)),
        nprocs=args.world_size,
        join=True,
    )

    print("Merging shards...")
    data = torch.load(args.data, map_location="cpu", weights_only=False)
    n = data.num_nodes
    final = np.zeros((n, 80), dtype=np.float32)
    final_resnet = np.zeros((n, 512), dtype=np.float32)
    final_fused = np.zeros((n, 256), dtype=np.float32)
    seen = np.zeros(n, dtype=bool)
    for r in range(args.world_size):
        sh_path = out_dir / f"shard_{r}.npz"
        sh = np.load(sh_path)
        final[sh["indices"]] = sh["preds"]
        final_resnet[sh["indices"]] = sh["feats_resnet"]
        final_fused[sh["indices"]] = sh["feats_fused"]
        seen[sh["indices"]] = True

    missing_n = int((~seen).sum())
    if missing_n:
        raise RuntimeError(f"{missing_n} spots missing predictions after merge")

    with open(args.cell_types, "rb") as f:
        cell_types = pickle.load(f)
    if len(cell_types) != 80:
        raise RuntimeError(f"cell_types pickle has {len(cell_types)} entries, expected 80")

    df = pd.DataFrame(final, columns=cell_types)
    df.insert(0, "spot_id", list(data.spot_id))
    df.insert(1, "X", data.pos[:, 0].numpy().astype(int))
    df.insert(2, "Y", data.pos[:, 1].numpy().astype(int))
    csv_path = out_dir / "predictions.csv"
    npy_path = out_dir / "predictions.npy"
    feat_resnet_path = out_dir / "features_resnet.npy"
    feat_fused_path = out_dir / "features_fused.npy"
    df.to_csv(csv_path, index=False)
    np.save(npy_path, final)
    np.save(feat_resnet_path, final_resnet)
    np.save(feat_fused_path, final_fused)

    for r in range(args.world_size):
        (out_dir / f"shard_{r}.npz").unlink()

    print(f"Saved:")
    print(f"  {csv_path}  ({n} spots × 80 cell types)")
    print(f"  {npy_path}")
    print(f"  {feat_resnet_path}  ({n} spots × 512 resnet18 features)")
    print(f"  {feat_fused_path}  ({n} spots × 256 pre-fused-head features)")


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    main()
