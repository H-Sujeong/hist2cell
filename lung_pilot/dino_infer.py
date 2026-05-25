"""DINOv2 ViT-B/14 inference on lung_pilot 224-graph patches.

읽기: PyG Data(.pt) 의 data.x [N, 3, 224, 224] (ImageNet-norm).
출력: features_dinov2.npy [N, 768] (CLS token of dinov2_vitb14).

multi-GPU: nn.DataParallel (4 GPU 환경에서 자동 batch split).
xformers 가 .venv 에 없어도 vanilla attention fallback 으로 동작 (DINOv2 코드의 warn).

DINOv2 repo + 가중치는 외부 위치:
  /home/sjhong/dinov2                       (facebookresearch/dinov2 clone)
  /home/sjhong/dinov2_vitb14_pretrain.pth   (official .pth, 330 MB)

Usage:
    .venv/bin/python lung_pilot/dino_infer.py \\
        --data   lung_pilot/graph_output/224/<slide>.pt \\
        --output lung_pilot/dino_output/<slide>
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

_REPO_ROOT = Path(__file__).resolve().parent.parent
DINOV2_REPO_DEFAULT = str(_REPO_ROOT / "external" / "dinov2")
DINOV2_WEIGHTS_DEFAULT = "/home/sjhong/dinov2_vitb14_pretrain.pth"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data", required=True, help="path to prepared PyG .pt")
    ap.add_argument("--output", required=True, help="output directory")
    ap.add_argument("--weights", default=DINOV2_WEIGHTS_DEFAULT,
                    help="DINOv2 ViT-B/14 pretrain .pth")
    ap.add_argument("--dinov2-repo", default=DINOV2_REPO_DEFAULT,
                    help="local facebookresearch/dinov2 clone (for import)")
    ap.add_argument("--batch-size", type=int, default=256,
                    help="total batch (split across GPUs by DataParallel)")
    args = ap.parse_args()

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA required.")

    sys.path.insert(0, args.dinov2_repo)
    from hubconf import dinov2_vitb14  # noqa: E402

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"loading data: {args.data}")
    data = torch.load(args.data, map_location="cpu", weights_only=False)
    x_all = data.x  # [N, 3, 224, 224]
    n = x_all.size(0)
    print(f"  N={n}, x shape={tuple(x_all.shape)}, dtype={x_all.dtype}")

    print(f"loading dinov2 weights: {args.weights}")
    model = dinov2_vitb14(pretrained=False)
    sd = torch.load(args.weights, map_location="cpu", weights_only=False)
    missing, unexpected = model.load_state_dict(sd, strict=False)
    if missing:
        print(f"  WARN missing keys: {len(missing)} (first 3: {missing[:3]})")
    if unexpected:
        print(f"  WARN unexpected keys: {len(unexpected)} (first 3: {unexpected[:3]})")
    model = model.eval().cuda()

    n_gpu = torch.cuda.device_count()
    if n_gpu > 1:
        model = nn.DataParallel(model, device_ids=list(range(n_gpu)))
        print(f"  DataParallel on {n_gpu} GPUs (effective batch = {args.batch_size})")

    feats = np.zeros((n, 768), dtype=np.float32)
    t0 = time.time()
    with torch.no_grad():
        for i in range(0, n, args.batch_size):
            batch = x_all[i:i + args.batch_size].cuda(non_blocking=True)
            f = model(batch)  # [B, 768] (CLS token)
            feats[i:i + args.batch_size] = f.cpu().numpy()
            done = min(i + args.batch_size, n)
            if (i // args.batch_size) % 5 == 0 or done == n:
                rate = done / (time.time() - t0 + 1e-6)
                print(f"  {done}/{n}  ({rate:.1f} spots/s)", flush=True)
    elapsed = time.time() - t0
    print(f"done in {elapsed:.1f}s ({n / elapsed:.1f} spots/s)")

    out_npy = out_dir / "features_dinov2.npy"
    np.save(out_npy, feats)
    print(f"saved: {out_npy}  ({n} × 768)")


if __name__ == "__main__":
    main()
