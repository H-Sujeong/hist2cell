"""
TCGA-LUAD WSI tiling driver.

WSI_tile_sampling_framework 의 WSITileSampler 를 그대로 사용하되:
  - tile_processing.py CLI 의 positional-arg 오정렬 버그를 피하려고
    WSITileSampler 를 명시적 kwargs 로 생성하고 서브메소드를 직접 호출한다.
  - 썸네일은 framework 의 matplotlib savefig 경로(작게 리스케일됨)를 쓰지 않고
    openslide thumbnail 배열을 그대로 저장한다 (Masks/Overlays 와 해상도 일치).

사용:
  python run_tiling_tcga_luad.py                              # tile_size 224 (기본)
  python run_tiling_tcga_luad.py --tile-size 112 --output ... # 임의 크기

산출물 (per slide), --output 디렉터리 안:
  <name>.h5                     full-res tile 좌표 + metadata (framework HDF5)
  Thumbnails/<basename>.png     슬라이드 썸네일 (실제 배열 해상도)
  Masks/<name>_tissue_mask.png  조직(foreground) mask
  Overlays/<name>_tiles.png     썸네일 위 tile 위치 overlay (QA)
"""
import argparse
import os
import sys
import time

import cv2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from TileSampling import WSITileSampler  # noqa: E402

WSI_DIR = "/mnt/fileserver/NAS2_pathology/Pathology_project/TCGA-LUAD/wsi"
WSI_FILES = [
    "TCGA-05-4245-01A-01-BS1.41d3cf23-4e36-4e42-9e08-adfea139f37e.svs",
    "TCGA-05-4245-01A-01-TS1.bf71c76b-e802-4a7a-b6c3-c5f46212fab0.svs",
    "TCGA-05-4390-01A-01-BS1.38f2a7ef-442a-4fa6-acad-6e5d567bdcfd.svs",
]
DEFAULT_OUT = "/home/sjhong/hist2cell/lung_pilot/tilitng_output/TCGA-LUAD"


def short_name(fname):
    """TCGA barcode 부분만 (UUID 앞)."""
    return os.path.basename(fname).split(".")[0]


def run(tile_size, out, overlap=0, min_tiles=5):
    for sub in ("", "Thumbnails", "Masks", "Overlays"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    # save_thumb=False — framework 의 matplotlib savefig 썸네일을 막고,
    # 아래에서 thumbnail 배열을 직접 저장한다.
    sampler = WSITileSampler(
        root=os.path.join(WSI_DIR, WSI_FILES[0]),
        output_dir=out,
        endswith="svs",
        save_thumb=False,
        tile_size=tile_size,
        overlap=overlap,
        min_tiles=min_tiles,
        is_normalized=False,
    )

    rows = []
    for fn in WSI_FILES:
        path = os.path.join(WSI_DIR, fn)
        name = short_name(fn)
        print(f"\n==================== {name}  (tile_size={tile_size}) ====================")
        t0 = time.time()

        # 1) load WSI + thumbnail 배열
        thumbnail, scaler = sampler.load_wsi(path)

        # 2) 썸네일을 실제 배열 해상도로 직접 저장 (matplotlib 우회)
        thumb_path = os.path.join(out, "Thumbnails", os.path.splitext(fn)[0] + ".png")
        cv2.imwrite(thumb_path, cv2.cvtColor(thumbnail, cv2.COLOR_RGB2BGR))

        # 3) foreground (tissue) mask
        _, mask = sampler.compute_foreground_mask(thumbnail)

        # 4) tile 좌표 sampling (full-res top-left x,y)
        coords = sampler.sample_tiles(mask, scaler)

        # 5) framework HDF5 저장 (coords + metadata)
        sampler.save_hdf5(name, coords)

        # 6) 조직 mask PNG 저장
        cv2.imwrite(os.path.join(out, "Masks", f"{name}_tissue_mask.png"), mask)

        # 7) tile overlay QA 이미지
        overlay = thumbnail.copy()
        ts_thumb = max(1, round(tile_size / scaler))
        for (x, y) in coords:
            xs, ys = round(x / scaler), round(y / scaler)
            cv2.rectangle(overlay, (xs, ys), (xs + ts_thumb, ys + ts_thumb),
                          (255, 0, 0), 1)
        cv2.imwrite(os.path.join(out, "Overlays", f"{name}_tiles.png"),
                    cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR))

        dt = time.time() - t0
        tissue_frac = float((mask > 0).mean())
        print(f"  tiles={len(coords)}  thumb={thumbnail.shape}  scaler={scaler}  "
              f"tissue={tissue_frac:.1%}  {dt:.1f}s")
        rows.append((name, len(coords), thumbnail.shape, scaler, tissue_frac, dt))

    print(f"\n==================== SUMMARY (tile_size={tile_size}) ====================")
    for name, n, shp, sc, tf, dt in rows:
        print(f"  {name:26s}  tiles={n:7d}  scaler={sc:>3}  "
              f"thumb={shp[1]}x{shp[0]}  tissue={tf:5.1%}  {dt:6.1f}s")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile-size", type=int, default=224,
                    help="tile 한 변 (level-0 px). stride = tile_size (overlap 0)")
    ap.add_argument("--output", default=DEFAULT_OUT, help="출력 디렉터리")
    ap.add_argument("--min-tiles", type=int, default=5,
                    help="contour cluster 당 최소 tile 수 (미만이면 제거)")
    args = ap.parse_args()
    run(args.tile_size, args.output, min_tiles=args.min_tiles)


if __name__ == "__main__":
    main()
