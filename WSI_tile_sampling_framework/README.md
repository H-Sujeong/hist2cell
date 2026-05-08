# Pathology WSI Tile Sampling System

A comprehensive pipeline for automated tile sampling from whole-slide images (WSI) in digital pathology applications, specifically designed for TCGA dataset processing.

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![OpenSlide](https://img.shields.io/badge/OpenSlide-supported-green.svg)](https://openslide.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🔬 Overview

This system automatically identifies tissue regions in pathology images and extracts meaningful tile coordinates for downstream deep learning applications. It combines advanced computer vision techniques with efficient processing pipelines to handle large-scale pathology datasets.

### Procedure Overview Figure
![Image](/Docs/tiling_overview.png)

### Key Features

- **🎯 Intelligent Tissue Detection**: Multi-modal foreground masking using YUV and HSV color spaces
- **⚡ High Performance**: Multiprocessing support for scalable tile sampling
- **🔧 Flexible Configuration**: Organ-specific parameter tuning (BRCA, ACC, Prostate, etc.)
- **📊 Quality Assurance**: Built-in visualization tools for result validation
- **💾 Efficient Storage**: HDF5 format for coordinate storage with metadata preservation
- **🖼️ Multiple Formats**: Support for SVS and other OpenSlide-compatible formats


## 📁 Repository Structure

```
├── ForegroundMasking.py      # Tissue detection and artifact removal
├── TileSampling.py          # Main WSI processing orchestration
├── Tiling.py                # Tile coordinate extraction and validation
├── tile_processing.py       # Command-line interface
├── tile_processed_visualize.py  # Visualization and QA tools
├── templates/               # Put stain normalization templates (.json)
├── README.md                # This file
└── examples/                # Usage examples and tutorials
```

## 🛠️ Core Components

![Image](/Docs/tiling_processing.png)

### WSITileSampler
Main orchestration class providing:
- WSI loading with automatic thumbnail generation
- Performance timing and monitoring
- Exception handling for problematic slides
- Integration of masking and tiling components
- save tile information to HDF5

### ForegroundMasker
Implements sophisticated tissue detection algorithms:
- **YUV Color Filtering**: Separate red, green, and blue tissue detection
- **Fat/Artifact Removal**: Intelligent filtering of non-tissue regions
- **Morphological Cleanup**: Noise reduction and hole filling
- **Background Removal**: HSV-based background elimination

### TileSampler
Advanced tile extraction with:
- **Boundary-based clustering**: Groups tiles by tissue regions
- **Multi-processing**: Parallel tile validation
- **Quality filtering**: Ensures meaningful tissue content
- **Coordinate management**: Full-resolution coordinate mapping

## 🚀 Quick Start

### Prerequisites

```bash
# Core dependencies
pip install openslide-python opencv-python numpy scikit-image
pip install h5py matplotlib pillow tqdm natsort

# System requirements
# - OpenSlide library installed on system
# - Python 3.7+
# - Sufficient RAM for WSI processing
```

### Basic Usage

```bash
# Process WSI files
python tile_processing.py \
    --root /path/to/wsi/files \
    --output_dir /path/to/output/coords \
    --endswith svs \
    --save_thumb True \
    --tile_size 224 \
    --min_tiles 5

# Visualize processed results
python tile_processed_visualize.py --organ ACC
```

> Thumbnails are automatically saved in the `output/Thumbnail/` directory when `--save_thumb True` is used.
>
> The parameter `--downscale_factor` is internally inferred from WSI metadata.

### Python API

`WSITileSampler` can process either:
- an entire directory containing multiple WSI files, or  
- a single WSI file directly.

```python
from TileSampling import WSITileSampler

# --- Case 1: Process all WSI files in a directory ---
sampler = WSITileSampler(
    root_dir='/path/to/wsi/files/',
    output_dir='/path/to/output',
    endswith='svs',
    save_thumb=True,
    tile_size=224,
    min_tiles=5,
    is_normalized=False
)

sampler.run()  # Automatically finds and processes all .svs files


# --- Case 2: Process a single WSI file ---
sampler = WSITileSampler(
    root_dir='/path/to/wsi/files/image1.svs',
    output_dir='/path/to/output',
    endswith='svs',
    save_thumb=True,
    tile_size=224,
    min_tiles=5,
    is_normalized=False
)

sampler.run()  # Processes only the specified WSI file
```

> `run()` automatically detects whether `root_dir` is a directory or a single file,
> and handles both cases accordingly.
>
> Thumbnails are saved in `output/Thumbnail/` when `save_thumb=True`.



## ⚙️ Configuration

### Organ-Specific Parameters

| Parameter | Description | BRCA | ACC | Prostate |
|-----------|-------------|------|-----|----------|
| `min_tiles` | Minimum tiles per boundary region | 50 | 5 | Variable |
| `tile_size` | Tile size in pixels | 256 | 256 | 256 |
| `overlap` | Tile overlap ratio | 0 | 0 | 0 |
| `green_thresh` | Fat region filtering threshold | 200 | 200 | Tunable |

### Command Line Arguments

Run the tile processing script with configurable options:

```bash
python tile_processing.py --help
```

#### Options

| Argument              | Type | Default | Description                                                                        |
| --------------------- | ---- | ------- | ---------------------------------------------------------------------------------- |
| `--root`              | str  | —       | Path to the directory containing Whole Slide Images (`.svs`) or a single WSI file. |
| `--output_dir`        | str  | —       | Directory where output HDF5 coordinate files will be saved.                        |
| `--endswith`, `--end` | str  | `"svs"` | File extension of WSI files.                                                       |
| `--max_depth`         | int  | `0`     | Maximum directory depth to search under the root directory.                        |
| `--save_thumb`        | bool | `False` | Whether to save thumbnails of WSIs (saved under `output/Thumbnail/`).              |
| `--tile_size`         | int  | `224`   | Patch size used for tiling WSIs.                                                   |
| `--min_tiles`         | int  | `5`     | Minimum number of tiles required in a valid region during post-filtering.          |
| `--is_normalized`     | bool | `False` | Whether to apply stain normalization (uses BRCA default template).                 |
| `--thumbnail_level`   | int  | `-1`    | Thumbnail resolution level (if available).                                         |
| `--downscale_factor`  | int  | `32`    | Downscaling factor used to generate thumbnails when none are provided.             |

#### Example

```bash
python tile_processing.py \
  --root /data/WSI \
  --output_dir /data/output \
  --endswith svs \
  --save_thumb True \
  --tile_size 224 \
  --min_tiles 5 \
  --downscale_factor 32
```

## 📊 Output Format

```bash
# Case 1) Directory input (process all .svs under /path/to/wsi/files)
/path/to/output
├── DB-00000-H.h5
├── DB-00001-A1.h5
├── DB-00002-K.h5
└── Thumbnail/
    ├── DB-00000-H.png
    ├── DB-00001-A1.png
    └── DB-00002-K.png
```

```bash
# Case 2) Single file input (process only image1.svs)
/path/to/output
├── image1.h5
└── Thumbnail/
    └── image1.png
```

> `--save_thumb True` 사용 시 썸네일이 `output/Thumbnail/`에 생성됩니다.

### HDF5 Coordinate Files
```python
# Structure of output .h5 files
{
    'tile_coords': [[x1, y1], [x2, y2], ...],  # Full-resolution coordinates
    'metadata': {
        'tile_size': 256,
        'overlap': 0,
        'total_tiles': 1847
    }
}
```

### Visualization Outputs
- **Thumbnails**: Low-resolution WSI overviews
- **Patch Overlays**: Red bounding boxes showing selected tile positions
- **Quality Assurance**: Visual validation of tissue detection

## 🔧 Advanced Usage

### Custom Tissue Detection

```python
from ForegroundMasking import ForegroundMasker

# Initialize with custom parameters
masker = ForegroundMasker(
    kernel_size=(5, 5),
    min_hole_area=1000,
    min_area=800,
    hole_area_threshold=1500
)

# Apply custom filtering
mask = masker.get_foreground(image)
```

### Batch Processing
```python
# Process multiple organs
organs = ['BRCA', 'ACC', 'LUAD']
for organ in organs:
    sampler = WSITileSampler(
        root_dir=f'/data/TCGA-{organ}/wsi',
        output_dir=f'/output/TCGA-{organ}/coords'
    )
    results = sampler.process_images(image_list)
```

## 📈 Performance Optimization

### Memory Management
- Process thumbnails instead of full-resolution images
- Stream processing for very large WSIs
- Efficient HDF5 storage format

### CPU Utilization
- Multiprocessing for tile validation
- Optimized OpenCV operations
- Parallel contour processing

### Recommended System Specs
- **RAM**: 16GB+ for large WSI collections
- **CPU**: Multi-core processor (8+ cores recommended)
- **Storage**: SSD for WSI files and outputs
- **GPU**: Optional for future morphological acceleration

## 🧪 Quality Assurance

### Validation Pipeline
1. **Visual Inspection**: Automatic thumbnail generation with tile overlays
2. **Coordinate Validation**: Boundary checking and edge case handling
3. **Metadata Verification**: Preservation of magnification and resolution info
4. **Statistical Analysis**: Tile count and distribution metrics

### Troubleshooting

| Issue | Solution |
|-------|----------|
| Low tile count | Adjust `min_tiles` parameter |
| Over-segmentation | Increase `green_thresh` for fat filtering |
| Missing tissue regions | Lower morphological kernel sizes |
| Memory errors | Reduce `downscale_factor` or process smaller batches |

## 🤝 Contributing

We welcome contributions! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Development Setup
```bash
# Clone repository
git clone https://github.com/your-org/pathology-tile-sampling.git
cd pathology-tile-sampling

# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
python -m pytest tests/
```

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@software{pathology_tile_sampling,
  title={Pathology WSI Tile Sampling System},
  author={Sungmin Lee},
  year={2025},
  url={https://github.com/CocoSungMin/Pathology-WSI-Tile-Sampling-System}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📋 Changelog

### v1.0.0 (2025-02-22)
- Initial release
- Multi-organ support (BRCA, ACC, Prostate)
- YUV/HSV-based tissue detection
- Multiprocessing tile sampling
- HDF5 coordinate storage
- Visualization tools


### v1.1.0 (2025-04-10)
- Added parameter `end` for file extension name  
- Added parameter `thumbnail_level`  
- Added parameter `downscaling_factor`  
- Generates thumbnail when not provided.

### v2.0.0 (2025-10-14)
- WSITileSampler refactor  
  - Integrated `save()` function into WSITileSampler  
  - Removed `original_magnification` (fixed as 0)  
  - Removed `thumbnail_dir`, added `save_thumb` (thumbnails saved in `output/Thumbnail/`)  
  - Renamed parameter `end` → `endswith`  
  - Added parameters `tile_size`, `min_tiles` (organ-dependent, requires further adjustment)  
  - Modified internal TileSampler `tiler` instance; now initialized within `load_wsi()` in `process_image()` (not returned)
  - `downscale_factor` automatically inferred from WSI metadata  
- TileSampler update  
  - Inherits parameters from WSITileSampler  
  - Collects tile coordinates from original WSI dimensions (not mask)  
  - Added pixel filtering (retains tiles with ≥20% valid tissue pixels)
  - Removed unused functions

## 🏥 Medical Disclaimer

This software is for research purposes only and is not intended for clinical diagnosis or treatment decisions. Always consult with qualified medical professionals for clinical applications.

---

**Made with ❤️ for the digital pathology community**
