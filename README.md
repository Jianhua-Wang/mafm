# mafm


[![pypi](https://img.shields.io/pypi/v/mafm.svg)](https://pypi.org/project/mafm/)
[![python](https://img.shields.io/pypi/pyversions/mafm.svg)](https://pypi.org/project/mafm/)
[![Build Status](https://github.com/Jianhua-Wang/mafm/actions/workflows/dev.yml/badge.svg)](https://github.com/Jianhua-Wang/mafm/actions/workflows/dev.yml)
[![codecov](https://codecov.io/gh/Jianhua-Wang/mafm/branch/main/graphs/badge.svg)](https://codecov.io/github/Jianhua-Wang/mafm)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



Multi-ancestry fine-mapping pipeline with interactive web visualization.


* Documentation: <https://Jianhua-Wang.github.io/mafm>
* GitHub: <https://github.com/Jianhua-Wang/mafm>
* PyPI: <https://pypi.org/project/mafm/>
* Free software: MIT


## Features

- **Multi-ancestry fine-mapping**: Support for multiple fine-mapping tools (SuSiE, FINEMAP, etc.)
- **Meta-analysis capabilities**: Combine results across populations and cohorts
- **Quality control**: Built-in QC metrics and visualizations
- **Interactive web interface**: Explore results through a modern web dashboard
- **Command-line interface**: Easy-to-use CLI for all operations

## Installation

### Basic Installation
```bash
pip install mafm
```

### With Web Visualization
To use the interactive web interface, install with web dependencies:
```bash
pip install mafm[web]
```

## Quick Start

### Command Line Usage

```bash
# Run complete fine-mapping pipeline
mafm pipeline input_loci.txt output_dir/

# Launch web visualization interface
mafm web /path/to/results --port 8080

# View specific loci files
mafm web /path/to/data \
  --allmeta-loci data/allmeta_loci.txt \
  --popumeta-loci data/popumeta_loci.txt \
  --nometa-loci data/nometa_loci.txt
```

### Web Interface

The web interface provides:
- **Home page**: Overview of all loci with interactive filtering
- **Locus pages**: Detailed views with LocusZoom-style plots
- **Quality control**: Comprehensive QC metrics and visualizations
- **Multi-tool comparison**: Compare results across different fine-mapping methods

Access the web interface at `http://localhost:8080` after running `mafm web`.

## Documentation

For detailed documentation, see <https://Jianhua-Wang.github.io/mafm>

## Web Visualization

The web module (`mafm.web`) provides interactive visualization of fine-mapping results. See [mafm/web/README.md](mafm/web/README.md) for detailed usage instructions.
