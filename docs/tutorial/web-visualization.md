# Web Visualization Tutorial

MAFM provides an interactive web interface for exploring fine-mapping results. This tutorial covers how to install, configure, and use the web visualization features.

## Prerequisites

Before using the web interface, ensure you have:

1. **Completed fine-mapping analysis** with MAFM
2. **Generated results** in the expected directory structure
3. **Loci information files** (at least one of: allmeta, popumeta, or nometa)

## Installation

### Install Web Dependencies

The web visualization requires additional dependencies that are not included in the base MAFM installation:

```bash
# Install MAFM with web dependencies
pip install mafm[web]

# Or install dependencies separately
pip install dash dash-bootstrap-components dash-mantine-components plotly
```

### Verify Installation

Check that the web command is available:

```bash
mafm web --help
```

## Quick Start

### Basic Usage

The simplest way to start the web interface:

```bash
# Navigate to your MAFM results directory
cd /path/to/your/mafm/results

# Launch web interface (uses current directory)
mafm web
```

This will:
1. Look for default loci files in standard locations
2. Process data for web visualization (if needed)
3. Start the web server on `http://localhost:8080`

### With Custom Data Directory

```bash
# Specify a different data directory
mafm web /path/to/mafm/results --port 8080
```

## Command Options

### Required Arguments

- `data_dir`: Base directory containing fine-mapping data (default: current directory)

### Optional Parameters

#### Data Configuration
- `--webdata-dir, -w`: Directory for processed web data (default: "webdata")
- `--allmeta-loci, -a`: Path to allmeta loci info file
- `--popumeta-loci, -p`: Path to popumeta loci info file
- `--nometa-loci, -n`: Path to nometa loci info file
- `--force-regenerate, -f`: Force regeneration of web data

#### Server Configuration
- `--port`: Port to run web server on (default: 8080)
- `--host`: Host to bind web server to (default: "0.0.0.0")
- `--debug`: Run in debug mode

#### Processing Configuration
- `--threads, -t`: Number of threads for data processing (default: 10)

## Data Structure

### Expected Directory Layout

Your MAFM results should follow this structure:

```
mafm_results/
├── data/
│   └── real/
│       ├── credset/           # Fine-mapping results
│       │   ├── allmeta/
│       │   ├── popumeta/
│       │   └── nometa/
│       └── qc/                # Quality control results
│           ├── allmeta/
│           ├── popumeta/
│           └── nometa/
├── loci_files/
│   ├── all_meta_loci_sig.txt  # Allmeta loci info
│   ├── loci_info_sig.txt      # Popumeta loci info
│   └── all_loci_list_sig.txt  # Nometa loci info
└── webdata/                   # Generated web data (created automatically)
```

### Loci Information Files

Each loci file should contain columns like:
- `locus_id`: Unique identifier for each locus
- `chr`: Chromosome
- `start`: Start position
- `end`: End position
- `prefix`: File prefix for sumstats and LD data
- `popu`: Population code
- `cohort`: Cohort identifier
- `sample_size`: Sample size

## Usage Examples

### Example 1: Basic Launch

```bash
# Simple launch from results directory
cd /path/to/mafm/results
mafm web
```

### Example 2: Custom Configuration

```bash
# Specify custom loci files and settings
mafm web /path/to/data \
  --allmeta-loci /path/to/allmeta_loci.txt \
  --popumeta-loci /path/to/popumeta_loci.txt \
  --nometa-loci /path/to/nometa_loci.txt \
  --webdata-dir custom_webdata \
  --port 8080 \
  --threads 20
```

### Example 3: Force Data Regeneration

```bash
# Force regeneration of web data
mafm web /path/to/data --force-regenerate
```

### Example 4: Debug Mode

```bash
# Run in debug mode for development
mafm web /path/to/data --debug --port 8081
```

## Web Interface Features

### Home Page

The main dashboard provides:

#### Interactive Filtering
- **Meta-analysis Method**: Choose from allmeta, popumeta, or nometa
- **Fine-mapping Tool**: Select specific tools (SuSiE, FINEMAP, etc.)
- **Plot Type**: View different metrics (credible sets, PIP statistics)

#### Visualizations
- **Bar Plot**: Interactive plots showing statistics per locus
- **Summary Table**: Sortable table with locus information and clickable links

#### Key Metrics
- Number of credible sets per locus
- Credible set sizes
- Number of SNPs with high PIP values
- SNP counts by p-value thresholds

### Locus Pages

Detailed views for individual loci include:

#### Association Plots
- **LocusZoom Plot**: Manhattan plot with LD coloring
- **Regional Association**: SNP-level fine-mapping results

#### Quality Control
- **QC Metrics Table**: Lambda, DENTIST-S, MAF correlation
- **LD Diagnostics**: Fourth moment, decay patterns
- **Expected vs Observed**: Z-score distributions

#### Fine-mapping Results
- **Multiple Tool Comparison**: Compare results across methods
- **Credible Set Visualization**: Highlight credible variants
- **PIP Distributions**: Posterior inclusion probabilities

## Data Processing

### Automatic Processing

When you run `mafm web`, it automatically:

1. **Checks for existing web data** in the webdata directory
2. **Processes raw results** if web data doesn't exist or is outdated
3. **Generates summary files** for the web interface
4. **Optimizes data format** for fast web access

### Manual Processing

You can also process data separately:

```python
from mafm.web.export import export_for_web

export_for_web(
    data_base_dir="/path/to/mafm/results",
    webdata_dir="webdata",
    allmeta_loci_file="allmeta_loci.txt",
    popumeta_loci_file="popumeta_loci.txt",
    nometa_loci_file="nometa_loci.txt",
    threads=10
)
```

## Advanced Configuration

### Custom Web Data Directory

```bash
# Use custom directory for processed web data
mafm web /path/to/data --webdata-dir /custom/path/webdata
```

### Multiple Meta-analysis Types

You can include results from different meta-analysis approaches:

```bash
mafm web /path/to/data \
  --allmeta-loci allmeta_results.txt \      # All-ancestry meta
  --popumeta-loci popumeta_results.txt \    # Population-specific meta
  --nometa-loci nometa_results.txt          # No meta-analysis
```

### Performance Tuning

For large datasets:

```bash
# Increase processing threads
mafm web /path/to/data --threads 30

# Use SSD storage for webdata
mafm web /path/to/data --webdata-dir /fast/ssd/webdata
```

## Troubleshooting

### Common Issues

#### Missing Dependencies
```bash
# Error: "Web dependencies not found"
pip install mafm[web]
```

#### No Data Found
```bash
# Error: "No loci files found"
# Solution: Specify loci files explicitly
mafm web /path/to/data --allmeta-loci /path/to/loci.txt
```

#### Port Already in Use
```bash
# Error: "Port 8080 is already in use"
# Solution: Use different port
mafm web /path/to/data --port 8081
```

#### Out of Memory
```bash
# For large datasets, reduce threads
mafm web /path/to/data --threads 5
```

### Debug Mode

Enable debug mode for troubleshooting:

```bash
mafm web /path/to/data --debug
```

This provides:
- Detailed error messages
- Live code reloading
- Enhanced logging

### Log Files

Check MAFM logs for processing issues:
- Processing errors during data export
- File format issues
- Memory usage problems

## Best Practices

### Data Organization

1. **Consistent naming**: Use consistent file naming conventions
2. **Complete data**: Ensure all required files are present
3. **Valid formats**: Check that input files have correct formats

### Performance

1. **SSD storage**: Use SSD for webdata directory
2. **Adequate memory**: Ensure sufficient RAM for large datasets
3. **Network access**: Use appropriate host/port settings

### Security

1. **Internal networks**: Bind to localhost for local use only
2. **Firewall rules**: Configure appropriate network access
3. **Data sensitivity**: Consider data privacy when sharing URLs

## Integration with Workflows

### Pipeline Integration

Integrate web visualization into your analysis pipeline:

```bash
#!/bin/bash
# Complete analysis pipeline

# 1. Run fine-mapping
mafm pipeline input_loci.txt results/

# 2. Launch web interface
mafm web results/ --port 8080
```

### Batch Processing

For multiple datasets:

```bash
# Process multiple result sets
for dataset in dataset1 dataset2 dataset3; do
    mafm web results/${dataset} --webdata-dir webdata/${dataset} --port 808${i}
    ((i++))
done
```

## Next Steps

- Explore the **home page** to get an overview of your results
- Click on **individual loci** for detailed views
- Use **filtering options** to focus on specific analyses
- Export **plots and tables** for presentations
- Integrate with your **analysis workflow**

For more advanced usage, see:
- [Advanced Tutorial](advanced.md)
- [API Documentation](../api.md)
- [Contributing Guide](../contributing.md) 