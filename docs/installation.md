# Installation

## Requirements

MAFM requires Python 3.9 or higher. The base installation includes all dependencies needed for fine-mapping analysis.

## Installation Options

### Standard Installation

To install the base MAFM package, run this command in your terminal:

```console
$ pip install mafm
```

This is the preferred method to install MAFM, as it will always install the most recent stable release.

### Web Visualization

To use the interactive web visualization features, install with web dependencies:

```console
$ pip install mafm[web]
```

This installs additional packages required for the web interface:
- `dash` - Web application framework
- `dash-bootstrap-components` - Bootstrap components for Dash
- `dash-mantine-components` - Mantine components for Dash  
- `plotly` - Interactive plotting library

### Development Installation

For development or to get the latest features:

```console
$ git clone https://github.com/Jianhua-Wang/mafm.git
$ cd mafm
$ pip install -e .[web]
```

## Verify Installation

Check that MAFM is installed correctly:

```console
$ mafm --help
```

For web visualization specifically:

```console
$ mafm web --help
```

If you see the help output, the installation was successful!

## Alternative Installation Methods

### Using conda

If you prefer conda:

```console
$ conda install -c conda-forge mafm
```

Note: Web dependencies may need to be installed separately with conda.

### From source

The source for MAFM can be downloaded from the [Github repo][].

You can either clone the public repository:

```console
$ git clone git://github.com/Jianhua-Wang/mafm
```

Or download the [tarball][]:

```console
$ curl -OJL https://github.com/Jianhua-Wang/mafm/tarball/master
```

Once you have a copy of the source, you can install it with:

```console
$ pip install .
```

For development with web features:

```console
$ pip install -e .[web]
```

## Troubleshooting

### Missing Dependencies

If you get import errors when using web features:

```console
$ pip install mafm[web]
```

### Permission Issues

If you get permission errors, try installing with `--user`:

```console
$ pip install --user mafm[web]
```

### Python Version

Ensure you're using Python 3.9 or higher:

```console
$ python --version
```

If you don't have [pip][] installed, this [Python installation guide][] can guide you through the process.

## Next Steps

After installation:

1. Check out the [Quick Start Guide](tutorial/quick-start.md)
2. Try the [Web Visualization Tutorial](tutorial/web-visualization.md)
3. See [Usage Examples](usage.md)

  [pip]: https://pip.pypa.io
  [Python installation guide]: http://docs.python-guide.org/en/latest/starting/installation/
  [Github repo]: https://github.com/Jianhua-Wang/mafm
  [tarball]: https://github.com/Jianhua-Wang/mafm/tarball/master
