#!/usr/bin/env python3
"""
Example script demonstrating MAFM web visualization usage.

This script shows how to:
1. Process fine-mapping results for web visualization
2. Launch the web interface programmatically
3. Handle different data configurations
"""

import os
import sys
from pathlib import Path

# Add the project root to Python path for imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

try:
    from mafm.web.app import run_app
    from mafm.web.export import export_for_web
except ImportError:
    print("Web dependencies not found. Please install with:")
    print("pip install mafm[web]")
    sys.exit(1)


def example_data_processing():
    """Example of processing data for web visualization."""
    # Example data directory structure
    data_dir = "/path/to/your/mafm/results"
    webdata_dir = "webdata"

    # Example loci files
    allmeta_loci = "data/real/meta/all/all_meta_loci_sig.txt"
    popumeta_loci = "data/real/meta/ancestry/loci_info_sig.txt"
    nometa_loci = "data/real/all_loci_list_sig.txt"

    print("Processing data for web visualization...")

    # Process data for web visualization
    try:
        export_for_web(
            data_base_dir=data_dir,
            webdata_dir=webdata_dir,
            allmeta_loci_file=allmeta_loci,
            popumeta_loci_file=popumeta_loci,
            nometa_loci_file=nometa_loci,
            threads=10,
        )
        print(f"Data processed successfully. Output in {webdata_dir}/")
    except Exception as e:
        print(f"Error processing data: {e}")
        return False

    return True


def example_web_launch():
    """Example of launching the web interface."""
    webdata_dir = "webdata"

    print("Launching web interface...")
    print("Access the interface at: http://localhost:8080")
    print("Press Ctrl+C to stop the server.")

    try:
        run_app(webdata_dir=webdata_dir, debug=True, port=8080, host="0.0.0.0")
    except KeyboardInterrupt:
        print("\nWeb server stopped.")
    except Exception as e:
        print(f"Error running web app: {e}")


def main():
    """Main example function."""
    print("MAFM Web Visualization Example")
    print("=" * 40)

    # Check if we should process data
    if len(sys.argv) > 1 and sys.argv[1] == "--process-data":
        success = example_data_processing()
        if not success:
            return

    # Check if webdata exists
    if not os.path.exists("webdata/all_loci_info.txt"):
        print("No processed web data found.")
        print("Run with --process-data flag to process data first, or")
        print("use the CLI command: mafm web /path/to/data")
        return

    # Launch web interface
    example_web_launch()


if __name__ == "__main__":
    main()
