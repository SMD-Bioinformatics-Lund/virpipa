#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib.ticker import NullFormatter
import sys
import os

def create_rug_jitter_plot(filename, custom_title=None):
    """
    Create a rug plot with a jitter plot of data from the given file.
    The plot title is either provided or derived from the filename.
    """
    # Use custom title if provided, otherwise extract from filename
    if custom_title:
        title = custom_title
    else:
        title = os.path.splitext(os.path.basename(filename))[0]

    # Read data from file
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse data
    data = [float(line.strip()) for line in lines]

    # Create figure and primary axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Create the KDE plot with rug on the primary axis
    sns.kdeplot(data, ax=ax1, bw_adjust=0.5, fill=True, alpha=0.4, color='#3498db', linewidth=2)
    sns.rugplot(data, ax=ax1, color='black', alpha=0.2, height=0.05)

    # Set up the primary axis
    ax1.set_xlabel('Value')
    ax1.set_ylabel('Density')
    ax1.set_xlim(0.4, 1)
    ax1.set_title(title, fontsize=14, fontweight='bold')

    # Create secondary axis for jitter plot
    ax2 = ax1.twinx()

    # Generate jittered y-positions for the scatter plot (between 0.9 and 0.99)
    jitter_y = np.random.uniform(0.9, 0.99, len(data))

    # Create the jitter plot on the secondary axis
    ax2.scatter(data, jitter_y, color='black', alpha=0.3, s=3, zorder=5)

    # Set up the secondary axis
    ax2.set_ylabel('Jitter')
    ax2.set_ylim(0, 1)
    ax2.yaxis.set_major_formatter(NullFormatter())  # Hide the tick labels
    ax2.grid(False)  # Turn off grid for the secondary axis

    # Improve the layout
    plt.tight_layout()

    # Save and show the plot
    # We use the 'title' variable for the filename to keep it consistent with the display
    clean_filename = "".join([c for c in title if c.isalnum() or c in (' ', '_', '-')]).strip().replace(' ', '_')
    output_filename = f"{clean_filename}_rug_kde_plot.png"
    
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_filename}")

if __name__ == "__main__":
    # Check if filename was provided as argument
    if len(sys.argv) < 2:
        print("Usage: python script.py <data_filename> [optional_title]")
        sys.exit(1)

    # Get filename from command line argument
    filename = sys.argv[1]
    
    # Check for an optional title in the second argument
    provided_title = sys.argv[2] if len(sys.argv) > 2 else None

    # Check if file exists
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    create_rug_jitter_plot(filename, provided_title)
