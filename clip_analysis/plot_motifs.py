import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

import matplotlib.pyplot as plt
from itertools import izip
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import os
import argparse
import pybedtools
import parsers as p
### do this to avoid making tons of dots in svg files:


rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Default REGIONS to look at;
REGIONS = ['CDS', '3utr', '5utr', 'intron']


# Plotting functions #

def plot(
        in_file,
        ax=None, title=""
):
    if ax is None:
        ax = plt.gca()

    ax.set_title(title)
    ax.legend()
    return ax

def make_plot(in_file, out_file):
    fig, ax = plt.subplots()
    plot(
        in_file,
        title="title",
        ax=ax
    )
    fig.savefig(out_file)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--in_file",
        required=True,
    )
    parser.add_argument(
        "--REGIONS",
        help="Specify which REGIONS to plot (default: {})".format(
            REGIONS
        ),
        nargs='+',
        default=None,
        required=False
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )

    # Process arguments
    args = parser.parse_args()
    in_file = args.in_file
    out_file = args.out_file

    # main func
    make_plot(
        in_file, out_file
    )


if __name__ == "__main__":
    main()
