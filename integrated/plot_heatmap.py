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

def make_plot(
        in_files, in_labels, l2fc, padj, out_file, out_tsv,
        gene_label_limit=100, fig_size=(15, 10)
):
    # Run once to get a subset of significant genes
    df = p.read_and_join_deseq2(
        in_files, in_labels, l2fc, padj
    )
    names_to_keep = df.index

    # Run again to get all fold changes of these genes
    df = p.read_and_join_deseq2(
        in_files, in_labels, 0, 1
    )
    df = df.ix[names_to_keep]

    # Set a limit to display gene names, too many genes
    # results in muddled plot.
    if df.shape[0] > gene_label_limit:
        set_yticklabels = False
    else:
        set_yticklabels = True

    # Save tabbed heatmap values
    df.to_csv(out_tsv, sep='\t')

    # Plot the heatmap.
    cg = sns.clustermap(
        df, xticklabels=True, yticklabels=set_yticklabels, figsize=fig_size,
        metric='euclidean',
        method='average', edgecolors='white', linewidths=0.001,
        cmap='bwr', vmin=-6, vmax=6
    )
    plt.title('Log2 FC')
    ticks = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    cg.savefig(out_file)
    return 0


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--in_files",
        required=True,
        nargs='+'
    )
    parser.add_argument(
        "--l2fc",
        required=False,
        default=1.5
    )
    parser.add_argument(
        "--padj",
        required=False,
        default=0.05
    )
    parser.add_argument(
        "--in_labels",
        help="Specify labels (by default will "
             "use the file basenames)",
        nargs='+',
        default=None,
        required=False
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )
    parser.add_argument(
        "--out_tsv",
        required=False,
        default=None
    )
    # Process arguments
    args = parser.parse_args()
    in_files = args.in_files
    in_labels = args.in_labels
    out_file = args.out_file
    l2fc = float(args.l2fc)
    padj = float(args.padj)
    out_tsv = args.out_tsv

    if in_labels is None:
        in_labels = [os.path.basename(fn) for fn in in_files]

    if out_tsv is None:
        out_tsv = out_file + '.tsv'

    # main func
    make_plot(
        in_files, in_labels, l2fc, padj, out_file, out_tsv
    )


if __name__ == "__main__":
    main()
