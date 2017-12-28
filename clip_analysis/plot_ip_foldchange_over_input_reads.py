import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import parsers as p

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Set default color scheme for stuff
PALETTE = sns.color_palette("hls", 8)
REGIONS = ['CDS','3utr','5utr']

def plot(
    l2fcwithpval_enr, inp_reads_by_loc,
    ax=None,
    regions=REGIONS,
    colors=PALETTE,
    alpha=0.3, title="Fold change over input reads"
):
    """
    Plots the region-based analysis of genes enriched in IP over
     reads in size-matched input. This is the same figure used in Figure2b
     of Ann's IMP paper.

    Parameters
    ----------
    ip_l2fc : basestring
    inp_reads_by_loc : basestring
    field_list : dict
        dictionary of {region: color} to plot
    Returns
    -------
    :param regions:
    :param l2fcwithpval_enr:
    :param inp_reads_by_loc: 
    :param alpha: 
    :param ax: 

    """
    i = 0  # color counter
    df = p.scatter_matrix(l2fcwithpval_enr, inp_reads_by_loc)
    if ax is None:
        ax = plt.gca()
    for region in regions:
        if region in df.columns:
            ax.scatter(
                np.log2(df[region] + 1),
                df["{} l2fc".format(region)],
                color=colors[i],
                alpha=alpha
            )
        else:
            print("region {} not found in dataframe matrix".format(
                region
            ))
        i += 1
    ax.set_title("Region-based analysis of genes enriched")
    ax.set_xlabel("Reads in SMInput (log2)")
    ax.set_ylabel("Fold Enrichment (log2)")
    ax.set_title(title)
    ax.legend()
    return ax



def make_plot(
    l2fcwithpval_enr, inp_reads_by_loc,
    out_file, regions=REGIONS,
    alpha=0.3
):

    fig, ax = plt.subplots()
    plot(
        l2fcwithpval_enr, inp_reads_by_loc,
        regions=regions,
        alpha=alpha,
        ax=ax
    )
    fig.savefig(out_file)

def main():
    parser = ArgumentParser()

    parser.add_argument(
        "--l2fcwithpval_enr",
        required=True
    )
    parser.add_argument(
        "--inp_reads_by_loc",
        required=True
    )
    parser.add_argument(
        "--out_file",
        required=True
    )

    args = parser.parse_args()
    l2fcwithpval_enr = args.l2fcwithpval_enr
    inp_reads_by_loc = args.inp_reads_by_loc
    out_file = args.out_file

    make_plot(
        l2fcwithpval_enr, inp_reads_by_loc,
        out_file, regions=REGIONS,
        alpha=0.3
    )

if __name__ == "__main__":
    main()