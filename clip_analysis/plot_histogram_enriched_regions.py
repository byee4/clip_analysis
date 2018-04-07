import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import parsers as p

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Default REGIONS to look at;
REGIONS = ['CDS', '3utr', '5utr']


def plot(
        l2fc_pval_enr,
        regions=REGIONS,
        xlabel='eCLIP log2 fold-enrichment',
        ylabel='fraction of REGIONS in bin',
        xlim=(-10, 10),
        ax=None,
        title="Enriched REGIONS"
):
    """
    Main plotting function.

    Plots a histogram of enriched REGIONS
    """
    if ax is None:
        ax = plt.gca()

    df = p.read_l2fcwithpval_enr(l2fc_pval_enr)
    bins = np.arange(0, 100 + 1, 1)
    for region in regions:
        n, bins = np.histogram(
            df['{} l2fc'.format(region)].dropna(),
            bins=100, range=xlim
        )
        pdf = [float(b) / sum(n) for b in n]  # sum all histogram bins to 1
        ax.bar(range(100), pdf, label=region, alpha=0.4)

        ## For some reason, getting the PDF doesn't work this way... ugh
        # sns.distplot(
        #     df['{} l2fc'.format(region)].dropna(),
        #     ax=ax, norm_hist=False,
        #     label=region # , hist_kws={'density':True}
        # )

    # set the x, ylabel
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(np.arange(0, 100 + 1, 10.0))  # set ticks every 10

    ax.set_xticklabels(bins[::10])  # label every 10
    ax.set_title(title)
    ax.legend()
    return ax


def make_plot(
    l2fc_pval_enr, out_file,
    regions, xlabel='eCLIP log2 fold-enrichment',
    ylabel='fraction of REGIONS in bin', xlim=(-10, 10)
):

    fig, ax = plt.subplots()
    plot(
        l2fc_pval_enr=l2fc_pval_enr,
        regions=regions,
        xlabel=xlabel,
        ylabel=ylabel,
        xlim=xlim,
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
        "--out_file",
        required=True
    )
    parser.add_argument(
        "--regions",
        help="Specify which REGIONS to plot (default: {})".format(
            REGIONS
        ),
        nargs='+',
        default=REGIONS,
        required=False
    )

    args = parser.parse_args()
    l2fcwithpval_enr = args.l2fcwithpval_enr
    out_file = args.out_file
    regions = args.regions

    make_plot(
        l2fcwithpval_enr, out_file,
        regions=regions, xlabel='eCLIP log2 fold-enrichment',
        ylabel='fraction of REGIONS in bin', xlim=(-10, 10)
    )

if __name__ == "__main__":
    main()