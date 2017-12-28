import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

import matplotlib.image as mpimg
from argparse import ArgumentParser
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import parsers as p
import os
from collections import OrderedDict

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Set default color scheme for stuff
palette = sns.color_palette("hls", 8)

REGIONS = OrderedDict()
REGIONS["all"] = "All"
REGIONS["cds"] = "CDS"
REGIONS["three_prime_utrs"] = "3' UTR"
REGIONS["five_prime_utrs"] = "5' UTR"
REGIONS["proxintron500"] = "Proximal\nIntron"
REGIONS["distintron500"] = "Distal\nIntron"

def build_common_motifs(motif_grid, homer_location=None, regions=REGIONS):
    """

    Find a the top 8 common motifs in each region as determined by homer
    and places their images on the motif grid

    motif_grid - a grid section to plot the most common motifs onto
    homer_location - the location on the file system where homer is
    (for getting the data, should factor out)

    """


    if homer_location is None or regions is None:
        raise NotImplementedError(
            "Pickle file doesn't have data to generate this figure")

    for i, region in enumerate(regions.keys()):

        # make a gridspec for the top 8 motifs
        gs_homer_motifs = gridspec.GridSpecFromSubplotSpec(8, 1, subplot_spec=(
        motif_grid[i]))

        # for each top 8 motifs
        for j, gs in enumerate(gs_homer_motifs):

            # get it from where homer stored the output
            motif_logo = "motif" + str(j + 1) + ".logo.png"
            motif_pvalue = 'motif' + str(j + 1) + '.motif'
            motif_logo_file = os.path.join(homer_location, region,
                                           "homerResults", motif_logo)
            motif_pvalue_file = os.path.join(homer_location, region,
                                             "homerResults", motif_pvalue)
            if os.path.exists(motif_logo_file) and os.path.exists(
                    motif_pvalue_file):
                motif = mpimg.imread(motif_logo_file)
                pvalue = float(
                    open(motif_pvalue_file).readline().strip().split(",")[
                        -1].split(":")[-1])

                # print title only once
                if j == 0:
                    ax = plt.subplot(gs, frameon=False,
                                     xticks=[],
                                     yticks=[],
                                     title=regions[
                                               region] + "\n" + '{:.2e}'.format(
                                         pvalue))
                else:
                    ax = plt.subplot(gs, frameon=False, xticks=[], yticks=[],
                                     title='{:.2e}'.format(pvalue))
                ax.imshow(motif)
            else:
                print "no motif %s" % motif_logo_file



def make_plot(
    pickle1, pickle2, region, highlights, label, out_file
):
    fig, ax = plt.subplots()
    rep1_scores = p.read_kmer_enrichment_from_pickle(pickle1, region)
    rep2_scores = p.read_kmer_enrichment_from_pickle(pickle2, region)
    build_common_motifs(
        rep1_scores, rep2_scores, highlights,
        label=label, color=palette[4], highlight_color=palette[5],
        ax=ax
    )
    fig.savefig(out_file)

def main():
    parser = ArgumentParser()

    parser.add_argument(
        "--pickle1",
        required=True
    )
    parser.add_argument(
        "--pickle2",
        required=True
    )
    parser.add_argument(
        "--out_file",
        required=True
    )
    parser.add_argument(
        "--region",
        required=False,
        default='all'
    )
    parser.add_argument(
        "--highlights",
        required=False,
        nargs='+',
        default=[]
    )
    parser.add_argument(
        "--label",
        required=False,
        default='all 6mers'
    )
    args = parser.parse_args()
    pickle1 = args.pickle1
    pickle2 = args.pickle2
    region = args.region
    label = args.label
    highlights = args.highlights
    out_file = args.out_file

    make_plot(
        pickle1, pickle2, region, highlights, label, out_file
    )

if __name__ == "__main__":
    main()