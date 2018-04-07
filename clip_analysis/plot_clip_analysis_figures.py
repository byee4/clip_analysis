import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import parsers as p
import plot_compare_rbp_enriched_regions
import plot_histogram_enriched_regions
import plot_kmer_enrichment
import plot_region_distribution
import plot_ip_foldchange_over_input_reads
import plot_repetitive_elements_bar
import seaborn as sns

### do this to avoid making tons of dots in svg files:
from matplotlib import rc

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Set default color scheme for stuff
palette = sns.color_palette("hls", 8)


def plot_all(
    l2fcwithpval_enr_r1, l2fcwithpval_enr_r2,
    inp_reads_by_loc_r1, inp_reads_by_loc_r2,
    pickle_r1, pickle_r2,
    rep_element_parsed_r1_ip, rep_element_parsed_r1_input,
    rep_element_parsed_r2_ip, rep_element_parsed_r2_input,
    motifs, regions, out_file, annotated_files, annotation_source,
    suptitle
):

    nrows = 5
    ncols = 3

    full_grid = gridspec.GridSpec(
        nrows, ncols,
        height_ratios=[1 for i in range(nrows)],
        width_ratios=[1 for i in range(ncols)],
        hspace=0.5, wspace=2
    )
    fig = plt.figure(figsize=(20, 20))

    map_rows = []

    for row in range(0, nrows):
        map_rows.append(gridspec.GridSpecFromSubplotSpec(
            1, ncols,
            subplot_spec=full_grid[row, :])
        )

    if l2fcwithpval_enr_r1 is not None:
        plot_histogram_enriched_regions.plot(
            l2fcwithpval_enr_r1, ax=plt.subplot(map_rows[0][0]),
            title='Rep 1 enriched REGIONS', regions=regions
        )

    if l2fcwithpval_enr_r2 is not None:
        plot_histogram_enriched_regions.plot(
            l2fcwithpval_enr_r2, ax=plt.subplot(map_rows[0][1]),
            title='Rep 2 enriched REGIONS', regions=regions
        )

    if l2fcwithpval_enr_r1 is not None and inp_reads_by_loc_r1 is not None:
        plot_ip_foldchange_over_input_reads.plot(
            l2fcwithpval_enr_r1, inp_reads_by_loc_r1,
            ax=plt.subplot(map_rows[1][0]),
            title='Rep 1 fold change over input reads',
            regions=regions
        )

    if l2fcwithpval_enr_r2 is not None and inp_reads_by_loc_r2 is not None:
        plot_ip_foldchange_over_input_reads.plot(
            l2fcwithpval_enr_r2, inp_reads_by_loc_r2,
            ax=plt.subplot(map_rows[1][1]),
            title='Rep 2 fold change over input reads',
            regions=regions
        )

    if l2fcwithpval_enr_r1 is not None and l2fcwithpval_enr_r2 is not None:
        plot_compare_rbp_enriched_regions.plot(
            l2fcwithpval_enr_r1, l2fcwithpval_enr_r2,
            ax=plt.subplot(map_rows[2][0]),
            title='Per-gene enrichment correlations',
            regions=regions
        )

    if annotated_files is not None:
        counts = p.get_counts(annotated_files, src=annotation_source)
        plot_region_distribution.plot(
            counts, ax=plt.subplot(map_rows[2][1]),
            title='Fraction of Peaks among RBPs'
        )

    if pickle_r1 is not None and pickle_r2 is not None:
        zscores_all_r1 = p.read_kmer_enrichment_from_pickle(
            pickle_r1, 'all', col_name='Rep1'
        )
        zscores_all_r2 = p.read_kmer_enrichment_from_pickle(
            pickle_r2, 'all', col_name='Rep2'
        )
        plot_kmer_enrichment.plot_zscores(
            zscores_all_r1, zscores_all_r2, label='all 6mers',
            ax=plt.subplot(map_rows[3][0]),
            highlights=motifs
        )

    if rep_element_parsed_r1_ip is not None and rep_element_parsed_r1_input is not None:
        plot_repetitive_elements_bar.plot(
            rep_element_parsed_r1_ip, rep_element_parsed_r1_input,
            ax=plt.subplot(map_rows[4][0])
        )
    if rep_element_parsed_r2_ip is not None and rep_element_parsed_r2_input is not None:
        plot_repetitive_elements_bar.plot(
            rep_element_parsed_r2_ip, rep_element_parsed_r2_input,
            ax=plt.subplot(map_rows[4][1])
        )
    try:
        plt.tight_layout(pad=5)
    except ValueError:
        pass
    fig.suptitle(suptitle)
    fig.savefig(out_file)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--l2fcwithpval_enr_r1",
        required=False,
        default=None
    )
    parser.add_argument(
        "--l2fcwithpval_enr_r2",
        required=False,
        default=None
    )
    parser.add_argument(
        "--inp_reads_by_loc_r1",
        required=False,
        default=None
    )
    parser.add_argument(
        "--inp_reads_by_loc_r2",
        required=False,
        default=None
    )
    parser.add_argument(
        "--pickle_r1",
        required=False,
        default=None
    )
    parser.add_argument(
        "--pickle_r2",
        required=False,
        default=None
    )
    parser.add_argument(
        "--rep_element_parsed_r1_ip",
        required=False,
        default=None
    )
    parser.add_argument(
        "--rep_element_parsed_r1_input",
        required=False,
        default=None
    )
    parser.add_argument(
        "--rep_element_parsed_r2_ip",
        required=False,
        default=None
    )
    parser.add_argument(
        "--rep_element_parsed_r2_input",
        required=False,
        default=None
    )
    parser.add_argument(
        "--annotated_files",
        required=False,
        nargs='+',
        default=None
    )
    parser.add_argument(
        "--annotation_script",
        required=False,
        default='brian'
    )
    parser.add_argument(
        "--motifs",
        required=False,
        nargs='+',
        default=[]
    )
    parser.add_argument(
        "--suptitle",
        required=False,
        default=''
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )
    parser.add_argument(
        "--regions",
        required=False,
        nargs='+',
        default=['CDS','3utr','5utr','intron']
    )
    # Process arguments
    args = parser.parse_args()
    l2fcwithpval_enr_r1 = args.l2fcwithpval_enr_r1
    l2fcwithpval_enr_r2 = args.l2fcwithpval_enr_r2
    inp_reads_by_loc_r1 = args.inp_reads_by_loc_r1
    inp_reads_by_loc_r2 = args.inp_reads_by_loc_r2
    rep_element_parsed_r1_ip = args.rep_element_parsed_r1_ip
    rep_element_parsed_r1_input = args.rep_element_parsed_r1_input
    rep_element_parsed_r2_ip = args.rep_element_parsed_r2_ip
    rep_element_parsed_r2_input = args.rep_element_parsed_r2_input
    out_file = args.out_file
    annotated_files = args.annotated_files
    annotated_source = args.annotation_script
    pickle_r1 = args.pickle_r1
    pickle_r2 = args.pickle_r2
    suptitle = args.suptitle
    motifs = args.motifs
    regions = args.regions

    # main func
    plot_all(
        l2fcwithpval_enr_r1=l2fcwithpval_enr_r1,
        l2fcwithpval_enr_r2=l2fcwithpval_enr_r2,
        inp_reads_by_loc_r1=inp_reads_by_loc_r1,
        inp_reads_by_loc_r2=inp_reads_by_loc_r2,
        pickle_r1=pickle_r1,
        pickle_r2=pickle_r2,
        rep_element_parsed_r1_ip=rep_element_parsed_r1_ip,
        rep_element_parsed_r1_input=rep_element_parsed_r1_input,
        rep_element_parsed_r2_ip=rep_element_parsed_r2_ip,
        rep_element_parsed_r2_input=rep_element_parsed_r2_input,
        motifs=motifs,
        regions=regions,
        out_file=out_file,
        annotated_files=annotated_files,
        annotation_source=annotated_source,
        suptitle=suptitle
    )


if __name__ == "__main__":
    main()
