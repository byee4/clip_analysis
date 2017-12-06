import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import argparse
import parsers as p
import plot_compare_rbp_enriched_regions
import plot_histogram_enriched_regions
import plot_motifs
import plot_region_distribution
import plot_ip_foldchange_over_input_reads
import seaborn as sns

### do this to avoid making tons of dots in svg files:
from matplotlib import rc

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Set default color scheme for stuff
palette = sns.color_palette("hls", 8)


def plot_all(l2fcwithpval_enr_r1, l2fcwithpval_enr_r2, inp_reads_by_loc_r1,
             inp_reads_by_loc_r2, out_file, annotated_files):
    nrows = 3
    ncols = 2

    full_grid = gridspec.GridSpec(
        nrows, ncols,
        height_ratios=[1 for i in range(nrows)],
        width_ratios=[1 for i in range(ncols)],
        hspace=0.5, wspace=3
    )
    fig = plt.figure(figsize=(15, 25))

    map_rows = []

    for row in range(0, nrows):
        map_rows.append(gridspec.GridSpecFromSubplotSpec(
            1, ncols,
            subplot_spec=full_grid[row, :])
        )

    plot_histogram_enriched_regions.plot(
        l2fcwithpval_enr_r1, ax=plt.subplot(map_rows[0][0])
    )
    plot_histogram_enriched_regions.plot(
        l2fcwithpval_enr_r2, ax=plt.subplot(map_rows[0][1])
    )

    plot_ip_foldchange_over_input_reads.plot(
        l2fcwithpval_enr_r1, inp_reads_by_loc_r1,
        ax=plt.subplot(map_rows[1][0])
    )
    plot_ip_foldchange_over_input_reads.plot(
        l2fcwithpval_enr_r2, inp_reads_by_loc_r2,
        ax=plt.subplot(map_rows[1][1])
    )
    plot_compare_rbp_enriched_regions.plot(
        l2fcwithpval_enr_r1, l2fcwithpval_enr_r2,
        ax=plt.subplot(map_rows[2][0])
    )

    counts = p.get_counts(annotated_files)
    plot_region_distribution.plot(counts, ax=plt.subplot(map_rows[2][1]))
    fig.savefig(out_file)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--l2fcwithpval_enr_r1",
        required=True
    )
    parser.add_argument(
        "--l2fcwithpval_enr_r2",
        required=False,
        default=''
    )
    parser.add_argument(
        "--inp_reads_by_loc_r1",
        required=True
    )
    parser.add_argument(
        "--inp_reads_by_loc_r2",
        required=False,
        default=''
    )
    parser.add_argument(
        "--annotated_files",
        required=True,
        nargs='+'
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )

    # Process arguments
    args = parser.parse_args()
    l2fcwithpval_enr_r1 = args.l2fcwithpval_enr_r1
    l2fcwithpval_enr_r2 = args.l2fcwithpval_enr_r2
    inp_reads_by_loc_r1 = args.inp_reads_by_loc_r1
    inp_reads_by_loc_r2 = args.inp_reads_by_loc_r2
    out_file = args.out_file
    annotated_files = args.annotated_files

    # main func
    plot_all(
        l2fcwithpval_enr_r1, l2fcwithpval_enr_r2,
        inp_reads_by_loc_r1, inp_reads_by_loc_r2,
        out_file, annotated_files
    )


if __name__ == "__main__":
    main()
