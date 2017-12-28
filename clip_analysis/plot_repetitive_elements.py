#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib.patches as patches
# from clip_analysis import parsers as p
import parsers as p
### do this to avoid making tons of dots in svg files:


rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Default REGIONS to look at;
NUM_ELEMENTS = 10

# Plotting functions #

def plot(ip_parsed, input_parsed, ax=None,
         title="Top {} unique/repetitive elements mapped",
         num_elements=10, read_threshold=10):
    if ax is None:
        ax = plt.gca()
    ax2 = plt.twinx()

    blues = sns.color_palette("Blues", 10)
    df = p.return_l2fc_entropy_from_parsed(
        ip_parsed, input_parsed
    ).sort_values(
        by=['Information_content', 'Fold_enrichment'],
        ascending=False
    )
    width = 0.25
    ind = np.arange(num_elements)
    dx = df.iloc[:num_elements]
    colors = []  # colors for plotting fold enrichment
    for read_num in dx['IP_read_num']:
        if read_num < read_threshold:
            colors.append(blues[1])
        else:
            colors.append(blues[9])
    ax.bar(ind, dx['Fold_enrichment'], width, color=colors)
    ax.set_ylabel('Fold enrichment')
    ax2.bar(ind + width, dx['Information_content'], width, color='red')
    ax2.set_ylabel('Information content')
    ax.set_xticks(range(num_elements))

    ax.set_title(title.format(num_elements))
    ticklabels = ax.set_xticklabels(dx.index, rotation='vertical')
    less = patches.Patch(color=blues[0],
                         label='less than {} IP reads'.format(read_threshold))
    greater = patches.Patch(color=blues[9],
                            label='Fold enrichment')
    information = patches.Patch(color='red',
                            label='Information content')
    leg = ax.legend(
        loc=1,
        handles=[less, greater, information],
        borderaxespad=0.
    )
    plt.tight_layout()
    return ax

def make_plot(ip_parsed, input_parsed, out_file, title, num_elements):
    fig, ax = plt.subplots()
    plot(
        ip_parsed=ip_parsed,
        input_parsed=input_parsed,
        title=title,
        ax=ax,
        num_elements=num_elements
    )
    fig.savefig(out_file)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ip_parsed",
        required=True,
    )
    parser.add_argument(
        "--input_parsed",
        required=True,
    )
    parser.add_argument(
        "--num_elements",
        help="Specify how many of the top N elements to plot"
             " (default: {})".format(
            NUM_ELEMENTS
        ),
        default=NUM_ELEMENTS,
        type=int,
        required=False
    )
    parser.add_argument(
        "--title",
        required=False,
        default="Top 10 unique/repetitive elements mapped"
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )

    # Process arguments
    args = parser.parse_args()
    ip_parsed = args.ip_parsed
    input_parsed = args.input_parsed
    out_file = args.out_file
    num_elements = args.num_elements
    title = args.title

    # main func
    make_plot(
        ip_parsed, input_parsed, out_file, title, num_elements
    )


if __name__ == "__main__":
    main()
