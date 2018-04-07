#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import parsers as p

### do this to avoid making tons of dots in svg files:

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

def get_start(current_width, previous_width, previous_step):
    """
    Returns the appropriate start offset based on the previous "pie" slice.
    """
    circle = np.pi * 2  # 1 radian = 6.28
    if previous_width is None:
        previous_width = 0
    stepper = circle * 0.05  # okay honestly don't really know why this needs to be 0.05
    return previous_step + stepper * (current_width * 10 + previous_width * 10)


def sunburst(widths, heights, labels, colors, ax):
    """
    Make a sunburst pie chart (pie chart with different heights for each pie)
    """
    try:
        assert float(sum(widths)) == 1.0  # make sure all slices equal 1
    except AssertionError:
        print(
            "Pie slices do not equal 1! Pie chart may be incomplete! {}".format(
                sum(widths)
            )
        )

    previous_step = 0
    previous_width = None

    starts = []
    for w in range(len(widths)):
        starts.append(get_start(widths[w], previous_width, previous_step))
        previous_step = get_start(widths[w], previous_width, previous_step)
        previous_width = widths[w]

    circle = np.pi * 2  # 1 radian = 6.28
    widths = [circle * x for x in widths]

    # let's make a hard cutoff of bar heights, setting at 14 for now...
    for h in range(len(heights)):
        if heights[h] > 14:
            heights[h] = 14

    bars = ax.bar(starts, heights, widths,
                  bottom=0.0, clip_on=False)  # change bottom to 1 for 'donut' look

    # Use custom colors and opacity
    i = 0
    for r, bar in zip(heights, bars):
        bar.set_alpha(0.5)
        bar.set_color(colors[labels[i]])
        bar.set_linewidth(0)
        bar.set_label("{}".format(labels[i]))
        i += 1

    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_rlabel_position(0)

    # plt.setp(ax.get_yticklabels(), visible=False)  # keep or get rid of yticks?
    # set_legend(ax, ax)
    ax.set_ylim(0, 10)
    ax.legend(bbox_to_anchor=(1.55, 1))
    return ax


def split_others(df, fold_change_cutoff, soft_cutoff, hard_cutoff, focus=[]):
    """
    Given a fold change and hard/soft read cutoffs, separate
    any element that passes filtering requirements from those
    that don't. Any element that does not pass either cutoff
    will be binned into an 'others' category.

    default 'passing' requirements:
    1) if ip_rpr >= soft_cutoff and fold_change >= 1
    2) if ip_rpr >= hard_cutoff and fold_change >= 4

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of a *parsed file (from repetitive element pipeline)
        Usually split by RBP (contains elements from only one RBP at a time)
    fold_change_cutoff : float
        fold change required in conjunction with hard cutoff to show
    soft_cutoff : float
        show element if it passes this threshold and does not have
        depletion
    hard_cutoff : float
        bin into 'others' if fails this cutoff no matter the fold change
    focus: list
        list of elements to plot regardless of fold or cutoff thresholds.
    Returns
    -------
    binned : pandas.DataFrame
    """
    pass1 = df[
        (df['IP_clip_rpr'] >= soft_cutoff) & \
        (df['Fold_enrichment'] >= 1)
        ]
    pass2 = df[
        (df['IP_clip_rpr'] >= hard_cutoff) & \
        (df['Fold_enrichment'] >= fold_change_cutoff)
        ]
    kept_elements = set(list(pass1.index) + list(pass2.index) + list(focus))
    # print('elements to plot: {}'.format(kept_elements))
    passed = pd.Index(kept_elements)
    failed = df.index.difference(passed)
    return df.loc[passed], df.loc[failed]


def get_pie_values(df, fold_change_cutoff, soft_cutoff,
                   hard_cutoff, focus):
    """

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of a *parsed file (from repetitive element pipeline)
    fold_change_cutoff : float
        fold change required in conjunction with hard cutoff to show
    soft_cutoff : float
        show element if it passes this threshold and does not have
        depletion
    hard_cutoff : float
        bin into 'others' if fails this cutoff no matter the fold change
    focus : list
        list of elements to plot regardless of fold or cutoff thresholds.
    Returns
    -------
    merged : pandas.DataFrame
        dataframe of a filtered set of elements that show enrichment
        beyond fold change cutoffs and % cutoffs. 
    """
    dx = df[['IP_clip_rpr', 'Fold_enrichment']]
    dpass, dfail = split_others(
        dx, fold_change_cutoff, soft_cutoff, hard_cutoff, focus
    )
    others = pd.DataFrame(
        columns=['IP_clip_rpr', 'Fold_enrichment'],
        index=['others']
    )
    others.loc['others']['IP_clip_rpr'] = dfail['IP_clip_rpr'].sum()
    others.loc['others']['Fold_enrichment'] = dfail['Fold_enrichment'].mean()
    merged = pd.concat([dpass, others])
    # because the ip_circular_... doesn't always add up to 1
    if merged['IP_clip_rpr'].sum(axis=0) < 1:
        merged.loc['others', 'IP_clip_rpr'] += (
            1 - merged['IP_clip_rpr'].sum(axis=0)
        )
    return merged


# Plotting functions #

def plot(
    ip_parsed, input_parsed, color_file,
    title,
    fold_change_cutoff, soft_cutoff, hard_cutoff, focus, ax=None,
):
    if ax is None:
        ax = plt.subplot(111, projection='polar')

    df = p.return_l2fc_entropy_from_parsed(
        ip_parsed, input_parsed
    )
    color_dict = make_colordict(set(df.index), color_file)
    vals = get_pie_values(df, fold_change_cutoff, soft_cutoff, hard_cutoff, focus=focus)
    sunburst(vals['IP_clip_rpr'], vals['Fold_enrichment'], vals.index, color_dict,
             ax)
    plt.tight_layout()
    return ax


def make_plot(ip_parsed, input_parsed, color_file, out_file, title,
              fold_change_cutoff, soft_cutoff, hard_cutoff, focus):
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection='polar')
    plot(
        ip_parsed=ip_parsed,
        input_parsed=input_parsed,
        color_file=color_file,
        title=title,
        ax=ax,
        fold_change_cutoff=fold_change_cutoff,
        soft_cutoff=soft_cutoff,
        hard_cutoff=hard_cutoff,
        focus=focus
    )
    fig.savefig(out_file)


def make_colordict(elements, color_fn):
    colors = {}
    with open(color_fn, 'r') as f:
        colors[
            'others'] = f.readline().rstrip()  # don't like black, it's the first color so... skip
        for element in elements:
            colors[element] = f.readline().rstrip()

    return colors


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
        "--fold_change_cutoff",
        required=False,
        default=4,
        type=int,
    )
    parser.add_argument(
        "--soft_cutoff",
        required=False,
        default=0.04,
        type=float,
    )
    parser.add_argument(
        "--hard_cutoff",
        required=False,
        default=0.01,
        type=float
    )
    parser.add_argument(
        "--focus",
        required=False,
        default=[],
        nargs='+'
    )
    parser.add_argument(
        "--color_file",
        help='line delimited file containing hex colors',
        default=None,
        type=str,
        required=False
    )
    parser.add_argument(
        "--title",
        required=False,
        default="Sunburst pie chart"
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )

    # Process arguments
    args = parser.parse_args()
    ip_parsed = args.ip_parsed
    input_parsed = args.input_parsed
    fold_change_cutoff = args.fold_change_cutoff
    soft_cutoff = args.soft_cutoff
    hard_cutoff = args.hard_cutoff
    focus = args.focus
    out_file = args.out_file
    color_file = args.color_file
    title = args.title


    # main func
    make_plot(
        ip_parsed, input_parsed, color_file, out_file, title,
        fold_change_cutoff, soft_cutoff, hard_cutoff, focus
    )


if __name__ == "__main__":
    main()
