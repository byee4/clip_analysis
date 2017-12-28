#!/usr/bin/env python

"""

Filters an input normalized peak BED6 file (log10 pvalue as the 'name' column,
log2 foldchange as the 'score' column).

"""
import argparse
import os
from clip_analysis import parsers
import pybedtools


def filter(peaks_file, ko_peaks_file, output_file, l10p, l2fc):

    df = parsers.filter_input_norm(peaks_file, l10p, l2fc)
    if ko_peaks_file is not None:
        print("removing {} from {}".format(
            os.path.basename(ko_peaks_file), os.path.basename(peaks_file))
        )
        wt_peaks = pybedtools.BedTool.from_dataframe(df)
        ko_peaks = pybedtools.BedTool(ko_peaks_file)
        df = parsers.remove_peaks(wt_peaks, ko_peaks)

    df.to_csv(
        output_file, sep='\t', header=False, index=False
    )

def main():
    """
    Main program.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
    )
    parser.add_argument(
        "--output",
        required=True,
    )
    parser.add_argument(
        "--l10p",
        required=False,
        type=int,
        default=3
    )
    parser.add_argument(
        "--l2fc",
        required=False,
        type=int,
        default=3
    )
    parser.add_argument(
        "--ko_peaks",
        required=False,
        default=None
    )
    args = parser.parse_args()

    l10p = args.l10p
    l2fc = args.l2fc
    i = args.input
    o = args.output
    ko_peaks = args.ko_peaks

    filter(i, ko_peaks, o, l10p, l2fc)

if __name__ == "__main__":
    main()
