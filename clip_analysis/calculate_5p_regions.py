'''
Created on April 2nd, 2018

Intersects forward-stranded reads (seCLIP, peCLIP.r2) with significant peaks
and saves the 5' regions (+/- 10nt) of each read for more granular motif analysis.

@author: brianyee (bay001@ucsd.edu)
'''

import numpy as np
import pandas as pd
import os
import argparse
import pysam
import subprocess
import pybedtools
from subprocess import check_call

# Gabe's entropy functions #


def return_chrom_sizes_dict(chrom_sizes_file):
    """
    Returns chrom_sizes[chrom] = int(sizes) dict.
    """
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            chrom, sizes = line.strip().split('\t')
            chrom_sizes[chrom] = int(sizes)
    return chrom_sizes


def bedtools_intersect_bam_with_peaks(bam_file, peaks_file, output_intersected_bam):
    """
    bedtools intersect -abam $bam -s -b $bed > $readsinpeaksbam
    """
    # bam = pybedtools.BedTool(bam_file)
    # peaks = pybedtools.BedTool(peaks_file)
    # reads_in_peaks_bam = bam.intersect(peaks, s=True).saveas(output_intersected_bam)
    
    priming_call = "bedtools intersect -abam {} -s -b {}".format(
        bam_file, peaks_file
    )
    with open(output_intersected_bam, 'w') as f:
        check_call(priming_call, shell=True, stdout=f)  # TODO: convert to pybedtools func

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam",
        required=True,
    )
    parser.add_argument(
        "--chrom_sizes",
        required=True,
    )
    parser.add_argument(
        "--peaks",
        required=True,
    )
    parser.add_argument(
        "--output_5p_bed",
        required=True,
    )
    parser.add_argument(
        "--output_intersected_bam",
        required=True,
    )
    
    args = parser.parse_args()
    bam = args.bam
    chrom_sizes = args.chrom_sizes
    peaks = args.peaks
    output_5p_bed = args.output_5p_bed
    output_intersected_bam = args.output_intersected_bam

    chrom_sizes_dict = return_chrom_sizes_dict(chrom_sizes)
    bedtools_intersect_bam_with_peaks(bam, peaks, output_intersected_bam)
    
if __name__ == "__main__":
    main()
