'''
Created on April 2nd, 2018

Calculates entropy score for peak file

@author: brianyee (bay001@ucsd.edu)
@author: gabriel pratt (gpratt@ucsd.edu)
'''

import numpy as np
import pandas as pd
import os
import argparse
import functools

# Gabe's entropy functions #

annotated_bedtool_header = ['chrom', 'start', "stop", "name", "score",
                            "strand", "annotation", "gene_id"]
full_header = ["chrom", "start", "stop", "full_name", "ip_reads",
               "input_reads", "p_val", "chisq", "test_type",
               "enrichment", "log10_p_val", "log2_fold_change"]


# def get_full_from_annotated(fn):
#     stripped_fn = ".".join(fn.split(".")[:-3])
#     return stripped_fn + ".full.compressed2.bed.full"


def calculate_entropy(row, total_ip_reads, total_input_reads):
    """
    Calculates the entropy score given

    Parameters
    ----------
    row
    total_ip_reads
    total_input_reads

    Returns
    -------

    """
    p_ip = float(row.ip_reads) / total_ip_reads
    p_input = float(row.input_reads) / total_input_reads

    return p_ip * np.log2(p_ip / p_input)


def get_entropy_from_annotated(fn):
    fn = os.path.basename(fn)
    stripped_fn = ".".join(fn.split(".")[:-3])
    stripped_fn = os.path.join(out_dir, stripped_fn)
    return stripped_fn + ".full.compressed2.bed.full.entropy.bed"


def sum_entropy(filtered_peaks, original_peaks, out_file=None):
    if out_file is None:
        out_file = filtered_peaks + ".entropy.bed"

    entropy = pd.read_table(get_entropy_from_annotated(original_peaks))
    filtered_peaks = pd.read_table(filtered_peaks,
                                   names=annotated_bedtool_header)

    merged_peaks = pd.merge(filtered_peaks, entropy,
                            left_on=['chrom', 'start', 'stop'],
                            right_on=['chrom', 'start', 'stop'])
    merged_peaks.to_csv(out_file,
                        sep="\t", header=None, index=False)
    return merged_peaks.entropy.sum()


def sum_entropy_row(row):
    # Sadly the majority of the time in this operation is opening the files, can't make it faster :(
    try:
        entropy = sum_entropy(row.filtered_moderate, row['input_norm'])
    except Exception as e:
        print e, "something went wrong"
        entropy = 0
    return entropy

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--full_file",
        required=True,
    )
    parser.add_argument(
        "--output_file",
        required=True,
    )

    args = parser.parse_args()
    full_file = args.full_file
    output_file = args.output_file

    full_fn = full_file
    out_fn = output_file
    if os.path.exists(out_fn):
        pass

    ip_reads = row['CLIP_counts']
    input_reads = row['INPUT_counts']

    read_counts = pd.read_table(full_fn, names=full_header)
    if len(read_counts) != 0:
        tool = functools.partial(calculate_entropy,
                                 total_ip_reads=ip_reads,
                                 total_input_reads=input_reads)
        read_counts['entropy'] = read_counts.apply(tool, axis=1)
        read_counts.to_csv(out_fn, sep="\t", index=False, header=True)


if __name__ == "__main__":
    main()
