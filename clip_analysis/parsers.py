import cPickle as pickle
import numpy as np
import os
import pandas as pd
import pybedtools
from collections import defaultdict, OrderedDict

# from clip_analysis_legacy import kmerdiff


def split_single_cols(df, col, sep='|'):
    """
    Splits a df['col'] into two separated by 'sep'
    ie. -0.9201|0.00000 -> -0.9201  0.00000
    """
    df["{} l2fc".format('-'.join(col.split(sep)[1:]))], \
    df["{} l10p".format('-'.join(col.split(sep)[1:]))] = zip(
        *df[col].map(lambda x: x.split(sep))
    )
    return df


def split_l2fcwithpval_enr(df, discard = True):
    """
    Splits a dataframe into its l2fc and log10 pvalue
    ie. -0.9201|0.00000 -> -0.9201  0.00000
    
    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of l2fcwithpval_enr file
    discard : bool
        if True, discard the original column.
        if False, keep the original column.

    Returns
    -------
    df : pandas.DataFrame
    """

    for col in df.columns:
        df = split_single_cols(df, col)
        if discard:
            del df[col]
    return df


def read_l2fcwithpval_enr(fn, reps=1):
    """
    Reads in a *l2fcwithpval_enr.csv file and returns a dataframe.

    :param fn: basestring
    :param encode: boolean
        if True, parse using 2IP/1Input style
        if False, parse assuming 1IP/1Input
    :param reps: int
        number of replicates that the "combined.csv.l2fcwithpval_enr.csv" has
    :return: 
    """
    df = pd.read_table(
        fn,
        index_col=0,
    )
    if reps == 1:
        df = split_l2fcwithpval_enr(df)
        df = df.replace('NaN', np.nan)
        df = df.apply(pd.to_numeric)
        return df
    else:
        dfs = []
        for rep in range(1, reps+1):
            substr = '_{:02d}'.format(rep)
            dfx = df[[c for c in df.columns if substr in c]]
            dfx = split_l2fcwithpval_enr(dfx)
            dfx = dfx.replace('NaN', np.nan)
            dfx = dfx.apply(pd.to_numeric)
            dfs.append(dfx)
        return dfs

def filter_l2fcwithpval_enr(l2fcwithpval_enr, region, l10p, l2fc):
    """
    Reads in a l2fcwithpval_enr file and returns just those that pass
    l10p and l2fc cutoffs. Must also specify a region within the l2fcwithpval_enr file.

    Parameters
    ----------
    l2fcwithpval_enr
    region
    l10p
    l2fc

    Returns
    -------
    df : pandas.DataFrame
        table containing the subset of the original l2fcwithpval_enr file.

    """
    df = read_l2fcwithpval_enr(l2fcwithpval_enr)
    df = df.filter(regex=("{}\s".format(region)))
    df = df[(df['{} l10p'.format(region)] >= l10p) & (df['{} l2fc'.format(region)] >= l2fc)]
    return df

def get_avg_l2fcwithpval_enr(
    l2fcwithpval_enr1, l2fcwithpval_enr2, region,
    l10p=-np.log10(0.05), l2fc=0
):
    """
    Returns the average l2 fold from two l2fc files for a given region

    Parameters
    ----------
    l2fcwithpval_enr1
    l2fcwithpval_enr2
    region
    l10p
    l2fc

    Returns
    -------
    df : pandas.DataFrame
        table containing the average l2 fold

    """
    merged_r1 = filter_l2fcwithpval_enr(l2fcwithpval_enr1, region, l10p, l2fc)[['{} l2fc'.format(region)]]
    merged_r2 = filter_l2fcwithpval_enr(l2fcwithpval_enr2, region, l10p, l2fc)[['{} l2fc'.format(region)]]
    merged = pd.merge(merged_r1, merged_r2, how='outer', left_index=True, right_index=True)
    merged.fillna(0, inplace=True)
    merged.index.names = ['Gene']
    merged = pd.DataFrame(merged.mean(axis=1))
    merged.columns = [region]
    return merged

def get_avg_l2fcwithpval_enr_regions(
    l2fcwithpval_enr1, l2fcwithpval_enr2, regions,
    l10p=-np.log10(0.05), l2fc=0
):
    merged = get_avg_l2fcwithpval_enr(
        l2fcwithpval_enr1, l2fcwithpval_enr2, regions[0], l10p, l2fc
    )
    for region in regions[1:]:
        merged = pd.merge(merged, get_avg_l2fcwithpval_enr(
            l2fcwithpval_enr1, l2fcwithpval_enr2, region, l10p, l2fc
        ), how='outer', left_index=True, right_index=True)
    merged.fillna(0, inplace=True)
    return merged


def scatter_matrix(ip_l2fc, inp_reads_by_loc):
    """
    Inner joins the ip l2fc and input reads by loc files
    and builds a dataframe containing: input reads by location, 
    ip 
    Parameters
    ----------
    ip_l2fc : basestring
        "IP*_ReadsByLoc_combined.csv.l2fcwithpval_enr.csv"
    inp_reads_by_loc : basestring
        "INPUT*reads_by_loc.csv"

    Returns
    -------
    x : pandas.DataFrame
        intersection of fold-enrichment and input reads covering each gene.
    """

    plot_x = pd.read_table(
        inp_reads_by_loc,
        index_col=0,
    )

    plot_y = read_l2fcwithpval_enr(ip_l2fc)

    x = pd.merge(plot_x, plot_y, how='inner', left_index=True,
                 right_index=True)
    return x

def filter_input_norm(file_name, l10p, l2fc, as_bedtool=False):
    """
    Filters an "input norm"-formatted file given
    log2 fold change and log10 pvalue thresholds.
    See data/input_norm_bed.bed file for an example.
    
    :param file_name: basestring
        filename of the input normalized file.
    :param l2fc: float
        log2 fold change cutoff
    :param l10p: float
        -log10 pvalue cutoff
    :param col_names: list
        column names of the bed or bedlike file.
    :param as_bedtool: bool
        if True, return as bedtool. 
        if False, return dataframe.
    :return: 
    """
    try:
        col_names = ['chrom','start','end','l10p','l2fc','strand']
        df = pd.read_table(file_name, names=col_names)
        dff = filter_df(df, l10p, l2fc)
        print("before: {}, after: {}".format(df.shape[0], dff.shape[0]))
        if as_bedtool:
            return pybedtools.BedTool.from_dataframe(dff)
        return dff

    except Exception as e:
        print(e)
        return 1

def filter_df(df, l10p, l2fc):
    """
    Returns a dataframe filtered using l10p and l2fc values.
    Assumes the columns are aptly named 'l2fc' and 'l10p'.
    
    :param df: 
    :param l2fc: 
    :param l10p: 
    :return: 
    """
    return df[(df['l2fc']>=float(l2fc)) & (df['l10p']>=l10p)]

def return_region_eric(row):
    """
    Given a row of a inputnormed bedfile, return region
    Row must be in the same format as a line in Eric's
    *.annotated file.

    """
    try:
        if row['annotation'] == 'intergenic':
            return 'intergenic'
        region = row['annotation'].split('|')[0]

        return region
    except Exception as e:
        print(e)
        print(row)

def read_annotated_file(fn, headers=None, src='brian'):
    """
    Reads in an annotated bedfile from either eric or me.
    Ensures that 'region' column is defined and set.
    Returns dataframe
    
    :param fn: basestring
        filename
    :param headers: list
        if src isn't from eric or me, you have to 
        explicitly set the header columns.
    :param src: basestring
        either 'brian' or 'eric' depending on whose script you use.
    :return df: pandas.DataFrame
    """
    if src == 'brian':
        headers=[
            'chrom','start','end','l10p','l2fc','strand','geneid',
            'genename','region','alloverlapping'
        ]
        df = pd.read_table(fn, names=headers)
    elif src == 'eric':
        headers=[
            'chrom','start','end','l10p','l2fc',
            'strand','annotation','geneid'
        ]
        df = pd.read_table(fn, names=headers)
        df['region'] = df.apply(return_region_eric, axis=1)
    elif src == 'eric_v2':
        pass # TODO: get eric's most up to date version
    else:
        assert 'region' in headers
        df = pd.read_table(fn, names=headers)
    return df

def get_region_counts(fn, headers, src):
    """
    Returns a dataframe of REGIONS and the number of peaks
    associated with each region.
    
    :param fn: 
    :return: 
    """
    return pd.DataFrame(
        read_annotated_file(fn, headers, src)['region'].value_counts()
    )

def get_counts_df(dfs):
    """
    Returns the counts of each region in the annotation file.
    
    :param dfs: dict
        dict of dataframes and their associated region value counts.
        Keys = filenames, values = dataframes belonging to that filename.
    :param headers: list
        optional if either 'brian' or 'eric' are src. Otherwise, 
        this MUST contain a 'region' category. 
    :param src: basestring
        either 'brian' or 'eric' to denote the annotation structure.
    :return merged: pandas.DataFrame
        table containing the sum of counts for each region in the annotated file.
    """
    columns = dfs.keys()
    assert len(columns) > 0
    annotated_df = dfs[columns[0]]
    merged = annotated_df
    # merged.columns = [os.path.basename(annotated_file)]
    for key in columns[1:]:
        merged = pd.merge(merged, dfs[key], how='left', left_index=True,
                          right_index=True)
    return merged

def get_counts(fns, headers=None, src='brian', basename=True):
    dfs = OrderedDict()
    for fn in fns:
        if basename:
            key = os.path.basename(fn)
        else:
            key = fn
        dfs[key] = get_region_counts(fn, headers, src)
        dfs[key].columns = [key]
    return get_counts_df(dfs)


def remove_peaks(wt_peaks, ko_peaks, as_bedtool=False):
    """
    Given a bed file, remove any peak contained with a 'peaks_to_remove'
    file.

    Parameters
    ----------
    wt_peaks : pybedtools.BedTool
    ko_peaks : pybedtools.BedTool
    as_bedtool : bool

    Returns
    -------

    """
    # ko = pybedtools.BedTool(ko_peaks)
    # wt = pybedtools.BedTool(wt_peaks)


    bedtool = wt_peaks.intersect(ko_peaks, v=True)
    print('number of peaks before: {}, after: {}'.format(len(wt_peaks), len(bedtool)))

    if as_bedtool:
        return bedtool
    df = bedtool.to_dataframe()
    return df


def read_parsed(fn):
    """
    Reads Eric's parsed file from the repetitive element pipeline.

    Parameters
    ----------
    fn : basestring
        the *.parsed file

    Returns
    -------
    total_df : pandas.DataFrame
        dataframe of total reads per unique/repetitive element family.
    element_df : pandas.DataFrame
        dataframe of unique repetitive/unique elements that each unique read
        mapped to.
    total_reads : int
    total_genomic_reads : int
    total_usable_reads : int
    total_repfamily_reads : int
    """
    df = pd.read_table(fn, names=[
        'total_or_element', 'element', 'read_num',
        'clip_rpr', 'annotation', 'gene'
    ])
    total_reads = df[
        (df['total_or_element'] == '#READINFO') & (df['element'] == 'AllReads')
        ]['read_num'].values[0]
    total_genomic_reads = df[
        (df['total_or_element'] == '#READINFO') & (
        df['element'] == 'GenomicReads')
        ]['read_num'].values[0]
    total_usable_reads = df[
        (df['total_or_element'] == '#READINFO') & (
        df['element'] == 'UsableReads')
        ]['read_num'].values[0]
    total_repfamily_reads = df[
        (df['total_or_element'] == '#READINFO') & (
        df['element'] == 'RepFamilyReads')
        ]['read_num'].values[0]

    total_df = df[df['total_or_element'] == 'TOTAL'][
        ['element', 'read_num', 'clip_rpr']
    ]
    element_df = df[df['total_or_element'] == 'ELEMENT'][
        ['element', 'read_num', 'clip_rpr']
    ]
    return total_df, element_df, \
           total_reads, total_genomic_reads, \
           total_usable_reads, total_repfamily_reads


def return_l2fc_entropy_from_parsed(ip_parsed, input_parsed, nopipes=True):
    """
    From 2 parsed rep element pipeline outputs (ip and input),
    compute fold change and information content. Usually fold changes of > 3+
    and information content of 0.1? can be considered enriched.

    Parameters
    ----------
    ip_parsed : str
        filename of the ip parsed string
    input_parsed : str
        filename of the input parsed string
    nopipes : bool
        if True, return just the uniquely mapped rep family mappings
        if False, return all unique and nonunique
    Returns
    -------
    merged : Pandas.DataFrame
        table consisting of fold enrichment and information content params
    """
    total_ip, _, _, _, _, _ = read_parsed(ip_parsed)
    total_input, _, _, _, total_input_usable_reads, _ = read_parsed(
        input_parsed)
    # a pipe indicates read totals mapping to more than one element/rep family.
    if nopipes:
        total_ip = total_ip[total_ip['element'].str.contains('\|') == False]
        total_input = total_input[
            total_input['element'].str.contains('\|') == False]
    # index columns by their element
    total_ip.set_index('element', inplace=True)
    total_input.set_index('element', inplace=True)
    # rename the IP and input columns separately
    total_ip.columns = ["IP_{}".format(c) for c in total_ip.columns]
    total_input.columns = ["Input_{}".format(c) for c in total_input.columns]
    # merge the two on element id
    merged = pd.merge(total_ip, total_input, how='left', left_index=True,
                      right_index=True)
    # deal with missing values
    merged['Input_read_num'].fillna(
        1, inplace=True
    )  # Pseudocount all missing values
    merged['Input_clip_rpr'].fillna(
        merged['Input_read_num'] / (total_input_usable_reads), inplace=True)
    # calculate fold enrichment and information content
    merged['Fold_enrichment'] = merged['IP_clip_rpr'].div(
        merged['Input_clip_rpr'])
    merged['Information_content'] = merged['IP_clip_rpr'] * np.log2(
        merged['IP_clip_rpr'].div(merged['Input_clip_rpr']))

    return merged

# LEGACY functions to handle some of the old gabe and eric stuff #

def read_kmer_enrichment_from_pickle(
        pickle_file, region='all', k=6, col_name='zscore delta'
):
    """
    Reads in a pickle file from gabe's clip_analysis script and returns a
    dataframe containing kmers and their enriched z-scores

    :param pickle_file: basestring
        pickle filename output from gabe's clip_analysis script.
    :param region: basestring
        one of:
        'all', 'three_prime_utrs', 'five_prime_utrs', 'distintron500',
        'cds', 'proxintron500'
    :return df: pandas.DataFrame
    """
    loaded = pickle.load(open(pickle_file, 'rb'))
    try:
        df = pd.DataFrame(loaded['kmer_results'][region][k]).T
    except KeyError:
        df = pd.DataFrame(loaded['kmer_results'][region][str(k)]).T
    df.columns = ['fg', 'bg', col_name]
    return df[[col_name]]

def extract_kmer_enrichment_from_pickle_and_save(pickle_file, out_file):
    df = pd.read_table(pickle_file)
    df.to_csv(out_file, sep='\t', index=True, headers=True)

def get_top_enriched_kmers(merged, top=5):
    average_zscore = merged.mean(axis=1)
    return list(
        average_zscore.sort_values(ascending=False).iloc[:top].index
    )

def bed6_to_bed8(interval):
    """
    Basically appends the start/stop fields to 'thickStart/thickStop' fields
    Turns BED6 into BED8 (formerly called: make_clipper_ish)
    (Helps with plot_clip_analysis_figures.py in CLIPper). 
    
    Parameters
    ----------
    interval : pybedtools.Interval

    Returns
    -------

    """
    interval.name = interval[7]
    interval[6] = interval.start
    interval[7] = interval.stop

    return interval


def filter_data(interval, l2fc, pval):
    """
    col4 is -log10 p-val
    col5 is -log2 fold enrichment

    Expects the standard input norm file format.

    Parameters
    ----------
    interval : pybedtools.Interval
    l2fc : float
    pval : float

    Returns
    -------

    """

    return (float(interval[4]) >= pval) and (float(interval[3]) >= l2fc)


def split_annotated_file_into_region_beds_and_save(
    annotated, headers=None, src='brian',
    split_beds_directory=os.getcwd()
):
    if not os.path.exists(split_beds_directory):
        os.mkdir(split_beds_directory)

    df = read_annotated_file(annotated, headers, src)
    regions = set(df['region'])
    for region in regions:
        region_df = df[df['region']==region][[
            'chrom','start','end','l10p','l2fc','strand'
        ]]
        region_df.to_csv(
            os.path.join(
                split_beds_directory,
                os.path.basename(annotated) + ".{}.bed".format(
                    region
                )
            ),
            sep='\t',
            header=False,
            index=False
        )

### This is modified from Emily's code ###

def get_shuffled_region(bed_file_tool, bed_file_to_shuffle_within):
    """
    Makes a shuffled file
    :param bed_file: Bed file of peaks to shuffle
    :param bed_file_to_shuffle_within: Bed file of REGIONS to include in the shuffling.
        This is typicalls a file that matches the genic region of peaks (exons, 3'UTRs, etc.)
    :return: a shuffled bedtool object
    """

    # bed_file_tool = pybedtools.BedTool(bed_file)
    shuffled = bed_file_tool.shuffle(genome="hg19",
                                     incl=bed_file_to_shuffle_within,
                                     noOverlapping=True)
    return shuffled
