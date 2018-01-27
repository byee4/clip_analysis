import pybedtools
import gffutils
import pandas as pd
import numpy as np
import cPickle as pickle
import os
from collections import defaultdict


def filter_deseq2_df(df, l2fc, padj, filter_unknown_contigs=False):
    """
    Returns a dataframe that is filtered for l2fc and padj.
    """
    genes = df.index
    kept = []

    for gene in genes:
        if filter_unknown_contigs:
            if 'ENS' in gene or 'Uni' in gene:
                kept.append(gene)
        else:
            kept.append(gene)
    df = df.ix[kept]
    return df[(np.abs(df['log2FoldChange']) >= l2fc) & (df['padj'] <= padj)]


def get_sig_l2fc_from_deseq2_df(df, l2fc, padj, label=None):
    """
    Reads CSV file from a directory and returns filtered fold change only.

    diffexp_dir: directory where the *.diffexp.txt is found. MUST just contain a single file like this.
    label: label for dataframe column
    l2fc: l2fc filter
    padj: padj filter

    """

    df = filter_deseq2_df(df, l2fc, padj)
    # print(df.shape)
    df = df[['log2FoldChange']]
    df.columns = [label]
    return df


def get_sig_l2fc_from_deseq2_fn(fn, l2fc, padj, label):
    """
    Reads a file and returns just the log2foldchange part
    for any significant gene
    
    Parameters
    ----------
    fn
    l2fc
    padj
    label

    Returns
    -------

    """
    df = read_deseq2_file(fn)
    df = get_sig_l2fc_from_deseq2_df(df, l2fc, padj, label)
    return df

def read_deseq2_file(fn):
    """
    Reads a DESeq2 output file and returns a dataframe.
    
    Parameters
    ----------
    fn

    Returns
    -------

    """
    df = pd.read_table(
        fn, sep=',', index_col=0
    )
    df = df.replace('NA',np.nan)
    df = df.astype(float)
    return df

def read_and_join_deseq2(lst, labels, l2fc, padj):
    """
    Using a comparisons-dict key set, read in and merge filtered log2 fold change
    values. If any gene passes filter params for at least one condition (key), 
    then it will be reported in the final merged dataframe. Returns a merged
    dataframe containing 'significant' genes in at least one condition.
    """
    if labels is None:
        labels = [os.path.basename(l) for l in lst]

    i = 1  # i counter for l in list
    # use the first DESeq2 file to set the index initially. This index will be appended to in the loop below.
    # get_and_read_csv() will filter based on the l2fc and padj values.
    merged = get_sig_l2fc_from_deseq2_fn(lst[0], l2fc, padj, labels[0])

    # outer join the next DESeq2 files. get_and_read_csv() will filter based on the l2fc and padj values.
    for l in lst[1:]:
        df = get_sig_l2fc_from_deseq2_fn(l, l2fc, padj, labels[i])
        # only merge if there are signicant events to merge.
        # This means empty dataframes after filtering won't show up in the final figure.
        if df.shape[0] > 0:
            merged = pd.merge(merged, df, how='outer', left_index=True, right_index=True)

        i += 1
    merged.fillna(0, inplace=True)  # mask all of the insignificant values with zero.
    return merged


def gene_id_to_name(db):
    '''
    Returns a dictionary containing a gene_id:name translation
    Note: may be different if the 'gene_id' or 'gene_name' 
    keys are not in the source GTF file
    (taken from gscripts.region_helpers)
    '''
    genes = db.features_of_type('gene')
    gene_name_dict = {}
    for gene in genes:
        gene_id = gene.attributes['gene_id'][0] if type(gene.attributes['gene_id']) == list else gene.attributes[
            'gene_id']
        try:
            gene_name_dict[gene_id] = gene.attributes['gene_name'][0]
        except KeyError:
            print(gene.attributes.keys())
            print("Warning. Key not found for {}".format(gene))
            return 1
    return gene_name_dict


def convert_gene_id_to_name(df, id2name_dict):
    """
    Using an {id:name} dictionary, replace the dataframe index of gene ids with 
    corresponding names.
    
    Parameters
    ----------
    df
    id2name_dict

    Returns
    -------

    """
    assert 'gene_name' not in df.columns
    df['gene_name'] = df.index.to_series().map(id2name_dict)
    df.set_index('gene_name', inplace=True)
    return df

def get_featuredb(db_file):
    return gffutils.FeatureDB(db_file)