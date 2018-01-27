import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

from fastcluster import linkage
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from scipy.cluster import hierarchy
from scipy.spatial import distance
from Clusters import get_cluster_classes

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import argparse
import parsers as p

### do this to avoid making tons of dots in svg files:
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

### Default REGIONS to look at;
REGIONS = ['CDS', '3utr', '5utr', 'intron']
MAX_CLUSTERS = 15
CLUSTER_COLORS = [str(x) for x in sns.color_palette('hls', MAX_CLUSTERS).as_hex()]

# Clustering functions



def get_row_linkage(df, metric='euclidean', method='average'):
    hierarchy.set_link_color_palette(
        CLUSTER_COLORS
    )
    row_linkage = hierarchy.linkage(
        distance.pdist(
            df, metric=metric
        ),
        method=method,
        metric=metric,
        optimal_ordering=True,
    )
    return row_linkage

def get_dendrogram(df, row_linkage, color_threshold, out_dendrogram):
    fig, ax = plt.subplots(figsize=(75,15)) # make this really big so we can see individual genes

    den = dendrogram(
        row_linkage, labels=df.index,
        color_threshold=color_threshold,
        no_plot=False
    )

    fig.savefig(out_dendrogram)
    return den

def get_sidebar_colors(df, leaf_height, out_dendrogram):
    row_colors = pd.DataFrame(index=df.index)
    row_linkage = get_row_linkage(df)
    den = get_dendrogram(df, row_linkage, leaf_height, out_dendrogram)
    row_colors['Cluster'] = 'White' # initialize with some variable
    clusters = get_cluster_classes(den)
    if len(clusters.keys()) > MAX_CLUSTERS:
        print("Warning: more clusters ({}) than there are colors ({}). "
              "Please use higher leaf height to bin genes into larger groups!").format(len(clusters.keys()), MAX_CLUSTERS)

    for key in clusters.keys():
        for gene in get_cluster_classes(den)[key]:
            row_colors.loc[gene]['Cluster'] = key
    return row_colors

def set_heatmap_minmax(df):
    """
    Gets equivalent pos and neg boundaries for heatmap
    based on dataframe values
    
    Parameters
    ----------
    df

    Returns
    -------

    """
    if abs(min(df.min())) > abs(max(df.max())):
        vmin = int(min(df.min())) - 1
        vmax = -vmin
    else:
        vmax = int(max(df.max())) + 1
        vmin = -vmax
    return vmin, vmax

def write_clusters_vertically(clusters, out_cluster):
    """
    Writing the clusters as columns and genes vertically makes it 
    easier to copy/paste into GO enrichment from excel. 
    
    Parameters
    ----------
    clusters
    out_cluster

    Returns
    -------

    """
    merged = pd.DataFrame(index=[0])

    for color in set(clusters['Cluster']):
        dx = pd.DataFrame(clusters[clusters['Cluster'] == color].index)
        dx.columns = [color]
        merged = pd.merge(merged, dx, how='outer', left_index=True, right_index=True)

    merged.to_csv(out_cluster, sep='\t', index=False, header=True)

# Plotting functions
def make_plot(
        in_files, in_labels, id2name_dict,
        l2fc, padj, leaf_height,
        out_file, out_tsv, out_cluster, out_dendrogram,
        gene_label_limit=100, fig_size=(15, 10)
):
    # Run once to get a subset of significant genes
    df = p.read_and_join_deseq2(
        in_files, in_labels, l2fc, padj
    )
    names_to_keep = df.index

    # Run again to get all fold changes of these genes
    df = p.read_and_join_deseq2(
        in_files, in_labels, 0, 1
    )
    df = df.ix[names_to_keep]
    print(df.shape)
    # Set a limit to display gene names, too many genes
    # results in muddled plot.
    if df.shape[0] > gene_label_limit:
        set_yticklabels = False
    else:
        set_yticklabels = True

    # set the vmin and vmax (centers the heatmap)
    vmin, vmax = set_heatmap_minmax(df)

    # get the clusters
    clusters = get_sidebar_colors(df, leaf_height, out_dendrogram)

    # convert gene_ids to names in cluster and df indices:
    if id2name_dict is not None:
        df = p.convert_gene_id_to_name(df, id2name_dict)
        clusters = p.convert_gene_id_to_name(clusters, id2name_dict)

    # Save tabbed heatmap values
    df.to_csv(out_tsv, sep='\t')
    write_clusters_vertically(clusters, out_cluster)

    # Plot the heatmap.
    cg = sns.clustermap(
        df, xticklabels=True, yticklabels=set_yticklabels, figsize=fig_size,
        metric='euclidean', row_colors=clusters,
        method='average', edgecolors='white', linewidths=0.000,
        cmap='bwr', vmin=vmin, vmax=vmax
    )
    plt.title('Log2 FC')
    ticks = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    cg.savefig(out_file)
    return 0


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--in_files",
        required=True,
        nargs='+'
    )
    parser.add_argument(
        "--l2fc",
        required=False,
        default=1.5
    )
    parser.add_argument(
        "--padj",
        required=False,
        default=0.05
    )
    parser.add_argument(
        "--in_labels",
        help="Specify labels (by default will "
             "use the file basenames)",
        nargs='+',
        default=None,
        required=False
    )
    parser.add_argument(
        "--out_file",
        required=True,
    )
    parser.add_argument(
        "--out_tsv",
        required=False,
        default=None
    )
    parser.add_argument(
        "--out_cluster",
        required=False,
        default=None
    )
    parser.add_argument(
        "--out_dendrogram",
        required=False,
        default=None
    )
    parser.add_argument(
        "--gtfdb",
        required=False,
        default=None
    )
    parser.add_argument(
        "--leaf_height",
        help="height at which to cut the dendrogram to form"
             " clusters (default: 5)",
        default=5
    )
    # Process arguments
    args = parser.parse_args()
    in_files = args.in_files
    in_labels = args.in_labels
    out_file = args.out_file
    l2fc = float(args.l2fc)
    padj = float(args.padj)
    out_tsv = args.out_tsv
    out_cluster = args.out_cluster
    gtfdb = args.gtfdb
    leaf_height = float(args.leaf_height)
    out_dendrogram = args.out_dendrogram

    if in_labels is None:
        in_labels = [os.path.basename(fn) for fn in in_files]

    if out_tsv is None:
        out_tsv = out_file + '.log2fc.tsv'
    if out_cluster is None:
        out_cluster = out_file + '.clusters.tsv'
    if out_dendrogram is None:
        out_dendrogram = out_file + '.dendrogram.svg'

    if gtfdb is not None:
        id2name_dict = p.gene_id_to_name(
            p.get_featuredb(gtfdb)
        )
    else:
        id2name_dict = None

    # main func
    make_plot(
        in_files, in_labels, id2name_dict, l2fc, padj, leaf_height, out_file, out_tsv, out_cluster, out_dendrogram
    )


if __name__ == "__main__":
    main()
