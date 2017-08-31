# CLIPSEQ

## Python3:
```
conda create -n clip_analysis python=3.5
source activate py3_clip_analysis
conda install bedtools pybedtools numpy pandas seaborn future
```

## Python2:
```
conda create -n clip_analysis python=2.7
source activate py2_clip_analysis
conda install bedtools pybedtools numpy pandas seaborn future
```

Usage:

```
import os
import glob
from clip_analysis import plotting_methods as pm

### plot the region distributions given eric's annotated input norm files

wd = '/path/to/input_norm_annotated_bedfiles/'
annotated_files = glob.glob(os.path.join(wd,'*.annotated'))
distribution_png = '/path/to/distribution.png'
l10p = 3
l2fc = 3

pm.filter_and_plot_region_distribution(
    annotated_files, distribution_png, l10p, l2fc, format='eric'
)

### plot Figure 2b from Ann's IMP paper

palette = sns.color_palette("hls", 8)

field_list = {
    'noncoding_exon':palette[0],
    '3utr':palette[1],
    '5utr':palette[2],
    'intron':palette[3],
    'noncoding_intron':palette[4],
    'CDS':palette[5]
}

input_reads = '/path/to/input.reads_by_loc.csv'
ip_l2fc = '/path/to/l2fcwithpval_enr.csv'
scatter_png = '/path/to/scatter.png'

pm.plot_ip_foldchange_over_input_reads(
    ip_l2fc, input_reads, scatter_png, field_list
)
```