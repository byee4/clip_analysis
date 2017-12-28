# CLIP ANALYSIS

makes some useful clip analysis figures for QC and other things.

# Requirements:
- HOMER
- EMBOSS
- bedtools
- python 2.7
- matplotlib
- seaborn
- numpy
- pandas
- pybedtools

(or run the run_environments.sh script - don't forget to ```source activate clip_analysis```)

# GATK -> eCLIP 0.1.x
- made with table markdown generator (https://www.tablesgenerator.com/markdown_tables#)

|                                     | GATK                                                        | eCLIP 0.1.7                                                                                |
|-------------------------------------|-------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| Demuxed + adapter trimmed reads     | ```204.01_RBFOX2.A01.r*.fqTrTr.fqgz```                      | ```RBFOX2-204-CLIP_S1_R*.A01_204_01_RBFOX2.adapterTrim.round2.fastq.gz```                  |
| Repetitive element filtered reads   | ```204.01_RBFOX2.A01.r-.fqTrTrU*.fq```                      | ```RBFOX2-204-CLIP_S1_R1.A01_204_01_RBFOX2.adapterTrim.round2.rep.bamUnmapped.out.mate*``` |
| Unique genome aligned reads         | ```204.01_RBFOX2.A01.r-.fqTrTrU-SoMaSo.bam```               | ```RBFOX2-204-CLIP_S1_R1.A01_204_01_RBFOX2.adapterTrim.round2.rmRep.bam```                 |
| PCR duplicate removed aligned reads | ```204.01_RBFOX2.A01.r-.fqTrTrU-SoMaSoCpSo.bam```           | ```RBFOX2-204-CLIP_S1_R1.A01_204_01_RBFOX2.adapterTrim.round2.rmRep.rmDup.sorted.bam```    |
| Barcode merged alignments           | ```204.01_RBFOX2.---.r-.fqTrTrU-SoMaSoCpSoMeV2.bam```       | ```204_01_RBFOX2.merged.r2.bam```                                                          |
| CLIPper peaks                       | ```204.01_RBFOX2.---.r-.fqTrTrU-SoMaSoCpSoMeV2Cl.bed```     | ```204_01_RBFOX2.merged.r2.peaks.bed```                                                    |
| Input-normalized peaks              | ```204.01_RBFOX2.---.r-.fqTrTrU-SoMaSoCoSoMeV2ClNpCo.bed``` | ```204_01.basedon_204_01.peaks.l2inputnormnew.bed.compressed.bed```                        |

# Usage:
```
clip_analysis \
--l2fcwithpval_enr_r1 inputs/204_01_RBFOX2.l2fc_significant_regioncalls.txt \  # output from region_based_enrichment (l2fcWithPvalEnr)
--l2fcwithpval_enr_r2 inputs/204_02_RBFOX2.l2fc_significant_regioncalls.txt \
--inp_reads_by_loc_r1 inputs/204_01_RBFOX2.input.broadfeaturecounts.txt \  # output from region_based_enrichment (inputBroadFeatures)
--inp_reads_by_loc_r2 inputs/204_02_RBFOX2.input.broadfeaturecounts.txt \
--pickle_r1 inputs/204_01_RBFOX2.pickle \  # output from clip_analysis_legacy (--pickle FILE.pickle)
--pickle_r2 inputs/204_02_RBFOX2.pickle \
--annotated_files \
inputs/204_01_RBFOX2 \  # output from annotator (*.annotated)
inputs/204_02_RBFOX2 \
--rep_element_parsed_r1_ip inputs/204_01_RBFOX2.combined_w_uniquemap.rmDup.sam.parsed \  # output from rep_element_pipeline (parsed file)
--rep_element_parsed_r1_input inputs/204_INPUT_RBFOX2.combined_w_uniquemap.rmDup.sam.parsed \
--rep_element_parsed_r2_ip inputs/204_02_RBFOX2.combined_w_uniquemap.rmDup.sam.parsed \
--rep_element_parsed_r2_input inputs/204_INPUT_RBFOX2.combined_w_uniquemap.rmDup.sam.parsed \
--regions CDS intron 3utr \  # specify regions of interest (default: CDS intron 3utr 5utr)
--out_file outputs/all_figures.png  # output image file
```

(or see examples/run_clip_analysis.sh script)