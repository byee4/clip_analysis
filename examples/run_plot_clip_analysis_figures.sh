#!/bin/bash

python ../clip_analysis/plot_clip_analysis_figures.py \
--l2fcwithpval_enr_r1 data/rep1.l2fcwithpval_enr.txt \
--l2fcwithpval_enr_r2 data/rep2.l2fcwithpval_enr.txt \
--inp_reads_by_loc_r1 data/rep1.input_reads_by_loc.txt \
--inp_reads_by_loc_r2 data/rep2.input_reads_by_loc.txt \
--annotated_files data/rep1.annotated.bed data/rep2.annotated.bed data/some_reference.annotated.bed \
--out_file all_figures.png
