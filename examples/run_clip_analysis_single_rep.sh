#!/bin/bash

clip_analysis \
--l2fcwithpval_enr_r1 inputs/204_01_RBFOX2.l2fc_significant_regioncalls.txt \
--inp_reads_by_loc_r1 inputs/204_01_RBFOX2.input.broadfeaturecounts.txt \
--pickle_r1 inputs/204_01_RBFOX2.pickle \
--annotated_files \
inputs/204_01_RBFOX2 \
--rep_element_parsed_r1_ip inputs/204_01_RBFOX2.combined_w_uniquemap.rmDup.sam.parsed \
--rep_element_parsed_r1_input inputs/204_INPUT_RBFOX2.combined_w_uniquemap.rmDup.sam.parsed \
--regions CDS intron 3utr \
--out_file outputs/all_figures_single_rep.png
