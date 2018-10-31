#!/usr/bin/env cwltool

doc: |
  The main workflow that produces two reproducible peaks via IDR given
  two eCLIP samples (1 input, 1 IP each).

cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement


inputs:


  ### REQUIRED IDR PEAKS INPUTS ###


  rep1_clip_bam_file:
    type: File
  rep1_input_bam_file:
    type: File
  rep1_peaks_bed_file:
    type: File

  rep2_clip_bam_file:
    type: File
  rep2_input_bam_file:
    type: File
  rep2_peaks_bed_file:
    type: File

  merged_peaks_bed:
    type: string
  merged_peaks_custombed:
    type: string


  ### ANNOTATOR INPUTS ###

  species:
    type: File
  gtfdb_file:
    type: File
  rep1_peak_file:
    type: File

  rep2_peak_file:
    type: File

  ### REQUIRED REPEAT MAPPING INPUTS ###


  rep1_barcode1r1FastqGz:
    type: File
  rep1_barcode1r2FastqGz:
    type: File
  rep1_barcode1rmRepBam:
    type: File

  rep1_barcode2r1FastqGz:
    type: File
  rep1_barcode2r2FastqGz:
    type: File
  rep1_barcode2rmRepBam:
    type: File

  rep1_barcode1Inputr1FastqGz:
    type: File
  rep1_barcode1Inputr2FastqGz:
    type: File
  rep1_barcode1InputrmRepBam:
    type: File

  rep2_barcode1r1FastqGz:
    type: File
  rep2_barcode1r2FastqGz:
    type: File
  rep2_barcode1rmRepBam:
    type: File

  rep2_barcode2r1FastqGz:
    type: File
  rep2_barcode2r2FastqGz:
    type: File
  rep2_barcode2rmRepBam:
    type: File

  rep2_barcode1Inputr1FastqGz:
    type: File
  rep2_barcode1Inputr2FastqGz:
    type: File
  rep2_barcode1InputrmRepBam:
    type: File

  bowtie2_db:
    type: Directory
  fileListFile1:
    type: File
  fileListFile2:
    type: File

  gencodeGTF:
    type: File
  gencodeTableBrowser:
    type: File
  repMaskBEDFile:
    type: File


  ### REGION LEVEL ENRICHMENT INPUTS (USE IDR INPUTS) ###


  ### MOTIF ANALYSIS INPUTS ###

  ### SVG ###

  svg_output:
    type: string

outputs:

  merged_peaks_file:
    type: File
    outputSource: step_merge_peaks/merged_peaks_bed_file
  annotated_rep1_peaks_file:
    type: File
    outputSource: step_rep1_clip_analysis/output_annotated
  annotated_rep2_peaks_file:
    type: File
    outputSource: step_rep2_clip_analysis/output_annotated
  annotated_merged_peaks_file:
    type: File
    outputSource: step_annotate_idr_peaks/output_annotated

  rep1_ip_parsed_file:
    type: File
    outputSource: step_rep1_clip_analysis/output_repeat_ip_parsed_file
  rep2_ip_parsed_file:
    type: File
    outputSource: step_rep2_clip_analysis/output_repeat_ip_parsed_file

  rep1_input_parsed_file:
    type: File
    outputSource: step_rep1_clip_analysis/output_repeat_input_parsed_file
  rep2_input_parsed_file:
    type: File
    outputSource: step_rep2_clip_analysis/output_repeat_input_parsed_file

  rep1_l2fc_with_pval_enr_file:
    type: File
    outputSource: step_rep1_clip_analysis/output_l2fc_with_pval_enr_file
  rep2_l2fc_with_pval_enr_file:
    type: File
    outputSource: step_rep2_clip_analysis/output_l2fc_with_pval_enr_file

  rep1_pickle_file:
    type: File
    outputSource: step_rep1_clip_analysis/output_pickle_file
  rep2_pickle_file:
    type: File
    outputSource: step_rep2_clip_analysis/output_pickle_file

  svg_output_file:
    type: File
    outputSource: step_plot_clip_analysis_figures/output_file




steps:

  step_merge_peaks:
    run: wf_get_reproducible_eclip_peaks.cwl
    in:
      rep1_clip_bam_file: rep1_clip_bam_file
      rep2_clip_bam_file: rep2_clip_bam_file
      rep1_input_bam_file: rep1_input_bam_file
      rep2_input_bam_file: rep2_input_bam_file
      rep1_peaks_bed_file: rep1_peaks_bed_file
      rep2_peaks_bed_file: rep2_peaks_bed_file
      merged_peaks_bed: merged_peaks_bed
      merged_peaks_custombed: merged_peaks_custombed
    out:
      - rep1_clip_read_num
      - rep2_clip_read_num
      - rep1_input_read_num
      - rep2_input_read_num
      - rep1_input_normed_bed
      - rep2_input_normed_bed
      - rep1_entropy_full
      - rep2_entropy_full
      - idr_output
      - idr_output_bed
      - rep1_idr_output_input_normed_bed
      - rep2_idr_output_input_normed_bed
      - rep1_idr_output_input_normed_full
      - rep2_idr_output_input_normed_full
      - rep1_reproducing_peaks_full
      - rep2_reproducing_peaks_full
      - merged_peaks_bed_file
      - merged_peaks_custombed_file
      - reproducing_peaks_count

  step_annotate_idr_peaks:
    run: annotate.cwl
    in:
      input: step_merged_peaks/merged_peaks_bed_file
      gtfdb: gtf_db
    out:
      - output_annotated

  ### Replicate 1 ###

  step_rep1_clip_analysis:
    run: wf_clip_analysis.cwl
    in:
      species: species
      peak_file: rep1_peak_file
      gtfdb_file: gtfdb_file
      barcode1r1FastqGz: rep1_barcode1r1FastqGz
      barcode1r2FastqGz: rep1_barcode1r2FastqGz
      barcode1rmRepBam: rep1_barcode1rmRepBam
      barcode2r1FastqGz: rep1_barcode2r1FastqGz
      barcode2r2FastqGz: rep1_barcode2r2FastqGz
      barcode2rmRepBam: rep1_barcode2rmRepBam
      barcode1Inputr1FastqGz: rep1_barcode1Inputr1FastqGz
      barcode1Inputr2FastqGz: rep1_barcode1Inputr2FastqGz
      barcode1InputrmRepBam: rep1_barcode1InputrmRepBam
      bowtie2_db: bowtie2_db
      fileListFile1: fileListFile1
      fileListFile2: fileListFile2
      gencodeGTF: gencodeGTF
      gencodeTableBrowser: gencodeTableBrowser
      repMaskBEDFile: repMaskBEDFile
      clip_bam_file: rep1_clip_bam_file
      input_bam_file: rep1_input_bam_file
    out:
      - output_repeat_ip_parsed_file
      - output_repeat_input_parsed_file
      - output_ip_concatenated_pre_rmDup_sam_file
      - output_input_concatenated_pre_rmDup_sam_file
      - output_ip_concatenated_rmDup_sam_file
      - output_input_concatenated_rmDup_sam_file
      - output_clip_broad_feature_counts_file
      - output_input_broad_feature_counts_file
      - output_l2fc_with_pval_enr_file
      - output_pickle_file
      - output_homer_dir
      - output_homer_svg_file
      - output_annotated

  ### Replicate 2 ###


  step_rep2_clip_analysis:
    run: wf_clip_analysis.cwl
    in:
      species: species
      peak_file: rep2_peak_file
      gtfdb_file: gtfdb_file
      barcode1r1FastqGz: rep2_barcode1r1FastqGz
      barcode1r2FastqGz: rep2_barcode1r2FastqGz
      barcode1rmRepBam: rep2_barcode1rmRepBam
      barcode2r1FastqGz: rep2_barcode2r1FastqGz
      barcode2r2FastqGz: rep2_barcode2r2FastqGz
      barcode2rmRepBam: rep2_barcode2rmRepBam
      barcode1Inputr1FastqGz: rep2_barcode1Inputr1FastqGz
      barcode1Inputr2FastqGz: rep2_barcode1Inputr2FastqGz
      barcode1InputrmRepBam: rep2_barcode1InputrmRepBam
      bowtie2_db: bowtie2_db
      fileListFile1: fileListFile1
      fileListFile2: fileListFile2
      gencodeGTF: gencodeGTF
      gencodeTableBrowser: gencodeTableBrowser
      repMaskBEDFile: repMaskBEDFile
      clip_bam_file: rep2_clip_bam_file
      input_bam_file: rep2_input_bam_file
    out:
      - output_repeat_ip_parsed_file
      - output_repeat_input_parsed_file
      - output_ip_concatenated_pre_rmDup_sam_file
      - output_input_concatenated_pre_rmDup_sam_file
      - output_ip_concatenated_rmDup_sam_file
      - output_input_concatenated_rmDup_sam_file
      - output_clip_broad_feature_counts_file
      - output_input_broad_feature_counts_file
      - output_l2fc_with_pval_enr_file
      - output_pickle_file
      - output_homer_dir
      - output_homer_svg_file
      - output_annotated

  ### PLOT STUFF ###

  step_plot_clip_analysis_figures:
    run: plot_clip_analysis_figures.cwl
    in:
      l2fcwithpval_enr_r1: step_rep1_clip_analysis/output_l2fc_with_pval_enr_file
      l2fcwithpval_enr_r2: step_rep2_clip_analysis/output_l2fc_with_pval_enr_file
      inp_reads_by_loc_r1: step_rep1_clip_analysis/output_input_broad_feature_counts_file
      inp_reads_by_loc_r2: step_rep2_clip_analysis/output_input_broad_feature_counts_file
      pickle_r1: step_rep1_clip_analysis/output_pickle_file
      pickle_r2: step_rep2_clip_analysis/output_pickle_file
      rep_element_parsed_r1_ip: step_rep1_clip_analysis/output_repeat_ip_parsed_file
      rep_element_parsed_r2_ip: step_rep2_clip_analysis/output_repeat_ip_parsed_file
      rep_element_parsed_r1_input: step_rep1_clip_analysis/output_repeat_input_parsed_file
      rep_element_parsed_r2_input: step_rep2_clip_analysis/output_repeat_input_parsed_file
      annotated_files: [
        step_annotate_idr_peaks/output_annotated,
        step_rep1_clip_analysis/output_annotated,
        step_rep2_clip_analysis/output_annotated
      ]
      output_file: svg_output

    out:
      - output_file