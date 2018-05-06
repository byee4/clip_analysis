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

  # REQUIRED IDR PEAKS INPUTS

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

  # REQUIRED REPEAT MAPPING INPUTS
  barcode1r1FastqGz:
    type: File
  barcode1r2FastqGz:
    type: File
  barcode1rmRepBam:
    type: File

  barcode2r1FastqGz:
    type: File
  barcode2r2FastqGz:
    type: File
  barcode2rmRepBam:
    type: File

  barcode1Inputr1FastqGz:
    type: File
  barcode1Inputr2FastqGz:
    type: File
  barcode1InputrmRepBam:
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

  # ANNOTATOR DB
  dataset:
    type: string

  gtf_db:
    type: File

  # FINAL IDR OUTPUTS

  merged_peaks_bed:
    type: string
  merged_peaks_custombed:
    type: string


outputs:


  rep1_clip_read_num:
    type: File
    outputSource: rep1_input_norm_and_entropy/clip_read_num


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
      input: merged_peaks_bed_file
      gtfdb: gtf_db

  step_annotate_rep1_peaks:
    run: annotate.cwl
    in:
      input: rep2_peaks_bed_file
      gtfdb: gtf_db

  step_annotate_rep2_peaks:
    run: annotate.cwl
    in:
      input: rep2_peaks_bed_file
      gtfdb: gtf_db

  step_repeat_mapping:
    run: wf_ecliprepmap_pe.cwl
    in:
      barcode1r1FastqGz:
      barcode1r2FastqGz:
      barcode1rmRepBam:
      barcode2r1FastqGz:
      barcode2r2FastqGz:
      barcode2rmRepBam:
      barcode1Inputr1FastqGz:
      barcode1Inputr2FastqGz:
      barcode1InputrmRepBam:
      bowtie2_db:
      fileListFile1:
      fileListFile2:
      gencodeGTF:
      gencodeTableBrowser:
      repMaskBEDFile:
    out:
      - output_barcode1_concatenated_pre_rmDup_sam_file
      - output_barcode2_concatenated_pre_rmDup_sam_file
      - output_ip_concatenated_pre_rmDup_sam_file
      - output_input_concatenated_pre_rmDup_sam_file
      - output_barcode1_concatenated_rmDup_sam_file
      - output_barcode2_concatenated_rmDup_sam_file
      - output_ip_concatenated_rmDup_sam_file
      - output_input_concatenated_rmDup_sam_file
      - output_ip_parsed
      - output_input_parsed

  step_regionlevel_enrichment:
    run: wf_region_based_enrichment.cwl
    in:
      clipBamFile: rep1_clip_bam_file
      inputBamFile: rep1_input_bam_file
      gencodeGTFFile: gencodeGTF
      gencodeTableBrowserFile: gencodeTableBrowser
    out:
      - clipBroadFeatureCountsFile
      - inputBroadFeatureCountsFile
      - combinedOutputFile
      - l2fcWithPvalEnrFile
      - l2fcFile