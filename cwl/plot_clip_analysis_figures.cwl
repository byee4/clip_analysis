#!/usr/bin/env cwltool

### doc: "clipper cwl tool (https://github.com/yeolab/clipper)" ###

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 8
    coresMax: 16
    ramMin: 32000
    #tmpdirMin: 10000
    #outdirMin: 10000

baseCommand: [plot_clip_analysis_figures.py]

inputs:

  l2fcwithpval_enr_r1:
    type: File
    inputBinding:
      position: 1
      prefix: --l2fcwithpval_enr_r1
    doc: "rep 1 region-level enrichment output file with pval and fold change"

  l2fcwithpval_enr_r2:
    type: File
    inputBinding:
      position: 1
      prefix: --l2fcwithpval_enr_r2
    doc: "rep 2 region-level enrichment output file with pval and fold change"

  inp_reads_by_loc_r1:
    type: File
    inputBinding:
      position: 1
      prefix: --inp_reads_by_loc_r1
    doc: "rep1 broad feature counts (readsByLoc) file"

  inp_reads_by_loc_r2:
    type: File
    inputBinding:
      position: 1
      prefix: --inp_reads_by_loc_r2
    doc: "rep2 broad feature counts (readsByLoc) file"

  pickle_r1:
    type: File
    inputBinding:
      position: 1
      prefix: --pickle_r1
    doc: "rep1 motif analysis/kmer enrichment output pickle"

  pickle_r2:
    type: File
    inputBinding:
      position: 1
      prefix: --pickle_r2
    doc: "rep2 motif analysis/kmer enrichment output pickle"

  rep_element_parsed_r1_ip:
    type: File
    inputBinding:
      position: 1
      prefix: --rep_element_parsed_r1_ip
    doc: "rep1 parsed repeat element mapping file"

  rep_element_parsed_r1_input:
    type: File
    inputBinding:
      position: 1
      prefix: --rep_element_parsed_r1_input
    doc: "rep1 parsed repeat element mapping file"

  rep_element_parsed_r2_ip:
    type: File
    inputBinding:
      position: 1
      prefix: --rep_element_parsed_r2_ip
    doc: "rep2 parsed repeat element mapping file"

  rep_element_parsed_r2_input:
    type: File
    inputBinding:
      position: 1
      prefix: --rep_element_parsed_r2_input
    doc: "rep2 parsed repeat element mapping file"

  annotated_files:
    type: File[]
    inputBinding:
      position: 1
      prefix: --pickle_r1
    doc: ".annotated files from annotator script."

  out_file:
    type: string
    inputBinding:
      position: 1
      prefix: --out_file

outputs:

  output_file:
    type: File
    outputBinding:
      glob: inputs.out_file


doc: |
  This tool plots outputs from clip_analysis scripts/workflows.
