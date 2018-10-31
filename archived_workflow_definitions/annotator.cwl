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

baseCommand: [annotator]

inputs:

  species:
    type: string
    inputBinding:
      position: 0
      prefix: --species
    doc: "species: one of ce10 ce11 dm3 hg19 GRCh38 mm9 mm10"

  gtfdb_file:
    type: File
    inputBinding:
      position: 1
      prefix: --gtfdb
    doc: "sqlite database file created using gffutils createdb()"

  peak_file:
    type: File
    inputBinding:
      position: 1
      prefix: --bed
    doc: "bed6 file describing peak regions"

  out_file:
    type: string
    default: ""
    inputBinding:
      position: 10
      prefix: --outfile
      valueFrom: |
        ${
          if (inputs.out_file == "") {
            return inputs.peak_file.nameroot + ".annotated";
          }
          else {
            return inputs.out_file;
          }
        }

outputs:

  output_annotated:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.out_file == "") {
            return inputs.peak_file.nameroot + ".annotated";
          }
          else {
            return inputs.out_file;
          }
        }


doc: |
  This tool annotates a BED6 file using a region priority.
