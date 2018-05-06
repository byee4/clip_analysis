#!/bin/bash

python plot_repetitive_elements_sunburst.py \
--ip_parsed inputs/204_01_RBFOX2.combined_w_uniquemap.rmDup.sam.parsed \
--input_parsed inputs/204_INPUT_RBFOX2.combined_w_uniquemap.rmDup.sam.parsed \
--color_file inputs/color_list_269.lines \
--title "sunburst plot for RBFOX2 204 01" \
--out_file outputs/sunburst_example.png
