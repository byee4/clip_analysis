#!/bin/bash

conda create -y -n clip_analysis python=2.7;

source activate clip_analysis;

conda install -y matplotlib seaborn numpy pandas;
conda install -y bedtools pybedtools;
conda install -y mock;
