#!/usr/bin/env bash

#note that the attributes.txt file needs to be created beforehand

#read in variables
source common_variables.sh

PolarityCalculation -f $RiboMiner_dir/LSU_attributes.txt -c $RiboMiner_annot/longest.transcripts.info.txt -o $RiboMiner_dir/LSU_polarity -n 50
PolarityCalculation -f $RiboMiner_dir/D30_attributes.txt -c $RiboMiner_annot/longest.transcripts.info.txt -o $RiboMiner_dir/D30_polarity -n 50


#-n 50 is default for at least 50 counts
