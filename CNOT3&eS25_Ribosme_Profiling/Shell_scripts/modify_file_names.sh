#!/usr/bin/env bash


for file in $fastq_dir/*_S*_R1_001.fastq.gz; do
  base=$(basename "$file")
  newname=$(echo "$base" | sed -E 's/_S[0-9]+_R1_001//')
  mv "$file" "$fastq_dir/$newname"
done