#!/usr/bin/env bash
# 05-text-to-tsv.sh
# Kamil Slowikowski
# February 27, 2015
#
# Convert the text file to TSV.

(
echo -e "TF_name\tTF_target\tscore\ttarget_type"
perl -ne 'if (/^  \w/) { if(!/TF_name/) { s/^  //; s/ +/\t/g; print; } }' \
    regulated_genes.txt 
) > regulated_genes.tsv

