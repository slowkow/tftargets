#!/usr/bin/env bash
# 04-html-to-text.sh
# Kamil Slowikowski
# February 27, 2015
#
# Convert all TF targets HTML tables to text.

parallel 'pandoc -f html -t plain {}' ::: tables-html/*.html > regulated_genes.txt
