#!/usr/bin/env bash
# 02-tables.sh
# Kamil Slowikowski
# February 27, 2015
#
# Get the URLs for the TF target tables.

cat html/*.html \
    | grep -P -o 'system_data\S+_regulated_gene.html' \
    | sort -u > tables.txt
