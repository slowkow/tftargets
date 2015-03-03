#!/usr/bin/env bash
# download02.sh
# Kamil Slowikowski
# February 27, 2015
#
# Download the HTML for each TF target table.

BASE='http://itfp.biosino.org/itfp/'

mkdir -p tables-html

cat tables.txt | while read line; do
    out="tables-html/$(basename $line)"
    log="tables-html/$(basename $line .html).log"
    echo "$(date) -- $out"
    sleep 1
    wget -o $log -O $out "${BASE}${line}"
done

