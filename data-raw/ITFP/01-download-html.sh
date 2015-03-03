#!/usr/bin/env bash
# download.sh
# Kamil Slowikowski
# February 27, 2015
#
# TODO Write a description

URL='http://itfp.biosino.org/itfp/TFViewer/human_tf.jsp?table=human_tf&page='

mkdir -p html

for i in {1..206}; do
    echo "$(date) -- $i"
    sleep 1
    wget -o html/${i}.log -O html/${i}.html "${URL}${i}"
done
