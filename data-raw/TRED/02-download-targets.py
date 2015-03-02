#!/usr/bin/env python
"""
02-download-targets.py
Kamil Slowikowski
February 27, 2015

TODO Write description
"""

from bs4 import BeautifulSoup as bs
import urllib
import re
import time
import sys

URL = 'https://cb.utdallas.edu/cgi-bin/TRED/tred.cgi?process=searchTFGene&sel_type=factor_name&factor_organism=any&tx_search_terms={name}&target_organism=human&prom_quality=1&prom_quality=2&prom_quality=3&prom_quality=4&prom_quality=5&bind_quality=0&submit=SEARCH&start={start}'

# html = urllib.urlopen(URL.format(start=0, name="STAT3")).read()
# b = bs(html)
# rows = b.find_all('tr')
# rows = filter(lambda x: 'tred.cgi?process=geneInfo' in str(x), rows)

def main():
    lines = open('factors.tsv').readlines()
    fields = [x.split('\t') for x in lines]
    factors = [x[0] for x in fields if 'human' in x[2]]
    for i, factor in enumerate(factors, 1):
        log('{} / {} -- {}'.format(i, len(factors), factor))
        time.sleep(1)
        rows = get_rows(factor)
        print_rows(factor, rows)

def log(msg):
    sys.stderr.write('{} -- {}\n'.format(time.asctime(), msg))

def get_gid(row):
    m = re.search(r'gid=(\d+)', str(row))
    if m:
        return m.groups()[0]
    return None

def print_rows(name, rows):
    for row in rows:
        # Convert row to tsv text
        fid = get_gid(row)
        line = row.text.strip().replace('\n', '\t')
        if fid:
            line += '\t' + str(fid)
        print name + '\t' + line

def num_targets(rows):
    items = filter(lambda x: ' to ' in x.text and ' of ' in x.text, rows)
    if not items:
        return None
    m = re.search(r'(\d+) to (\d+) of (\d+)', items[0].text)
    if not m:
        return None
    return int(m.groups()[2])

def get_rows(name):
    html = urllib.urlopen(URL.format(start=0, name=name)).read()
    b = bs(html)
    retval = []
    rows = b.find_all('tr')
    # Find the number of targets for this TF.
    n_targets = num_targets(rows)
    if not n_targets:
        return retval
    # Extract table rows that have data.
    data_row = lambda x: 'tred.cgi?process=geneInfo' in str(x)
    data_rows = filter(data_row, rows)
    retval.extend(data_rows)
    # If there are more than 20 targets, we need to send more requests.
    if n_targets > 20:
        for i in xrange(20, n_targets, 20):
            html = urllib.urlopen(URL.format(start=i, name=name)).read()
            b = bs(html)
            rows = b.find_all('tr')
            data_rows = filter(data_row, rows)
            retval.extend(data_rows)
    return retval

if __name__ == '__main__':
    main()
