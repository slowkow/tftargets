#!/usr/bin/env python
"""
01-download-factors.py
Kamil Slowikowski
February 27, 2015

Download a table of all the transcription factors in the database.
"""

from bs4 import BeautifulSoup as bs
import urllib
import re
import time

URL = 'https://cb.utdallas.edu/cgi-bin/TRED/tred.cgi' \
    '?hitCount=179&process=dataBrowseResult&browseFactor=Browse+Factors'\
    '&start={start}&anotherPage=GO'

def main():
    print "\t".join(["Primary_Name", "All_Names", "Species", "FID"])
    for start in xrange(0, 473, 20):
        time.sleep(1)
        rows = get_rows(start)
        print_rows(rows)

def get_fid(row):
    m = re.search(r'fid=(\d+)', str(row))
    if m:
        return m.groups()[0]
    return None

def print_rows(rows):
    for row in rows:
        # Convert row to tsv text
        fid = get_fid(row)
        line = row.text.strip().replace('\n', '\t')
        if fid:
            line += '\t' + str(fid)
        print line

def get_rows(start):
    html = urllib.urlopen(URL.format(start=start)).read()
    b = bs(html)
    rows = b.find_all('tr')
    rows = filter(lambda x: 'tred.cgi?process=factorInfo' in str(x), rows)
    return rows

if __name__ == '__main__':
    main()
