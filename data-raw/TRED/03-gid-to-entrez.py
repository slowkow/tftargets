#!/usr/bin/env python
"""
03-gid-to-entrez.py
Kamil Slowikowski
March  4, 2015

TODO Write description
"""

from bs4 import BeautifulSoup as bs
import urllib
import re
import time
import sys
import os

URL = 'https://cb.utdallas.edu/cgi-bin/TRED/tred.cgi?' \
    'process=geneInfo&gid={gid}'

# html = urllib.urlopen(URL.format(gid='6850')).read()
# b = bs(html)
# rows = b.find_all('tr')
# rows = filter(lambda x: 'tred.cgi?process=geneInfo' in str(x), rows)

def main():
    # Get a list of gids
    lines = open('targets.tsv').readlines()
    fields = [x.strip().split('\t') for x in lines]
    gids = [x[7] for x in fields if 'human' in x[2]]
    gids = list(set(gids))

    # Save HTML
    if not os.path.exists('gid-html'):
        os.mkdir('gid-html')

    # Print Entrez Gene IDs
    print "gid\tentrez"
    for i, gid in enumerate(gids, 1):
        log('{} / {} -- {}'.format(i, len(gids), gid))
        time.sleep(1)
        html = urllib.urlopen(URL.format(gid=gid)).read()
        # Save the HTML
        with open('gid-html/{}.html'.format(gid), 'w') as out:
            out.write(html)
        # Find the Entrez Gene ID
        m = re.search(r'LocusLink\S+=(\d+)', html)
        entrez = ''
        if m:
            entrez = m.groups()[0]
        print "{}\t{}".format(gid, entrez)

def log(msg):
    sys.stderr.write('{} -- {}\n'.format(time.asctime(), msg))

if __name__ == '__main__':
    main()
