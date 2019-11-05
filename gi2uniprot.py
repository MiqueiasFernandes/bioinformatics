#!/usr/bin/python3

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# gi2uniprot.py - convert gi to uniprot id with uniprot api


#################################################

import sys
import urllib.parse
import urllib.request

#################################################
usage = "gi2uniprot.py gis.txt"
#################################################

try:
    script,gis=sys.argv
except: sys.exit("Correct usage is:\n"+usage)

sys.stdout = sys.stderr

def mapGI2Uniprot(gis):
    entries = set(gis)
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': 'P_GI',
    'to': 'ID',
    'format': 'tab',
    'query': ' '.join(entries)
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    res = {x: [] for x in entries}
    with urllib.request.urlopen(req) as f:
        res.update({x[0]:x[1] for x in 
                    [x.split('\t') for x in f.read().decode('utf-8').split('\n')] if x[0].isdigit()})
    return res

for gi,uniprot in mapGI2Uniprot([l.strip() for l in open(gis).readlines() if len(l) > 3 and l.strip().isdigit()]).items():
    print("%s\t%s" % (gi, uniprot))

# by mikeias.net
