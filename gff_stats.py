#!/usr/bin/python3

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# gff_stats.py - obtain basic gff3 statistics of features


#################################################

import sys

#################################################
usage = "gff_stats.py file.gff3"
#################################################

try:
    script,gff=sys.argv
except: sys.exit("Correct usage is:\n"+usage)

sys.stdout = sys.stderr

def getStats(file, field=3):
    gff = [l.strip().split('\t') for l in open(file).readlines() if l.count('\t') >7]
    parts = {t: [(int(y[4]) - int(y[3]) + 1) for y in gff if y[field-1] == t] for t in set([x[field-1] for x in gff])}
    print('%d types' % len(parts))
    print('type    \tqtd\tmean\tsum')
    for t in sorted(parts):
        print('%s\t%d\t%.1f\t%d' % (t.ljust(10, '.'), len(parts[t]), sum(parts[t]) / len(parts[t]), sum(parts[t])))

getStats(gff)

print('by mikeias.net')
