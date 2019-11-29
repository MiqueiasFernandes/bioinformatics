#!/usr/bin/python

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# marks_anchors.py - a script for mark sintenic regions surround specific CDS in anchors of https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)


#################################################

import sys

#################################################
usage = "marks_anchors.py anchors.simple genome.gff3 cds.ids"
#################################################

try: 
    script,anc,gff3,cds=sys.argv
except: sys.exit("Correct usage is:\n"+usage)

gff = [l.strip().split('\t') for l in open(gff3).readlines() if l.count('mRNA') > 0]
mrnas = {x[8].split('ID=')[1].split(';')[0]:(x[0],int(x[3]),int(x[4])) for x in gff}
anchors = [l.strip().split('\t') for l in open(anc).readlines()]
tps = [l.strip() for l in open(cds).readlines()]

with open(anc + '.marked', 'w') as o:
    for a in anchors:
        m1 = mrnas[a[0]]
        m2 = mrnas[a[1]]
        n_m = True
        c = [a]
        for t in tps:
            m = mrnas[t]
            s = m[0]
            i = m[1]
            f = m[2]
            if m1[0] == s and m1[1] <= i and m2[2] >= f:
                c.append(t)
                if n_m:
                    o.write('g*' + '\t'.join(a) + '\n')
                    n_m = False
        if n_m:
            o.write('\t'.join(a) + '\n')
            n_m = False
        if len(c) > 2:
            print('regiao sintenica comeca em: ' + c[0][0] + 
                  ' na goiaba (' + m1[0] + ') termina em ' + c[0][1] + 
                    ' e comeca no outro em ' + c[0][2] + ' a ' + c[0][3] + ' com: ' +  ', '.join(c[1:]))


