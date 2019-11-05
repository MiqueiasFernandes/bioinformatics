#!/usr/bin/python3

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# gff_remap.py - remap gff3 coords by agp from ALLMAPS


#################################################

import sys
import datetime
import re
from Bio import SeqIO

sys.path.append('/home/cluster/shared/pybio')
from SeqTools import rc
from Gff3 import GFF

#################################################
usage = "gff_remap.py agp_file old.fasta new.fasta old.gff3 new.gff3"
#################################################

try:
    script, agp, old_fasta, new_fasta, old_gff, new_gff=sys.argv
except: sys.exit("Correct usage is:\n"+usage)

def updateGFF(agp, fastaold, fastanew, gffold, gffnew, mapIDS={}):
    print(datetime.datetime.now())
    agp = [l.strip().split('\t') for l in open(agp).readlines() if l.count('\t') > 1]
    chr2scf = {x: [(z[5], int(z[1]), int(z[2]), z[8].replace('?', '+'), int(z[2]) - int(z[1]) + 1) 
                   for z in agp if z[0] == x and z[5].startswith('pg.scf')] for x in set([y[0] for y in agp])}
    fasta1 = SeqIO.to_dict(SeqIO.parse(fastaold, 'fasta'))
    fasta2 = SeqIO.to_dict(SeqIO.parse(fastanew, 'fasta'))
    
    ok, nok = 0 , 0
    

    somar = {}
    subt = {}

    for c in chr2scf:
        for s in chr2scf[c]:
            a = str(fasta2[c][int(s[1])-1:int(s[2])].seq)
            b = rc(str(fasta1[s[0]].seq), s[3])
            if a != b:
                print('ERR', c, s)
                nok += 1
            else: ## as sequencias são identicas
                if s[0] == c: ## possuem mesmo nome
                    if s[3] == '+': ## as sequencias possuem mesma orientação
                        ok += 1
                    else:
                        subt[s[0]] = [c, s[2] + 1, s[3]]
                        ok += 1
                else:  ## scaffold dentro de contig
                    if s[3] == '+':
                        somar[s[0]] = [c, s[1] - 1, s[3]]
                    else:
                        subt[s[0]] = [c, s[2] + 1, s[3]]
                    ok += 1
        
    print('OK: ', ok, 'ERR: ', nok)
    gff = [l.strip().split('\t') for l in open(gffold).readlines() if l.count('\t') > 7]

    for e in gff:
        if e[0] in somar:
            x = somar[e[0]]
            e[0] = x[0]
            e[3] = str(int(e[3]) + x[1])
            e[4] = str(int(e[4]) + x[1])
        elif e[0] in subt:
            x = subt[e[0]]
            e[0] = x[0]
            a = abs(x[1] - int(e[3]))
            b = abs(x[1] - int(e[4]))
            e[3] = str(min(a, b))
            e[4] = str(max(a, b))
            e[6] = '+' if e[6] == '-' else '-'
        if e[0] in mapIDS:
            e[0] = mapIDS[e[0]]
    print(datetime.datetime.now())
    print('persistindo novo gff ...')
    with open(gffnew, 'w') as o:
        o.write('##gff-version 3\n')
        o.write('\n'.join(['\t'.join(x) for x in sorted(gff,
                                                       key=lambda e: (1 if e[0].startswith('Chr') else 2, 
                                          int(re.sub("\D*", "", e[0])), int(e[3])))]) + '\n')
    
    if len(mapIDS) > 1:
        fastanew = fastanew + '.new.fasta'
        print('persistindo novo fasta ...')
        seqs = [fasta2[x] for x in fasta2 if not x in mapIDS]
        for k in mapIDS:
            a = fasta2[k]
            a.id = mapIDS[k]
            a.name = a.description = ''
            seqs.append(a)
        SeqIO.write(sorted(seqs, 
                           key=lambda e: (1 if e.id.startswith('Chr') else 2, 
                                          int(re.sub("\D*", "", e.id)))), fastanew, "fasta")
    print(datetime.datetime.now())
    print('validando ...')
    return validar(fastaold, gffold, fastanew, gffnew)
        
                    
                    
def validar(fastaold, gffold, fastanew, gffnew):
    from Bio import SeqIO
    from SeqTools import rc
    errs = []
    print('importando fastas ...')
    fasta1 = SeqIO.to_dict(SeqIO.parse(fastaold, 'fasta'))
    fasta2 = SeqIO.to_dict(SeqIO.parse(fastanew, 'fasta'))
    print('importando gffs ...')
    gff1 = [l.strip().split('\t') for l in open(gffold).readlines() if l.count('\t') > 7]
    gff2 = [l.strip().split('\t') for l in open(gffnew).readlines() if l.count('\t') > 7]
    print('analisando ...')
    print(datetime.datetime.now())
    prts1 = {x[8].split('ID=')[1].split(';')[0]:[x[0], int(x[3]), int(x[4]), x[6]] for x in gff1}
    prts2 = {x[8].split('ID=')[1].split(';')[0]:[x[0], int(x[3]), int(x[4]), x[6]] for x in gff2}
    pc = 0
    po = 0
    ids = set(list(prts1.keys()) + list(prts2.keys()))
    for i in ids:
        a = prts1[i]
        b = prts2[i]
        s1 = rc(str(fasta1[a[0]].seq[a[1]-1:a[2]]), a[3])
        s2 = rc(str(fasta2[b[0]].seq[b[1]-1:b[2]]), b[3])
        if s1 != s2:
            errs.append((a, b))
            break
        pc += 1
        t = int(pc * 100 / len(ids))
        if t != po:
            po = t
            if t % 10 == 0:
                print("%d%% ..." % t)
    if len(errs) < 1:
        print("Ok!")
    else:
        print("%d ERROS!" % len(errs))
    return errs


updateGFF(agp, old_fasta, new_fasta, old_gff, new_gff)
print(datetime.datetime.now())
print('by mikeias.net')
