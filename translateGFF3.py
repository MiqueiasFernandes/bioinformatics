#!/usr/bin/env python
# coding: utf-8

# Copyright(c) 2020 - Miquéias Fernandes <bio@mikeias.net>

import sys
import csv
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#################################################
usage = "python3 translateGFF3.py genome.fasta file.gff3"
#################################################

try: _, file_fasta,file_gff3=sys.argv
except: sys.exit("Correct usage is:\n"+usage)

out_file = file_gff3 + '-proteins.faa'

key_ID = 'ID'
key_Parent = 'Parent'


def translate(fast, mrna, gene, cds):
    return SeqRecord(Seq(''.join([str(
    fast[x[0]].seq[x[1]-1:x[2]] if x[3] else fast[x[0]].seq[x[1]-1:x[2]].reverse_complement()) for x in 
         sorted([(c[0], int(c[3]), int(c[4]), c[6] != '-') for c in cds], 
                key=lambda e: e[1] if e[3] else -e[2])])).translate(), mrna, description='gene=' + gene)


def process(fasta, mrna2cds, mrna2gene):
    cont = 0
    qtd = len(mrna2cds)
    res = []
    tp = 0
    for m, c in mrna2cds.items():
        res.append(translate(fasta, m, mrna2gene[m], c))
        cont += 1
        pc = int(cont / qtd * 100)
        if pc != tp:
            tp = pc
            print('[4/4] %d%%' %pc )
    return res


print('[1/4] opeining FASTA at ' + file_fasta)
fasta = SeqIO.to_dict(SeqIO.parse(file_fasta, 'fasta'))

print('[2/4] opeining GFF3 at ' + file_gff3)
gff3 = [x for x in list(csv.reader(open(file_gff3), delimiter='\t')) if len(x) == 9 and not x[0].startswith('#')]

print('[3/4] preprocessing %s genes ' % len([x for x in gff3 if x[2] == 'gene']))
mrnas = {x[8].split('%s=' % key_ID)[1].split(';')[0]: x for x in gff3 if x[2] == 'mRNA'}
mrna2gene = {m: x[8].split('%s=' % key_Parent)[1].split(';')[0] for m, x in mrnas.items()}
all_cds = [(x[8].split('%s=' % key_Parent)[1].split(';')[0].split(','), x) for x in gff3 if x[2] == 'CDS']
mrna2cds = {}
for x in all_cds:
    ms, c = x
    for m in ms: 
        if m in mrna2cds:
            mrna2cds[m].append(c)
        else:
            mrna2cds[m] = [c]

print('[4/4] translating %s mrnas... ' % len(mrna2cds))
print(SeqIO.write(process(fasta, mrna2cds, mrna2gene), out_file, 'fasta'))

print('all finished. sotored at: ' + out_file)

