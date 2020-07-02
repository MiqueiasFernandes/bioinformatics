#!/usr/bin/env python
# coding: utf-8

# Copyright(c) 2020 - Miquéias Fernandes <bio@mikeias.net>

import os
import sys
import csv
from Bio import SeqIO

#################################################
usage = "python3 maskLow.py file.fasta file.gff3 margin"
#################################################

try: _, file_fasta, file_gff, margin=sys.argv
except: sys.exit("Correct usage is:\n"+usage)
    
out_file = file_fasta + '.MASKED.fasta'
conv = {
    'A': '1', 'C': '2', 'T': '3', 'G': '4',
    'a': '1', 'c': '2', 't': '3', 'g': '4',
    '1': '1', '2': '2', '3': '3', '4': '4', 
    'N': 'N'
}
dec = {
    '1': 'A', 
    '2': 'C', 
    '3': 'T', 
    '4': 'G'
}
tmp_fasta = '__tmp__.fasta'
margin = int(margin)

print('Margin (%d)<-gene->(%d)' % (margin, margin))

print('[1/4] importing fasta %s ...' % file_fasta)
fasta = SeqIO.to_dict(SeqIO.parse(file_fasta, 'fasta'))

print('[2/4] importing GFF3 %s ...' % file_gff)
gff_genes = [ x for x in list(csv.reader(open(file_gff), delimiter='\t')) if len(x) > 2 and x[2] == 'gene']
n_genes = len(gff_genes)

print('[3/4] masking 0% ...')
seqs = set([x[0] for x in gff_genes])
total = len(seqs)

seq2gene = {seq: [] for seq in seqs}
for gene in gff_genes:
    seq2gene[gene[0]].append(gene)
    
cont =  0
for seq in seqs:
    scaffold = [x for x in str(fasta[seq].seq)]
    for gene in seq2gene[seq]:
        for i in range(max(0, int(gene[3])-(1+margin)), min(int(gene[4])+margin, len(scaffold))):
            scaffold[i] = conv[scaffold[i]]
    open(tmp_fasta, 'a').write('>' + seq + '\n' + ''.join([dec[b] if b.isdigit() else 'N' for b in scaffold]) +'\n')
    cont += 1
    if cont % 100 == 0:
        print('[3/4] masking %.1f%% ...' % (cont/total*100))

del fasta
del gff_genes

print('[4/4] persisting ...')
ss = SeqIO.write(SeqIO.parse(tmp_fasta, 'fasta'), out_file, 'fasta')
os.remove(tmp_fasta)
print('\nwrited %d seqs with %d genes masked in %s.' % (ss, n_genes, out_file))
