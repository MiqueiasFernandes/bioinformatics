#!/usr/bin/env python
# coding: utf-8

# Copyright(c) 2020 - Miquéias Fernandes <bio@mikeias.net>

from optparse import OptionParser
import csv
from Bio import SeqIO

usage = "usage: %prog -s specie cluster_file.txt"
parser = OptionParser(usage)
parser.add_option("-f", "--fasta", dest="fasta", help="Filter proteins from gene/mrna list")
parser.add_option("-s", "--spec", dest="spec", help="Specie name in onrthovenn cluster file")
parser.add_option("-g", "--gff", dest="gff", help="gff file to get gene structure", metavar="FILE")
(options, args) = parser.parse_args()

if len(args) < 1:
    parser.error('arg file.txt is obrigatory!')

FILTRO = not options.fasta is None
GENE = not options.gff is None

gff = []
mrna2gene = {}

if GENE:
    gff = [x for x in list(csv.reader(open(options.gff).readlines(), delimiter='\t')) if len(x) > 8 and x[2] == 'mRNA']
    mrna2gene = {m[8].split('ID=')[1].split(';')[0]: m[8].split('Parent=')[1].split(';')[0] for m in gff}

if FILTRO:
    names = [l.strip() for l in open(args[0]).readlines() if len(l.strip()) > 0]
    gs, ms = [], []
    
    for k, v in mrna2gene.items():
        if v in names:
            gs.append(v)
            names.append(k)
            ms.append(k)
        if k in names:
            ms.append(k)
            
    names = set(names)
    proteins = [v for k, v in SeqIO.to_dict(SeqIO.parse(options.fasta, 'fasta')).items() if k in names]
    SeqIO.write(proteins, options.fasta + '-filtered1.faa', 'fasta')
    for p in proteins:
        p.id = 'x' + p.id
    SeqIO.write(proteins, options.fasta + '-filtered2.faa', 'fasta')
    print('%d proteins writed of %d genes to %d seqs in fasta' % (len(set(ms)), len(set(gs)), len(proteins)))
    exit()

if options.spec is None:
    parser.error('[-s] Specie arg is mandatory!')    

data = list(csv.reader(open(args[0]).readlines(), delimiter='\t'))[1:]
parsed = {}

for l in data:
    go = l[3]
    if not go in parsed:
        parsed[go] = []
    mrnas = [x.split('|')[1] for x in l[4].split(';') if x.startswith(options.spec + '|')]
    parsed[go] = list(set(mrnas  + parsed[go]))
    


tba = []
tbg = []
tbp = []
for k, v in parsed.items():
    if not '; ' in k:
        continue
    go = k.split('; ')[0]
    anot = k.split('; ')[1]
    t = anot.split(':')[0]
    d = anot[2:]
    mrnas = set(v)
    genes = []
    if GENE:
        genes = set([mrna2gene[m] for m in mrnas])
        tbg.extend([(g, t, go, d) for g in genes])
    tbp.extend([(p, t, go, d) for p in mrnas])
    tba.append((t, go, d, str(len(genes)), str(len(mrnas))))

open('annotattions.csv', 'w').writelines(['\t'.join(l)+'\n' for l in tba])
open('annotattions_by_proteins.csv', 'w').writelines(['\t'.join(l)+'\n' for l in tbp])
if GENE:
    open('annotattions_by_genes.csv', 'w').writelines(['\t'.join(l)+'\n' for l in tbg])

print('finished.')
