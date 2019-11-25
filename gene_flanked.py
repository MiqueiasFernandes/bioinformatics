#!/usr/bin/python

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# gene_flanked.py - a script for get gene sequenvece of gff with flanked region

#################################################

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#################################################
usage = "gene_flanked.py genome.fasta genome.gff3 [flanked region, ex 5000] proteins.or.gene.ids.txt"
#################################################

try:
    script,genome,gff,region,ids=sys.argv
except: sys.exit("Correct usage is:\n"+usage)

region = int(region)

print('load id file %s' % ids)
proteins = [l.strip() for l in open(ids).readlines()]
print('load genome file %s' % genome)
guava_fa = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
print('load gff file %s' % gff)
gff3 = [l.strip().split('\t') for l in open(gff).readlines() if l.count('\t') == 8]
genes = {x[8].split('ID=')[1]: x for x in [x for x in gff3 if x[2] == 'gene']}
outf =  "genes.flanked." + str(region) + ".fa"
print("%d genes flanked at +- %d in %s" % (SeqIO.write([
    SeqRecord(guava_fa[genes[g][0]][
        max(0,int(int(genes[g][3])-region)):
            min(int(genes[g][4])+region,len(guava_fa[genes[g][0]]))].seq, id=g, description='')
    for g in set([x.split('.')[0] if x.count('.') > 0 else x for x in proteins])], outf, "fasta"), region, outf))

print('terminado com sucesso!\nby mikeias.net')


