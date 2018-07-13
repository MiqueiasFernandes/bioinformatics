#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 03/07/2018

import os, sys, csv
from sortscaffolds import Gene, Scaffold, Data, getTabelas

__author__ = "Miquéias Fernandes"
__copyright__ = "Copyright 2018 mikeias.net"
__license__ = "MIT"
__version__ = "1.1"


# 0'qseqid ' \'ID=Eucgr.A00001.v2.0;Name=Eucgr.A00001'
# 1'sseqid ' \'pg.scf.311'
# 2'pident ' \'85.774'
# 3'length ' \'2144'
# 4'mismatch ' \'210'
# 5'gapopen ' \'48'
# 6'qstart ' \'1',
# 7'qend ' \'2093'
# 8'sstart ' \'136631'
# 9'send ' \'134532'
# 0'evalue ' \'0.0',
# 1'bitscore''2182'


def getLine(data, genes):
    l = ''
    return l


def process(file, fasta, gff):

    print('importando as informações ...')
    ret = getTabelas(fasta, gff)

    dict = {}
    for g in ret[1]:
        dict[g.nome] = g


    fo = open(file + '.genes.gff', "w")
    print('processando as informações ...')
    with open(file, 'r') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            gene = dict[row[0]]


    for gene in genes:
        fo.write(">" + gene.nome + '\n')
        fo.write(gene.getSeq() + '\n')
    fo.close()



if __name__ == "__main__":
    sys.argv.extend([
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/homology_ecalyptus/blast_genes_ecalyptus_csv',
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/homology_ecalyptus/egrand.fa',
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/homology_ecalyptus/egrand.gff'])

    if len(sys.argv) < 4:
        print("use: ./csv_to_gff.py blast.csv genome.fa original.gff")
    else:
        process(sys.argv[1], sys.argv[2], sys.argv[3])




