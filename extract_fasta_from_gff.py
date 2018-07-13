#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 03/07/2018

import os, sys
from sortscaffolds import Gene, Scaffold, Data, getTabelas

__author__ = "Miquéias Fernandes"
__copyright__ = "Copyright 2018 mikeias.net"
__license__ = "MIT"
__version__ = "1.1"


def process(fasta, gff):

    print('importando as informações ...')
    ret = getTabelas(fasta, gff)
    genes = ret[1]

    fo = open(fasta + '.genes.fa', "w")

    for gene in genes:
        fo.write(">" + gene.nome + '\n')
        fo.write(gene.getSeq() + '\n')
    fo.close()


if __name__ == "__main__":
    sys.argv.extend(['/home/mfernandes/Documentos/relatorio_mestrado/julho/homology_ecalyptus/egrand.3.fa',
                     '/home/mfernandes/Documentos/relatorio_mestrado/julho/homology_ecalyptus/egrand.gff'])

    if len(sys.argv) < 3:
        print("use: ./extract_fasta_from_gff.py scaffolds.fasta genes.gff")
    else:
        process(sys.argv[1], sys.argv[2])
