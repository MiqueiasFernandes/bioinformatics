#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 08/06/2018

import sys
import operator
from threading import Thread

__author__ = "Miquéias Fernandes"
__copyright__ = "Copyright 2018 mikeias.net"
__license__ = "MIT"
__version__ = "1.1"


class Gene:
    def __init__(self, data, scaffold):
        self.nome = data[8][:-1]
        self.scaffold = scaffold
        self.start = int(data[3])
        self.stop = int(data[4])
        self.len = self.stop - self.start
        scaffold.addGene(self)

    def sobrepoe(self, gene):
        tipo = 0
        if self.stop >= gene.start >= self.start:
            tipo += 2
        if self.start <= gene.stop <= self.stop:
            tipo += 5
        print('comp ' + str(self) + ' com ' + str(gene) + ' = ' + str(tipo))
        return tipo

    def __str__(self):
        return 'Gene {' + self.nome + '}' + ' [ ' + str(self.start) + ' - ' + str(self.stop) + ' ] (' + str(self.len) + ')'


class Scaffold:
    def __init__(self, nome):
        self.nome = nome
        self.genes = []
        self.sobp = []

    def addGene(self, gene):
        self.genes.append(gene)

    def processar(self):
        for i in range(len(self.genes)):
            for j in range(1, len(self.genes)):
                tipo = self.genes[i].sobrepoe(self.genes[j])
                if tipo > 0:
                    self.sobp.append(Sobrepoem(self, self.genes[i], self.genes[j], tipo))
        self.sobp = sorted(self.sobp, key=operator.attrgetter('score'))


class Sobrepoem:
    def __init__(self, scaffold, gene1, gene2, tipo):
        self.scaffold = scaffold
        self.gene1 = gene1
        self.gene2 = gene2
        self.tipo = tipo
        self.score = min(gene1.start, gene2.start)

    def __str__(self):
        if self.tipo == 2:
            print("o gene " + self.gene2 + " inicia dentro do gene " + self.gene1)
        elif self.tipo == 5:
            print("o gene " + self.gene2 + " TERMINA dentro do gene " + self.gene1)
        elif self.tipo == 7:
            print("o gene " + self.gene2 + " ESTÁ CONTIDO NO gene " + self.gene1)
        else:
            print("erro de processamento entre gene " + self.gene1 + " e " + self.gene2)


class Computar(Thread):

    terminaram = 0

    def __init__(self, scaffolds):
        Thread.__init__(self)
        self.scaffolds = scaffolds

    def run(self):
        for scaffold in self.scaffolds:
            scaffold.processar()
        Computar.terminaram += 1
        if Computar.terminaram % 10 == 0:
            print(str(Computar.terminaram) + ' threads terminaram')

# # start gene g1
# pg.scf.589      AUGUSTUS        gene    95      610     0.16    +       .       g1
# pg.scf.589      AUGUSTUS        transcript      95      610     0.16    +       .       g1.t1
# pg.scf.589      AUGUSTUS        start_codon     95      97      .       +       0       transcript_id "g1.t1"; gene_id "g1";
# pg.scf.589      AUGUSTUS        CDS     95      610     0.16    +       0       transcript_id "g1.t1"; gene_id "g1";
# pg.scf.589      AUGUSTUS        stop_codon      608     610     .       +       0       transcript_id "g1.t1"; gene_id "g1";
# # protein sequence = [MSWSPELHRMFVNAVQQLGIDSARPQHILEIMKVEGLTKGNVSSHLQKYRKRVKEHNAAQNQENQGGTDPSMNVESPR
# # HHKFHRLAWGEPLPIINQASSRAINPSGVNQSIPQNHPNALQEYRPQGQYIQTPWFEHQPVARIHPMISASAQPENQTRNYPEHFHFSNPRME]
# # end gene g1

THREADS = 20


def getscaffold(scfs, nome):
    for s in scfs:
        if s.nome == nome:
            return s
    s = Scaffold(nome)
    scfs.append(s)
    return s


def process(file):
    genes = []
    scaffolds = []
    print('importando informações ... ')
    emGene = False
    for line in file:
        if emGene:
            if line.find("\tAUGUSTUS\tgene\t") >= 0:
                data = line.split("\t")
                if len(data) > 5:
                    scf = getscaffold(scaffolds, data[0])
                    genes.append(Gene(data, scf))
                    emGene = False
        elif line.startswith("# start gene "):
            emGene = True
    file.close()
    print(str(len(scaffolds)) + ' scaffolds importados ... ')
    print(str(len(genes)) + ' genes importados ... ')
    print('verificando sobreposições ... ')

    ranges = []
    mul = round(len(scaffolds)/THREADS)
    for k in range(THREADS):
        if k == THREADS-1:
            ranges.append((k * mul, len(scaffolds) - 1))
        else:
            ranges.append((k*mul, (((k+1)*mul)-1)))

    ts = []
    for r in ranges:
        t = Computar(scaffolds[r[0]:r[1]])
        ts.append(t)
        t.start()
        if len(ts) == len(ranges):
            print('as ' + str(THREADS) + ' threads foram startadas ...')

    for t in ts:
        t.join()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("use: ./gene_over.py genes.gff")
    else:
        process(open(sys.argv[1], "r"))
