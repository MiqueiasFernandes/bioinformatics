#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A constant-space parser for the GeneOntology OBO v1.2 & v1.4 format

https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html

Version 1.1: Python3 ready & --verbose CLI option
"""
from __future__ import with_statement
from collections import defaultdict
from py2neo import Graph, Node, Relationship
from threading import Thread
import sys

__author__    = "Uli Köhler"
__copyright__ = "Copyright 2013 Uli Köhler"
__license__   = "Apache v2.0"
__version__   = "1.1"

def processGOTerm(goTerm):
    """
    In an object representing a GO term, replace single-element lists with
    their only member.
    Returns the modified object as a dictionary.
    """
    ret = dict(goTerm) #Input is a defaultdict, might express unexpected behaviour
    for key, value in ret.items():
        if len(value) == 1:
            ret[key] = value[0]
    return ret

def parseGOOBO(filename):
    """
    Parses a Gene Ontology dump in OBO v1.2 format.
    Yields each
    Keyword arguments:
        filename: The filename to read
    """
    with open(filename, "r") as infile:
        currentGOTerm = None
        for line in infile:
            line = line.strip()
            if not line: continue #Skip empty
            if line == "[Term]":
                if currentGOTerm: yield processGOTerm(currentGOTerm)
                currentGOTerm = defaultdict(list)
            elif line == "[Typedef]":
                #Skip [Typedef sections]
                currentGOTerm = None
            else: #Not [Term]
                #Only process if we're inside a [Term] environment
                if currentGOTerm is None: continue
                key, sep, val = line.partition(":")
                currentGOTerm[key].append(val.strip())
        #Add last term
        if currentGOTerm is not None:
            yield processGOTerm(currentGOTerm)



class GO:
    def __init__(self, dict):
        self.dict = dict
        self.namespace = dict['namespace']
        self.id = dict['id']
        self.name = dict['name']
        self.definition = dict['def']
        self.rel = []
        self.node = Node(self.namespace, id=self.id, name=self.name, definition=self.definition)
        if 'relationship' in goTerm.keys():
            relc = goTerm['relationship']
            if type(goTerm['relationship']) is str:
                relc = [goTerm['relationship']]
            for r in relc:
                self.rel.append(Relacionamento(r, self))

    def relacionar(self, gos):
        count = 0
        for go in gos:
            for rel in self.rel:
                if rel.destino_id == go.id:
                    count += 1
                    rel.destino = go
        for rel in self.rel:
            if rel.destino is None:
                print('WARN: Não foi encontrado TERMO com id ' + rel.destino_id)
        return (len(self.rel), count)



class Relacionamento:
    def __init__(self, string, origem):
        self.origem  = origem
        self.origem_id = origem.id
        self.destino = None
        spl = string.split(' ')
        self.destino_id = spl[1]
        self.tipo = spl[0]

    def getNoOrigem(self):
        return self.origem.node

    def getNoDestino(self):
        return self.destino.node

    def ok(self):
        return not (self.origem is None or self.destino is None)


class Relacionar(Thread):

    def __init__(self, todos, gos, nome, ini, end):
        Thread.__init__(self)
        self.todos = todos
        self.gos = gos
        self.procs = 0
        self.total = 0
        self.nome = nome
        self.ini = ini
        self.end = end

    def rep(self):
        if self.nome > 1:
            return (' ' * ((self.nome - 1) * 10))[:self.nome]
        return ''

    def run(self):
        print('thread ' + str(self.nome) + ' rodando de ' + str(self.ini) + ' a ' + str(self.end) + ' com ' + str(len(self.gos)))
        for go in self.gos:
            ret = go.relacionar(self.todos)
            self.total += ret[0]
            self.procs += ret[1]
        print('thread ' + str(self.nome) + ' terminou com ' + str(self.procs) + ' terms processados')

THREADS = 100

if __name__ == "__main__":

    print('Conectando ao banco de dados ...')
    g = Graph(password="123")
    tx = g.begin()
    gos = []

    print('Importando termos ...')
    termCounter = 0
    for goTerm in parseGOOBO('/home/mfernandes/Documentos/neo4j/go-basic.obo'):
        gos.append(GO(goTerm))
        termCounter += 1
    print(str(termCounter) + ' termos importados com sucesso ...')

    print('Relacionando termos ...')
    relCounter = 0
    totalrel = 0

    ranges = []
    mul = round(termCounter/THREADS)
    for k in range(THREADS):
        if k == THREADS-1:
            ranges.append((k * mul, termCounter - 1))
        else:
            ranges.append((k*mul, (((k+1)*mul)-1)))

    ts = []
    for r in ranges:
        t = Relacionar(gos, gos[r[0]:r[1]], len(ts)+1, r[0], r[1])
        ts.append(t)
        t.start()

    for t in ts:
        t.join()

    for t in ts:
        relCounter += t.procs
        totalrel += t.total

    if totalrel - relCounter > 0:
        print(str(totalrel) + ' total de relacionamentos ...')
        print(str(totalrel - relCounter) + ' relacionamentos falharam ser criados ...')
    print(str(relCounter) + ' relacionamentos criados com sucesso ...')

    print('Persistindo termos e relacionamentos ...')
    termCounter2 = 0
    relCounter2 = 0
    for go in gos:
        tx.create(go.node)
        termCounter2 += 1
        for rel in go.rel:
            if rel.ok():
                r = Relationship(rel.getNoOrigem(), rel.tipo, rel.getNoDestino())
                tx.create(r)
                relCounter2 += 1

    print('total ' + str(termCounter - termCounter2) + ' termos falharam ao persistir ...')
    print('total ' + str(relCounter - relCounter2) + ' relacionamentos falharam ao persistir ...')

    print('commit ...')
    tx.commit()

    print('by mikeias.net')