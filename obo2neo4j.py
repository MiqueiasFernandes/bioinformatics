#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 08/06/2018

from __future__ import with_statement
from collections import defaultdict
from py2neo import Graph, Node, Relationship
from py2neo.packages.httpstream import http
from threading import Thread
import sys
f = open('/var/www/html/saida/res.txt', 'w')
sys.stdout = f


#Thanks Uli Köhler https://techoverflow.net/2013/11/18/a-geneontology-obo-v1-4-parser-in-python/
#and Thanks https://medium.com/labcodes/graph-databases-discutindo-o-relacionamento-dos-seus-dados-com-python-79688b557eec

__author__ = "Miquéias Fernandes"
__copyright__ = "Copyright 2018 mikeias.net"
__license__ = "MIT"
__version__ = "1.1"

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
        self.node = Node(self.namespace,
                         id=self.id,
                         name=dict['name'],
                         definition=dict['def'])

        self.handleStrList('synonym')
        self.handleStrList('subset')
        self.handleStrList('comment')
        self.handleStrList('xref')
        if 'is_obsolete' in dict.keys() and 'true' == dict['is_obsolete']:
            self.node['is_obsolete'] = 'true'

        self.rel = []
        self.handleRelacionamento('relationship', False)
        self.handleRelacionamento('is_a')
        self.handleRelacionamento('alt_id')
        self.handleRelacionamento('consider')
        self.handleRelacionamento('replaced_by')

    def tolist(self, key):
        relc = []
        if key in self.dict.keys():
            relc = self.dict[key]
            if type(self.dict[key]) is str:
                relc = [self.dict[key]]
        return relc

    def handleStrList(self, key):
        if key in self.dict.keys():
            self.node[key] = "; ".join(self.tolist(key))

    def handleRelacionamento(self, key, tipo_eq=True):
        for r in self.tolist(key):
            self.rel.append(Relacionamento(r, self, key if tipo_eq else None))

    def relacionar(self, gos):
        count = 0
        nenc = 0
        for go in gos:
            for rel in self.rel:
                if rel.destino_id == go.id:
                    count += 1
                    rel.destino = go
        for rel in self.rel:
            if rel.destino is None:
                nenc += 1
                # print('WARN: Não foi encontrado TERMO com id ' + rel.destino_id)
        return (len(self.rel), count, nenc)



class Relacionamento:
    def __init__(self, string, origem, tipo):
        self.origem  = origem
        self.origem_id = origem.id
        self.destino = None
        spl = string.split(' ')
        if tipo is None:
            self.destino_id = spl[1]
            self.tipo = spl[0]
        else:
            self.destino_id = spl[0]
            self.tipo = tipo

    def getNoOrigem(self):
        return self.origem.node

    def getNoDestino(self):
        return self.destino.node

    def ok(self):
        return not (self.origem is None or self.destino is None)


class Relacionar(Thread):

    terminaram = 0

    def __init__(self, todos, gos, nome, ini, end):
        Thread.__init__(self)
        self.todos = todos
        self.gos = gos
        self.procs = 0
        self.total = 0
        self.nenc = 0
        self.nome = nome
        self.ini = ini
        self.end = end

    def rep(self):
        if self.nome > 1:
            return (' ' * ((self.nome - 1) * 10))[:self.nome]
        return ''

    def run(self):
        for go in self.gos:
            ret = go.relacionar(self.todos)
            self.total += ret[0]
            self.procs += ret[1]
            self.nenc += ret[2]
        Relacionar.terminaram += 1
        if Relacionar.terminaram % 10 == 0:
            print(str(Relacionar.terminaram) + ' threads terminaram')

THREADS = 100

if __name__ == "__main__":

    gos = []

    print('Importando termos ...')
    termCounter = 0
    for goTerm in parseGOOBO(sys.argv[1]):
        gos.append(GO(goTerm))
        termCounter += 1
    print(str(termCounter) + ' termos importados com sucesso ...')

    print('Relacionando termos ...')
    relCounter = 0
    totalrel = 0
    relnenc = 0

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
        if len(ts) == len(ranges):
            print('as ' + str(THREADS) + ' threads foram startadas ...')

    for t in ts:
        t.join()

    for t in ts:
        relCounter += t.procs
        totalrel += t.total
        relnenc += t.nenc

    if totalrel - relCounter > 0:
        print(str(totalrel) + ' total de relacionamentos ...')
        print(str(relnenc) + ' relacionamentos não encontrados ...')
        print(str(totalrel - relCounter - relnenc) + ' relacionamentos falharam ser criados ...')
    print(str(relCounter) + ' relacionamentos criados com sucesso ...')

    print('Conectando ao banco de dados ...')
    http.socket_timeout = 9999
    g = Graph(password="123")

   # g = Graph("http://neo4j:123@localhost:7474/db/data/", bolt=False)
    tx = g.begin()

    print('Persistindo termos e relacionamentos ...')
    termCounter2 = 0
    relCounter2 = 0
    for go in gos:
        tx.create(go.node)
        termCounter2 += 1
        for rel in go.rel:
            if rel.ok():
                r = Relationship(rel.getNoOrigem(), rel.tipo, rel.getNoDestino())
                try:
                    tx.create(r)
                    relCounter2 += 1
                except:
                    print('impossivel criar relacionamento ' + str(r))

    print('total ' + str(termCounter - termCounter2) + ' termos falharam ao persistir ...')
    print('total ' + str(relCounter - relCounter2) + ' relacionamentos falharam ao persistir ...')

    print('commit ...')
    tx.commit()

    print('by mikeias.net')




