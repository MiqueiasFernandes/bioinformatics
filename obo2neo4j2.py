#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 08/06/2018

from __future__ import with_statement
from collections import defaultdict
from py2neo import Graph, Node, Relationship
import sys
import time
start_time = time.time()
from py2neo.packages.httpstream import http
http.socket.timeout = 99999
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

    def prep(self, gos, tx):
        falt = []
        keys = gos.keys()
        cria = 0
        for r in self.rel:
            if r.destino_id in keys:
                rel = Relationship(r.getNoOrigem(), r.tipo, gos[r.destino_id].node)
                try:
                    tx.create(rel)
                    # print('-', end='')
                    # sys.stdout.flush()
                    cria += 1
                except:
                    print('\nimpossivel criar relacionamento ' + str(rel))
            else:
                falt.append(r)
        return (falt, cria)


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

class Relacionamento:
    def __init__(self, string, origem, tipo):
        self.origem = origem
        self.origem_id = origem.id
        spl = string.split(' ')
        if tipo is None:
            self.destino_id = spl[1]
            self.tipo = spl[0]
        else:
            self.destino_id = spl[0]
            self.tipo = tipo

    def getNoOrigem(self):
        return self.origem.node


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('usage: python obo2neo4j.py file.gff')
        exit(-1)


    print('Conectando ao banco de dados ...')
    #http.socket_timeout = 9999
    g = Graph(password="123")
    tx = g.begin()

    gos = {}
    falt = {}

    print('Importando termos ...')
    termCounter = 0
    relCounter = 0
    for goTerm in parseGOOBO(sys.argv[1]):
        print(' tim ' + str(time.time() - start_time))
        go = GO(goTerm)
        gos[go.id] = go
        tx.create(go.node)
        ret = go.prep(gos, tx)
        for r in ret[0]:
            if r.destino_id in falt.keys():
                falt[r.destino_id].append(r)
            else:
                falt[r.destino_id] = [r]

        if go.id in falt.keys():
            for r in falt[go.id]:
                try:
                    r = Relationship(r.getNoOrigem(), r.tipo, go.node)
                    tx.create(r)
                    # print('*', end='')
                    # sys.stdout.flush()
                    relCounter += 1
                except:
                    print('impossivel criar relacionamento ' + str(r))
            falt[go.id] = None
        relCounter += ret[1]
        termCounter += 1
        if termCounter > 0 and relCounter > 0 and (termCounter % 1000 == 0 or relCounter % 1000 == 0):
            print(' T = ' + str(termCounter) + ' R = ' + str(relCounter))
            sys.stdout.flush()

    print('\n' + str(termCounter) + ' termos importados com sucesso ...')
    print(str(relCounter) + ' relacionamentos importados com sucesso ...')

    print('commit ...')
    tx.commit()
    print("--- %s horas ---" % ((time.time() - start_time)/3600))
    print('by mikeias.net')





