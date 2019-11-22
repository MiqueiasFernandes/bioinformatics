#!/usr/bin/python

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# seqping2gff.py - a script for fix gff from by pipeline seqping


#################################################

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

#################################################
usage = "eqping2gff.py genoma.fasta original.gff3 filtrado.gff3 PREFIX"
#################################################

try:
    script,genoma,gff3,filtro,prefix=sys.argv
    maps = sys.argv[4:]
except: sys.exit("Correct usage is:\n"+usage)

gff1 = gff3
gffp = gff3
gff2 = gffp + '.introns.gff3'
red=False

print('[1/14] open filtro ' + filtro)
filtrar_mrnas = list(set([m[8].split('ID=')[1].split(';')[0] for m in [l.strip().split('\t') for l in open(filtro).readlines() if l.count('\t') == 8] if m[2] == 'mRNA']))

print('[2/14] open gff ' + gff3)
gff = [l.strip().split('\t') for l in open(gff3).readlines() if l.count('\t') == 8]

print('[3/14] filtrar %d mrnas do gff' % (len([x for x in gff if x[2] == 'mRNA']) - len(filtrar_mrnas)))
n_tem_introns = len([x for x in gff if x[2] == 'intron']) < 1
if n_tem_introns:
    mrna2exon = {x[8].split('ID=')[1].split(';')[0]:[[], x] for x in [x for x in gff if x[2] == 'mRNA']}
    for e in [x for x in gff if x[2] == 'exon']:
        for p in [x for x in e[8].split('Parent=')[1].split(';')[0].split(',') if x in filtrar_mrnas]:
            mrna2exon[p][0].append((int(e[3]), int(e[4])))

print('[4/14] criando introns')
introns = []
if n_tem_introns:
    for m, d in mrna2exon.items():
        exons = d[0]
        l = len(exons)
        if l > 1:
            es = sorted(exons, key=lambda e: e[0])
            cont = 1
            for i in range(0, l-1):
                introns.append((m, es[i][1] + 1, es[i+1][0] - 1, d[1], cont))
                cont += 1

print('[5/14] adicionando %d introns' % len(introns))

gff.extend([[
    i[3][0], 
    i[3][1], 
    'intron', 
    str(min(i[1:3])), 
    str(max(i[1:3])), 
    '.', 
    i[3][6], 
    '.', 
    'ID=%s:intron%d;Parent=%s' % (i[0], i[-1], i[0])] for i in introns])


class Entry:
    def __init__(self, raw):
        self.raw = raw
        self.start = int(raw[3])
        self.end = int(raw[4])
        self.type = raw[2]
        self.key = '\t'.join(raw[3:8])
        self.oldID = raw[8].split('ID=')[1].split(';')[0]
        self.parents = raw[8].split('Parent=')[1].split(';')[0].split(',') if 'Parent=' in raw[8] else []
        self.newParents = []
        
    def parse(self, redundant=False, parent=''):
        return self.raw[:-1] + ['ID=' + self.newID + ((';Parent=' + ','.join([parent] if redundant else self.newParents))
                                                      if len(self.newParents) > 0 else '')]
              
class Gene(Entry):
    def __init__(self, raw):
        Entry.__init__(self, raw)
        self.mrnas = []
        
    def nomear(self, prefix, name, alph=range(1, 1000), sep='.'):
        self.newID = prefix + str(name)
        s = sorted(self.mrnas, key=lambda e: e.start)
        entryes_uniq = {}
        self.nmrnas = {}
        for i in range(len(s)):
            s[i].nomear(self.newID, alph[i], sep=sep)
            self.nmrnas[s[i].newID] = s[i]
            for t, es in s[i].entryes.items():
                if not t in entryes_uniq:
                    entryes_uniq[t] = {}
                for e in sorted(es, key=lambda e: e.start):
                    if not e.key in entryes_uniq[t]:
                        entryes_uniq[t][e.key] = [self.newID + sep + t + sep + str(len(entryes_uniq[t])+1),s[i].newID,[]]
                    entryes_uniq[t][e.key][2].append(s[i].newID)
        for t, es in entryes_uniq.items():
            for r, e in es.items():
                n = Entry(self.raw[:2] + [t] + r.split('\t') + ['ID=' + e[0] + ';Parent=' + ','.join(e[2])])
                n.newID = n.oldID
                n.newParents = n.parents
                self.nmrnas[e[1]].childs.append(n)
                for p in e[2]:
                    self.nmrnas[p].childs_redundantes.append(n)
        
    def addMRNA(self, m):
        m.gene = self
        self.mrnas.append(m)
           
    def parse(self, redundant=False):
        t = [Entry.parse(self)]
        for x in [m.parse(redundant) for m in sorted(self.mrnas, key=lambda e: e.start)]:
            t.extend(x)
        return t

class mRNA(Entry):
    def __init__(self, raw):
        Entry.__init__(self, raw)
        self.entryes = {}
        self.childs = []
        self.childs_redundantes = []
        
    def addEntry(self, entry):
        if entry.type in self.entryes:
            self.entryes[entry.type].append(entry.key)
        else:
            self.entryes[entry.type] = [entry.key]
        
    def nomear(self, geneID, name, sep='.'):
        self.newID = geneID + sep + str(name)
        self.newParents = [geneID]
                
    def parse(self, redundant=False):
        return [Entry.parse(self)] + [x.parse(redundant, self.newID) for x in 
                                      sorted((self.childs if not redundant else self.childs_redundantes), 
                                             key=lambda e: (e.start, (
                                                 1 if 'UTR' in e.type else (
                                                     2 if e.type == 'CDS' else 3))))]


print('[6/14] corrgindo nomes')

genes = {}
mrnas = {}

for x in sorted(gff, key=lambda y: (int(y[3]), 1 if y[2] == 'gene' else (2 if y[2] == 'mRNA' else 3))):
    if x[2] == 'gene':
        g = Gene(x)
        genes[g.oldID] = g
    elif x[2] == 'mRNA':
        m = mRNA(x)
        if m.oldID in filtrar_mrnas:
            genes[m.parents[0]].addMRNA(m)
            mrnas[m.oldID] = m
    else:
        e = Entry(x)
        for p in e.parents:
            if p in mrnas:
                m = mrnas[p]
                if not e.type in  m.entryes:
                    m.entryes[e.type] = []
                m.entryes[e.type].append(e)

print('[7/14] ordenando')

gs = sorted([g for g in genes.values() if len(g.mrnas) > 0], key=lambda e: (e.raw[0], e.start))
for i in range(len(gs)):
    gs[i].nomear(prefix, str(i+1).rjust(5, '0'))

print('[8/14] persistindo em ' + gff2)

with open(gff2, 'w') as o:
    o.write('##gff-version 3\n' + '\n'.join(['\n'.join([ '\t'.join(e) for e in g.parse(redundant=red)]) for g in gs]) + '\n')

print('[9/14] importando genoma em ' + genoma)

seqs_dict = SeqIO.to_dict(SeqIO.parse(genoma, "fasta"))

print('[10/14] persistindo ' + gffp + '.genes.fna')
SeqIO.write([SeqRecord((seqs_dict[g.raw[0]][g.start-1:g.end] if g.raw[6] == '+' else seqs_dict[g.raw[0]][g.start-1:g.end].reverse_complement()
                        ).seq,
                       id=g.newID,
                       description="contig=" + g.raw[0]) for g in gs], gffp + ".genes.fna", "fasta")

print('[11/14] persistindo ' + gffp + '.exons.fna')
exons = []
for m in mrnas.values():
    for e in [z for z in m.childs_redundantes if z.type == 'exon']:
        exons.append(SeqRecord((seqs_dict[e.raw[0]][e.start-1:e.end] if e.raw[6] == '+' else seqs_dict[e.raw[0]][e.start-1:e.end].reverse_complement() ).seq,
                               id=e.newID,
                               description="mRNA=" + m.newID))
SeqIO.write(sorted(exons, key=lambda e: e.id), gffp + ".exons.fna", "fasta")

print('[12/14] persistindo ' + gffp + '.cds.fna')

def getTrincas(codons):
    return codons[:3*int(len(codons)/3)]


cds = sorted([[m, Seq(getTrincas(''.join([
    str(x.seq if m.raw[6] == '+' else x.reverse_complement().seq) for x in [
        seqs_dict[m.raw[0]][c.start-1:c.end] for c in sorted([
            z for z in m.childs_redundantes if z.type == 'CDS'], key=lambda e: (e.start if e.raw[6] == '+' else -e.end))]
])), generic_dna)] for m in mrnas.values()], key=lambda x: x[0].newID)

SeqIO.write([SeqRecord(c[1], id=c[0].newID, description="gene="+c[0].gene.newID) for c in cds], gffp+ ".cds.fna", "fasta")

print('[13/14] persistindo ' + gffp + '.proteins.faa')
SeqIO.write([SeqRecord(c[1].translate(), id=c[0].newID, description="gene="+c[0].gene.newID) for c in cds], gffp+".proteins.faa", "fasta")

print('[14/14] validando')

def getFeatures(file):
    gff = [l.strip().split('\t') for l in open(file).readlines() if l.count('\t') == 8]
    comp = {x: set() for x in set([x[2] for x in gff])}
    ids = {}
    for e in sorted(gff, key=lambda x: (1 if x[2] == 'gene' else (2 if x[2] == 'mRNA' else 3))):
        ids[e[8].split('ID=')[1].split(';')[0]] = e[0]+'-'+e[2]+'-'+e[3]+'-'+e[4]+'-'+e[6]
        if e[8].count('Parent=') > 0:
            for p in e[8].split('Parent=')[1].split(';')[0].split(','):
                if len(p) > 0:
                    comp[e[2]].add(e[0]+'-'+ids[p]+'-'+e[3]+'-'+e[4]+'-'+e[6]+'-'+e[7])
        else:
            comp[e[2]].add(e[0]+'-'+e[3]+'-'+e[4]+'-'+e[6]+'-'+e[7])
    return comp


a = getFeatures(gff1)
b = getFeatures(gff2)

for e, es in b.items():
    if e in a:
        d = len(a[e] - es)
        n = len(es - a[e])
        if d > 0:
            print('new gff missed %d %s' % (d, e))
        if n > 0:
            print('new gff has %d new %s' % (n, e))
        if d + n < 1:
            print('%s [%d] ok!' % (e, len(es)))
    else:
        print('new gff has %d new %s' % (len(es), e))



print('terminado com sucesso!\nby mikeias.net')

