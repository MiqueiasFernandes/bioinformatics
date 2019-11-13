#!/usr/bin/python3

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# gff2html.py - a script for visualization of gff features


#################################################

import sys
from Bio import SeqIO

#################################################
usage = "gff2html.py genome.fasta genome.gff3 anotattion.tsv"
#################################################

try: 
    scr,gen,gff,anot=sys.argv
except: sys.exit("Correct usage is:\n"+usage)


from Bio import SeqIO

class Gene:
    def __init__(self, name, data, anot, exons, seqs):
        self.name = name
        self.seq = data[0]
        self.ini = int(data[1])
        self.end = int(data[2])
        self.strand = data[3] == '+'
        self.interpro = anot[0] if len(anot) > 0 else ''
        self.go = anot[1].split(',') if len(anot) > 1 else []
        self.title = anot[2] if len(anot) > 2 else ' UNKNOWN'
        self.exons = exons
        self.parts = [['e', e[0], e[1]] for e in self.exons]
        self.parts.extend([['i', x[1]+1, min([y[0] for y in self.exons if y[0] > x[1]])-1]  for x in self.exons if len([y for y in self.exons if y[0] > x[1]+1]) > 0])
        if min([min(x[1:]) for x in self.parts]) > self.ini:
            self.parts.append(['u', self.ini, min([min(x[1:]) for x in self.parts])])
        if max([max(x[1:]) for x in self.parts]) < self.end:
            self.parts.append(['u', max([max(x[1:]) for x in self.parts]), sef.end])
        seq = seqs[self.seq].seq
        self.hasN = seq[self.ini-1:self.end].count('N') > 0
        out = str(seq[self.ini -1:min([min(e) for e in self.exons])-1] if self.strand else seq[max([max(e) for e in self.exons])+1:self.end+1].reverse_complement())
        for p in sorted(self.parts, key=lambda x: (x[1] if self.strand else -x[2])):
            e = p[1:]
            out += ('<b class="%s">%s</b>' % (p[0], str(seq[e[0]-1:e[1]] if self.strand else seq[e[0]-1:e[1]].reverse_complement())))
        out += str(seq[max([max(e) for e in self.exons])+1: self.end+1] if self.strand else seq[self.ini-1:min([min(e) for e in self.exons])-1].reverse_complement())
        self.html = out
             
    def asHTML(self):
        out = "<h1><b>%s</b> : %s</h1>" %  (self.name, self.title)
        out += ('<h3>Uniprot: <a href="https://www.uniprot.org/uniprot/%s" target="_blank">%s</a></h3>' % (self.interpro, self.interpro) ) if len(self.interpro) > 1 else ''
        out += ("<h3>GO: " + ', '.join(['<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank">%s</a>' % (g, g) 
                             for g in self.go]) + "</h3><br>") if len(self.go) > 0 else ''
        out += self.html
        return out
        
        
def buildHTML(genes):
    out = '<html><head><style>body {overflow-wrap: break-word; margin: 50px; font-family: monospace} .e {background: lightgreen; -webkit-print-color-adjust: exact}</style></head><body>'
    tot = len(genes)
    writ = [x.asHTML() for x in genes.values() if not x.hasN]
    out += "<hr>".join(writ)
    out += "</body></html>"
    with open('result.html', 'w') as o:
        o.write(out)
    if tot > len(writ):
        print('%d/%d writed' % (len(writ), tot))

        
        
        
        
        
print('[1/5] load annotations of ' + anot) 
anotacao =  {x[0]: x[1:] for x in 
             [l.strip().split('\t') for l in open(anot).readlines() if len(l) > 3]
            if len(x) > 1}

print('[2/5] load gff3 ' + gff)
gff3 = [l.strip().split('\t') for l in open(gff).readlines() if not l.startswith('#') and l.count('\t') == 8 ]
gs = {x[8].split('Name=')[1].split(';')[0]: (x[0],x[3],x[4],x[6]) for x in [g for g in gff3 if g[2] == 'gene']}
des = {}
for e in [g for g in gff3 if g[2] == 'exon']:
    g = '.'.join(e[8].split('ID=')[1].split('.')[:2])
    if g in des:
        des[g].append(e)
    else:
        des[g] = [e]

genes = {g: gs[g] for g in anotacao}
exons =  {g: [(int(x[3]), int(x[4])) for x in des[g]] for g in anotacao}

print('[3/5] load genome of ' + gen)  
seqs = SeqIO.to_dict(SeqIO.parse(gen, "fasta"))

print('[4/5] parsing...')  
gns = {x: Gene(x, genes[x], anotacao[x] if x in anotacao else [], exons[x], seqs) for x in genes}

print('[5/5] persist...') 
buildHTML(gns)

print('stored at result.html')
print('by mikeias.net')
