#!/usr/bin/python

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# parse_allmaps.py - a script for simplify use of https://github.com/tanghaibao/jcvi/wiki/ALLMAPS with blast identity


#################################################

import os
import sys
sys.path.append('/home/cluster/shared/pybio')
from Blast import Blastn
from SeqTools import *
from Bio import SeqIO

#################################################
usage = "parse_allmaps.py genome.fasta primers_seq.txt filter.txt map1.tsv map2.tsv map3.tsv ..."
#################################################

try: 
    script,genome,primers,filtro=sys.argv[:4]
    maps = sys.argv[4:]
except: sys.exit("Correct usage is:\n"+usage)
    
sys.stdout = sys.stderr
    
class Mapa:
    def __init__(self, assembly, primers_seq, mapa_files=[], filter_file=None, run=True):
        self.mapa_files = mapa_files
        self.assembly = assembly
        self.filter_file = filter_file
        self.genome = SeqIO.to_dict(SeqIO.parse(assembly, 'fasta'))
        self.glen = sum([len(x) for x in self.genome.values()])
        self.fasta_primers = assembly + '.primers.query.fa'
        self.raw_m = {m: [l.strip().split('\t') for l in open(m).readlines() if len(l) > 2] for m in self.mapa_files}
        self.gls = {}
        self.marks = set()
        self.not_f = []
        self.marker2scf = {}
        self.map2gl = {}
        self.map2mark = {}
        for m, raw in self.raw_m.items():
            gl = None
            bf = []
            self.map2gl[m] = set()
            self.map2mark[m] = set()
            for l in raw:
                if l[0].startswith('#'):
                    if len(bf) > 0:
                        self.gls[gl] = [(float(x[0]), x[1]) for x in bf]
                    gl = l[0][1:]
                    bf = []
                    self.map2gl[m].add(gl)
                elif not gl is None and len(l) == 2:
                    bf.append(l)
                    self.marks.add(l[1])
                    self.map2mark[m].add(l[1])
            if len(bf) > 0:
                self.gls[gl] = [(float(x[0]), x[1]) for x in bf]
                self.map2gl[m].add(gl)
        for m in self.mapa_files:
            print('MAP %s: GL %d MK %d' % (m, len(self.map2gl[m]), len(self.map2mark[m])))
        print('GL importado: %d' % len(self.gls))
        print('markers importado: %d' % len(self.marks))
        self.primers2query(mapas=self.mapa_files, primers=primers_seq, query=self.fasta_primers)
        if run:
            self.run()
        
    def primers2query(self, primers, mapas, query):
        ls = [l.strip().split('\t') for l in open(primers).readlines() if not l.startswith('#')]
        mps = {x[0]: x[1:] for x in ls if x[0].startswith('mPg')}
        me = {x[0][:-1]: x[1] for x in ls if x[0].startswith('Me')}
        em = {x[0][:-1]: x[1] for x in ls if x[0].startswith('Em')}
        mks = set()
        for m in mapas:
            mks = mks.union(
                set([
                    l.strip().split('\t')[1] for l in open(m).readlines() if not l.startswith('#') and l.count('\t') == 1
                ]))

        fast = []

        not_f = set()
        ok = 0
        for m in mks:
            if m.startswith('mPg'):
                if m in mps:
                    ok += 1
                    fast.append(">%s\n%sNNNNN%s\n" % (m, mps[m][0], mps[m][1]))
                    fast.append(">R_%s\n%sNNNNN%s\n" % (m, mps[m][0], rc(mps[m][1])))
                else:
                    not_f.add(m)
            elif m.startswith('Me'):
                s = m.split('E')
                if s[0] in me and 'E' + s[1] in em:
                    ok += 1
                    fast.append(">%s\n%sNNNNN%s\n" % (m, me[s[0]],em['E' + s[1]]))
                    fast.append(">R_%s\n%sNNNNN%s\n" % (m, me[s[0]], rc(em['E' + s[1]])))
                else:
                    not_f.add(m)

        if len(not_f) > 0:
            print('%d primers dos mapas nao foram encontrados: %s' % (len(not_f), str(not_f)))
        self.not_f = list(not_f)
        print("%d primers enviados para %s" % (ok, query))
        with open(query, 'w') as o:
            o.writelines(fast)
   
    def alns(self, evalue, cov, blast=None, show=True):
        self.blast = Blastn(db=self.assembly, query=self.fasta_primers) if blast is None else blast
        b = self.blast.run(evalue=evalue, cov=cov, short=True, seq=False) if blast is None else None
        self.anch = sum([x[0].slen for x in self.blast.subject2query.values()])
        if show:
            print('total size aln: %d Mpb (%.2f%%) [%d scfs]' % (
                int(self.anch / 10**6), self.anch*100/self.glen, len(self.blast.subject2query)))
        tmp = {x[0]: [(y.sseqid, y.sstart) for y in x[1]] for x in self.blast.query2subject.items()}
        self.marker2scf = {x: set([s for s in tmp[x] + (tmp['R_' + x] if ('R_' + x) in tmp else []) ]) 
                           for x in set([y[2:] if y.startswith('R_') else y for y in tmp])}
        return self.marker2scf
             
    def filter_alns(self, rules={'evalue': ['<', 10], 'qcovhsp': ['>', 40]}, persist=False, show=True):
        b = self.blast
        old = b.alignments
        old_a = self.anch
        b.parseBlast(parsed=b.fiter(rules=rules))
        marker2scf = self.alns(evalue=0, cov=0, blast=b, show=show)
        if show:
            print('reduzira %.2f%% do anterior' % (100-(self.anch*100/old_a)))
        if not persist:
            b.parseBlast(parsed=old)
            self.alns(evalue=0, cov=0, blast=b, show=show)
        return marker2scf
            
    def toLinks(self, marker2scf={}, file=None, pref='new', forMAP=None, show=True):
        if len(marker2scf) < 1 and len(self.marker2scf) > 0:
            marker2scf = self.marker2scf
        if not forMAP is None and show:
            print('tratando mapa %s ...' % forMAP)
        marker2scf = marker2scf if forMAP is None else {x: marker2scf[x] for x in self.map2mark[forMAP] if x in marker2scf}
        mrks = self.marks if forMAP is None else self.map2mark[forMAP]
        marks_not_in_map = [x for x in marker2scf if not x in mrks]
        marks_not_in_blast = [x for x in mrks if not x in marker2scf and not x in self.not_f]
        if len(marks_not_in_map) > 0:
            print("%d marks of blast not in map: %s" % (len(marks_not_in_map), str(marks_not_in_map)))
        if len(marks_not_in_blast) > 0:
            print("%d marks of map not blasted: %s" % (len(marks_not_in_blast), str(marks_not_in_blast)))
        outfile = pref + ('.'  if forMAP is None else '.' + forMAP) + '.out.map.csv' if file is None else file
        with open(outfile, 'w') as o:
            for gl in (self.gls if forMAP is None else self.map2gl[forMAP]):
                for mark in self.gls[gl]:
                    if mark[1] in marker2scf:
                        us = []
                        for scf in marker2scf[mark[1]]:
                            if not scf[0] in us:
                                us.append(scf[0])
                                o.write('%s,%d,%s,%f\n' % (scf[0], scf[1], gl, mark[0]))
        if show:
            print('mapa parseado salvo em %s ' % outfile)
        return outfile
    
    
    def run(self, evalue="100", cov=30, file=None):
        self.alns(evalue, cov)
        if not self.filter_file is None:
            self.filter_alns(rules={x[0]: [x[1], int(x[2]) if x[2].isdigit() else float(x[2])] for x in 
                              [l.strip().split('\t') for l in open(self.filter_file).readlines() if l.count('\t') == 2]
                             }, 
                             persist=True)
        nf = " ".join([self.toLinks(forMAP=self.mapa_files[m], file=('m%d.csv'%(m+1))) for m in range(len(self.mapa_files))])
        self.outfile = 'merged.bed'
        print('bed salvo em %s ' % self.outfile)
        os.system("python -m jcvi.assembly.allmaps merge " + nf +  " -o " +  self.outfile)
        self.filters = self.storeFilters()      
        
    def storeFilters(self, threshold=70):
        q = threshold
        atr = []
        tt = len(self.blast.alignments)
        for e in [l for l in range(int(self.blast.stats['bad evalue']), int(self.blast.stats['best evalue']), -1)]:
            for c in [l for l in range(int(self.blast.stats['bad coverage']), int(self.blast.stats['best coverage']), 1)]:
                f = self.blast.fiter(rules={'evalue': ['<', e], 'qcovhsp': ['>', c]})
                g = int(sum([z[1] for z in set([(y.sseqid, y.slen) for y in f])]) *100/ self.glen)
                if g > q:
                    atr.append([e, c, g, len(f)])
                else:
                    break

        d = {}
        for x in sorted(atr, key=lambda e: (-e[2], e[0], e[3], -e[1])):
            if not int(x[2]) in d:
                d[int(x[2])] = x
                
        fs = []
        for k, v in sorted(d.items(), key=lambda e:(-e[0])):
            print('[rule for %d%%] => evalue < %d & qcovhsp > %d => alns %d/%d' % (k, v[0], v[1], v[3], tt))
            nf = " ".join([self.toLinks(show=False,forMAP=self.mapa_files[m], marker2scf=self.filter_alns(show=False,rules={'evalue': ['<', v[0]], 'qcovhsp': ['>', v[1]]}), file=('f' + str(k) + 'm' + str(m+1))) for m in range(len(self.mapa_files))])
            o = "filtered" + str(k) + ".bed"
            os.system("python -m jcvi.assembly.allmaps merge " + nf +  " -o " + o)
            fs.append(o)
            print('bed salvo em %s ' % o)
            f = 'filtro.%d.txt' % k
            with open(f, 'w') as o:
                o.write('evalue\t<\t%d\nqcovhsp\t>\t%d\n' % (v[0], v[1])) 
        return fs
        
def handleMaps(genome, maps, primers, filtro=None):
    print('parsing maps ...')
    mapa = Mapa(assembly=genome, mapa_files=maps, primers_seq=primers, filter_file=filtro)
    print('run allmaps ...')
    cmd = "python -m jcvi.assembly.allmaps path %s " + genome
    for f in [mapa.outfile] + mapa.filters:
        print("rodando %s" % f)
        os.system("cut -f4 %s | cut -d- -f1 | uniq | sed 's/$/ 1/' > weights.txt" % f)
        os.system(cmd % f)
        print(cmd % f)
        os.system("for i in chr*.pdf ; do mv $i %s.$i ; done" % f.split('.')[0])
    print("paste -d@ <(seq -w 1 $(grep -c \> filtered78.fasta) | sed 's/^/>scaffold/' ) <(sed 's/>.*/>/' filtered78.fasta | tr -d \\n | tr \> \\n | tr -s N N | awk '{ print length, $0 }' | sort -nrs | cut -d' ' -f2- ) | tr @ \\n | fold > final.scaffolds.fasta")
    print("for f in `ls */*.sizes` ; do ;[ "$(head -11 $f | grep -c chr)"  == "11" ] && echo $f && paste <(paste -d. <(head $f -n11 | cut -f2 | rev | cut -c7- | rev )  <(head $f -n11 | cut -f2 | rev | cut -c4-6 | rev )) | sort -nrk1 ; done  | paste - - - - - - - - - - - - | sort -nk2")

handleMaps(genome=genome, maps=maps, primers=primers, filtro=filtro)


print('by mikeias.net')
