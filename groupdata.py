#!/usr/bin/env python
# coding: utf-8


## LICENCE MIT
## REV 07/19
## Join data of annotations into one table for relational databases
## Usage: groupdata.py proteins.fasta genome.gff3 blast.tsv eggnog_gos.tsv eggnog_ort.tsv interpro.tsv OUT_FILE.tsv
## www.mikeias.net
## bio@mikeias.net


import sys
from datetime import datetime
from Bio import SeqIO


class Anotattion:
    def __init__(self, raw, scaffold, seq):
        blast = ['' for i in range(47)]
        for r in range(len(raw)):
            blast[r] = raw[r]
        d = blast[2].split("|")
        self.mrna = d[0]
        self.gene = d[1]
        self.reviwed = blast[17] == 'reviewed'
        self.gis = [x for x in blast[3].replace(' ', '').split(';') if len(x) > 2]
        dic = {
            'InterPro': set([x for x in blast[23].replace(' ', '').split(';') if len(x) > 2]),
            'Pfam': set([x for x in blast[24].replace(' ', '').split(';') if len(x) > 2]),
            'KEGG': set([x for x in blast[25].replace(' ', '').split(';') if len(x) > 2]),
            'GeneID': set([x for x in blast[26].replace(' ', '').split(';') if len(x) > 2]),
            'RefSeq': set([x for x in blast[30].replace(' ', '').split(';') if len(x) > 2]),
            'PANTHER': set([x for x in blast[31].replace(' ', '').split(';') if len(x) > 2]),
            'BioCyc': set([x for x in blast[32].replace(' ', '').split(';') if len(x) > 2]),
            'BRENDA': set([x for x in blast[33].replace(' ', '').split(';') if len(x) > 2]),
            'Reactome': set([x for x in blast[34].replace(' ', '').split(';') if len(x) > 2]),
            'SABIO-RK': set([x for x in blast[35].replace(' ', '').split(';') if len(x) > 2]),
            'SignaLink': set([x for x in blast[36].replace(' ', '').split(';') if len(x) > 2]),
            'SIGNOR': set([x for x in blast[37].replace(' ', '').split(';') if len(x) > 2]),
            'UniPathway': set([x for x in blast[38].replace(' ', '').split(';') if len(x) > 2]),
            'eggNOG': set([x for x in blast[39].replace(' ', '').split(';') if len(x) > 2]),
            'OrthoDB': set([x for x in blast[40].replace(' ', '').split(';') if len(x) > 2]),
            'PDB': set([x for x in blast[41].replace(' ', '').split(';') if len(x) > 2]),
            'DBJ': set(),
            'EMBL': set(),
            'GENBANK': set(),
            'SWISS-Prot': set()
        }
        db = ''
        mp = {'dbj': 'DBJ', 'emb': 'EMBL', 'gb': 'GENBANK', 'gi': 'GeneID', 'ref': 'RefSeq', 'sp': 'SWISS-Prot'}
        for k in blast[4].split('|'):
            if len(k) < 1:
                continue
            if k in mp:
                db = mp[k]
            elif db != '':
                dic[db].add(k)
                db = ''
        self.dbs = dic
        self.domains = {}
        self.protein_names = []
        if len(blast[5].strip()) > 2:
            self.protein_names.append(blast[5].strip())
        if len(blast[18].strip()) > 2:
            self.protein_names.append(blast[18].strip())
        self.evalue = blast[6]
        self.coverage = blast[14]
        self.uniprot_entry = blast[15]
        self.uniprot_entry_name = blast[16]
        self.gene_names = blast[19]
        self.organism_blast = blast[20]
        self.gos = {k.split('[')[1].split(']')[0]:k.split('[')[0].strip() for k in blast[22].split(';')} if blast[22].count(";") > 0 else ( {blast[22].split('[')[1].split(']')[0]:blast[22].split('[')[0].strip() } if len(blast[22]) > 5 else {})
        self.organisms_ID = blast[28].replace(" ", '').split(";")
        self.protein_repet = blast[42]
        self.proteins_families = blast[43]
        self.proteins_motif = blast[44]
        self.proteins_domain1 = blast[45]
        self.proteins_domain2 = blast[46]
        self.ortologs = []
        self.kegg_ortologs = []
        self.BiGG_Reactions = []
        self.OGs = []
        self.COG = ''
        self.prot_anot = ''
        self.protein_seq = seq
        self.protein_len = str(len(seq))
        self.pathways = {}
        self.prize = []
        self.score = 0
        self.scaffold = scaffold
        self.name = ''
    
    def calc_prize(self):
        if self.reviwed:
            self.prize.append("R")
        if len(self.pathways) > 0:
            self.prize.append("P")
        if len(self.ortologs) > 0 and len(self.kegg_ortologs) > 0 and len(self.OGs) > 0 and len(self.COG) > 0:
            self.prize.append("O")
        if len(self.gos) > 0:
            self.prize.append("G")
        if len(self.domains) > 0:
            self.prize.append("D")
        if self.coverage == '100' and self.evalue == '0.0':
            self.prize.append("M")
        if (len(self.gene_names) > 3) and (len("".join(self.protein_names)) > 4) and (len(self.prot_anot) > 3):
            self.prize.append("A")
        self.score = sum([
            7 if 'R' in self.prize else 0,
            6 if 'P' in self.prize else 0,
            5 if 'O' in self.prize else 0,
            4 if 'G' in self.prize else 0,
            3 if 'D' in self.prize else 0,
            2 if 'M' in self.prize else 0,
            1 if 'A' in self.prize else 0
        ])
        return self.score
    
    def mesclar(self, o):
        self.gis.extend([g for g in o.gis if not g in self.gis])
        for k in self.dbs:
            self.dbs[k] = self.dbs[k].union(o.dbs[k])
        self.organisms_ID.extend([i for i in o.organisms_ID if not i in self.organisms_ID])
        self.protein_names.extend([p.strip() for p in o.protein_names if not p in self.protein_names and len(p.strip()) > 2])
        for k in o.gos:
            if not k in self.gos:
                self.gos[k] = o.gos[k]
        self.protein_repet += "; " + o.protein_repet if len(o.protein_repet) > 1 else ""
        self.proteins_families += "; " + o.proteins_families if len(o.proteins_families) > 1 else ""
        self.proteins_motif += "; " + o.proteins_motif if len(o.proteins_motif) > 1 else ""
        self.proteins_domain1 += "; " + o.proteins_domain1 if len(o.proteins_domain1) > 1 else ""
        self.proteins_domain2 += "; " + o.proteins_domain2 if len(o.proteins_domain2) > 1 else ""
        if self.uniprot_entry != o.uniprot_entry:
            print("CONFLITO uniprot_entry %s => %s & %s" % (self.mrna, self.uniprot_entry, o.uniprot_entry))
        if self.uniprot_entry_name != o.uniprot_entry_name:
            print("CONFLITO uniprot_entry_name %s => %s & %s" % (self.mrna, self.uniprot_entry_name, o.uniprot_entry_name))
        if self.gene_names != o.gene_names:
            print("CONFLITO gene_names %s => %s & %s" % (self.mrna, self.gene_names, o.gene_names))
        if self.organism_blast != o.organism_blast:
            print("CONFLITO organism_blast %s => %s & %s" % (self.mrna, self.organism_blast, o.organism_blast))
        if self.organisms_ID != o.organisms_ID:
            print("CONFLITO organisms_ID %s => %s & %s" % (self.mrna, self.organisms_ID, o.organisms_ID))
            
    def setOrt(self, orts):
        if orts and len(orts) > 1:
            self.ortologs = [x for x in orts.replace(' ', '').split(',') if len(x) > 2]
            
    def toFile(self):
        for k in [k for k in self.dbs if len(self.dbs[k]) < 1]:
            del self.dbs[k]
        for k in [k for k in self.gos if len(k) < 4]:
            del self.gos[k]
        return "\t".join([x.replace("\t", " ") for x in 
                     [
                         self.scaffold,
                         self.gene,
                         self.mrna,
                         self.name,
                         self.protein_len,
                         self.evalue,
                         self.coverage,
                         str(self.score),
                         "R" if self.reviwed else "N",
                         ",".join(self.gis),
                         str(self.dbs),
                         str(self.domains),
                         ",".join(self.protein_names),
                         self.gene_names,
                         self.uniprot_entry,
                         self.uniprot_entry_name,
                         self.organism_blast,
                         ','.join(self.gos),
                         ','.join(self.organisms_ID),
                         self.protein_repet,
                         self.proteins_families,
                         self.proteins_motif,
                         self.proteins_domain1,
                         self.proteins_domain2,
                         ','.join(self.ortologs),
                         ','.join(self.kegg_ortologs),
                         ','.join(self.BiGG_Reactions),
                         ','.join(self.OGs),
                         self.COG,
                         self.prot_anot,
                         str(self.pathways),
                         ','.join(self.prize)
                     ]
                     ])    
        
def mesclar(mrnas):
    m = mrnas[0]
    if len(mrnas) < 2: ## R1: so tendo uma anotação persiste ela
        return m
    e = float(m.evalue)
    for m1 in mrnas[1:]:
        if m.evalue != m1.evalue: ## R2: tendo evalues diferentes persiste a de evalue mais proximo de zero
            e1 = float(m1.evalue)
            if e1 < e:
                m = m1
                e = e1
            continue
        ## R3: anotações com evalue igual são mescladas
        m.mesclar(m1)
    return m
        
def importDATA(proteins, gff3, blast, eggnog_gos, eggnog_ort, interpro_tsv, out_file):
    seq_dict = SeqIO.to_dict(SeqIO.parse(proteins, "fasta"))
    gff = {x[8].split(';')[0].split('=')[1]: x[0] for x in [l.strip().split('\t') for l in open(gff3).readlines() if l.count('\tmRNA\t') > 0]}
    anotations = [Anotattion(a, gff[a[2].split('|')[0]], str(seq_dict[a[2]].seq)) for a in [l.strip().split('\t') for l in open(blast).readlines() if not l.startswith("#")]]
    mrnas = {}
    print('importando mrnas ...')
    for mrna in set([x.mrna for x in anotations]):
        ms = [x for x in anotations if x.mrna == mrna]
        mrna_rev = [x for x in ms if x.reviwed]
        mrnas[mrna] = mesclar(mrna_rev if len(mrna_rev) > 0 else ms)
        if (len(mrnas) % 1000 == 0):
            print("%d mrnas importados ..." % len(mrnas))
            anotations = [x for x in anotations if not x.mrna in mrnas]
    print("%d mrnas importados ..." % len(mrnas))
    print("# genes: ", len(set([mrnas[a].gene for a in mrnas])))
    print("# genes rev: ", len(set([mrnas[a].gene for a in mrnas if mrnas[a].reviwed])))
    print("# genes evalue 0: ", len(set([mrnas[a].gene for a in mrnas if mrnas[a].evalue == '0.0'])))
    print("# genes evalue 0 rev: ", len(set([mrnas[a].gene for a in mrnas if mrnas[a].evalue == '0.0' and mrnas[a].reviwed])))
    print("# genes N rev: ", len(set([mrnas[a].gene for a in mrnas if not mrnas[a].reviwed])))
    print("# proteins: ", len(mrnas))
    print("# proteins rev: ", len([a for a in mrnas if mrnas[a].reviwed]))
    print("# proteins evalue 0: ", len([a for a in mrnas if mrnas[a].evalue == '0.0']))
    print("# proteins evalue 0 rev: ", len([a for a in mrnas if mrnas[a].evalue == '0.0' and mrnas[a].reviwed]))
    print("# proteins N rev: ", len([a for a in mrnas if not mrnas[a].reviwed]))
    print("importando ortogos ...")
    eggnog = [l.strip().split('\t') for l in open(eggnog_gos).readlines() if not l.startswith("#")]
    for egg in eggnog:
        if len(egg) < 5:
            continue
        m = egg[0].split("|")[0]
        if not m in mrnas:
            mrnas[m] = Anotattion(['', '', egg[0]], gff[egg[0].split('|')[0]], str(seq_dict[egg[0]].seq))
        mrna = mrnas[m]
        if len(egg[4]) > 2:
            mrna.gene_names += ';' + egg[4]
        if len(egg) > 5 and len(egg[5]) > 3:
            for g in egg[5].replace(" ", '').split(","):
                if not g in mrna.gos:
                    mrna.gos[g] = ""
        if len(egg) > 6 and len(egg[6]) > 3:
            mrna.kegg_ortologs = egg[6].replace(" ", '').split(",")
        if len(egg) > 7 and len(egg[7]) > 3:
            mrna.BiGG_Reactions = egg[7].replace(" ", '').split(",")
        if len(egg) > 9 and len(egg[9]) > 3:
            mrna.OGs = egg[9].replace(" ", '').split(",")
        if len(egg) > 11 and len(egg[11]) > 1:
            mrna.COG = egg[11].strip()
        if len(egg) > 12 and len(egg[12]) > 2:
            mrna.prot_anot = egg[12].strip()
    orts = [l.strip().split('\t') for l in open(eggnog_ort).readlines() if not l.startswith("#")]
    for o in [x for x in orts if len(x) > 1]:
        mrnas[o[0].split("|")[0]].setOrt(o[1])
    print("importando interpro ...")
    interpro = [l.strip().split('\t') for l in open(interpro_tsv).readlines() if not l.startswith("#")]
    for x in interpro:
        if len(x) > 0:
            m = x[0].split('|')[0]
            if not m in mrnas:
                mrnas[m] = Anotattion(['', '', x[0]], gff[x[0].split('|')[0]], str(seq_dict[x[0]].seq))
            mrna = mrnas[m]
            if len(x) > 3:
                banco = x[3]
                if not banco in mrna.domains:
                    mrna.domains[banco] = {}
                pp = (x[6], x[7])
                if not pp in mrna.domains[banco]:
                    mrna.domains[banco][pp] = []
                if len(x[4]) > 1:
                    mrna.domains[banco][pp].append([x[4], x[5]])
                if len(x) > 11:
                    mrna.domains[banco][pp].append([x[11], x[12]])
                    if len(x) > 13:
                        for g in x[13].replace(" ", '').split("|"):
                            if not g in mrna.gos:
                                mrna.gos[g] = ""
                        if len(x) > 14 and len(x[14]) > 4:
                            for pathway in x[14].split('|'):
                                y = pathway.split(':')
                                db = y[0].strip()
                                if not db in mrna.pathways:
                                    mrna.pathways[db] = []
                                mrna.pathways[db].append(y[1].strip())
    
    ps = [mrnas[m].calc_prize() for m in mrnas]
    for i in range(max([m.score for m in mrnas.values()]) + 1):
        print('%d mRNAS with prize %d' % (ps.count(i), i))
    
    orde = []
    ll = "ABCDEFGHIJKLMNOPQRSTUVXWYZ"
    for i in range(max([m.score for m in mrnas.values()]) + 1, -1, -1):
        cont = 1
        v = sorted([m for m in mrnas.values() if m.score == i], key=lambda e: -int(e.protein_len))
        if len(v) > 1:
            l = ll[0]
            ll = ll[1:]
            for mrna in v:
                mrna.name = l + str(cont)
                cont += 1
                orde.append(mrna)
    with open(out_file, 'w') as f:
        f.write('\n'.join([x.toFile() for x in orde]) + '\n')
    return mrnas, orde
    


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: ./groupdata.py proteins.fasta genome.gff3 blast.tsv eggnog_gos.tsv eggnog_ort.tsv interpro.tsv OUT_FILE.tsv")
    else:
        print("import data (%s) ...\n\n" % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        a = importDATA(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
        print("tsv stored in %s (%s) ...\n\nterminado com sucesso\nby mikeias.net" % (sys.argv[7], datetime.now().strftime("%d/%m/%Y %H:%M:%S")))



