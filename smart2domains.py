#!/usr/bin/env python
# coding: utf-8

## LICENCE MIT
## REV 08/19
## Extract fasta domains from file get`ed of SMART database
## Usage: smart2domains.py smart.file proteins.fasta output.fasta
## www.mikeias.net
## bio@mikeias.net

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord


class SMART:
    def __init__(self, raw, fasta):
        self.raw = raw
        self.user_protein_id = raw[0].split('=')[1].strip()
        self.smart_protein_id = raw[1].split('=')[1].strip()
        self.features_found = int(raw[2].split('=')[1].strip())
        self.features = []
        for i in range(self.features_found):
            self.features.append({y.split('=')[0].strip():y.split('=')[1].strip() for y in [x for x in raw[3 + (i*7):3 + ((i+1)*7)] if x.count('=') > 0]})
        self.export = [Feature(self, f, fasta) for f in self.features] 
        
    def __repr__(self):
        return "UPID: %s | SPID: %s | %d | %s" % (self.user_protein_id, self.smart_protein_id, self.features_found, str(self.features))
       
class Feature:
    def __init__(self, smart, feature, fasta):
        self.smart = smart
        self.feature = feature
        self.seq = SeqRecord(
            Seq(str(fasta[smart.user_protein_id][int(feature['START'])-1:int(feature['END'])].seq), generic_protein),
            id='|'.join([smart.user_protein_id, feature['DOMAIN'], feature['EVALUE']]), 
            description=smart.smart_protein_id+ " => " + "; ".join(["%s: %s" % (x,feature[x]) for x in ('STATUS', 'TYPE', 'START', 'END')]))

def parseFile(file, fasta):
    print('loading fasta ...')
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    print(str(len(record_dict)) + " seqs loaded!")
    print('loading smart ...')
    smarts = []
    seqs = []
    with open(file) as f:
        for l in [l.strip() for l in f.readlines()]:
            if l.startswith('-- SMART RESULTS TEXTFORMAT --'):
                tmp = []
            elif l.startswith('-- FINISHED --'):
                s = SMART(tmp, record_dict)
                smarts.append(s)
                print("%s => %d domains" % (s.user_protein_id, len(s.export)))
                seqs.extend(s.export)
            else:
                tmp.append(l)
    print(str(len(smarts)) + " SMART records loaded!")
    return smarts, [x.seq for x in seqs]

def exportFasta(seqs, output):
    print('exporting fasta ...')
    SeqIO.write(seqs, output, "fasta")
        
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: smart2domains.py smart.file proteins.fasta output.fasta")
    else:
        exportFasta(parseFile(sys.argv[1], sys.argv[2])[1], sys.argv[3])
        print("Terminado com sucesso!\nby mikeias.net")

