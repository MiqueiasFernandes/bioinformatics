#!/usr/bin/python3

import sys
from BCBio import GFF
from Bio import SeqIO

def process(in_file, file_pos):
    
    limit_info = dict(gff_type = ["gene"])
    in_handle = open(in_file)
    gff3 = GFF.parse(in_handle, limit_info=limit_info) 
    genes = {}
    for g in list(gff3):
        genes[g.id] = list(g.features)
    
    in_handle.close()
    
    for line in open(file_pos):
        spl = line.split(' ')
        if spl[1] in genes:
            for g in genes[spl[1]]:
                start = int(spl[2])
                end = start + int(spl[3])
                if start in g.location or end in g.location:
                    print(g.id, end="")
                    break
        print("," + spl[0])


if __name__ == "__main__":
    if len(sys.argv) < 3:
        ### file.pos.of.seqs
        ### id seq   start   len
        ### 01 Chr01 1490    2
        ### 02 Chro5 3535645 200
        ###...
        print("usage: ./subpart_from_gff3.py file.gene.gff3 file.pos.of.seqs")
    else:
        process(sys.argv[1], sys.argv[2])



