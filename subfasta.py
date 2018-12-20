#!/usr/bin/python3

from Bio import SeqIO
import sys



def process(fasta, filtro, out):
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    sequences = []
    ids = []
    for line in open(filtro):
        if not line in ids:
            sequences.append(record_dict[line.strip()])
            ids.append(line)
    
    SeqIO.write(sequences, out, "fasta")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: sub_fasta.py file.fasta filtros outfile")
    else:
        process(sys.argv[1], sys.argv[2], sys.argv[3])
        print("Terminado com sucesso!")



