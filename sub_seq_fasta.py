#!/usr/bin/python3

from Bio import SeqIO
import sys


def process(file, sequence, pos_init, length):
    record_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    seq = record_dict[sequence]
    if not seq is None:
        if pos_init is None:
            return str(seq.seq)
        if length is None:
            return str(seq.seq[int(pos_init):])
        start = int(pos_init)
        end = start + int(length)
        return str(seq.seq[start:min(end, len(seq.seq)-1)])
    return ""


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: sub_seq_fasta.py file.fasta sequence [pos_init opcional] [length opcional]")
    else:
        print(process(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]))





