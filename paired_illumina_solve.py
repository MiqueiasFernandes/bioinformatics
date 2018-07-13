

from Bio import SeqIO
import sys


def getIds(fasta):
    ids = []
    for record in SeqIO.parse(fasta, "fastq"):
        ids.append(record.description)
        if len(ids) % 100000 == 0:
            print("cont " + str(len(ids)))
    return ids



if __name__ == "__main__":

    sys.argv.extend([
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/base/illumina1.fa',
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/base/illumina1.fa'
    ])

    if len(sys.argv) < 3:
        print("usage: ./paired_illumina_solve illumina1.fa illuma2.fa")
    else:
        ill1 = sys.argv[1]
        ill2 = sys.argv[2]

        id1 = getIds(ill1)
        id2 = getIds(ill2)



        print('terminado com sucesso...')




