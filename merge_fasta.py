

from Bio import SeqIO
from shutil import copyfile
import sys
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import threading


class Seq (threading.Thread):
    def __init__(self, princ, base, dirs, out):
        threading.Thread.__init__(self)
        self.princ = princ
        self.dirbase = base
        self.arqbase = base + '/' + princ
        self.dirs = dirs
        self.out = out
        self.final = []
        self.eliminados = -1

    def run(self):
        print('a processar ' + self.princ + ' ...')
        self.process()

    def getIds(self, fasta):
        ids = []
        ill = self.isIllumina() and fasta == self.arqbase
        if ill:
            print("tratado como illumina ...")
        # try:
        for record in SeqIO.parse(fasta, "fastq" if ill else "fasta"):
            ids.append(record.description)
        # except:
        #     print('erro em ' + fasta + ' modo ' + ('ILL' if ill else 'FASTA'))
        return ids

    def merge(self, ids):
         qtd = len(ids)
         enc = {}
         for listIds in ids:
             for id in listIds:
                 i = id.split(' ')[0]
                 if i in enc.keys():
                     enc[i] += 1
                 else:
                     enc[i] = 1
         final = []
         elim = 0
         for id in enc.keys():
            if enc[id] == qtd:
                final.append(id)
            else:
                elim += 1
         self.final = final
         self.eliminados = elim
         print('projeto: ' + self.princ + ': ' + str(elim) + ' eliminados ...')

    def isIllumina(self):
        return self.princ.startswith('illumina')

    def process(self, ):
        ids = []
        dirs = []
        dirs.append(self.dirbase)
        dirs.extend(self.dirs)
        for dir in dirs:
            print('tratando file ' + dir + "/" + self.princ)
            ids.append(self.getIds(dir + "/" + self.princ))
        self.merge(ids)

    def persist(self, pareado=None):
        #self.process()
        if self.eliminados > 0:
            out_handle = open(self.out + '/' + self.princ, 'w')
            if self.isIllumina():
                add = []
                for id in self.final:
                    str = id.split(" ")[0] + " 2" + id.split(" ")[1][1:]  ## preservar paridade
                    if str in pareado.final:
                        add.extend([id, str])
                with open(self.arqbase) as in_handle:
                    for title, seq, qual in FastqGeneralIterator(in_handle):
                        if title in add:
                            out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            else:
                with open(self.arqbase) as in_handle:
                    for title, seq in SimpleFastaParser(in_handle):
                        if title.split(' ')[0] in self.final:
                            out_handle.write(">%s\n%s\n" % (title, seq))
            out_handle.close()
        else:
            copyfile(self.arqbase, self.out + '/' + self.princ)





if __name__ == "__main__":


    sys.argv.extend([
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/base',
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/arab' ,
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/arroz',
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/arrjp',
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/clr1' ,
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/clr2' ,
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/sorgo' ,
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/uva',
        '/home/mfernandes/Documentos/relatorio_mestrado/julho/removerseq/out'
    ])

    if len(sys.argv) < 4:
        print("usage: ./merge_fasta.py dir_base dir_rem1 dir_rem2 dir_rem3 ... output_dir [Ñ USE /]")
    else:
        base = sys.argv[1]
        outdir = sys.argv[-1]
        dirs = sys.argv[2:-1]

        proj = []

        for arq in os.listdir(base):
            s = Seq(arq, base, dirs, outdir)
            proj.append(s)
            s.start()
            print('projeto ' + s.princ + ' criado ...')


        for s in proj:
            s.join()

        for s in proj:
            print('a persistir ' + s.princ + ' ...')
            s.persist()

        print('terminado com sucesso...')



