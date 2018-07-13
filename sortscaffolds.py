
## Autor: Miquéias Fernandes 25/05/18
## www.mikeias.net

from Bio import SeqIO
import operator
import sys
import re

class Gene:
    def __init__(self, nome, scaffold, offset):
        self.nome = nome
        self.scaffold = scaffold
        self.offset = offset
        self.score = int(offset)
        self.id = -1
        scaffold.genes.append(self)
        self.gffs = []
        self.protein = ''
        self.pgID = None
        self.data = None
        self.sobrepoe = None

    def addData(self, line, l, gene=False):
        dt = Data(line, l)
        self.gffs.append(dt)
        if gene:
            self.data = dt

    def getSeq(self):
        return str(self.scaffold.sequence)[self.data.start:self.data.end]


class Scaffold:
    def __init__(self, nome, id, seq):
        self.nome = nome
        self.size = id
        self.id = id
        self.score = id
        self.genes = []
        self.genesByID = {}
        self.sequence = seq

    def update(self):
        for g in self.genes:
            self.genesByID[g.id] = g


class Data:
    def __init__(self, cst, linha):
        self.raw = cst.rstrip()
        self.split = self.raw.split('\t')
        self.seqid = self.split[0]
        self.source = self.split[1]
        self.feature = self.split[2]
        try:
            self.start = int(self.split[3])
            self.end = int(self.split[4])
            self.score = None
            if not self.split[5] == '.':
                self.score = float(self.split[5])
        except:
            print(str(linha) + ' => error na linha ' + cst)
            exit(-1)
        self.strand = self.split[6]
        self.frame = self.split[7]
        self.attributes = self.split[8]


def sortTable(tabela, rev=False):
    sorted_x = sorted(tabela, key=operator.attrgetter('score'), reverse=rev)
    cont = 0
    for k in sorted_x:
        cont += 1
        k.id = cont

# PADRAO PG_________G___________
#           SCAFFOLD   OFFSET


def getTabelas(fasta, gff, cloroplast=[], verificarsobreposicao=False, verbose=False, storefasta=False, min=0):
    with open(fasta, "rU") as handle:
        scaffolds = []
        menores = []
        rm = 0
        sums = 0
        maxi = 0
        cort = 0
        print("Importando scaffolds do fasta ...")
        for record in SeqIO.parse(handle, "fasta"):
            if not record.name in cloroplast:
                if len(record.seq) >= min:
                    scaffolds.append(
                        Scaffold(record.name, len(record.seq), record.seq)
                    )
                else:
                    cort += 1
                    menores.append(record.name)
            else:
                rm += 1
                sums += len(record.seq)
                maxi = max(maxi, len(record.seq))
        print(str(rm) + " scaffold de cloroplasto removidos ...")
        if verbose:
            print(str(sums) + " len total removidos ...")
            print(str(sums/rm) + " media size removidos ...")
            print(str(maxi) + " maior scaffold removido ...")
        if cort > 0:
            print(str(cort) + " scaffolds menores que " + str(min) + " removidos ...")
        if rm > 0:
            fw = open(fasta + ".without-ncrna.fa", "w")
            for scf in scaffolds:
                fw.write(">" + scf.nome + "\n" + str(scf.sequence) + "\n")
            fw.close()

        print("Ordenando scaffolds do fasta ...")
        sortTable(scaffolds, True)
        scf_dict = {}
        for scf in scaffolds:
            scf_dict[scf.nome] = scf

        print("Importando informações de genes ...")
        file = open(gff, "r")
        emData = False
        genes = []
        gene = None
        l = 0
        skip = 0
        for line in file:
            l += 1
            if "\tgene\t" in line:
                emData = True
                dt = line.split("\t")
                if len(dt) > 5:
                    nscf = dt[0]
                    pos = dt[3]
                    nome = dt[8][:-1]
                    scaffold = scf_dict[nscf] if nscf in scf_dict.keys() else None
                    if scaffold is None:
                        if nscf in cloroplast or nscf in menores:
                            skip += 1
                            gene = None
                            emData = False
                            continue
                        else:
                            print('gene ' + nome + ' não tem scaffold!')
                            exit(-1)
                    gene = Gene(nome, scaffold, pos)
                    gene.addData(line, l, True)
                    genes.append(gene)
            elif emData:
                if line.startswith("# protein sequence = ["):
                    gene.protein += re.sub(r'\W', "", line[22:])
                elif line.startswith("# end gene g"):
                    gene = None
                    emData = False
                elif line.startswith("#"):
                    gene.protein += re.sub(r'\W', "", line[2:])
                else:
                    gene.addData(line, l)
        print(str(skip) + " genes sem scaffolds pulados ...")
        print(str(len(genes)) + " genes a ordenar ...")
        k = 0
        for scaffold in scaffolds:
            sortTable(scaffold.genes)
            scaffold.update()
            if len(scaffold.genes) > 1 and verificarsobreposicao:
                pt = []
                if verbose:
                    print('Verificando sobreposição em scaffold ' + scaffold.nome)
                cont = 0
                for i in range(1, len(scaffold.genes)):
                    g1 = scaffold.genesByID[i]
                    if verbose:
                        print(str(g1.id), end=' ')
                    g2 = scaffold.genesByID[i+1]
                    if g2.data.start <= g1.data.end:
                        g2.sobrepoe = g1
                        if verbose:
                            print('<= ')
                        cont += 1
                        k += 1
                        pt.append('WARN: sobreposição entre os genes ' + g1.nome + ' e ' + g2.nome + ' em ' + str(g1.data.end - g2.data.start))
                if verbose:
                    for p in pt:
                        print(p)
                    print('\nscaffold ' + scaffold.nome + ' possui ' + str(cont) + ' sobreposicoes.')
        if verificarsobreposicao:
            print('foram encontradas ' + str(k) + ' sobreposicoes.')

        print("Criando ids dos genes...")
        for gene in genes:
            gene.pgID = "PG" + str(gene.scaffold.id) + "G" + str(gene.id)

        return scaffolds, genes


if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("use: ./geneID scaffolds.fasta genes.gff cloroplast_scaffold_list.txt")
        exit(-1)
    else:
        clr = []
        for line in open(sys.argv[3]):
            clr.append(line[:-1])
        print(str(len(clr)) + " cloroplast nuc importados...")
        ot = sys.argv[4:]
        min = int(ot[ot.index('-m') + 1]) if '-m' in ot else 0
        ret = getTabelas(sys.argv[1], sys.argv[2], clr, '-s' in ot, '-v' in ot, '-f' in ot, min)
        print("Salvando ids dos genes...")
        fo = open("output-genes-ids.csv", "w")
        fo.write("gene,id\n");
        for gene in ret[1]:
            fo.write(str(gene.nome) + "," + gene.pgID + "\n")
        fo.close()
        print("terminado com sucesso!")
        print("by mikeias.net")
