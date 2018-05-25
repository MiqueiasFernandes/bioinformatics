
## Autor: Miquéias Fernandes 25/05/18
## www.mikeias.net

from Bio import SeqIO
import operator
import sys

class Gene:
    def __init__(self, nome, scaffold, offset):
        self.nome = nome
        self.scaffold = scaffold
        self.score = int(offset)
        self.id = -1
        scaffold.genes.append(self)

class Scaffold:
    def __init__(self, nome, id, tabela):
        self.nome = nome
        self.id = id
        self.score = (id + (contScaffold(id, tabela) * 0.00001))
        self.genes = []

def contScaffold(id, tb):
    cont = 0
    for scf in tb:
        if scf.id == id:
            cont += 1
    return(cont)

def getSCF(nome, tb):
    for scf in tb:
        if scf.nome == nome:
            return(scf)

def sortTable(table):
    sorted_x = sorted(tabela, key=operator.attrgetter('score'))
    cont = 1
    for k in sorted_x:
        k.id = cont
        cont += 1

# PADRAO PG_________G___________
#           SCAFFOLD   OFFSET

if len(sys.argv) < 3:
    print("use: ./geneID scaffolds.fasta genes.gff")
else:
    with open(sys.argv[1], "rU") as handle:
        tabela = []

        print("Ordenando scaffolds no fasta ...")
        for record in SeqIO.parse(handle, "fasta"):
            tabela.append(
                Scaffold(record.name, len(record.seq), tabela)
            )

        sortTable(tabela)

        print("Salvando ids scaffolds ...")
        fo = open("output-scaffolds-ids.csv", "w")
        fo.write("scaffold,id\n");
        for scf in tabela:
            fo.write(str(scf.nome) + "," + str(scf.id) + "\n")
        fo.close()

        print("Importando informações de genes ...")
        file = open(sys.argv[2], "r")
        k = 0
        emGene = False
        tabela2 = []
        for line in file:
            if (emGene):
                if (line.find("AUGUSTUS	gene") >= 0):
                    dt = line.split("\t")
                    if (len(dt) > 5):
                        scf = dt[0]
                        pos = dt[3]
                        nome = dt[8][:-1]
                        scaffold = getSCF(scf, tabela)
                        tabela2.append(Gene(nome, scaffold, pos))
                        k += 1
                    emGene = False
            elif (line.startswith("# start gene ")):
                emGene = True
        print(str(k) + " genes a orndenar ...")
        for scaffold in tabela:
            sortTable(scaffold.genes)
        print("Salvando ids dos genes...")

        fo = open("output-genes-ids.csv", "w")
        fo.write("gene,id\n");
        for gene in tabela2:
            fo.write(str(gene.nome) + "," +
                     "PG" + str(gene.scaffold.id) + "G" + str(gene.id) +
                     "\n")
        fo.close()
        print("terminado com sucesso!")
        print("by mikeias.net")