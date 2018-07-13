#!/usr/bin/env python
## Autor: Miquéias Fernandes 01/06/18
## www.mikeias.net

## thanks https://github.com/PyMySQL/PyMySQL

import sys
from Bio import Entrez
import pymysql.cursors



Entrez.email = 'contatomiqueiasfernandes@hotmail.com'

class Gene:
    def __init__(self, code, scaffold):
        self.code = code
        self.scaffold = scaffold
        handle = Entrez.efetch(db="nuccore", id=code, rettype="gb", retmode="xml")
        self.gene = Entrez.read(handle)
        self.organism = self.gene[0]["GBSeq_organism"]
        self.definition = self.gene[0]["GBSeq_definition"]
        self.mol_type = self.gene[0]["GBSeq_moltype"]
        for feature in self.gene[0]["GBSeq_feature-table"]:
            if feature['GBFeature_key'] == 'gene':
                for qual in feature['GBFeature_quals']:
                    if qual['GBQualifier_name'] == 'note':
                        self.note = qual['GBQualifier_value']
            if feature['GBFeature_key'] == 'CDS':
                for qual in feature['GBFeature_quals']:
                    if qual['GBQualifier_name'] == 'product':
                        self.product = qual['GBQualifier_value']
                    elif qual['GBQualifier_name'] == 'protein_id':
                        self.protein_id = qual['GBQualifier_value']
                    elif qual['GBQualifier_name'] == 'translation':
                        self.translation = qual['GBQualifier_value']

    def hasorganism(self):
        return ('organism' in self.__dict__.keys())
    def hasdefinition(self):
        return ('definition' in self.__dict__.keys())
    def hasmol_type(self):
        return ('mol_type' in self.__dict__.keys())
    def hasproduct(self):
        return ('product' in self.__dict__.keys())
    def hasprotein_id(self):
        return ('protein_id' in self.__dict__.keys())
    def hastranslation(self):
        return ('translation' in self.__dict__.keys())
    def hasnote(self):
        return ('note' in self.__dict__.keys())

    def fields(self):
        return (self.code +
                 ("," + self.organism if self.hasorganism() else "") +
                 ("," + self.definition if self.hasdefinition() else "") +
                 ("," + self.mol_type if self.hasmol_type() else "") +
                 ("," + self.product if self.hasproduct() else "") +
                 ("," + self.protein_id if self.hasprotein_id() else "") +
                 ("," + self.translation if self.hastranslation() else "") +
                 ("," + self.note if self.hasnote() else "") )

    def __str__(self):
        return ("{code: " + self.code +
                 (",\n organism: " + self.organism if self.hasorganism() else "") +
                 (",\n definition: " + self.definition if self.hasdefinition() else "") +
                 (",\n mol_type: " + self.mol_type if self.hasmol_type() else "") +
                 (",\n product: " + self.product if self.hasproduct() else "") +
                 (",\n protein_id: " + self.protein_id if self.hasprotein_id() else "") +
                 (",\n translation: " + self.translation if self.hastranslation() else "") +
                 (",\n note: " + self.note if self.hasnote() else "") + "}" )

class Hint:
    def __init__(self, scaffold, gene):
        self.scaffold = scaffold
        self.gene = gene

def search():
    fo = open("gene_desc.csv", "w")
    fo.write("code,organism,definition,mol_type,product,protein_id,translation,note\n");
    fo.close()
    print("Buscando sequencias no NCBI ... ")
    for gene_id in ["XM_022518206.1", "XM_014803286.1", "XM_014803286.1"]:
        gene = Gene(gene_id)
        fo = open("gene_desc.cst", "a")
        fo.write(gene.fields() + "\n")
        fo.close()
        print(".", end='')

def getHints(gff):
    file = open(gff, "r")
    k = 0
    emHint = False
    tabela = []
    for line in file:
        if emGene:
            if line.find("\tAUGUSTUS\tgene\t") >= 0:
                dt = line.split("\t")
                if len(dt) > 5:
                    scf = dt[0]
                    pos = dt[3]
                    nome = dt[8][:-1]
                    scaffold = getscf(scf, tabela)
                    tabela2.append(Gene(nome, scaffold, pos))
                    k += 1
                    emGene = False
        elif line.startswith("# start gene "):
            emGene = True

def insert(connection, table, nome, code, note):
    try:
        with connection.cursor() as cursor:
            # Create a new record
            sql = "INSERT INTO `" + table + "` (`nome`, `code`, `note`, `scaffold_id`, `source_id`, `feature_id`) VALUES (%s, %s, %s, %s, %s, %s)"
            cursor.execute(sql, (nome, code, note, "1", "1", "1"))
            connection.commit()
    except:
        print("ERRO: impossivel inserir tupla " + str((nome, code, note)))


def main(gff, database, table, password):
    try:
        connection = pymysql.connect(host='localhost',
                                     user='root',
                                     password=password,
                                     db=database,
                                     charset='utf8mb4',
                                     cursorclass=pymysql.cursors.DictCursor)
    except:
        print("ERRO: Impossivel conectar a base de dados.")
    else:
        print("conectado com sucesso ... ")
        with connection.cursor() as cursor:
            # Read a single record
            sql = "SELECT `id`, `nome`, `organism`, `note` FROM `" + table + "` WHERE `code`=%s"
            cursor.execute(sql, ('XM_022518206.1',))
            result = cursor.fetchone()
            print(result)
            insert(connection, table, "gen09g35", "XM_022518207.1", "note")


if __name__ == "__main__":
    if not len(sys.argv) == 5:
        print("usage: ./geneimport gff database table password")
    else:
        main(*sys.argv[1:])





