#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 08/06/2018

import os, sys
from sortscaffolds import Gene, Scaffold, Data, getTabelas
import pymysql.cursors

__author__ = "Miquéias Fernandes"
__copyright__ = "Copyright 2018 mikeias.net"
__license__ = "MIT"
__version__ = "1.1"


MAX_SIZE_SEQ_SCAFFOLD = 1000000  # up to 1MB store on disk


class ScaffoldMySql:

    SQL_INSERT = "INSERT INTO `scaffold`(`name`, `score`, `jhi_size`, `sequencetype`, `jhi_sequence`) VALUES (%s, %s, %s, %s, %s)"
    SQL_SELECT = "SELECT `id` FROM `scaffold` WHERE `name`=%s"
    SQL_SELECT_IDS = "SELECT `name`,`id` FROM `scaffold` WHERE `id`>0"

    def __init__(self, scf: Scaffold):
        self.id = -1
        self.name = scf.nome #String,
        self.score = scf.score #Integer,
        self.size = scf.size #Integer,
        self.sequence_type = scf.size > MAX_SIZE_SEQ_SCAFFOLD #Boolean, ## true if path
        self.sequence = str(scf.sequence)  # TextBlob

    def getSql(self):
        if self.sequence_type:
            return self.name, self.score, self.size, self.sequence_type, 'scaffolds/' + self.name + '.fa'
        return self.name, self.score, self.size, self.sequence_type, self.sequence


class GeneMySql:
    SQL_INSERT = "INSERT INTO `gene`(`nome`, `scaffold_id`) VALUES (%s, %s)"
    SQL_SELECT = "SELECT `id` FROM `gene` WHERE `nome`=%s"
    SQL_SELECT_IDS = "SELECT `nome`,`id` FROM `gene` WHERE `id`>0"

    def __init__(self, gene: Gene):
        self.nome = gene.pgID  # String,
        self.ortologo_eucalipto = None # String,
        self.ortologo_arabidopsis = None # String
        self.scaffold = gene.scaffold.nome
        self.id = -1

    def getSql(self, scaffold):
        return self.nome, scaffold


class FeatureMySql:
    SQL_INSERT = "INSERT INTO `feature`(`name`) VALUES (%s)"
    SQL_SELECT = "SELECT `id` FROM `feature` WHERE `name`=%s"
    SQL_SELECT_IDS = "SELECT `name`,`id` FROM `feature` WHERE `id`>0"

    def __init__(self, data: Data):
        self.name = data.feature  # String
        self.id = -1

    def getSql(self):
        return self.name


class GffMySql:
    SQL_INSERT = "INSERT INTO `gff`(`jhi_start`,`jhi_stop`,`jhi_sequence`,`score`,`strand`,`frame`,`attribute`,`feature_id`,`gene_id`) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)"

    def __init__(self, data: Data, gene: Gene):
        self.start = data.start # Integer,
        self.stop = data.end # Integer,
        self.sequence = str(gene.scaffold.sequence[data.start:data.end]) # TextBlob,  ## VER RANGES ################################
        self.score = data.score # Float,
        self.strand = data.strand # String,
        self.frame = data.frame # String,
        self.attribute = data.attributes # String,
        self.feature = data.feature
        self.gene = gene.pgID
        self.id = -1

    def getSql(self, f, g):
        return self.start, self.stop, self.sequence, self.score, self.strand, self.frame, self.attribute, f, g


# class BestMatchMySql:
#     def __init__(self, bestm):
#     self.code =  if 'code' in bestm.keys(): bestm['code'] else: None  #String,
#     self.maxScore #Integer,
#     self.totalScore #Integer,
#     self.queryCover #Integer,
#     self.eValue #Float,
#     self.ident #Integer,
#     self.offset #Integer,
#     self.organism #String,
#     self.definition #String,
#     self.mol_type #String,
#     self.product #String,
#     self.protein_id #String,
#     self.translation #String,
#     self.note #String


def insert(connection, sql, data):
    try:
        with connection.cursor() as cursor:
            cursor.execute(sql, data)
            connection.commit()
    except:
        print("ERRO: impossivel inserir tupla " + str(data) + ' query: ' + sql)


def select(connection, sql, data, one=True):
    try:
        with connection.cursor() as cursor:
            cursor.execute(sql, data)
            if one:
                result = cursor.fetchone()
            else:
                result = cursor.fetchall()
            return result
    except:
        print("ERRO: impossivel consultar tupla " + str(data) + ' query: ' + sql)


def process(fasta, gff, skip_scaffolds=False, skip_genes=False, skip_features=False, skip_gff=False):

    print('importando as informações ...')
    ret = getTabelas(fasta, gff)

    scaffolds = ret[0]
    scf_dict = {}

    genes = ret[1]
    gene_dict = {}

    feat_dict = {}

    tb_scaffold = []
    tb_gene = []
    tb_feature = []
    tb_gff = []

    print('mapeando as informações ...')
    for sc in scaffolds:
        s = ScaffoldMySql(sc)
        tb_scaffold.append(s)
        scf_dict[s.name] = s

    feat = []
    for gn in genes:
        g = GeneMySql(gn)
        tb_gene.append(g)
        gene_dict[g.nome] = g
        for dt in gn.gffs:
            tb_gff.append(GffMySql(dt, gn))
            if dt.feature not in feat:
                f = FeatureMySql(dt)
                tb_feature.append(f)
                feat_dict[f.name] = f
                feat.append(dt.feature)

    print(str(len(tb_scaffold)) + ' scaffolds importados ...')
    print(str(len(tb_gene)) + ' genes importados ...')
    print(str(len(tb_feature)) + ' features importadas ...')
    print(str(len(tb_gff)) + ' gff row importados ...')

    print('conectando ao banco de dados ...')
    try:
        connection = pymysql.connect(host='localhost',
                                     user='root',
                                     password='123',
                                     db='guava',
                                     charset='utf8mb4',
                                     cursorclass=pymysql.cursors.DictCursor)
    except:
        print("ERRO: Impossivel conectar a base de dados.")
    else:
        print("conectado com sucesso ... ")

        if skip_scaffolds:
            print('skip scaffolds ...')
            ret = select(connection, ScaffoldMySql.SQL_SELECT_IDS, None, False)
            for s in ret:
                scf_dict[s['name']].id = s['id']
        else:
            print('persistindo scaffolds')
            dirpath = False
            id = 0
            for scaffold in tb_scaffold:
                if not skip_scaffolds:
                    if scaffold.sequence_type:
                        try:
                            if not dirpath:
                                os.mkdir('scaffolds')
                                dirpath = True
                            fo = open("scaffolds/" + scaffold.name + '.fa', "w")
                            fo.write(">" + scaffold.name + '\n')
                            fo.write(scaffold.sequence)
                            fo.close()
                        except:
                            print('impossivel salvar scaffold local ' + scaffold.name)
                            exit(-1)
                    insert(connection, ScaffoldMySql.SQL_INSERT, scaffold.getSql())
                    scaffold.id = select(connection, ScaffoldMySql.SQL_SELECT, scaffold.name)['id']
                print('.', end='')
                sys.stdout.flush()
                id += 1
            print('\n ' + str(len(tb_scaffold) - id) + ' scaffolds falharam ser persistidos ....')

        if skip_genes:
            print('skip genes ...')
            ret = select(connection, GeneMySql.SQL_SELECT_IDS, None, False)
            for g in ret:
                gene_dict[g['nome']].id = g['id']
        else:
            print('persistindo genes ...')
            id = 0
            for gene in tb_gene:
                s = scf_dict[gene.scaffold].id
                if s > 0:
                    insert(connection, GeneMySql.SQL_INSERT, gene.getSql(s))
                    print('.', end='')
                    sys.stdout.flush()
                    id += 1
                    gene.id = select(connection, GeneMySql.SQL_SELECT, gene.nome)['id']
                    gene_dict[gene.nome] = gene
                else:
                    print('Gene ' + gene.nome + ' não pode ser inserido por erro no scaffold ' + gene.scaffold)
            print('\n ' + str(len(tb_gene) - id) + ' genes falharam ser persistidos ....')

        if skip_features:
            print('skip features ...')
            ret = select(connection, FeatureMySql.SQL_SELECT_IDS, None, False)
            for f in ret:
                feat_dict[f['name']].id = f['id']
        else:
            print('persistindo features ...')
            id = 0
            for feature in tb_feature:
                insert(connection, FeatureMySql.SQL_INSERT, feature.getSql())
                print('.', end='')
                sys.stdout.flush()
                id += 1
                feature.id = select(connection, FeatureMySql.SQL_SELECT, feature.name)['id']
                feat_dict[feature.name] = feature
            print('\n ' + str(len(tb_feature) - id) + ' features falharam ser persistidos ....')

        if not skip_gff:
            print('persistindo gffs ...')
            id = 0
            for gff in tb_gff:
                f = feat_dict[gff.feature].id
                g = gene_dict[gff.gene].id
                if f > 0 and g > 0:
                    insert(connection, GffMySql.SQL_INSERT, gff.getSql(f, g))
                    print('.', end='')
                    sys.stdout.flush()
                    id += 1
                    gff.id = id
                else:
                    print('falhou inserir gff porque gene => ' + str(g) + ' ou feature => ' + str(f) + ' falhou.')
            print('\n ' + str(len(tb_gff) - id) + ' gff falharam ser persistidos ....')
    finally:
        connection.close()










if __name__ == "__main__":

    sys.argv.extend(['/home/mfernandes/Documentos/genome.v0.1/pguajava.genome.ufes.v0.1.fasta', '/home/mfernandes/Documentos/genome.v0.1/pguajava.v0.1.augustus.abinitio.gff'])

    if len(sys.argv) < 3:
        print("use: ./geneID scaffolds.fasta genes.gff")
    else:
        process(sys.argv[1], sys.argv[2])
