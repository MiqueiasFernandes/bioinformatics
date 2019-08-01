#!/usr/bin/env python
# coding: utf-8


## LICENCE MIT
## REV 07/19
## Persist data of bioinformatic analisys on MySQL data base
## Usage: loaddb.py congif_file.txt
## www.mikeias.net
## bio@mikeias.net


import sys
from datetime import datetime
from Bio import SeqIO
import ast
import mysql.connector

MAX_ROW_INSERT_PER_COMMIT = 100000



def connect(host, user, passwd, db):
    return mysql.connector.connect( host=host,user=user,passwd=passwd,database=db)

def select(db, table):
    mycursor = db.cursor()
    mycursor.execute("SELECT * FROM " + table)
    return mycursor.fetchall()

def insert(db, table, values, debug=False, file=None):
    try:
        mycursor = db.cursor()
        cont = 0
        for thead in set([','.join(list(x.data.keys())) for x in values]):
            columns = thead.split(',')
            subvals = [x for x in values if ','.join(list(x.data.keys())) == thead]
            sql = "INSERT INTO " + table + " (" + ", ".join(columns) + ") VALUES (%(" + ")s, %(".join(columns) + ")s)"
            for i in range(0, len(subvals), MAX_ROW_INSERT_PER_COMMIT):
                ate = min(len(subvals), i+MAX_ROW_INSERT_PER_COMMIT)
                for record in subvals[i:ate]:
                    mycursor.execute(sql, record.data)
                    record.setSqlId(mycursor.lastrowid)
                    cont += 1
                    if not file is None:
                        file.write((sql % record.data) + ';\n')
                print("%d records of `%s` inserted (%.0f%%)." % (cont, table.upper(), (100*cont/len(values))))
        db.commit()
        return cont
    except Exception as e:
        print('ERRO:')
        print(e)
        print(sql + ' => ' + str(record.data))
        
    
def select(db, table):
    mycursor = db.cursor()
    mycursor.execute("SELECT * FROM " + table)
    return mycursor.fetchall()



class Record:
    def __init__(self, data):
        self.data = data
        self.sqlid = None
        self.updateFields = {}
        
    def setSqlId(self, _sqlid):
        self.sqlid = _sqlid
        for field in self.updateFields:
            for record in self.updateFields[field]:
                record.data[field] = _sqlid
    
    def addFieldUpdate(self, field, record):
        if not field in self.updateFields:
            self.updateFields[field] = []
        self.updateFields[field].append(record)



def loadDB(fastaGenome, fastaProteins, gffGenome, tsvRRNA, tsvTRNA, tsvAnnotation, tsvRepeatsMISA, gffRepeats, gos_data):
    ##load fastas
    seq_dict_genome = SeqIO.to_dict(SeqIO.parse(fastaGenome, "fasta"))
    seq_dict_proteins = SeqIO.to_dict(SeqIO.parse(fastaProteins, "fasta"))
    ##load GFFS
    genomicGFF = [l.strip().split('\t') for l in open(gffGenome).readlines() if not l.startswith('#') and l.count('\t') == 8]
    repeatsGFF = [l.strip().split('\t') for l in open(gffRepeats).readlines() if not l.startswith('#') and l.count('\t') == 8]
    ##load TSVs
    repeatsTSV = [l.strip().split('\t') for l in open(tsvRepeatsMISA).readlines() if not l.startswith('#')]
    rrnaTSV = [l.strip().split('\t') for l in open(tsvRRNA).readlines() if not l.startswith('#')]
    trnaTSV = [l.strip().split('\t') for l in open(tsvTRNA).readlines() if not l.startswith('#')]
    annotationTSV = [l.strip().split('\t') for l in open(tsvAnnotation).readlines() if not l.startswith('#')]
    gosTSV = {x[0]: x[1:] for x in [l.strip().split('\t') for l in open(gos_data).readlines() if not l.startswith('#')]}
    
    sequences = {}
    scaffolds = {}
    
    rrnas = {}
    trnas = {}
    repetitiveels = {}
    
    annotations = {}
    geneontologys = {}
    pathways = {}
    domains = {}
    dblinks = {}
    annotation_gene_ontology = {}
    annotation_dblink = {}
    domain_dblink = {}
    annotation_pathways = {}
    transcripts = {}
    
    genes = {}
    mrnas = {}
    five_prime_utrs = {}
    exnos = {}
    introns = {}
    cds = {}
    three_prime_utrs = {}
    
    mrna_cds = {}
    mrna_exon = {}
    mrna_five_primeutr = {}
    mrna_intron = {}
    mrna_three_primeutr = {}
    
    
    ## parse scaffolds
    print('parse scaffolds: ' + fastaGenome)
    for seq in seq_dict_genome:
        scaffold = str(seq_dict_genome[seq].seq)
        length = str(len(scaffold))
        sequences[seq] = Record({
            'source': 'MaSuRCA',
            'start' : '1',
            'end': length,
            'len': length,
            'strand': '+',
            'frame': '.',
            'attributes': 'ID=' + seq + (';Type=scaffold;' if 'N' in scaffold.upper() else ';Type=contig;'),
            'sequence': scaffold,
            'type': 'NUCLEOTIDE' 
        })
        
        scaffolds[seq] = Record({
            'name': seq
        })
        
        sequences[seq].addFieldUpdate('sequence_id', scaffolds[seq]) 
    
    ## parse rRNA
    print('parse rRNA: ' + tsvRRNA)
    cont = {'8s': 1, '28s': 1, '18s': 1}
    for rrna in rrnaTSV:
        mtype = rrna[8].split('_')[0]
        name = 'rRNA' + rrna[8].split('_')[0] + str(cont[mtype])
        cont[mtype] += 1
        sequences[name] = Record({
            'source': rrna[1],
            'start' : rrna[3],
            'end': rrna[4],
            'len': str(int(rrna[4]) - int(rrna[3]) + 1),
            'strand': rrna[6],
            'frame': rrna[7],
            'attributes': 'ID=' + name + ';MoleculeType=' + mtype + ';',
            'sequence': sequences[rrna[0]].data['sequence'][int(rrna[3])-1:int(rrna[4])],
            'type': 'NUCLEOTIDE' 
        })
        
        rrnas[name] = Record({
            'name': name,
            'mulecule': 'M' + mtype.upper(),
            'note': None
        })
        
        sequences[name].addFieldUpdate('sequence_id', rrnas[name])
        scaffolds[rrna[0]].addFieldUpdate('scaffold_id', rrnas[name])
        
        
    ##parse tRNA
    print('parse tRNA: ' + tsvTRNA)
    cont = {}
    for trna in trnaTSV:
        anticodon = trna[3]
        if not anticodon in cont:
            cont[anticodon] = 1
        name = 'tRNA.%s%d' % (anticodon, cont[anticodon])
        cont[anticodon] += 1
        note = trna[7] if len(trna) > 7 else False
        sequences[name] = Record({
            'source': trna[1],
            'start' : trna[4],
            'end': trna[5],
            'len': str(int(trna[5]) - int(trna[4]) + 1),
            'strand': trna[6],
            'frame': '.',
            'attributes': 'ID=' + name + ';Type=' + trna[2] + ';' + (('Note=' + note + ';') if note and len(note) < 10 else ''),
            'sequence': sequences[trna[0]].data['sequence'][int(trna[4])-1:int(trna[5])],
            'type': 'NUCLEOTIDE' 
        })
        
        trnas[name] = Record({
            'name': name,
            'peptide': trna[2],
            'anticodon': anticodon,
            'asc_art':  note if note and len(note) > 10 else '',
            'note':  note if note and len(note) < 10 else ''
        })
    
        sequences[name].addFieldUpdate('sequence_id', trnas[name])
        scaffolds[trna[0]].addFieldUpdate('scaffold_id', trnas[name])

    ## parse repetitive els
    print('parse SSR: ' + tsvRepeatsMISA)
    cont = 1
    for repet in repeatsTSV:
        name = 'SSR' + str(cont)
        sequences[name] = Record({
            'source': 'misa',
            'start' : repet[5],
            'end': repet[6],
            'len': repet[4],
            'strand': '+',
            'frame': '.',
            'attributes': 'ID=' + name + ';Type=' + repet[2] + ';SSR=' + repet[3] + ';SSR_NUM=' + repet[1] + ';',
            'sequence': sequences[repet[0]].data['sequence'][int(repet[5])-1:int(repet[6])],
            'type': 'NUCLEOTIDE' 
        })
        
        repetitiveels[name] = Record({
            'name': name,
            'is_complete': 1,
            'is_chimera': 0,
            'repeat_class': 'SSR',
            'super_type': repet[2].upper().replace('*', '1'),
            'sub_type': repet[3],
            'note': '# SSR' + repet[1]
        })
        
        sequences[name].addFieldUpdate('sequence_id', repetitiveels[name])
        scaffolds[repet[0]].addFieldUpdate('scaffold_id', repetitiveels[name])
        cont += 1
    
    cont = {}
    print('parse Tranposons: ' + gffRepeats)
    for repet in repeatsGFF:
        tipo = repet[8].split('=')[1].split('_')[2]
        classe = tipo.split('-')[0]
        superTipo = '-'.join(tipo.split('-')[1:])
        if not classe in cont:
            cont[classe] = 1
        name = classe + str(cont[classe])
        cont[classe] += 1
        sequences[name] = Record({
            'source': 'REPET',
            'start' : repet[3],
            'end': repet[4],
            'len': str(int(repet[4]) - int(repet[3]) + 1),
            'strand': repet[6],
            'frame': repet[7],
            'attributes': repet[8],
            'sequence': sequences[repet[0]].data['sequence'][int(repet[3])-1:int(repet[4])],
            'type': 'NUCLEOTIDE' 
        })
        
        repetitiveels[name] = Record({
            'name': name,
            'is_complete': 1 if repet[2] == 'match' else 0,
            'is_chimera': max(tipo.count('-chim'), 1),
            'repeat_class': 'DNA' if classe.startswith('D') else ('RNA' if classe.startswith('R') else 'OTHER'),
            'super_type': classe.upper(),
            'sub_type': superTipo,
            'note': 'Score=' + repet[5]
        })
        
        sequences[name].addFieldUpdate('sequence_id', repetitiveels[name])
        scaffolds[repet[0]].addFieldUpdate('scaffold_id', repetitiveels[name])
 
        
    ## parse annotations
    print('parse annotations: ' + tsvAnnotation)
    for annotation in annotationTSV:
        protein = annotation[2] + '|' + annotation[1]
        protein_seq = str(seq_dict_proteins[protein].seq)
        annotations[protein] = Record({
            'name': annotation[3],
            'evalue': annotation[5],
            'coverage': annotation[6],
            'score': annotation[7],
            'reviwed': '1' if annotation[8] == 'R' else '0', 
            'gis': annotation[9],
            'uniprot': annotation[14] + ',' + annotation[15],
            'protein': annotation[12],
            'gene': annotation[13],
            'organism': annotation[16],
            'organisms': annotation[18],
            'repet': annotation[19],
            'families': annotation[20],
            'motif': annotation[21],
            'coments': annotation[22],
            'features': annotation[23],
            'ortologs': annotation[24],
            'kegg': annotation[25],
            'bigg': annotation[26],
            'ogs': annotation[27],
            'cog': annotation[28],
            'anotattion': annotation[29],
            'prize': annotation[31] if len(annotation) > 31 else None,
            'observation': None
        })
        
        for go in annotation[17].split(','):
            if len(go) < 4:
                continue
            name = protein + go
            if not go in geneontologys:
                geneontologys[go] = Record({
                    'name': go,
                    'description': gosTSV[go][0] if go in gosTSV else None,
                    'long_description': gosTSV[go][2] if go in gosTSV else None,
                    'aspect': gosTSV[go][1].upper() if go in gosTSV else None
                }) 
            annotation_gene_ontology[name] = Record({})
            annotations[protein].addFieldUpdate('annotation_id', annotation_gene_ontology[name])
            geneontologys[go].addFieldUpdate('go_id', annotation_gene_ontology[name])
        
        pws = ast.literal_eval(annotation[30])
        for db in pws:
            for link in pws[db]:
                if len(link) < 2: 
                    continue
                p = db + link
                name = protein + p
                if not p in pathways:
                    pathways[p] = Record({
                        'source': db,
                        'entry': link
                    }) 
                annotation_pathways[name] = Record({})
                annotations[protein].addFieldUpdate('annotation_id', annotation_pathways[name])
                pathways[p].addFieldUpdate('pathway_id', annotation_pathways[name])

        
        links = ast.literal_eval(annotation[10])
        for db in links:
            for l in links[db]:
                dbRec = db+l
                dbl = protein + dbRec
                if not dbRec in dblinks:
                    dblinks[dbRec] = Record({
                        'dbname': db,
                        'entry': l
                    })
                annotation_dblink[dbl] = Record({})
                dblinks[dbRec].addFieldUpdate('dblink_id', annotation_dblink[dbl])
                annotations[protein].addFieldUpdate('annotation_id', annotation_dblink[dbl])
                
        dms = ast.literal_eval(annotation[11])
        for db in dms:
            for pos in dms[db]:
                ini = int(pos[0])
                fim = int(pos[1])
                name = protein + '|' + db + '|' + str(pos)
                sequences[name] = Record({
                    'source': 'interproscan',
                    'start' : pos[0],
                    'end': pos[1],
                    'len':  str(fim - ini + 1),
                    'strand': '+',
                    'frame': '.',
                    'attributes': None,
                    'sequence': protein_seq[ini-1:fim],
                    'type': 'PEPT' 
                })
    
                domains[name] = Record({
                        'name': name,
                        'description': ' '.join(set([x[1] for x in dms[db][pos]])), 
                        'note': None
                    })
                
                annotations[protein].addFieldUpdate('annotation_id', domains[name])
                sequences[name].addFieldUpdate('sequence_id', domains[name])
                
                for l in dms[db][pos]:
                    dbname = 'InterPro' if l[0].startswith('IPR') else db
                    dbRec = dbname+l[0]
                    dbl = name + dbRec
                    if not dbRec in dblinks:
                        dblinks[dbRec] = Record({
                            'dbname': dbname,
                            'entry': l[0], 
                            'description': l[1]
                        })
                    elif len(l[1]) > 2 and ((not 'description' in dblinks[dbRec].data) or (len(dblinks[dbRec].data['description']) < 2)):
                        dblinks[dbRec].data['description'] = l[1]
                    domain_dblink[dbl] = Record({})
                    dblinks[dbRec].addFieldUpdate('dblink_id', domain_dblink[dbl])
                    domains[name].addFieldUpdate('domain_id', domain_dblink[dbl])
          
    ## parse transcripts
    print('parse transcripts: ' + fastaProteins)
    for protein in seq_dict_proteins:
        sequence = str(seq_dict_proteins[protein].seq)
        length = str(len(sequence))
        sequences[protein] = Record({
            'source': 'biopython',
            'start' : '1',
            'end': length,
            'len': length,
            'strand': '+',
            'frame': '.',
            'attributes': None,
            'sequence': sequence,
            'type': 'PROT' 
        })
        
        transcripts[protein] = Record({
            'name': protein
        })
        sequences[protein].addFieldUpdate('sequence_id', transcripts[protein])
        if protein in annotations:
            annotations[protein].addFieldUpdate('annotation_id', transcripts[protein])
         
    ## parse genome GFF
    print('parse genomic gff: ' + gffGenome)
    def reverse_complement(x):
        r = ""
        for i in range(len(x)-1, -1, -1):
            r += "A" if x[i] == "T" else ("T" if x[i] == "A" else ("C" if x[i] == "G"  else ("G" if x[i] == "C"  else "N")))
        return r
    cdsdosmrnas = {}
    for scaffold in set([x[0] for x in genomicGFF]):
        seq = sequences[scaffold].data['sequence']
        for feature in sorted([x for x in genomicGFF if x[0] == scaffold], key=lambda e: (1 if e[2] == 'gene' else(2 if e[2] == 'mRNA' else 3))):
            dic = {k.split('=')[0]: k.split('=')[1].split(',') for k in feature[8].split(';') if len(k) > 3}
            name = dic['ID'][0]
            parents = dic['Parent'] if 'Parent' in dic else []
            ini = int(feature[3])
            end = int(feature[4])
            sequence = seq[ini-1:end] if feature[6] == '+' else reverse_complement(seq[ini-1:end].upper())
            sequences[name] = Record({
                'source': feature[1],
                'start' : feature[3],
                'end': feature[4],
                'len': str(len(sequence)),
                'strand': feature[6],
                'frame': feature[7],
                'attributes': feature[8],
                'sequence': sequence,
                'type': 'NUCLEOTIDE' 
            })
            tip = feature[2]
            
            if tip == 'gene':
                genes[name] = Record({
                    'name': name,
                    'anotation': '',
                    'gis': '',
                    'note': None
                })
                sequences[name].addFieldUpdate('sequence_id', genes[name])
                scaffolds[scaffold].addFieldUpdate('scaffold_id', genes[name])
            elif tip == 'mRNA':
                mrnas[name] = Record({
                    'name': name,
                    'cdscomplete': None,
                    'note': None
                })
                protein =   name + '|' + parents[0]
                gene = genes[parents[0]]
                sequences[name].addFieldUpdate('sequence_id', mrnas[name])
                if protein in transcripts:
                    transcripts[protein].addFieldUpdate('transcript_id', mrnas[name])
                gene.addFieldUpdate('gene_id', mrnas[name])
                if protein in annotations:
                    gene.data['anotation'] += '|' + annotations[protein].data['gene']
                    gene.data['gis'] = ','.join(set([ x for x in (gene.data['gis'] + ',' + annotations[protein].data['gis']).split(',') if len(x) > 3]))
            elif tip == 'exon':
                exnos[name] = Record({
                    'name': name
                })
                sequences[name].addFieldUpdate('sequence_id', exnos[name])
                for p in parents:
                    n = name + p
                    mrna_exon[n] = Record({ })
                    mrnas[p].addFieldUpdate('mrna_id', mrna_exon[n])
                    exnos[name].addFieldUpdate('exon_id', mrna_exon[n])
            elif tip == 'intron':
                introns[name] = Record({
                    'name': name
                })
                sequences[name].addFieldUpdate('sequence_id', introns[name])
                for p in parents:
                    n = name + p
                    mrna_intron[n] = Record({})
                    mrnas[p].addFieldUpdate('mrna_id', mrna_intron[n])
                    introns[name].addFieldUpdate('intron_id', mrna_intron[n])
            elif tip == 'CDS':
                cds[name] = Record({
                    'name': name
                })
                sequences[name].addFieldUpdate('sequence_id', cds[name])
                for p in parents:
                    n = name + p
                    mrna_cds[n] = Record({ })
                    mrnas[p].addFieldUpdate('mrna_id', mrna_cds[n])
                    cds[name].addFieldUpdate('codingsequence_id', mrna_cds[n])
                    if not p in cdsdosmrnas:
                        cdsdosmrnas[p] = {}
                    cdsdosmrnas[p][(ini, end)] = sequence
            elif tip == 'five_prime_UTR':
                five_prime_utrs[name] = Record({
                    'name': name
                })
                sequences[name].addFieldUpdate('sequence_id', five_prime_utrs[name])
                for p in parents:
                    n = name + p
                    mrna_five_primeutr[n] = Record({})
                    mrnas[p].addFieldUpdate('mrna_id', mrna_five_primeutr[n])
                    five_prime_utrs[name].addFieldUpdate('fiveprimeutr_id', mrna_five_primeutr[n])
            elif tip == 'three_prime_UTR':
                three_prime_utrs[name] = Record({
                    'name': name
                })
                sequences[name].addFieldUpdate('sequence_id', three_prime_utrs[name])
                for p in parents:
                    n = name + p
                    mrna_three_primeutr[n] = Record({})
                    mrnas[p].addFieldUpdate('mrna_id', mrna_three_primeutr[n])
                    three_prime_utrs[name].addFieldUpdate('threeprimeutr_id', mrna_three_primeutr[n])
            else:
                print("ERROR type %s unknown" % tip)
                
    ##atualizar cds mrna
    print('parse genomic cds for mrnas')
    for m in cdsdosmrnas:
        mrnas[m].data['cdscomplete'] = "".join([
            cdsdosmrnas[m][c] for c in 
            sorted(cdsdosmrnas[m], key=lambda e: -e[1] if sequences[m].data['strand'] == '-' else e[0])
        ])
    
    
    return {
        'sequence': sequences,
        'scaffold': scaffolds,
        'rrna': rrnas, 
        'trna': trnas, 
        'repetitive_element': repetitiveels, 
        'annotation': annotations,
        'gene_ontology': geneontologys,
        'domain': domains, 
        'db_link': dblinks,
        'pathway': pathways,
        'annotation_pathway': annotation_pathways,
        'annotation_go': annotation_gene_ontology,
        'annotation_dblink': annotation_dblink,
        'domain_dblink': domain_dblink,
        'transcript': transcripts, 
        'gene': genes, 
        'mrna': mrnas, 
        'five_prime_utr': five_prime_utrs, 
        'exon': exnos, 
        'intron': introns, 
        'coding_sequence': cds, 
        'three_prime_utr': three_prime_utrs,
        'mrna_codingsequence': mrna_cds, 
        'mrna_exon': mrna_exon, 
        'mrna_fiveprimeutr': mrna_five_primeutr,
        'mrna_intron': mrna_intron,
        'mrna_threeprimeutr': mrna_three_primeutr
    }, ['sequence','scaffold',
        'rrna', 'trna','repetitive_element', 
        'annotation','gene_ontology','domain', 'db_link', 'pathway',
        'annotation_go', 'annotation_dblink', 'domain_dblink', 'annotation_pathway',
        'transcript', 'gene', 'mrna','five_prime_utr', 'exon', 'intron','coding_sequence','three_prime_utr',
        'mrna_codingsequence', 'mrna_exon', 'mrna_fiveprimeutr','mrna_intron','mrna_threeprimeutr']



def persistDB(database, data, tables):
    def populate(db, table, values):
        print('populate table ' + table)
        return insert(db, 
               table, 
               values)
    cont = 0
    for table in tables:
        cont += populate(database, table, list(data[table].values()))
    print('%d records loaded' % cont)
    return data
    


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("use: ./loaddb.py file_config.txt")
    else:
        try:
            with open(sys.argv[1]) as f:
                config = {x[0].strip():x[1].strip() for x in [l.strip().split(':') for l in f.readlines() if not l.startswith('#')]}
                print("import data (%s) ...\n\n" % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
                r = loadDB(
                    fastaGenome = config['fastaGenome'],
                    fastaProteins = config['fastaProteins'], 
                    gffGenome = config['gffGenome'], 
                    tsvRRNA = config['tsvRRNA'], 
                    tsvTRNA = config['tsvTRNA'], 
                    tsvAnnotation = config['tsvAnnotation'], 
                    tsvRepeatsMISA = config['tsvRepeatsMISA'], 
                    gffRepeats = config['gffRepeats'], 
                    gos_data = config['gos_data']
                )
                print("persist data (%s) ...\n\n" % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
                d = persistDB(
                    database = connect(config['host'], config['user'], config['password'], config['database']), 
                    data = r[0], 
                    tables = r[1]
                )
            print("done (%s) ...\n\n" % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
            print('\nterminado com sucesso\nby mikeias.net')

        except FileNotFoundError:
            with open(sys.argv[1], 'w') as f:
                f.write('\n'.join([
                    '##',
                    '##',
                    'host: bioserver1', 
                    'user: jhipster', 
                    'password: jhipster',
                    'database: guavadb',
                    '##',
                    '##',
                    'fastaGenome: /home/cluster/shared/data/guava_final/guava.fa',
                    'fastaProteins: /home/cluster/shared/data/guava_final/guava.protein.faa', 
                    'gffGenome: /home/cluster/shared/data/guava_final/guava.gff3', 
                    'tsvRRNA: /home/cluster/shared/relatorio_mestrado/julho/guavadb/gffs/out.gff', 
                    'tsvTRNA: /home/cluster/shared/relatorio_mestrado/julho/trna/arag_trnase.tsv', 
                    'tsvAnnotation: /home/cluster/shared/relatorio_mestrado/julho/guavadb/annotation_table.tsv', 
                    'tsvRepeatsMISA: /home/cluster/shared/relatorio_mestrado/julho/guavadb/gffs/genome.fa.misa',
                    'gffRepeats: /home/cluster/shared/relatorio_mestrado/julho/guavadb/gffs/guava_repet.gff3',
                    'gos_data: /home/cluster/shared/relatorio_mestrado/julho/guavadb/gos_mapped.tsv'
                ]) + '\n')
            print('File ' + sys.argv[1] + ' was created for define data to populate db')



