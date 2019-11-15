#!/usr/bin/python3

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# gff2html.py - a script for visualization of gff features


#################################################

import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#################################################
usage = "dart4snp.py report_single_row.csv genome.fasta genome.introns.gff3 gene_ids"
#################################################

try: 
    scr,file_report,genome_fasta,introns_gff3,gene_ids=sys.argv
except: sys.exit("Correct usage is:\n"+usage)


## configurações
colunas_genotipos_a_partir = 31
linhas_a_partir = 7
threads = 4

print('[1/6] parsear arquivo report')
report = [l.strip().split(',') for l in open(file_report).readlines()]
linhas = report[linhas_a_partir:]
print('linhas uteis do arquivo report: %d' % len(linhas))

linhas_filtradas = [(l[0], l[3], l[colunas_genotipos_a_partir:]) for l in linhas if len(set(l[colunas_genotipos_a_partir:])) > 1]
print('alelos polimorficos %.1f%% (%d)' % (len(linhas_filtradas)*100/len(linhas), len(linhas_filtradas)))

print('\n\n[2/6] construir query')
with open('query.fa', 'w') as o:
    o.write('\n'.join('>%s\n%s' % (x[0], x[1]) for x in linhas_filtradas))

print('\n\n[3/6] construir subject')
genes = [l.strip() for l in open(gene_ids).readlines()]
gff3 = [l.strip().split('\t') for l in open(introns_gff3).readlines() if not l.startswith('#') and l.count('\t') == 8]
genoma = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

sequencias = []

for gene in set(genes):
    entry = [g for g in gff3 if g[2] == 'gene' and g[8].split('ID=')[1].split(';')[0] == gene][0]
    seq = genoma[entry[0]].seq[int(entry[3])-1:int(entry[4])]
    sequencias.append(SeqRecord(seq if entry[6] == '+' else seq.reverse_complement(), id=gene, description=""))

print('%d genes salvos em database.fasta' % SeqIO.write(sequencias, 'database.fasta', 'fasta'))

os.mkdir('blast')
os.link('database.fasta', 'blast/database.fasta')
os.system('makeblastdb -in blast/database.fasta -dbtype nucl')


print('\n\n[4/6] alinhar sequencias de SNP nos genes')

if os.system('blastn -db blast/database.fasta -query query.fa -evalue 1e-5 -out blast.out.csv -num_threads %d \
          -outfmt "10 qseqid sseqid qstart qend sstart send evalue qseq sseq"' % threads) != 0:
    raise Exception('Houve um erro ao realizar o alinhamento!')

alinhamento = [l.strip().split(',') for l in open('blast.out.csv').readlines()]
query2gene = {q: [s for s in alinhamento if s[0] == q] for q in set([x[0] for x in alinhamento])}
print('%d SNP alinhados' % len(query2gene))
print('%d genes alinhados' % len(set([x[1] for x in alinhamento])))


print('\n\n[5/6] localizar posição exata da ocorrencia do SNP')
cont = 1
res = []
for snp, genes in query2gene.items():
    pos_snp = int(snp.split('-')[1].split(':')[0])
    expect = ''.join(sorted(snp.split(':')[-1].split('>')))
    for gene in genes:
        q_a = int(gene[2]) -1
        q_b = int(gene[3]) -1
        snp_start = min(q_a, q_b)
        snp_end = max(q_a, q_b)
        if pos_snp < snp_start or pos_snp > snp_end:
            print('ERRO:', snp, gene[1], ' snp fora do alinhamento: em %d snp[%d-%d]' % (pos_snp, snp_start, snp_end))
        else:
            snp_start = int(gene[2]) - 1
            strand_snp = 1 if snp_start < int(gene[3]) else -1
            gene_start = int(gene[4]) - 1
            strand_gene = 1 if gene_start < int(gene[5]) else -1
            p_s = snp_start
            p_g = gene_start
            for i in range(len(gene[7])):
                if p_s == pos_snp:
                    break
                p_s += strand_snp if gene[7][i] != '-' else 0
                p_g += strand_gene if gene[8][i] != '-' else 0
            obtido = ''.join(sorted(set([gene[7][i], gene[8][i]])))
            if obtido in expect:
                g = [g for g in gff3 if g[2] == 'gene' and g[8].split('ID=')[1].split(';')[0] == gene[1]][0]
                pos_no_genoma = (int(g[4]) if g[6] == '-' else int(g[3])) + ((-1 if g[6] == '-' else 1) * p_g)
                ocorreu = False
                for mrna in [g for g in gff3 if g[2] == 'mRNA' and g[8].split('Parent=')[1].split(';')[0] == gene[1]]:
                    m = mrna[8].split('ID=')[1].split(';')[0]
                    cds = [(int(g[3]),int(g[4])) for g in gff3 if g[2] == 'CDS' and m in g[8].split('Parent=')[1].split(';')[0].split(',')]
                    ocorre_em = [x for x in cds if x[0] <= pos_no_genoma and x[1] >= pos_no_genoma]
                    if len(ocorre_em) > 0:
                        ocorreu = True
                        genoma_seq = [x for x in str(genoma[g[0]].seq)]
                        nuc_original = genoma_seq[pos_no_genoma-1]
                        genoma_seq[pos_no_genoma-1] = 'X'
                        gs = Seq(''.join(genoma_seq), generic_dna)
                        tmp = ''.join([str(gs[c[0]-1:c[1]] if g[6] == '+' else gs[c[0]-1:c[1]].reverse_complement()) for c in sorted(cds, key=lambda e: (e[0] if g[6] == '+' else -e[1]))])
                        trinca = [x for x in [tmp[i:i + 3] for i in range(0, len(tmp), 3)] if 'X' in x][0]
                        p1 = str(Seq(trinca.replace('X', expect[0]), generic_dna).translate())
                        p2 = str(Seq(trinca.replace('X', expect[1]), generic_dna).translate())
                        print(cont, snp, gene[1], m, gene[6], p1 + p2, 'SINONIMA' if p1 == p2 else 'MUTANTE')
                        res.append([str(cont), snp, gene[1], m, gene[6], p1 + p2, 'SINONIMA' if p1 == p2 else 'MUTANTE'])
                        cont += 1
                    else:
                        e = [g[2] for g in gff3 if g[2] != 'gene' and m in g[8].split('Parent=')[1].split(';')[0].split(',') and int(g[3]) <= pos_no_genoma and int(g[4]) >= pos_no_genoma]
                        if len(e) > 0:
                            ocorreu = True
                            print(cont, snp, gene[1], m, gene[6], e[0][0].upper(), 'SINONIMA')
                            res.append([str(cont), snp, gene[1], m, gene[6], e[0][0].upper(), 'SINONIMA'])
                            cont += 1
                if not ocorreu:
                    print('ERRO:', snp, gene[1], gene[6], pos_no_genoma, 'NAO OCORREU EM MRNA')
            else:
                print('ERRO:', i, p_s, p_g, expect, obtido, 'esperado não foi obtido')
            
print('\n\n[6/6] salvar resultados')
with open('resultados.csv', 'w') as o:
      o.write('\n'.join([';'.join(x) for x in res]) + '\n')

print('by mikeias.net')
