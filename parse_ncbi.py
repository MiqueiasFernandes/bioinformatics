#!/usr/bin/python3

import argparse

def clean_ptnas(fasta, fout):
    valid = []
    tmp = ''
    for l in open(fasta).readlines():
        l = l.strip()
        if l.startswith('>'):
            valid.append([l, False])
            tmp = ''
        else:
            tmp += l
            valid[-1][1] = len(tmp) > 1 and tmp[:-1].isalpha()
            
    export = set([x.split(' ')[0][1:] for x, v in valid if v])
    print('Qtd total', len(valid))
    print('Qtd valid', len(export))
    print('Qtd invalid', len([x for x, v in valid if not v]))
    
    with open(fout, 'w') as fd:
        w = False
        for lo in open(fasta).readlines():
            l = lo.strip()
            
            if l.startswith('>'):
                if '[protein_id=' in lo:
                    w = l.split(' ')[0][1:] in export
                    pid = lo.split('[protein_id=')[1].split(']')[0].strip()
                    lo = '>' + pid + ' ' + l[1:] + '\n'
                else:
                    w = False
            elif not l[-1].isalpha():
                lo = lo[:-1]
            
            if w:
                fd.write(lo)
                    

def geneFromGTF(gtf, prefix, genoma):
    genes = []
    mrna2gene, mrna2ptna = {}, {}
    
    get = lambda e, f: e.split(f)[1].split(';')[0].strip().replace('"', '')
    
    for l in open(gtf).readlines():
        if l.count('\t') < 8:
            continue
        seq, _, ft, ini, fim, _, fita, _, anot = l.split('\t')
        if ft == 'gene':
            genes.append([get(anot, 'gene_id '), seq, int(ini)-1, int(fim), fita == '+'])
        elif ft == 'transcript':
            mrna2gene[get(anot, 'transcript_id ')] = get(anot, 'gene_id ')
        elif ft == 'CDS' and 'protein_id ' in anot:
            mrna2ptna[get(anot, 'transcript_id ')] = get(anot, 'protein_id ')
    
    ## persist table gene2mrna2ptna
    with open(prefix+'.txt', 'w') as fd:
        for m, g in mrna2gene.items():
            if m in mrna2ptna:
                fd.write(f'{mrna2gene[m]},{m},{mrna2ptna[m]}\n')
    
    print(len(genes))
    ## persist gene seqs
    def strand(seq, fita):
        if fita:
            return seq
        
        complement = seq.replace('a', '1').replace('A', '5') 
        complement = complement.replace('t', '2').replace('T', '6') 
        complement = complement.replace('c', '3').replace('C', '7') 
        complement = complement.replace('g', '4').replace('G', '8')
        
        complement = complement.replace('1', 't').replace('5', 'T')
        complement = complement.replace('2', 'a').replace('6', 'A')
        complement = complement.replace('3', 'g').replace('7', 'G')
        complement = complement.replace('4', 'c').replace('8', 'C')
        
        return complement[::-1]
                       
    limits = {}
    for _, c, _, f, _ in genes:
        limits[c] = max(f,limits[c]) if c in limits else f
    wrs = []
    with open(prefix+'.fna', 'w') as fd:
        with open(genoma) as fr:
            id, seq = None, ''
            for l in fr.readlines():
                lo = l.strip()
                if lo[0] == '>':
                    id = lo[1:].split(' ')[0]
                    skip = not id in limits
                    seq = ''
                elif not skip:
                    seq += lo
                if skip:
                    continue
                if len(seq) >= limits[id]:
                    for g, c, s, e, f in genes:
                        if c == id and not g in wrs:
                            fd.write(f'>{g}\n{strand(seq[s:e], fita)}\n')
                            wrs.append(g)
                    skip = True

parser = argparse.ArgumentParser(description="""
Script to extract:
    1. cleaned proteome
    2. gene sequences
    3. table GENE > TRANSCRIPT > PROTEIN
using ncbi files.
""")

parser.add_argument('genome')
parser.add_argument('proteome')
parser.add_argument('gtf')
args = parser.parse_args()

geneFromGTF(args.gtf, 'clean', args.genome)
clean_ptnas(args.proteome, 'clean_ptnas.faa')                   
