#!/usr/bin/env python
# coding: utf-8

# Copyright(c) 2020 - Miquéias Fernandes <bio@mikeias.net>

import sys
import csv
import os


#################################################
usage = "python3 annotGenes.py file.gff3 file_diamond_SP.tsv file_diamond_NR.tsv file_eggnog.tsv file_interpro.tsv '[pref keyword1];keyowordNR2; ...'"
#################################################

try: _, file_gff,file_diamond_SP,file_diamond_NR,file_eggnog,file_interpro,prefs=sys.argv
except: sys.exit("Correct usage is:\n"+usage)

out_file = 'genes_anotados.csv'
UNIPROT_API_PARALLEL=100


def bestName(names, pref=['['], filterPs = ['[', ']', 'protein', '_', 'PREDIC', 'PUTAT', 'putat', '.', 'hypothetic']):
    names = list(names)
    if len(names) < 1:
        return ''
    if len(names) == 1:
        return names[0]
    txt = ' '.join(names)
    pScore = {p: txt.count(p) for p in set(txt.split()) if not p.isdigit() and not any([fp in p for fp in filterPs])}
    def pontuar(name):
        a = min([(pref.index(p) if p in name else len(pref)) for p in pref]) ### best pref < melhor
        b = -sum([pScore[p] if p in pScore else 0 for p in name.split()])### melhor concenso > melhor
        c = -len(name) ### maior texto
        return (a, b, c)
    return sorted(names, key=pontuar)[0]


def parseDiamondNR(file, mrna2gene, pref=['[']):
    ## Needs Seq Title in Last Column!
    dt = [ x for x in list(csv.reader(open(file), delimiter='\t')) if len(x) > 2]
    res = {}
    for hit in dt:
        g = mrna2gene[hit[0]]
        if g in res:
            res[g][0].add(hit[1])
            res[g][1].add(hit[-1])
        else:
            res[g] = [set([hit[1]]), set([hit[-1]])]
    for k, v in res.items():
        res[k][1] = bestName(res[k][1], pref)
    return res


def parseDiamond(file, mrna2gene):
    dt = [ x for x in list(csv.reader(open(file), delimiter='\t')) if len(x) > 2]
    res = {}
    for hit in dt:
        g = mrna2gene[hit[0]]
        if g in res:
            res[g].add(hit[1])
        else:
            res[g] = set([hit[1]])
    return res


def parseInterpro(file, mrna2gene):
    dt = [ x for x in list(csv.reader(open(file), delimiter='\t')) if len(x) > 2 and not x[0].startswith('#')]
    res = {}
    #0  Protein Accession (e.g. P51587)
    #1  Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    #2  Sequence Length (e.g. 3418)
    #3  Analysis (e.g. Pfam / PRINTS / Gene3D)
    #4  Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    #5  Signature Description (e.g. BRCA2 repeat profile)
    #6  Start location
    #7  Stop location
    #8  Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
    #9  Status - is the status of the match (T: true)
    #10 Date - is the date of the run
    #11 (InterProScan annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
    #12 (InterProScan annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
    #13 (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
    #14 (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
    def f(x):
        interesse = [3, 4, 5]
        return [[x[i] for i in interesse] + [x[13] if len(x) > 13 else '', x[14] if len(x) > 14 else ''], 
                ['interpro', x[11], x[12], '', ''] if len(x) > 12 else ['', '', '', '', '']]
    
    for hit in dt:
        g = mrna2gene[hit[0]]
        if g in res:
            res[g].extend(f(hit))
        else:
            res[g] = f(hit)
    return res

def parseEggNog(file, mrna2gene):
    dt = [ x for x in list(csv.reader(open(file), delimiter='\t')) if len(x) > 2 and not x[0].startswith('#')]
    res = {}
    #0  query_name
    #1  seed_eggNOG_ortholog
    #2  seed_ortholog_evalue
    #3  seed_ortholog_score
    #4  best_tax_level
    #5  Preferred_name
    #6  GOs
    #7  EC
    #8  KEGG_ko
    #9  KEGG_Pathway
    #10 KEGG_Module
    #11 KEGG_Reaction
    #12 KEGG_rclass
    #13 BRITE
    #14 KEGG_TC
    #15 CAZy
    #16 BiGG_Reaction
    def f(x):
        interesse = [5, 6, 8, 9, 16]
        return [x[i] for i in interesse]
    
    for hit in dt:
        g = mrna2gene[hit[0]]
        if g in res:
            res[g].append(f(hit))
        else:
            res[g] = [f(hit)]
    return res


def getUniprotReg(ids, limi=100, dbs_interesse=['GO','KEGG','PFAM','Pfam','Reactome','SMART','STRING','UniPathway']):
    baixados = {}
    total = len(ids)
    interesse = ['ID ', 'DE ', 'OS ', 'OC ', 'OX ', 'DR ', 'CC   -!- ']
    
    def parseReg(rows):
    
        def getByPrefix(p):
            l = len(p)
            return [ln[l:].strip() for ln in rows if ln.startswith(p)]
        tID = getByPrefix('ID ')
        idR = tID[0] if len(tID) > 0 else ''
        _ent = idR.split(" ")[0]
        _isRev = ' Reviewed;' in idR
        os =  getByPrefix('OS ')
        _os = os[0] if len(os) > 0 else ''
        ox =  getByPrefix('OX ')
        _ox = (ox[0].split("NCBI_TaxID=")[1].split(";")[0] if 'NCBI_TaxID=' in ox[0] else None) if len(ox) > 0 else None

        _oc = '; '.join(getByPrefix('OC ')).replace(";;", ';')

        names = getByPrefix('DE ')
        _name = (names[0].split('RecName: Full=')[1] if 'RecName: Full=' in names[0] else names[0]) if len(names) > 0 else ''

        _dbs = [(y[0], y[1], y[2]) for y in [x.split(';') for x in getByPrefix('DR ')] if y[0] in dbs_interesse] + [('Uniprot', _ent, 'Reviwed' if _isRev else 'Unreviwed')]
        _coments = '|'.join(getByPrefix('CC   -!- '))

        return _name, _os, _oc, _ox if not _ox is None else '', _dbs, _coments
    
    
    def fParse(lines):
        return parseReg([x for x in lines if any([x.startswith(y) for y in interesse])])
    
    while len(ids) > 0:
        atual = ids[:limi]
        fs = [x+'.txt' for x in atual]
        os.system(' '.join(['wget https://www.uniprot.org/uniprot/%s.txt 1>/dev/null 2>/dev/null &' % k for k in atual]) + ' wait')
        ids = ids[limi:]
        files = set(os.listdir()).intersection(fs)
        baixados.update({b.replace('.txt', ''): fParse([x.strip() for x in open(b).readlines()]) for b in files})
        for f in files:
            os.remove(f)
        print('faltam: %d...' % len(ids))
    print('%d de %d obtidos' % (len(baixados), total))
    return baixados


def mesclar(mult, king='Viridiplantae'):
    if len([x for x in mult if king in x[2]]) > 0:
        mult = [x for x in mult if king in x[2]]
    ret = list(mult[0])
    if len(mult) > 1:
        ret[0] = max([x[0] for x in mult])
        idsj = [x[1] for x in ret[4]]
        for o in mult[1:]:
            for dbE in o[4]:
                if not dbE[1] in idsj:
                    ret[4].append(dbE)
            for com in o[5]:
                if not com in ret[5]:
                    ret[5] + '|' + com
    return tuple(ret)


def addInterpro(gene_mesclado, ipro_data):
    for gene, datas in ipro_data.items():
        for data in datas:
            dbKs = [(data[0], data[1], data[2]), ('GO', data[3], '')]
            ## fix pathways
            if gene in gene_mesclado:
                
                try:
                    ks = [x[1] for x in gene_mesclado[gene][4]]
                except:
                    print(gene_mesclado[gene])
                for dbK in dbKs:
                    if not dbK in ks:
                        gene_mesclado[gene][4].extend([dbK])
            else:
                gene_mesclado[gene] = ('','','','',dbKs,'')
    return gene_mesclado


def addEggnog(gene_mesclado, eggnog_data):
    for gene, datas in eggnog_data.items():
        if gene in gene_mesclado:
            for data in datas:
                if len(data[0]) > 1:
                    gene_mesclado[gene] = list(gene_mesclado[gene])
                    gene_mesclado[gene][0] += '|' + data[0]
                if 'GO' in data[1]:
                    entryes = [x[1] for x in gene_mesclado[gene][4]]
                    for go in data[1].split(","):
                        if not go in entryes:
                            gene_mesclado[gene][4].append(('GO', go, ''))
                    if ':' in data[2]:
                        ko = data[2].split(":")[1]
                        if not ko in entryes:
                            gene_mesclado[gene][4].append(('KO', ko, ''))
                if len(data[3]) > 1:
                    gene_mesclado[gene][4].append(('KEGG_Pathway', data[3], ''))
                if len(data[4]) > 1:
                    gene_mesclado[gene][4].append(('BiGG_Reaction', data[4], ''))
        else:
            dbks = []
            if 'GO' in data[1]:
                dbks = [('GO', x, '') for x in data[1].split(",")]
            if len(data[3]) > 1:
                dbks.append(('KEGG_Pathway', data[3], ''))
            if len(data[4]) > 1:
                dbks.append(('BiGG_Reaction', data[4], ''))
            gene_mesclado[gene] = (data[0],'','','',dbks,'')
    return gene_mesclado


def regs2table(data, nrData):
    ret = []
    for gene, records in data.items():
        txt = [x for x in (records[0] + ' ' + records[5]).replace('|', ' ').replace('SubName: Full=', '').replace('Full=', '').split(' ') if not x in [
            'MISCELLANEOUS:',
            'uncharacterized',
            'protein.',
            'PROTEIN.',
            'Putative',
            'PUTATIVE',
            'PROTEIN-LIKE.',
            'hypothetical',
            'Hypothetical',
            'Uncharacterized'
        ]]
        
        old_name = records[0]
        nrN = ''
        if gene in nrData:
            nm = nrN = (nrData[gene][1] + ' ')
            if len(nm) > 1: 
                old_name  = nm + ' | ' + old_name
        nome = nrN + '[%s] %s' % (sorted(txt, key=lambda e: (-len(e), txt.count(e)))[0], ' :: '.join(set(old_name.split("|"))))
        
        ks = ','.join(set([x[0]+':'+x[1] for x in records[4] if len(x[1]) > 3])).replace('|', ',').split(",")
        gos = ','.join(set([x.replace(' ', '').replace("GO:GO:", 'GO:') for x in ks if x.startswith('GO:') and len(x) > 6]))
        dbs = ','.join(set([x for x in ks if not x.startswith('GO:')]))
        ret.append((
            gene.replace('\t', ' ').replace(';', ','), 
            nome.replace('\t', ' ').replace(';', ','), 
            records[1].replace('\t', ' ').replace(';', ','), records[2].replace('\t', ' ').replace(';', ','), 
            gos.replace('\t', ' ').replace(';', ','),
            ','.join(nrData[gene][0] if gene in nrData else []),
            dbs.replace('\t', ' ').replace(';', ','),
            (records[5] + ' | '.join(set([x[0]+':'+x[2] for x in records[4] if len(x[2]) > 3]))).replace('\t', ' ').replace(';', ',')
        ))
    for gene, records in nrData.items():
        if not gene in data:
            ret.append((
                gene.replace('\t', ' ').replace(';', ','), 
                records[1], 
                '-', 
                '-',
                ','.join(nrData[gene][0]),
                dbs.replace('\t', ' ').replace(';', ','),
                '-'
            )) 
    return ret



print('[1/7] parseando gff3 ...')
mrna2gene = {
    x.split("ID=")[1].split(";")[0]: x.split("Parent=")[1].split(";")[0] for x in
    [m[8] for m in list(csv.reader(open(file_gff), delimiter='\t')) if len(m) > 2 and m[2] == 'mRNA']
}

print('[2/7] parseando diamond NR ...')
diamond_NR = parseDiamondNR(file_diamond_NR, mrna2gene, prefs.split(';'))
genes_NR = set(diamond_NR)

print('[3/7] parseando diamond swiss prot ...')
diamond_SP = parseDiamond(file_diamond_SP, mrna2gene)
genes_SP = set(diamond_SP)
swps = list(set(','.join([','.join(v) for v in diamond_SP.values()]).split(',')))
regs = getUniprotReg(swps, limi=UNIPROT_API_PARALLEL)

print('[4/7] parseando EggNOG ...')
eggnog = parseEggNog(file_eggnog, mrna2gene)
genes_EN = set(eggnog)

print('[5/7] parseando InterproScan ...')
interpro = parseInterpro(file_interpro, mrna2gene)
genes_IP = set(interpro)

print('[6/7] mesclando dados ...')
g_mesc = {k: mesclar([regs[x] for x in v if x in regs]) for k, v in diamond_SP.items()} 
g_mesc_int = addInterpro(g_mesc, interpro)
g_mesc_int_egg = addEggnog(g_mesc_int, eggnog)
parsed_mescled_table = regs2table(g_mesc_int_egg, diamond_NR)

print('[7/7] persistindo ...')
open(out_file, 'w').writelines(['\t'.join(x) + '\n' for x in parsed_mescled_table])
open(out_file + '_agbase.txt', 'w').writelines(['\n'.join([(x[0] +'\t'+ y) for y in x[4].split(',')]) + '\n' for x in parsed_mescled_table if 'GO:' in x[4]])

all_G = genes_NR.union(genes_SP).union(genes_EN).union(genes_IP)
genes_N, genes_S, genes_E, genes_I = len(genes_NR), len(genes_SP), len(genes_EN), len(genes_IP)
total = len(all_G)

print('%d (%.1f%%) (%d so nele) genes anotados no NR' % (genes_N, genes_N/total*100, len(genes_NR.difference((genes_SP).union(genes_EN).union(genes_IP)))))
print('%d (%.1f%%) (%d so nele) genes anotados no SWPROT' % (genes_S, genes_S/total*100, len(genes_SP.difference(genes_NR.union(genes_EN).union(genes_IP)))))
print('%d (%.1f%%) (%d so nele) genes anotados no Eggnog' % (genes_E, genes_E/total*100, len(genes_EN.difference(genes_NR.union(genes_SP).union(genes_IP)))))
print('%d (%.1f%%) (%d so nele) genes anotados no InterproScan' % (genes_I, genes_I/total*100, len(genes_IP.difference(genes_NR.union(genes_SP).union(genes_EN)))))
print('TOTAL: %d genes anotados' % total)

print('\nall finished. sotored at: ' + out_file)

