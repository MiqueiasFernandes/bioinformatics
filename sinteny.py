#!/usr/bin/env python3


import os, threading, time, sys, re
from datetime import datetime


from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord


class Sequence:
    def __init__(self, rec):
        self.rec = rec
        self.mrnas = {}
        self.density = []

class ImportFasta(threading.Thread):
    
    def setData(self, prefix, seq_dict, sequences, gene_density=False):
        self.prefix = prefix
        self.seq_dict = seq_dict
        self.gene_density = gene_density
        self.sequences = sequences
        self.sequences_out = []
        
    def run(self):
        name = self.getName() + "_" + self.prefix
        
        fasta =  name  + ".fa"
        gff = name + ".gff"
        
        SeqIO.write([self.seq_dict[s] for s in self.sequences], fasta, "fasta")
        for s in self.sequences:
            os.system("grep -P '^%s\t.+' %s >> %s" % (s, self.prefix + '.gff3', gff))
        
        seq_dict = SeqIO.to_dict(SeqIO.parse(open(fasta), "fasta"))
        limit_info = dict(gff_type = ['gene', 'mRNA', 'CDS'], gff_id = self.sequences)
        with open(gff) as in_handle:
            n = 0
            c = 0
            ids = []
            for rec in GFF.parse(in_handle, limit_info=limit_info, base_dict=self.seq_dict):
                if not rec.id in self.sequences:
                    continue
                sequence = Sequence(rec)
                genes = []
                for gene in [g for g in rec.features if g.type == 'gene']:
                    genes.append((int(gene.location.start),int(gene.location.end)))
                    for mrna in [m for m in gene.sub_features if m.type == 'mRNA' and not m.id in ids]:
                        c += 1
                        ids.append(mrna.id)
                        sequence.mrnas[mrna.id] = [
                            Seq(''.join([str(cds.extract(rec).seq) for cds in sorted(mrna.sub_features, key=lambda e: e.location.start, reverse=mrna.strand != 1) if cds.type == 'CDS']), generic_dna), 
                            mrna.location ]
                
                if self.gene_density:
                    dic = {m: len([g for g in genes if g[0] <= m and g[1] >= m]) for m in set([len(rec)] + [g[0] for g in genes] + [g[1] for g in genes] + [1+genes[i][1] for i in range(len(genes)-1)] + [genes[i][0]-1 for i in range(1, len(genes))])}
                    ini = 1
                    den = 0 if min(dic) > 1 else dic[min(dic)]
                    out = []
                    for p in sorted(set(dic)):
                        if dic[p] != den:
                            sequence.density.append([ini, p-1, den])
                            den = dic[p]
                            ini = p
                    if ini < len(rec):
                        sequence.density.append([ini, len(rec), den])
                self.sequences_out.append(sequence)
        os.system("rm " + name + "*")


def parseParallel(prefix, density=False, n_threads=24):
    with open(prefix + '.fa') as in_seq_handle:
        seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        
        by_s = sorted(seq_dict, key=lambda e: -len(seq_dict[e]))
        m = int(len(by_s)/2)
        seq_maiores = by_s[:m]
        seq_maiores.reverse()
        seq_menores = by_s[m:]
        parts = {x: [] for x in range(n_threads)}
        p = 0
        while len(seq_maiores) > 0 or len(seq_menores) > 0:
            if len(seq_maiores) > 0:
                parts[p%n_threads].append(seq_maiores.pop())
            if len(seq_menores) > 0:    
                parts[p%n_threads].append(seq_menores.pop())
            p += 1
        

        threads = []
        for k in range(n_threads):
            impF = ImportFasta(name = "Thread-{}".format(k+1))
            impF.setData(prefix, seq_dict, parts[k], gene_density=density)
            threads.append(impF)
            impF.start()

        c = 0
        runn = n_threads
        while runn > 0:
            runn = len([t for t in threads if t.isAlive()])
            seqs = sum([len(t.sequences_out) for t in threads])
            print("[%d] %d threads running (%.1f%%)..." % (c, runn, (seqs/len(seq_dict)*100)))
            time.sleep(60)
            c += 1
    return threads


def export(data, color):
    prefix = data[0].prefix
    cont_seq = 1
    cont_mrna = 1
    karyotype = open("%s.kar" % prefix, "w")
    bed = open("%s.bed" % prefix, "w")
    chrMAP = open("%s.map" % prefix, "w")
    seqs_store = []
    den_store = []
    for t in data:
        for seq in t.sequences_out:
            idS = "CHR%s%d" % (prefix, cont_seq)
            karyotype.write("chr\t-\t%s\t%s\t0\t%d\t%s\n" % (idS, seq.rec.id, len(seq.rec), color))
            chrMAP.write("%s\t%s\t\t\n" % (seq.rec.id, idS))
            for m in seq.mrnas:
                mrna = seq.mrnas[m]
                idmrna = 'MRNA%s%d' % (prefix, cont_mrna)
                bed.write("%s\t%s\n" % (idS, '\t'.join(str(mrna[1])[1:-1].replace('](', ':%s:0:' % idmrna).split(':'))))
                seqs_store.append(SeqRecord(mrna[0], id=idmrna, name="", description=""))
                chrMAP.write("%s\t%s\t%s\t%s\n" % (seq.rec.id, idS, m, idmrna))
                cont_mrna += 1
            for d in [x for x in seq.density if x[2] > 0]:
                den_store.append("%s\t%s\t%s\t%s" % (idS, d[0], d[1], d[2]))
            cont_seq += 1
    SeqIO.write(seqs_store, "%s.cds" % prefix, "fasta")
    karyotype.close()
    bed.close()
    chrMAP.close()
    if len(den_store) > 0:
        geneDensity = open("%s.gene.density" % prefix, "w")
        geneDensity.write('\n'.join(den_store) + '\n')
        geneDensity.close()


def importKaryotypes(base, prefixs):
    karyotypes = {}
    maps = {}
    for p in prefixs:
        if (not (os.path.isfile('./' + p + '.bed') and os.path.isfile('./' + p + '.cds') and os.path.isfile('./' + p + '.kar'))):
            export(parseParallel(p, base != p), ['red', 'green', 'grey', 'black', 'white'][min(4, prefixs.index(p))])
        with open(p + '.kar') as kar:
            karyotypes[p] = {k.split('\t')[2]: k.strip() for k in kar.readlines()}
        with open(p + '.map') as map:
            maps[p] = {k.split('\t')[1]: k.split('\t')[0] for k in map.readlines()}
        print("%s importado ..." % p)
    return karyotypes, maps



def runSinteny(base, sp2, minspan=1, cscore=7, rewrite=False):
    sp1=base
    if rewrite:
        os.remove("%s.%s.last.filtered" % (sp1, sp2))
        print('rewrite links ...')
    cmd1 = 'lastdb %s %s.cds -P 24' % (sp2, sp2)
    cmd2 = 'python -m jcvi.compara.catalog ortholog %s %s --cscore=.%d' % (sp1, sp2, cscore)
    cmd3 = 'python -m jcvi.compara.synteny screen --minspan=%d --simple %s.%s.anchors %s.%s.anchors.new' % (minspan, sp1, sp2, sp1, sp2)
    print('[2:1] run sinteny %s X %s => %s ' % (sp1, sp2, cmd1))
    if not os.path.isfile('./' + sp2 + '.bck') :
        print(os.popen(cmd1).read())
    print('[2:2] run sinteny %s X %s => %s ' % (sp1, sp2, cmd2))
    print(os.popen(cmd2).read())
    print('[2:3] run sinteny %s X %s => %s ' % (sp1, sp2, cmd3))
    print(os.popen(cmd3).read())
    return [f.strip().split('\t') for f in open("%s.%s.anchors.simple" % (sp1, sp2)).readlines()]


def anchors2links(data, base, sp2, filtr=[], store=False, karyotypeBase=None):
    sp1 = base
    print('convert anchors to links %s X %s ...' % (sp1, sp2))
    cont = 0
    bed1 = {k[3]: [k[0], int(k[1]), int(k[2])] for k in [f.strip().split('\t') for f in open(sp1 + '.bed').readlines()]}
    bed2 = {k[3]: [k[0], int(k[1]), int(k[2])] for k in [f.strip().split('\t') for f in open(sp2 + '.bed').readlines()]}
    seqs = {}
    seqs2 = []
    seqs3 = {}
    seqs4 = {}
    remo = 0
    tsc = []
    for l in data:
        ma1 = bed1[l[0]]
        ma2 = bed1[l[1]]
        mb1 = bed2[l[2]]
        mb2 = bed2[l[3]]
        seq = bed2[l[2]][0]
        if seq in filtr:
            remo += 1
            continue
        if not seq in seqs:
            seqs[seq] = []
            seqs3[seq] = set()
            seqs4[seq] = [0, '']
        seqs3[seq].add(bed1[l[0]][0])
        if not karyotypeBase is None and not bed1[l[0]][0] in tsc:
            seqs4[seq][0] += int(karyotypeBase[bed1[l[0]][0]].split('\t')[5])
            tsc.append(bed1[l[0]][0])
        seqs[seq].append('%s\t%d\t%d\t%s\t%d\t%d' % (
                bed1[l[0]][0],
                min([ma1[1],ma1[2],ma2[1],ma2[2]]),
                max([ma1[1],ma1[2],ma2[1],ma2[2]]),
                seq,
                min([mb1[1],mb1[2],mb2[1],mb2[2]]),
                max([mb1[1],mb1[2],mb2[1],mb2[2]]))
               )
        seqs2.append(bed1[l[0]][0])
        cont += 1
    print('%d links removeds ...' % remo)
    print('%d links parsed ...' % cont)
    ls = {}
    us = []
    for k in seqs:
        if store:
            l = '%s.%s.%s.links' % (sp1, sp2, k)
            ls[l] = k
            with open(l, 'w') as f:
                f.write('\n'.join(seqs[k]))
        nus = [s.split('\t')[0] for s in seqs[k] if not s.split('\t')[0] in us]
        seqs4[k][1] = nus[int(len(nus)/2)].split('\t')[0]
        us.append(seqs4[k][1])
    if store:
        print('%d CHR LINKS stored [ %s X %s ] ...' % (len(seqs), sp1, sp2))
        
    return list(set(seqs.keys())), list(set(seqs2)), seqs, ls, seqs3, seqs4



def generateCircosConfig(base, spcs, links, rewrite=False, img='circos.png', meio="-", lab="100 Mbp", replaces={}, paleta="grays", min_g=1, max_g=3):
    
    def verifyAndWrite(file, txt):
        if (not os.path.isfile(file)) or rewrite:
            with open(file, 'w') as conf:
                conf.write(txt)
            print('%s file stored ...' % file)
        else:
            print('%s file existis, skiping ...' % file)
    
    verifyAndWrite('circos.conf', "\n".join([
                    "<<include colors_fonts_patterns.conf>>",
                    "<<include ideogram.conf>>",
                    "<<include ticks.conf>>",
                    "<<include colors.brewer.conf>>",
                    "<<include colors.conf>>",
                    "<image>",
                    "<<include etc/image.conf>>",
                    "file*   = %s" % img,
                    "</image>",
                    "",
                    "karyotype = %s.karyotype,%s.karyotype" % (base, ".karyotype,".join(spcs)),
                    "chromosomes_units = 1",
                    "spacing = 0",
                    "chromosomes_display_default = yes",
                    "chromosomes       = /CHR\w+1$/",
                    "chromosomes_scale =  /CHRguava*/:4",
                    "#chromosomes_radius = /CHRegrandis*/:0.8r;",
                    "<myblock>",
                    "paleta = '%s-'" % paleta,
                    "</myblock>",
                    '<links>',
                    "#show = no",
                    '',
                   'bezier_radius = 0r',
                   'thickness = 5p',
                   'radius = 0.98r',
                    '',
                    "\n".join(['\n<link>\n' + 
                               'file = ' + l + '\n' +
                               'color = eval(conf(myblock,paleta) . (counter(link)+1))\n' +
                               '</link>\n' for l in sorted(links, key=lambda e: -int(re.sub("[^0-9]*", "",links[e])))]),
                    '</links>',
                    "",
                    "<plots>",
                    "#show = no",
                    "<plot>",
                    "",
                    "type      = histogram",
                    "file      = %s.gene.density" % spcs[0],
                    "",
                    "r1        = 1.05r",
                    "r0        = 1.01r",
                    "min       = %d" % min_g,
                    "max       = %d" % max_g,
                    "fill_color = black",
                    "#stroke_type = outline",
                    "thickness   = 4",
                    "color       = vdgrey",
                    "extend_bin  = no",
                    "</plot>",
                    "<plot>",
                    "",
                    "type      = histogram",
                    "file      = %s.scale" % base,
                    "",
                    "r1        = 1.05r",
                    "r0        = 1r",
                    "max       = 1",
                    "min       = -1",
                    "fill_color = black",
                    "stroke_type = outline",
                    "thickness   = 1",
                    "color       = black",
                    "extend_bin  = yes",
                    "</plot>",
                    "</plots>",
                    "",
                    '<<include etc/housekeeping.conf>>',
                    'data_out_of_range* = trim'
    ]) + '\n')
    
    verifyAndWrite('ideogram.conf', "\n".join([
        '<ideogram>',
        '<spacing>',
        'default = 0.00007r',
        '</spacing>',
        '',
        '<<include ideogram.position.conf>>',
        '<<include ideogram.label.conf>>',
        '</ideogram>'
    ])) 
                                                        
    verifyAndWrite('ideogram.label.conf', "\n".join([
        'show_label       = yes',
        'label_font       = bold',
        'label_radius     = dims(ideogram,radius) + 0.056r',
        'label_with_tag   = yes',
        'label_size       = 3 ## 36',
        'label_parallel   = no ## yes',
        '#label_case       = lower',
        'label_format     = eval( sprintf("%s",var(label)) )',
        '#label_format     =  eval( %s var(label) eq "%s" ? "%s" :  var(chr) =~ /CHRe.*/ ? var(label) : "")' % (  " ".join(['var(label) eq "%s" ? "%s" :' % (k, replaces[k]) for k in replaces]),   meio, lab),
        '## ',
        '## ' + "## ".join("%s\t%s\t%s\n" % (k[0], k[1], replaces[k[1]]) for k  in sorted([[links[l], l.split('.')[2]] for l in links], key=lambda e: int(re.sub("[^0-9]*", "",e[0]))))
    ]))
    verifyAndWrite('ticks.conf', "\n".join([
        '',
        'radius           = 0.9r',
        'thickness        = 30p',
        'fill             = yes',
        'fill_color       = black',
        '#stroke_thickness = 2',
        '#stroke_color     = black',
        '',
        'show_ticks          = no',
        'show_tick_labels    = no',
        ''
    ]))
    verifyAndWrite('ideogram.position.conf', "\n".join([
        'radius           = 0.9r',
        'thickness        = 30p',
        'fill             = yes',
        'fill_color       = black',
        'stroke_thickness = 1',
        'stroke_color     = black'
    ]))
    verifyAndWrite('colors.conf', "\n".join([
        '',
        '<colors>',
        'red_faint = 255,0,0,0.50',
        'red_faint-1 = 0,0,255,0.80',
        'red_faint-2 = 0,255,0,0.100',
        'red_faint-3 = 0,255,255,0.120',
        'red_faint-4 = 255,0,0,0.60',
        'red_faint-5 = 255,0,255,0.70',
        'red_faint-6 = 255,255,0,0.90',
        'red_faint-7 = 255,100,0,0.55',
        'red_faint-8 = 255,0,100,0.65',
        'red_faint-9 = 100,0,100,0.75',
        'red_faint-10 = 0,255,100,0.85',
        'red_faint-11 = 255,125,50,0.95',
        '</colors>'
    ]))
    
    print('EXECUTE: circos --conf circos.conf')
    os.system('circos --conf circos.conf')
    print('IMAG: %s ... ' % img)
    

def processSinteny(specie_base, species, mod=3, min_links=50, rewrite=False, palet=0, img='circos.png', jcvi_args={}):
    print("[1] parse files to mcscan (%s) ..." % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    specie=specie_base
    s = [specie]
    s.extend(species)
    im = importKaryotypes(specie, s)
    ks = im[0]
    mps = im[1]
    tams = {k: sum([int(ks[k][l].split('\t')[5]) for l in ks[k]]) for k in ks}
    for k in tams:
        print("%s => %s bp" % (k, tams[k]))
    print("Total => %s bp" % sum([tams[k] for k in tams]))
    keep = {}
    links = {}
    paleta = "grays"
    palets = [
        'spectral-6-div',
        'brbg-6-div', 
        'piyg-6-div',
        'prgn-6-div',
        'puor-6-div',
        'rdbu-6-div',
        'rdgy-6-div',
        'rdylbu-6-div',
        'rdylgn-6-div',
        'paired-6-qual',
        'set3-6-qual'    
    ]
    print("\n\n[2] run mcscan (%s) ...\n\n" % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    for sp in species: 
        mins = jcvi_args[sp]['minspan'] if (sp in jcvi_args) and  ('minspan' in jcvi_args[sp]) else 1
        csco = jcvi_args[sp]['cscore'] if (sp in jcvi_args) and  ('cscore' in jcvi_args[sp]) else 99
        rw = jcvi_args[sp]['rewrite'] if (sp in jcvi_args) and  ('rewrite' in jcvi_args[sp]) else False
        dt = runSinteny(specie, sp, minspan=mins, cscore=csco, rewrite=rw)
        seqs = anchors2links(dt, specie, sp)[2]
        keep[sp] = anchors2links(dt, specie, sp, [k for k in seqs if len(seqs[k]) < min_links], True, karyotypeBase=ks[specie])
        paleta = palets[(len(keep)-1) if palet < 1 else (palet-1)].replace('6', str(len(keep[sp][3])))
        print('for %s use palet %s' % (sp, paleta))
        for arq_link in keep[sp][3]:
            links[arq_link] = mps[sp][keep[sp][3][arq_link]]
    print("\n\n[3] parse output mcsan to circus (%s) ...\n\n" % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    kpBase = set()
    print('salvando karyotypes ...')
    all_scfs = []
    replaces = {}
    for k in species:
        with open(k + '.karyotype', 'w') as f:
            tl = 0
            for j in keep[sp][5]:
                ts = (keep[sp][5][j][0] / tams[specie] * 100)
                tl += ts
                keep[sp][5][j].append(("%.1f%%" % ts).replace(".0", ""))
            print('total anc %.2f%%' % tl)
            replaces = {a: keep[sp][5][a][2] for a in keep[sp][5]}
            for l in sorted([ks[k][l] for l in ks[k] if l in keep[k][0]], key=lambda e: int(re.sub("[^0-9]*", "", e.split('\t')[3]))):
                f.write(l + '\n')
                x = [s for s in keep[sp][4][l.split('\t')[2]] if s not in all_scfs]
                x.extend(all_scfs)
                all_scfs = x
            kpBase = kpBase.union(keep[k][1])
    all_scfs.extend([s for s in ks[specie] if s not in all_scfs])
    with open(specie + '.anchored', 'w') as anc:
        anc.write('\n'.join(kpBase) + '\n')
    with open(specie + '.nanchored', 'w') as nanc:
        nanc.write('\n'.join([s for s in ks[specie] if s not in kpBase]) + '\n')
    print('%.1f%% (-%d) scaffolds anchored ...' % (len(kpBase)/len(ks[specie])*100, len(ks[specie])-len(kpBase)))
    with open(specie + '.karyotype', 'w') as f:
        if mod == 1:  ### view all
            print('mod: view all')
            scaf_sorted = ks[specie]  
        elif mod == 2:  ### view anchored
            print('mod: view anchored')
            scaf_sorted = [l for l in all_scfs if l in kpBase]  
        elif mod == 3:  ### view anchored first
            print('mod: view anchored first')
            tt = []
            c = 0
            for k in sorted([ks[specie][s].split('\t') for s in ks[specie] if s not in kpBase], key=lambda e: int(e[5])):
                tt.append("%s\t0\t%s\t%s" % (k[2], k[5], '1' if len(tt) < 15 else '0.2'))
                c += int(k[5])
                if c < 50000000:
                    m = k[2]
            for k in range(len(tt)-3, len(tt)):
                tt[k] = tt[k][:-3] + '1'
            with open(specie + '.scale', 'w') as sc:
                sc.write('\n'.join(tt) + '\n')
            tt = [t.split('\t')[0] for t in tt]
            scaf_sorted = [l for l in all_scfs if l in kpBase]  
            scaf_sorted.extend(tt)
            scaf_sorted.extend([l for l in all_scfs if not l in kpBase and not l in tt])
        else:  ### view anchored spaced
            print('mod: view anchored spaced')
            scaf_sorted = all_scfs 
        f.write('\n'.join([ks[specie][l] for l in scaf_sorted]) + '\n')
    print("\n\n[4] run circos (%s) ...\n\n" % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    generateCircosConfig(specie, species, links, rewrite, img, m, "%.1f%%" % (c/tams[specie]*100), replaces=replaces, paleta=paleta)
    print('\n\nterminado com sucesso (%s) ...\nby mikeias.net' % datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    
   

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("use: ./sinteny.py [-rewrite] [-palet=1] [-img=fig.png] -b specie_base specie2 specie3 ...")
    else:
        rewrite = '-rewrite' in sys.argv
        p = [x for x in sys.argv if x.startswith('-palet=')]
        palet= int(p[0].split('=')[1]) if len(p) > 0 else 1
        img= [x for x in sys.argv if x.startswith('-img=')]
        bindex = sys.argv.index('-b')
        processSinteny(sys.argv[bindex + 1], sys.argv[(bindex + 2):], rewrite=rewrite, palet=palet, img=(img[0].split('=')[1] if (len(img) > 0) else (sys.argv[bindex + 1] +'.png')))


