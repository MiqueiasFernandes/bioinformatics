#!/usr/bin/env python3

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
import os
import threading
import sys


def parseGFF(prefix, color='green'):
    with open(prefix + '.fa') as in_seq_handle:
        seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        limit_info = dict(gff_type = ['gene', 'mRNA', 'CDS'])
        with open(prefix + '.gff3') as in_handle:
            gff = GFF.parse(in_handle, limit_info=limit_info, base_dict=seq_dict)
            ids = []
            print('parsing [%s] files ...' % prefix)
            karyotype = {}
            bed = open(prefix + '.bed', 'w')
            mrnas = open(prefix + '.cds', 'w')
            kar = open(prefix + '.kar', 'w')
            n = 0
            c = 0
            for rec in gff:
                n += 1
                idS = "CHR%s%d" % (prefix, n)
                karyotype[idS] = ("chr\t-\t%s\t%d\t0\t%d\t%s\n" % (idS, n, len(rec), color))
                kar.write(karyotype[idS])
                for gene in [g for g in rec.features if g.type == 'gene']:
                    for mrna in [m for m in gene.sub_features if m.type == 'mRNA' and not m.id in ids]:
                        c += 1
                        ids.append(mrna.id)
                        idmrna = 'MRNA%s%d' % (prefix, c)
                        SeqIO.write(SeqRecord(Seq(''.join([str(cds.extract(rec).seq) for cds in sorted(mrna.sub_features, key=lambda e: e.location.start, reverse=mrna.strand != 1) if cds.type == 'CDS']), generic_dna), id=idmrna, name="", description=""), mrnas, "fasta")
                        bed.write("%s\t%s\n" % (idS, '\t'.join(str(mrna.location)[1:-1].replace('](', ':%s:0:' % idmrna).split(':'))))
            print('%d mrnas parsed in %d seqs' % (c, n))
            bed.close()
            mrnas.close()
            return karyotype
    return False


class Karyotype(threading.Thread):
    force = False
    
    def file(self, prefix, color=0):
        self.p = prefix
        self.c = color
        
    def run(self):
        if (not (os.path.isfile('./' + self.p + '.bed') and os.path.isfile('./' + self.p + '.cds'))) or Karyotype.force:
            r = parseGFF(self.p, color=['red', 'green', 'grey', 'black', 'white'][min(4, self.c)])
            if r:
                self.karyotype = r
            else:
                print('ERRO! fail on import %s' % self.p)
                exit(1)
        else:
            with open(self.p + '.kar') as kar:
                self.karyotype = {k.split('\t')[2]: k.strip() for k in kar.readlines()}
            print(self.p + ' [bed&cds] existis, skipping ...')
        print("%s importado ..." % self.p)


def importFiles(prefixs, force=False):
    Karyotype.force=force
    ts = []
    for p in prefixs:
        k = Karyotype()
        k.file(p, len(ts))
        ts.append(k)
        k.start()
        
    for t in ts:
        t.join()
        
    print('%d importados ...' % len(ts))
    return {k.p: k.karyotype for k in ts}


def runSinteny(base, sp2, minspan=1, cscore=7, rewrite=False):
    sp1=base
    if rewrite:
        os.remove("%s.%s.last.filtered" % (sp1, sp2))
        print('rewrite links ...')
    cmd1 = 'lastdb %s %s.cds -P 24' % (sp2, sp2)
    cmd2 = 'python -m jcvi.compara.catalog ortholog %s %s --cscore=.%d' % (sp1, sp2, cscore)
    cmd3 = 'python -m jcvi.compara.synteny screen --minspan=%d --simple %s.%s.anchors %s.%s.anchors.new' % (minspan, sp1, sp2, sp1, sp2)
    print('[1] run sinteny %s X %s => %s ' % (sp1, sp2, cmd1))
    if not os.path.isfile('./' + sp2 + '.bck') :
        print(os.popen(cmd1).read())
    print('[2] run sinteny %s X %s => %s ' % (sp1, sp2, cmd2))
    print(os.popen(cmd2).read())
    print('[3] run sinteny %s X %s => %s ' % (sp1, sp2, cmd3))
    print(os.popen(cmd3).read())
    return [f.strip().split('\t') for f in open("%s.%s.anchors.simple" % (sp1, sp2)).readlines()]


def anchors2links(data, base, sp2, filtr=[], store=False):
    sp1 = base
    print('convert anchors to links %s X %s ...' % (sp1, sp2))
    cont = 0
    bed1 = {k[3]: [k[0], int(k[1]), int(k[2])] for k in [f.strip().split('\t') for f in open(sp1 + '.bed').readlines()]}
    bed2 = {k[3]: [k[0], int(k[1]), int(k[2])] for k in [f.strip().split('\t') for f in open(sp2 + '.bed').readlines()]}
    seqs = {}
    seqs2 = []
    seqs3 = {}
    remo = 0
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
        seqs3[seq].add(bed1[l[0]][0])
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
    ls = []
    if store:
        for k in seqs:
            l = '%s.%s.%s.links' % (sp1, sp2, k)
            ls.append(l)
            with open(l, 'w') as f:
                f.write('\n'.join(seqs[k]))
        print('%d seqs stored [ %s X %s ] ...' % (len(seqs), sp1, sp2))
    return list(set(seqs.keys())), list(set(seqs2)), seqs, ls, seqs3



def generateCircosConfig(base, spcs, links, rewrite=False, img='circos.png', meio="-", lab="100 Mbp"):
    
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
                    '<links>',
                    "\n".join(['\n<link>\n' + 
                               'file = ' + l + '\n' +
                               'color = ' + links[l] + '\n' +
                               'bezier_radius = 0r\n' +
                               'thickness = 5p\n' +
                               'radius = 0.98r\n' +
                               '</link>\n' for l in links]),
                    '</links>',
                    "",
                    "<plots>",
                    "<plot>",
                    "",
                    "type      = histogram",
                    "file      = histogram.txt",
                    "",
                    "r1        = 1.043r",
                    "r0        = 0.995r",
                    "max       = 1",
                    "min       = -1",
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
        'label_size       = 36',
        'label_parallel   = yes',
        '#label_case       = lower',
        'label_format     =  eval( var(chr) eq "%s" ? "%s" :  var(chr) eq "CHRguava6" ? "" : (var(label) < 12 ? replace(replace(var(chr), "CHRegrandis", "Chr" ), "CHRguava", "") : "") )' % (meio,lab),
        ''
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
    specie=specie_base
    s = [specie]
    s.extend(species)
    ks = importFiles(s)
    keep = {}
    links = {}
    palets = [
        'spectral-6-div-',
        'brbg-6-div-', 
        'piyg-6-div-',
        'prgn-6-div-',
        'puor-6-div-',
        'rdbu-6-div-',
        'rdgy-6-div-',
        'rdylbu-6-div-',
        'rdylgn-6-div-',
        'paired-6-qual-',
        'set3-6-qual-'    
    ]
    for sp in species: 
        mins = jcvi_args[sp]['minspan'] if (sp in jcvi_args) and  ('minspan' in jcvi_args[sp]) else 1
        csco = jcvi_args[sp]['cscore'] if (sp in jcvi_args) and  ('cscore' in jcvi_args[sp]) else 99
        rw = jcvi_args[sp]['rewrite'] if (sp in jcvi_args) and  ('rewrite' in jcvi_args[sp]) else False
        dt = runSinteny(specie, sp, minspan=mins, cscore=csco, rewrite=rw)
        seqs = anchors2links(dt, specie, sp)[2]
        keep[sp] = anchors2links(dt, specie, sp, [k for k in seqs if len(seqs[k]) < min_links], True)
        p = palets[(len(keep)-1) if palet < 1 else (palet-1)].replace('6', str(len(keep[sp][3])))
        print('for %s use palet %s' % (sp, p[:-1]))
        c = 1
        for l in keep[sp][3]:
            links[l] = p + str(c)
            c += 1
    kpBase = set()
    print('salvando karyotypes ...')
    all_scfs = []
    for k in species:
        with open(k + '.karyotype', 'w') as f:
            for l in sorted([ks[k][l] for l in ks[k] if l in keep[k][0]], key=lambda e: int(e.split('\t')[3])):
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
                if c < 100000000:
                    tt.append("%s\t0\t%s\t%s" % (k[2], k[5], '1' if len(tt) < 15 else '0.2'))
                    c += int(k[5])
                else:
                    break
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
    
    generateCircosConfig(specie, species, links, rewrite, img, m, str(c)[:3]+ ' Mbp')
    print('terminado com sucesso ...\nby mikeias.net')
    
    
    

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

