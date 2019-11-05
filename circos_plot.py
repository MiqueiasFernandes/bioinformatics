#!/usr/bin/python3

# Copyright(c) 2019 - Mikeias Fernandes <bio@mikeias.net>

# circos_plot.py - simplify use of circos for plot whole genome view


#################################################

import sys
import os
from Bio import SeqIO

#################################################
usage = "circos_plot.py fasta gff3 links_file cuffdiff_file misa_file snp_file"
#################################################

try:
    script,fasta,gff,links_file,cuffdiff_file,misa_file,snp_file=sys.argv
except: sys.exit("Correct usage is:\n"+usage)


class Circos:
    def __init__(self, fasta, gffs, gff_keys, links_file, limit_size, 
                 cuffdiff_file, misa_file, snp_file,
                 out_dir, 
                 chrs='Chr'):
        self.fasta = fasta
        self.chrs = chrs
        self.out = out_dir
        self.gffs = {"GFF%d"%(x+1): [gffs[x], gff_keys[x]] for x in range(len(gffs))}
        self.links = links_file
        self.limit_size = limit_size
        self.cuffdiff_file = cuffdiff_file
        self.misa_file = misa_file
        self.snp = snp_file
    
    def prepare(self):
        print(" [1/5] importando fasta ...")
        self.fasta_dict = SeqIO.to_dict(SeqIO.parse(self.fasta, 'fasta'))
        
        print(" [2/5] preparando arquivos em %s ..." % self.out)
        os.mkdir(self.out)
        os.symlink(self.snp, self.out + '/snp.hist')
        os.symlink(self.fasta, self.out + '/genome.fa')
        for k, v in self.gffs.items():
            os.link(v[0], self.out + '/' + k)
        self.misa2gff3(self.misa_file, self.out + '/misa.gff3')
        self.geneateConfig(chrs=self.chrs)
        
        print(" [3/5] importando gffs  " + self.out)
        for k, v in self.gffs.items():
            self.gffdensity(self.out + '/' + k, v[1], self.out + '/' + k + '.bars', self.chrs)
        self.gffdensity(self.out + '/misa.gff3', 'SSR', self.out + '/misa.bars', self.chrs)
            
        print(" [4/5] importando links  " + self.out)
        self.minimize_links(self.links, self.out + '/links', self.limit_size)
        
        print(" [5/5] importando gene expression  " + self.out)
        for f in self.importGeneExp(self.cuffdiff_file, self.out, self.chrs):
            self.normWindow(f, f + ".norm")
        print("run circos --config circos.conf ... ")
        os.system("cd " + self.out + " && circos -conf circos.conf")
    
    def karyotype(self, out, chrs='Chr'):
        print("file" + out + ' ... OK!')
        with open(out, 'w') as o:
            for c in [c for c in self.fasta_dict if c.startswith(chrs)]:
                o.write('chr\t-\t%s\t%s\t0\t%d\tgreen\n' % (c, c, len(self.fasta_dict[c])))
        return out
    
    def geneateConfig(self, chrs='Chr', rewrite=False):
        self.karyotype(self.out + '/genome.karyotype', chrs)
        
        conf = """
<<include etc/colors_fonts_patterns.conf>>
<<include ideogram.conf>>
<<include ticks.conf>>
<<include colors.brewer.conf>>
<<include colors.conf>>
<image>
    <<include etc/image.conf>>
</image>
chromosomes_units           = 1000000
chromosomes_display_default = yes
chromosomes_color  = /.*/:piyg-3-div-3
karyotype = genome.karyotype
<links>
    <link>
        file          = links
        radius        = 0.2r
        color         = piyg-3-div-1
        bezier_radius = 0.1r
        thickness     = 5
        #ribbon = yes
        <rules>
            <rule>
                condition  = var(intrachr) && abs(var(pos1)-var(pos2)) > 50Kb
                color        = piyg-3-div-2
            </rule>
            <rule>
                condition  = var(intrachr) && abs(var(pos1)-var(pos2)) > 60Kb
                color        = piyg-3-div-3
            </rule>
        </rules>
    </link>
</links>
<plots>
    <plot>
        type      = scatter
        file      = misa.bars
        r0        = 0.21r
        r1        = 0.26r
        color     = piyg-4-div
    </plot>
    <plot>
        type      = heatmap
        file      = flor_folha.heatmap.norm
        r0        = 0.27r
        r1        = 0.30r
        color     = piyg-4-div
        stroke_thickness = 1
        stroke_color     = black
    </plot>
    <plot>
        type      = heatmap
        file      = flor_fruto.heatmap.norm
        r0        = 0.31r
        r1        = 0.34r
        color     = piyg-4-div
        stroke_thickness = 1
        stroke_color     = black
    </plot>
    <plot>
        type      = heatmap
        file      = folha_fruto.heatmap.norm
        r0        = 0.35r
        r1        = 0.38r
        color     = piyg-4-div
        stroke_thickness = 1
        stroke_color     = black
    </plot>
    <plot>
        type      = histogram
        file      = GFF1.bars
        r0        = 0.39r
        r1        = 0.45r
        stroke_type = outline
        thickness   = 4
        color       = lgrey
        fill_color = lgrey
        extend_bin  = yes
    </plot>
    <plot>
        show      = no
        type      = line
        file      = GFF2.bars
        r0        = 0.39r
        r1        = 0.45r
        stroke_type = outline
        thickness   = 4
        extend_bin  = yes
        color       = piyg-3-div-2
        <rules>
            <rule>
                condition    = var(value) > 50
                color        = piyg-3-div-3
            </rule>
            <rule>
                condition    = var(value) < 25
                color        = piyg-3-div-1
            </rule>
        </rules>
    </plot>
    <plot>
        type      = line
        file      = snp.hist
        r0        = 0.46r
        r1        = 0.52r
        stroke_type = outline
        thickness   = 4
        extend_bin  = yes
        color       = piyg-3-div-2
        <rules>
            <rule>
                condition    = var(value) > 200
                color        = piyg-3-div-1
            </rule>
            <rule>
                condition    = var(value) < 100
                color        = piyg-3-div-3
            </rule>
        </rules>
    </plot>
</plots>
<<include etc/housekeeping.conf>>
        """

        ideogram = """
<ideogram>
    <spacing>
        default = 0.0025r
        break   = 0.5r
    </spacing>
<<include ideogram.position.conf>>
<<include ideogram.label.conf>>
</ideogram>
        """
        
        label = """
show_label       = yes
label_font       = default
label_radius     = 0.63r
label_with_tag   = yes
label_size       = 36
label_parallel   = yes
#label_case       = lower
label_format     = eval(sprintf("%s", replace(var(label), "Chr", "LG") ))
    """
        
        position = """
radius           = 1.5r
thickness        = 30p
fill             = yes
stroke_thickness = 2
stroke_color     = black
    """
        ticks = """
show_ticks          = no
show_tick_labels    = no
<ticks>
    skip_first_label = no
    skip_last_label  = no
    radius           = dims(ideogram,radius_outer)
    tick_separation  = 2p
    label_separation = 5p
    multiplier       = 1e-6
    color            = black
    thickness        = 4p
    size             = 20p
    <tick>
        spacing        = 1u
        show_label     = no
        thickness      = 2p
        color          = dgrey
    </tick>
    <tick>
        spacing        = 5u
        show_label     = no
        thickness      = 3p
        color          = vdgrey
    </tick>
    <tick>
        spacing        = 10u
        show_label     = yes
        label_size     = 20p
        label_offset   = 10p
        format         = %d
        grid           = yes
        grid_color     = dgrey
        grid_thickness = 1p
        grid_start     = 0.5r
        grid_end       = 0.999r
    </tick>
</ticks>
    """

        colors = """
<colors>
 chrs = 254,158,218
</colors>"""

        def persist(file, var):
            if not rewrite and os.path.exists(file) and os.path.isfile(file):
                raise BaseException('ERROR: file ' + file + ' EXISTIS! call with REWRITE arg!')
            with open(file, 'w') as o:
                o.write(var)
            print('file %s ... OK' % file)

        persist(self.out + '/circos.conf', conf)
        persist(self.out + '/ideogram.conf', ideogram)
        persist(self.out + '/ideogram.label.conf', label)
        persist(self.out + '/ideogram.position.conf', position)
        persist(self.out + '/ticks.conf', ticks)
        persist(self.out + '/colors.conf', colors)
    
    def minimize_links(self, file, out, limitMIN=10):
        seg_dup = [l.strip().split('\t') for l in open(file).readlines() if not l.startswith('#')]
        ss = [(x[0], int(x[1]), int(x[2]), int(x[2]) - int(x[1]), x[3], x[4], x[5]) for x in seg_dup if len(x) == 6]
        links = [x for x in ss if x[3] > limitMIN * 1000]
        with open(out, 'w') as o:
            for l in links:
                o.write("%s\t%d\t%d\t%s\t%s\t%s\n" % (l[0], l[1], l[2], l[4], l[5], l[6]))
        print("%d writed in %s ...." % (len(links), out))
    
        
    def run(self):
        print('iniciando ...')
        self.prepare()
        

    def gffdensity(self, gff, keys, out, chrs='Chr', window=100000):
        print("importando %s ..." % keys)
        genes = [l.strip().split('\t') for l in open(gff).readlines() if l.count("\t" + keys + "\t") > 0]
        print("parseando %s ..." % keys)
        chr2genes = {x: [(int(y[3]),int(y[4])) for y in genes if y[0] == x] for x in set([x[0] for x in genes])}
        print("salvando em " + out + ' ...')
        with open(out, 'w') as o:
            for c in [c for c in self.fasta_dict if c.startswith(chrs) and c in chr2genes]:
                ranges = list(sorted(set(list(range(1, len(self.fasta_dict[c]), window)) + [len(self.fasta_dict[c])+1])))
                for i in range(1, len(ranges)):
                    o.write('%s\t%d\t%d\t%d\n' % (c, ranges[i-1], ranges[i]-1, len([x for x in chr2genes[c] if 
                          (x[0] >= ranges[i-1] and x[0] <= ranges[i]-1) or 
                          (x[1] >= ranges[i-1] and x[1] <= ranges[i]-1)])))
    
    def misa2gff3(self,file, out):
        k = [l.strip().split("\t") for l in open(file).readlines() if l.count('\t') > 0]
        if k[0] == 'ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend'.split("\t"):
            with open(out, 'w') as o:
                o.write('\n'.join(['\t'.join([
                    x[0], 
                    'misa', 
                    'SSR', 
                    x[5], 
                    x[6], 
                    '.', '.', '.', 
                    'ID=' + x[2] + '.' + x[0] + '.' + x[1]]) for x in k[1:]]) +'\n')
        else:
            raise BaseException('File ' + file + ' not of misa output!')

    
    def importGeneExp(self, file, out_dir='./', chrs=None):
        header = 'test_id\tgene_id\tgene\tlocus\tsample_1\tsample_2\tstatus\tvalue_1\tvalue_2\tlog2(fold_change)\ttest_stat\tp_value\tq_value\tsignificant'.split('\t')
        lines = [x.strip().split("\t") for x in open(file).readlines() if x.count("\t") > 3]
        if lines[0] != header:
            raise BaseException('Header NOT OK, verify if file is from CUFFDIFF %s' % file)

        ls = [(x[3],'_'.join(sorted([x[4],x[5]])),float(x[9]),float(x[12]), x[13] == 'yes') for x in lines[1:]]
        condicoes = {x: [y for y in ls if y[1] == x] for x in set([z[1] for z in ls])}
        files = []
        for k, v in condicoes.items():
            with open(out_dir + '/' + k + '.heatmap', 'w') as o:
                ls = [(x[0].replace(':', '\t').replace('-', '\t'), str(x[2])) for x in v if x[4] and (chrs is None or x[0].startswith(chrs))]
                o.write('\n'.join(['\t'.join(x) for x in ls]) + '\n')
                print('file %s ... OK' % (out_dir + '/' + k + '.heatmap'))
                files.append(out_dir + '/' + k + '.heatmap')
        return files
    
    def normWindow(self, file, out):
        fasta = self.fasta_dict
        ls = [l.strip().split("\t") for l in open(file).readlines() if l.count("\t") > 0]
        chr2exp = {x: [(int(z[1]), int(z[2]), abs(float(z[3]))) for z in ls if z[0] == x] for x in set([y[0] for y in ls])}
        window = 10 * 10000
        with open(out, 'w') as o:
            for c in [c for c in fasta if c in chr2exp]:
                ranges = list(sorted(set(list(range(1, len(fasta[c]), window)) + [len(fasta[c])+1])))
                for i in range(1, len(ranges)):
                    t = [x[2] for x in chr2exp[c] if 
                          (x[0] >= ranges[i-1] and x[0] <= ranges[i]-1) or 
                          (x[1] >= ranges[i-1] and x[1] <= ranges[i]-1)]
                    o.write('%s\t%d\t%d\t%f\n' % (c, ranges[i-1], ranges[i]-1, sum(t)/len(t) if len(t) > 0 else 0))
                    

circos = Circos(fasta=fasta, gffs=[gff,gff], gff_keys=['gene', 'CDS'], links_file=links_file, limit_size=40,cuffdiff_file=cuffdiff_file, misa_file=misa_file,snp_file=snp_file,out_dir='circos_out')
circos.run()

print('by mikeias.net')
