#!/usr/bin/env python
# coding: utf-8

# In[18]:


import os
from Bio import SeqIO
import numpy as np
from IPython.display import Image, display


class Circos:
    def __init__(self, files_fasta=None, min_size={}, chr_pattern={}, sliding_window=1000000, paleta='spectral-div', exact_density=False):
        
        self.allPaletas = {
            '-seq': ['blues','bugn','bupu','gnbu','greens','greys','oranges','orrd','pubu','pubugn','purd','purples','rdpu','reds','ylgn','ylgnbu','ylorbr','ylorrd'],
            '-div': ['brbg','piyg','prgn','puor','rdbu','rdgy','rdylbu','rdylgn','spectral'],
            '-qual': ['paired', 'set3']}
        
        if not any([paleta.endswith(x) for x in self.allPaletas]):
            raise BaseException('INVALID paleta {}. VISIT: http://circos.ca/documentation/tutorials/configuration/colors/'.format(paleta))
        
        self.paletaDefault = paleta
        self.fastas = [] if files_fasta is None else ([files_fasta] if type(files_fasta) == str else files_fasta)
        self.min_size = min_size
        self.chr_pattern = chr_pattern
        self.sliding_window = sliding_window
        self.exact_density = exact_density
        
        self.circos_conf = {
            'chromosomes_units': str(self.sliding_window),
            'chromosomes':  '',
            'chromosomes_display_default': 'yes',
            '<ideogram>>': {
                '<spacing>>': {
                    'default': '0.005r'
                },
                'radius': '0.90r',
                'thickness': '20p',
                'fill': 'yes',
                'stroke_color': ' dgrey',
                'stroke_thickness': '2p',
                'show_label': 'yes',
                'label_font': 'default',
                'label_radius': 'dims(image,radius) - 60p',
                'label_size': '30',
                'label_parallel': 'yes'
            },
            '<plots>>': {},
            'show_ticks': 'yes',
            'show_tick_labels': 'yes',
            '<ticks>>': {
                'radius':'1r',
                'color':'black',
                'thickness':'2p',
                'multiplier':'1e-6',
                'format':'%d',
                '<tick>0': {
                    'spacing': '0.5u',
                    'size': '6p',
                    'color':'lgrey'
                },
                '<tick>1': {
                    'spacing': '1u',
                    'size': '10p',
                    'color':'grey'
                },
                '<tick>2': {
                    'spacing': '5u',
                    'size': '12p'
                },
                '<tick>3': {
                    'spacing': '10u',
                    'size': '18p',
                    'show_label': 'yes',
                    'label_size': '15p',
                    'label_offset': '10p',
                    'thickness':'3p',
                    'format': '%d Mb'
                }
            }
        }
        
        self.plots = {}
        self.seqs = {}
        self.seqNames = {}
        self.allSeqs = []
        
        ## Gerar o karyotype
        for fasta in self.fastas:
            self.addFasta(fasta, auto=True)
        self.LastFasta = None
        self.LastPlot = None
            
    def addFasta(self, fasta, min_size=None, chr_pattern=None, auto=False):
        if not min_size is None:
            self.min_size[fasta] = min_size
        if not chr_pattern is None:
            self.chr_pattern[fasta] = chr_pattern
        parsed = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
        parsed = {key: parsed[key] for key in parsed if 
                  ((self.min_size is None) or (fasta not in self.min_size) or (len(parsed[key]) >= self.min_size[fasta])) and 
                  ((self.chr_pattern is None) or (fasta not in self.chr_pattern) or (key.startswith(self.chr_pattern[fasta])))
                 }
        self.seqs[fasta] = parsed
        self.seqNames[fasta] = [x for x in parsed]
        if not auto:
            self.fastas.append(fasta)
        self.allSeqs.extend([x for x in parsed])
        
        minW = self.sliding_window
        for k, x in parsed.items():
            minW = min(len(x), minW)
        
        if minW < self.sliding_window:
            print('Set NEW sliding_window TO: %d for fasta %s' % (minW, fasta))
        
        if 'karyotype' in self.circos_conf:
            self.circos_conf['karyotype'] = ';'.join([self.circos_conf['karyotype'], fasta + '.karyotype'])
        else:
            self.circos_conf['karyotype'] = fasta + '.karyotype'
        open(fasta + '.karyotype', 'w').writelines([
            'chr\t-\t{0}\t{0}\t0\t{1}\tgreen\t\n'.format(x, len(parsed[x])) 
                                                for x in parsed])
        self.LastFasta = fasta
        return self
    
    def showAll(self):
        self.circos_conf['chromosomes_display_default'] = 'yes'
        return self
    
    def showOnly(self, chrs=None):
        self.circos_conf.update({
            'chromosomes':  ';'.join([chrs] if type(chrs) == str else([self.seqNames[x][0] for x in self.fastas] if chrs is None else chrs)),
            'chromosomes_display_default': 'no',
        })
        return self
    
    def withHist(self, file, is_gff=True, feature='gene', r1=None, r2=None, qualitativo=False):
        if self.LastFasta is None:
            raise Exception('Error in call')
        if is_gff:
            self.with_density(file, self.LastFasta, feature=feature, r1=None, r2=None, label=file, qualitativo=qualitativo)
        else:
            self.addPlot(file, r1=r1, r2=r2)
        return self
    
    def ideogram_radius(self, radius=90):
        self.circos_conf['<ideogram>>']['radius'] = '{}r'.format(radius/100)
        return self
    
    def ideogram_thickness(self, tk=20):
        self.circos_conf['<ideogram>>']['thickness'] = '{}p'.format(tk)
        return self
    
    def ideogram_fill(self, fill=True):
        self.circos_conf['<ideogram>>']['fill'] = 'yes' if fill == True else 'no'
        return self
    
    def ideogram_stroke_color(self, c='dgrey'):
        self.circos_conf['<ideogram>>']['stroke_color'] = c
        return self
    
    def ideogram_stroke_thickness(self, s=2):
        self.circos_conf['<ideogram>>']['stroke_thickness'] = '{}p'.format(s)
        return self
    
    def ideogram_show_label(self, s=True):
        self.circos_conf['<ideogram>>']['show_label'] = 'yes' if s == True else 'no'
        return self
    
    def ideogram_label_size(self, s=30):
        self.circos_conf['<ideogram>>']['label_size'] = str(s)
        return self
    
    def setEspaco(self, e=0.005, breaks=False, inter=None):
        if not inter is None:
            self.circos_conf['<ideogram>>']['<spacing>>']['<pairwise %s>\n spacing ' % ' '.join(inter)] =  str(e) + ' \n</pairwise>\n'
            return self
        self.circos_conf['<ideogram>>']['<spacing>>']['break' if breaks else 'default'] = '{}r'.format(e)
        return self    
    
    def show_ticks(self, s=True):
        self.circos_conf['show_ticks'] = 'yes' if s == True else 'no'
        return self
    
    def show_tick_labels(self, s=True):
        self.circos_conf['show_tick_labels'] = 'yes' if s == True else 'no'
        return self
    
    def tick_format(self, f='%d Mb'):
        self.circos_conf['<ticks>>']['<tick>3']['format'] = f
        return self
    
    def addLinks(self, file, raio='0.3r'):
        if not '<links>>' in self.circos_conf:
            self.circos_conf['<links>>'] = {}
        label = '<link>' + str(len(self.circos_conf['<links>>']))
        plot = self.circos_conf['<links>>'][label] = {
            'file': file,
            'color': 'black_a5', 
            'radius': raio, 
            'bezier_radius': '0r', 
            'thickness': '1'
        }
        self.LastPlot = self.plots[label] = plot
        return self
    
    def configureLink(self, key, value=None, file=None):
        link = self.getLink(file)
        if not link is None:
            if value is None:
                if value in link[key]:
                    del link[key]
                    return self
            link[key] = value
            return self
        print('Link `{}` não existe!'.format(file))
        return self

    def run_circos(self):
        
        try: os.remove('circos.png')
        except: pass
        
        write_seq = [
            'karyotype',
            'chromosomes_units',
            'chromosomes',
            'chromosomes_display_default',
            '<ideogram>>',
            'show_ticks',
            '<ticks>>',
            '<plots>>',
            'show_tick_labels'
        ]
        
        lines = []    
        for conf in write_seq + [x for x in self.circos_conf if not x in write_seq]:
            line = conf
            if conf.startswith("<"): ## 1 level
                line = line[:-1]
                line += '\n'
                for sub_confs_K, sub_confs_V  in self.circos_conf[conf].items():
                    if sub_confs_K.startswith("<"):
                        line += '  ' + sub_confs_K[:-1] + '\n' + '\n'.join(['    {} = {}'.format(k, v) for k, v in sorted([x for x in sub_confs_V.items()], 
                                                                                                                          key=lambda e: (e[0].count('##R'), 
                                                                                 - int(e[0].split('<rules>\n##R')[1].split('#')[0]) if e[0].count('##R') > 0 else e[0]))]) + '\n  </' + sub_confs_K[1:-1] + '\n'
                        line = line.replace('</rules>\n    <rules>', '')
                    else:
                        line += '  {} = {}\n'.format(sub_confs_K, sub_confs_V)
                line += ('</' + conf[1:-1] + '\n')
            else:                    ## 1 level
                line += ' = '  + self.circos_conf[conf] + '\n'
            lines.append(line + '\n')
        open('circos.conf', 'w').writelines(lines + [
            '<image>\n',
            '<<include etc/image.conf>>\n',
            '</image>\n',
            '<<include etc/colors_fonts_patterns.conf>>\n',
            '<<include etc/housekeeping.conf>>\n'
        ])
        os.system("circos --conf circos.conf 1> circos.log 2> circos.err")
        return self
    
    def addPlot(self, file, tipo='histogram', r1=None, r2=None, label=None):
        ind = len(self.circos_conf['<plots>>'])+1
        plot = {
            'type'        : tipo,
            'file'        : file,
            'r1'          : '0.{}r'.format((str(ind)+'8') if r1 is None else r1),
            'r0'          : '0.{}r'.format(ind if r2 is None else r2),
            'stroke_type' : 'outline',
            'thickness'   : '4',
            'color'       : 'vdgrey',
            #'extend_bin'  : 'yes',
            'fill_color'  : 'grey'
        }
        self.circos_conf['<plots>>']['<plot>%s' % ind] = plot
        if label is None:
            label = file
        self.LastPlot = self.plots[label] = plot
        return self
    
    def getPlot(self, label):
        if label in self.plots.values():
            return label
        if label in self.plots:
            return self.plots[label]
        print('plot não existe!')
        
    def getLink(self, file=None):
        if file is None:
            if len(self.circos_conf['<links>>']) == 1:
                return list(self.circos_conf['<links>>'].values())[0]
            return None
        if not '<links>>' in self.circos_conf:
            return None
        for k, v in self.circos_conf['<links>>'].items():
            if v['file'] == file:
                return self.circos_conf['<links>>'][k]
            
    def getTrack(self, label=None):
        if label is None:
            return self.LastPlot
        l = self.getLink(label)
        if l is None:
            return self.getPlot(label)
        return l
    
    def trackTipe(self, label=None):
        if label is None:
            return 0 if self.LastPlot is None else (1 if self.LastPlot in list(self.circos_conf['<links>>'].values()) else 2)
        l = self.getLink(label)
        if l is None:
            return 2
        return 1

    def configPlots(self, labels, config=None, value=None, configs={}):
        for l in labels:
            self.configPlot(label=l, config=config, value=value, configs=configs)
        return self
    
    def configAllPlots(self, excepts = [], config=None, value=None, configs={}):
        labels = [x for x in circos.plots if not x in excepts]
        return self.configPlots(labels, config=config, value=value, configs=configs)
        
    def configPlot(self, label=None, config=None, value=None, configs={}):
        
        if not configs is None and len(configs) > 0:
            for k, v in configs.items():
                self.configPlot(label=label, config=k, value=v)
            return self
        
        plot = self.LastPlot if label is None else self.getPlot(label)
        if plot:
            self.LastPlot = plot
            if config is None:
                return self
            if value is None:
                if config in plot:
                    del plot[config]
            else:
                plot[config] = value
        else:
            self.LastPlot = None
        return self
    
    def orientarFora(self, label=None, ok=True):
        self.configPlot(self.LastPlot if label is None else label, config='orientation', value='out' if ok else 'in')
        return self
    
    def drawAs(self, tipo, label=None):
        return self.configPlot(self.LastPlot if label is None else label, config='type', value=tipo)
    
    def __plotLimits(self, plot):
        ss = [int(l.strip().split('\t')[3]) for l in open(plot['file']).readlines() if l.count('\t') > 2 and l.split('\t')[0] in self.allSeqs]
        d, a = min(ss), max(ss)
        pts = (a - d) / 4
        return d, a, [d, pts, pts * 2, a - pts, a], (a - d) / 8
    
    def drawAsLine(self, label=None):
        self.configPlot(self.LastPlot if label is None else label, config='type', value='line')
        plot = self.LastPlot
        if plot:
            lims = self.__plotLimits(plot)[3]
            plot['<axes>#'] = """
        <axis>
            color     = vlgreen
            thickness = 2
            position  = {}
        </axis>
        <axis>
            color     = lgreen
            thickness = 1
            position  = {}
        </axis>
        <axis>
            color     = lgrey
            thickness = 1
            position  = {}
        </axis>
        <axis>
            color     = grey
            thickness = 2
            position  = {}
        </axis>
        <axis>
            color     = lgrey
            thickness = 1
            position  = {}
        </axis>
        <axis>
            color     = vlred
            thickness = 1
            position  = {}
        </axis>
        <axis>
            color     = lred
            thickness = 2
            position  = {}
        </axis>
    </axes>        
            """.format(lims * 1, lims * 2, lims * 3, lims * 4, lims * 5, lims * 6, lims * 7)
        return self ##     vlG      lg        lgy       g          lgy      vlR       lR
    
    def drawAsHist(self, label=None):
        return self.configPlot(self.LastPlot if label is None else label, config='type', value='histogram')
    
    def drawAsScat(self, label=None):
        return self.configPlot(self.LastPlot if label is None else label, config='type', value='scatter')
    
    def drawAsHeat(self, label=None, div=9):
        return self.configPlot(self.LastPlot if label is None else label, configs={
            'type': 'heatmap',
            'color': self.paletaDefault.replace('-', '-%d-' % div),
            'fill_color': None
        })
    
    def setRaio(self, r0, r1=None, label=None, esp=None):
        plot = self.getTrack(label)
        if plot and not plot is None:
            self.LastPlot = plot
            plot['r0'] = '{}r'.format(r0/100)
            plot['r1'] = '{}r'.format((r1/100) if not r1 is None else (((r0+esp)/100) if not esp is None else ((r0+1)/100)))
        return self
                                      
    def tracksIn(self, tracks=None, r0=30, r1=90):
        tracks = [x for x in self.plots] if tracks is None else tracks
        ini = r0
        step = (r1-r0)/len(tracks)
        for t in tracks:
            fim = ini + step
            self.setRaio(label=t, r0=ini, r1=fim-1)
            ini = fim
        return self
    
    def hideAll(self, tracks=None):
        for t in [x for x in self.plots] if tracks is None else tracks:
            self.hide(t)
        return self
                                      
    def hide(self, label=None):
        plot = None
        if label:
            plot = self.getPlot(label)
        elif self.LastPlot:
            plot = self.getPlot(label)
        if plot:
            plot['show'] = 'no'
        return self
            
    def unHide(self, label=None, labels=None):
        if not labels is None:
            for l in labels:
                self.unHide(self, label=l)
            return self
        plot = None
        if label:
            plot = self.getPlot(label)
        elif self.LastPlot:
            plot = self.getPlot(label)
        if plot:
            plot['show'] = 'yes'
        return self
                
    def removePlot(self, file):
        rm = None
        for k, v in self.circos_conf['<plots>>'].items():
            if v['file'] == file:
                rm = k
                break
        if not rm is None:
            del self.circos_conf['<plots>>'][rm]
            print('ok!')
        return self
        
    def with_density(self, gff, spec, feature='gene', tipo='histogram', r1=None, r2=None, label=None, qualitativo=False):
        file = gff+'_'+feature+'.hist'
        self.__calc_density(gff, 
                          self.seqs[spec], 
                          file, 
                          feature=feature, 
                          seq_filters=[x for x in self.seqs[spec]],
                          qualitativo=qualitativo)
        self.addPlot(file, tipo=tipo, r1=r1, r2=r2, label=file if label is None else label)
    
    def __calc_density(self, gff, sizes, out, feature='gene', seq_filters=None, qualitativo=False):
        print('calc_density', gff)
        gff_parsed = [l.strip().split("\t") for l in open(gff).readlines() if l.count('\t'+feature+'\t') > 0]
        if seq_filters is None:
            seq_filters = list(set([x[0] for x in gff_parsed]))
        else:
            gff_parsed = [x for x in gff_parsed if x[0] in seq_filters]
            
        windows = {seq: {w: 0 for w in range(1, len(sizes[seq]) + self.sliding_window, self.sliding_window)} for seq in seq_filters}
        for g in gff_parsed:
            windows[g[0]][((int(g[4] if g[6] == '-' else g[3]) // self.sliding_window)*self.sliding_window)+1] += 1
        lines = []
        for seq, v in windows.items():
            for sliding_w, qtd in v.items():
                if self.exact_density and ((sliding_w + self.sliding_window) >  len(sizes[seq])):
                    qtd /= (len(sizes[seq])-sliding_w)/self.sliding_window
                if qualitativo:
                    if qtd > 0:
                        qtd = 1
                    else:
                        qtd = 0
                lines.append('{}\t{}\t{}\t{}\n'.format(seq, sliding_w, sliding_w + self.sliding_window - 1, qtd))
        open(out, 'w').writelines(lines)
           
    def show(self, width=300):
        display(Image(filename='circos.png', width=width))
        
    def plot(self):
        self.run_circos().show()
        
    def setRule(self, values={}, condition='1', label=None, key='R1'):
        rs = self.getTrack(label)
        if not rs is None:
            rs['<rules>\n##' + key + '##\n     <rule>\n            condition '] = condition + '\n            ' +             '\n            '.join(["{} = {}".format(k, v) for k, v in values.items()]) +             '\n        </rule>\n</rules>'
        return self
        
    def colorizeLinkBySize(self, a, b, p, track=None, paleta=None, inv=False):
        return self.colorize(a, b, p, "var(size1) >= {}", track=track, paleta=paleta, inv=inv)
    
    def colorizePlotByValue(self, a, b, p, track=None, paleta=None, inv=False):
        return self.colorize(a, b, p, "var(value) >= {}", track=track, paleta=paleta, inv=inv)
        
    def colorize(self, a, b, p, condicao, track=None, paleta=None, inv=False):
        paleta = (self.paletaDefault if paleta is None else paleta).replace('-', '-%d-') + '-%d'
        s = (b-a)//p
        cont = 1000
        for x in [x for x in np.arange(a, b, s)]: ##[x for x in range(a, b, s)]:
            self.setRule(values={
                'color': paleta % (p+1,abs(((p+2) if inv else 0) - (cont-999))),
                'z': str(cont)
            }, condition = condicao.format(x), label=track, key='R%d'%cont)
            cont += 1
        return self
    
    def autoColorize(self, parts=5, track=None, paleta=None, inv=False):
        tt = self.trackTipe(track)
        file = None
        if tt > 0:
            fName = self.getTrack(track)['file']
            file = [l.strip().split('\t') for l in open(fName).readlines() if l.count("\t") > 1]
        if tt == 1:
            ss = [abs(int(x[1])-int(x[2])) for x in file]
            #print('colorindo {} de {} a {} / {} partes como Link'.format(fName, min(ss), max(ss), parts))
            return self.colorizeLinkBySize(min(ss), max(ss), p=parts, track=track, paleta=paleta, inv=inv)
        if tt == 2:
            ss = [float(x[3]) for x in file]
            #print('colorindo {} de {} a {} / {} partes como Plot'.format(fName, min(ss), max(ss), parts))
            return self.colorizePlotByValue(min(ss), max(ss), p=parts, track=track, paleta=paleta, inv=inv)
        
    def testPaletas(self, allS=False):
        if allS:
            self.showAll()            
        self.plot()
        pd = self.paletaDefault.split('-')[0]
        td = '-' + self.paletaDefault.split('-')[1]
        os.system('rm -rf paletas && mkdir paletas')
        for t, v in self.allPaletas.items():
            for p in v:
                cmd = "rm -rf {2}{3} && mkdir {2}{3} && sed 's/{0}/{2}/g ; s/{1}/{3}/' circos.conf >  {2}{3}/circos.conf && cd {2}{3} && circos --conf circos.conf 1> circos.log 2> circos.err && cp circos.png ../paletas/{2}{3}.png &".format(
                    pd, td, p, t)
                os.system(cmd)
        print('on end, run: "zip -r paletas.zip paletas/" and get all images')
                
    def renameImg(self, new):
        os.system('cp circos.png %s.png' % new)
        os.system('cp circos.svg %s.svg' % new)
        
    def order(self, file='Pguajava.gff3_gene.hist', suff='f', add=False):
        os.system("bash -c 'cut -f4 {0} | sort -g | paste - - - -  > d && cat <(cut -f1 d) <(cut -f2 d | sort -rg) <(cut -f3 d) <(cut -f4 d | sort -gr) | paste - - - - - - - - - - - > d2 && paste <(cut -f-3 {0}) <(for i in `seq 1 11` ; do cat <(cut -f$i d2) ; done | grep -v ^$) > {0}_{1}_1 && rm d d2'".format(file, suff))
        os.system("bash -c 'paste <(cut -f-3 {0}) <(cut -f4 {0} | sort -g) > {0}_{1}_2'".format(file, suff))
        os.system("bash -c 'paste <(cut -f-3 {0}) <(cut -f4 {0} | sort -gr) > {0}_{1}_3'".format(file, suff))
        os.system("bash -c \"awk '{if ( \$4 > 1 ) print \$1,\$2,\$3,0 ; else print \$1,\$2,\$3,1}' %s > %s_%s_4\"" % (file, file, suff))
        os.system("bash -c \"awk '{if ( \$4 > 5 ) print \$1,\$2,\$3,0 ; else print \$1,\$2,\$3,1}' %s > %s_%s_5\"" % (file, file, suff))
        os.system("bash -c \"awk '{if ( \$4 > 10 ) print \$1,\$2,\$3,0 ; else print \$1,\$2,\$3,1}' %s > %s_%s_6\"" % (file, file, suff))
        if add:
            self.addPlot("{0}_{1}_1".format(file, suff))
            self.addPlot("{0}_{1}_2".format(file, suff)) 
            self.addPlot("{0}_{1}_3".format(file, suff)) 
            self.addPlot("{0}_{1}_4".format(file, suff)) 
            self.addPlot("{0}_{1}_5".format(file, suff)) 
            self.addPlot("{0}_{1}_6".format(file, suff)) 
        return self
        

# In[38]:


circosU = Circos(paleta='rdpu-seq', sliding_window=30000)\
.addFasta('Pguajava.fasta')\
.withHist('Pguajava.gff3')\
.addLinks('sedef.links.40k')\
.autoColorize()\
.configureLink('radius', '0.39r')\
.withHist("TE_DNA.gff", feature='Pguajava_REPET_TEs')\
.drawAsScat()\
.orientarFora(ok=False)\
.autoColorize()\
.withHist("TE_RNA.gff", feature='Pguajava_REPET_TEs')\
.drawAsScat()\
.autoColorize()\
.addPlot('to_circos_G')\
.drawAsLine()\
.addPlot('flor.coverage.bg')\
.drawAsHeat()\
.addPlot('folha.coverage.bg')\
.drawAsHeat()\
.addPlot('fruto.coverage.bg')\
.drawAsHeat()\
.configPlots([
    'folha.coverage.bg', 
    'flor.coverage.bg', 
    'fruto.coverage.bg'
], configs = {
    'scale_log_base': '.2'
})\
.ideogram_show_label(False)\
.show_ticks(False)\
.setEspaco(e=0.002)\
.setEspaco(e='500u',inter=['scaffold0200', 'scaffold0001'])\
.setRaio(label='Pguajava.gff3', r0=101, r1=108).configPlot('to_circos_G', configs= {
    'fill_color': None
}).setRaio(r0=91, r1=98).tracksIn(r0=70, r1=90, tracks= [
    'folha.coverage.bg', 
    'flor.coverage.bg', 
    'fruto.coverage.bg']).tracksIn(r0=40, r1=60, tracks= [
    "TE_DNA.gff", 
    "TE_RNA.gff",
]).tracksIn(r0=40, r1=70, tracks= [
    "TE_DNA.gff", 
    "TE_RNA.gff",
])\
.showOnly(['scaffold%s' % str(s).rjust(4, '0') for s in range(1,201)])\
.show_tick_labels(False)


# In[39]:


circosU.plot()




# In[12]:


circosU.testPaletas()

