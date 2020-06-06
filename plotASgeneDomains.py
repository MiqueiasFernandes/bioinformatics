#!/usr/bin/env python
# coding: utf-8

# Copyright(c) 2020 - Miquéias Fernandes <bio@mikeias.net>

from IPython.display import display, clear_output, SVG
import svgwrite
import csv
import os
from Bio import Seq, SeqIO
import cairosvg
from optparse import OptionParser

usage = "usage: %prog [-f] [-g] dir_mats1/ dir_mats2/ ..."
parser = OptionParser(usage)
parser.add_option("-e", "--events", dest="events", help="only plot this events Ids, comma separated")
parser.add_option("-x", "--xevents", dest="xevents", help="skip this exents Ids, comma separated")
parser.add_option("-f", "--fasta", dest="fasta", help="fasta for identify stop codons", metavar="FILE")
parser.add_option("-g", "--gff", dest="gff", help="gff file to get gene structure", metavar="FILE")
parser.add_option("-i", "--interpro", dest="interpro", help="InterproScan5 output tsv file with domains", metavar="FILE")
parser.add_option("-d", "--db", dest="dbs", help="InterproScan5 dbs to consider comma separated like: Pfam,SMRT")
parser.add_option("-r", "--relative", dest="relative", action="store_true", default=False, help="plot gene with relative sizes")
parser.add_option("-s", "--skip", dest="skip", action="store_true", default=False, help="skip genes with multiple AS events")
parser.add_option("-m", "--mrna", dest="mrna", type="int", default=50, help="mRNA height to plot default [%default]")
parser.add_option("-t", "--threshold", dest="fdr", type="float", default=.05, help="FDR threshold defalt [%default]")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="print details")
(options, dirs) = parser.parse_args()

if options.fasta is None:
    parser.error("options [-f] is obrigatory!")
if options.gff is None:
    parser.error("options [-g] is obrigatory!")
if len(dirs) < 1:
    parser.error("set at least one MATS dir!")

VERBOSE = options.verbose
    
class Draw:
    def __init__(self, file, width=600, height=200, marginX=10, marginY=10, pdf=True, png=True):
        
        self.file = file
        self.png, self.pdf = png, pdf
        self.width = width
        self.height = height
        self.marginX = marginX
        self.offSetX = marginX//2
        self.marginY = marginY
        self.offSetY = marginY//2
        self.boundingBox = ((self.offSetX, self.offSetY), (width-marginX, height-marginY))
        
        self.drawObj = svgwrite.Drawing(filename=file, height='{}px'.format(height), width='{}px'.format(width))
        self.drawObj.viewbox(0, 0, self.width, self.height)
        self.mainBox = self.drawObj.add(self.drawObj.g(id='mainBox'))
        
        self.PATTERNS_COLORS = ['green', 'blue', 'red', 'pink', 'gray', 'brown', 
                                'orange', 'yellow', 'black', 'white', 'purple', 'cyan']
        self.added_patterns = []
        self.boxes = []
        self.tops = []
        
    def updateTop(self):
        box = self.box()
        for top in self.tops:
            box.add(top)
    
    def store(self):
        self.updateTop()
        self.drawObj.save()
        if self.pdf:
            cairosvg.svg2pdf(url=self.file, write_to=self.file.replace('.svg', '.pdf'))
        if self.png:
            cairosvg.svg2png(url=self.file, write_to=self.file.replace('.svg', '.png'))
        
    def show(self, clear=True):
        self.updateTop()
        self.drawObj.save()
        if clear:
            clear_output(wait=True)
        display(SVG(filename=self.file))

    def listPatterns(self):
        return ' '.join([' '.join([c+'|', c+'-', c+'\\', c+'/', c+'0', c+'#']) for c in self.PATTERNS_COLORS]).split()
    
    def fillPattern(self, name, d=5):
        
        if not name in self.listPatterns():
            raise Exception('Options are: ' + ', '.join(self.listPatterns()))
        
        name = name.replace('\\', '@').replace('#', '1')
        if not name in self.added_patterns:
            c = name[:-1]
            o = name[-1]
            s = (3, 3) if (o in ['-', '|']) else (100, 100)
            pattern = self.drawObj.pattern(size=s, id='pattern-' + name, patternUnits="userSpaceOnUse")
            
            if (o in ['-', '|']):
                tLine = ((0, 2), (3, 2)) if o == '-' else ((2, 3), (2, 0))
                pattern.add(self.drawObj.line(start=tLine[0], end=tLine[1], stroke=c, stroke_width=0.5))
            else:
                if o in ['0', '1']:
                    ss = d
                    pd = ss/2
                    pd2 = ss/8
                    for i in range(100):
                        if i % ss == 0:
                            for j in range(100):
                                if j % ss == 0:
                                    if o == '0':
                                        pattern.add(self.drawObj.circle(center=(i+ pd, j+pd), r=pd-1, fill=c))
                                    if o == '1':
                                        pattern.add(self.drawObj.rect(insert=(i, j), size=(ss-(pd2*2), ss-(pd2*2)), fill=c))
                else:
                    a = 100 if o == '/' else 0
                    b = 0 if o == '/' else 100
                    for i in range(-100, 200, 4):
                        pattern.add(self.drawObj.line(start=(i, a), end=(i+100, b), stroke=c, stroke_width=1))
            
            self.drawObj.defs.add(pattern)
            self.added_patterns.append(name)
        return "url(#pattern-%s)" % name
    
    def box(self, name=None):
        name = ('box%s' % len(self.boxes)) if name is None else name
        return self.drawObj.add(self.drawObj.g(id=name))
    
    def square(self, x, y, w, h, c, box=None):
        box = self.mainBox if box is None else box
        box.add(self.drawObj.rect(insert=(x, y), size=(w, h), fill=c))
        
    def circle(self, x, y, r, c, box=None):
        box = self.mainBox if box is None else box
        box.add(self.drawObj.circle(center=(x, y), r=r, fill=c))
        
    def line(self, x1, y1, x2, y2=None, c='black', box=None, sw=2):
        box = self.mainBox if box is None else box
        y2 = y1 if y2 is None else y2
        box.add(self.drawObj.line(start=(x1, y1), end=(x2, y2), stroke_width=sw, stroke=c))
    
    def text(self, text, x, y, box=None, italic=False, family='Times', size=12, bold=False, 
             center=False, rg=False, down=False, c='black', diag=False, vcenter=False):
        box = self.mainBox if box is None else box
        anchor = 'end' if rg else 'start'
        anchor =  'middle' if center else anchor
        base = "hanging" if down else 'baseline'
        base = 'middle' if vcenter else base
        box.add(self.drawObj.text(
            text, insert=(x, y), 
            font_style="italic" if italic else 'normal', 
            font_family=family, 
            font_size=size, 
            text_anchor=anchor, 
            alignment_baseline=base,
            fill=c,
            transform="rotate(30 {},{})".format(x, y) if diag else "rotate(0)",
            font_weight= 'bold' if bold else 'normal'))
        
    def star(self, x, y, box=None, fill="gray", stroke="white", r=5, top=True):
        box = self.mainBox if box is None else box
        STAR8 = "0 8 L 8 12 L 7 2 L 14 -4 L 4 -6 L 0 -15 L -4 -6 L -14 -4 L -7 2 L -8 12 L 0 8"
        STAR5 = "0 5 L 5 8 L 4 1 L 9 -3 L 2 -4 L 0 -10 L -2 -4 L -9 -3 L -4 1 L -5 8 L 0 5"
        STAR = STAR8 if r == 8 else STAR5
        d = 'M ' + ' L '.join([str(int(p.split()[0])+x) + ' ' + str(int(p.split()[1])+y) for p in STAR.split(" L ")]) + ' z'
        _star = self.drawObj.path(d=d, fill=fill, stroke=stroke)
        box.add(_star)
        if top:
            self.tops.append(_star)

    def path(self, d, c, f=None, box=None, sw=2):
        box = self.mainBox if box is None else box
        box.add(self.drawObj.path(d=d, fill=c if f is None else f, stroke=c, stroke_width=sw))
        
    def seta(self, x1, x2, y, box=None, c='black', h=20, start=False, pS=4, baixo=True, top=True):
        box = self.mainBox if box is None else box
        m, n, o = (y+h) if baixo else (y-h), (y+pS) if baixo else (y-pS), pS if baixo else -pS
        s1 = self.drawObj.path(d='M {0} {1} C {0} {2} {3} {2} {3} {1}'.format(x1, y+o, m, x2), 
                          fill="none", stroke=c, stroke_width=pS//2)
        box.add(s1)
        xs = x1 if start else x2
        s2 = self.drawObj.path(
            d='M {0} {1} L {2} {1} L {3} {4} z'.format(
                xs-pS-(pS/8), n, xs+pS-(pS/8), xs-(pS/8), y), stroke=c, fill=c)
        box.add(s2)
        if top:
            self.tops.append(s1)
            self.tops.append(s2)
        
    def seta2(self, x, y, box=None, c='black', pS=4, top=True):
        box = self.mainBox if box is None else box
        s1 = self.drawObj.line(start=(x, y+pS), end=(x, y+(pS*2)), stroke_width=pS, stroke=c)
        s2 = self.drawObj.path(
            d='M {0} {1} L {2} {1} L {3} {4} z'.format(
                x-pS, y+pS, x+pS, x, y), stroke=c, fill=c)
        box.add(s1)
        box.add(s2)
        
        if top:
            self.tops.append(s1)
            self.tops.append(s2)


class Canvas:
    def __init__(self, x, y, w, h, draw, box='Box', parents={}):
        self.x, self.y, self.w, self.h = x, y, w, h
        self.bottom = y + h
        self.endX = x + w
        self.half_y = y + (h/2)
        self.half_x = x + (w/2)
        self.quarterTop = self.half_y - (h/4)
        self.quarterBottom = self.half_y + (h/4)
        self.draw = draw
        parents.update({box: self})
        self.parents = parents
        self.boxName = box
        self.box = draw.box(box)
        
    @staticmethod
    def create(name=None, x=0, y=0, w=600, h=400, box='Box'):
        name = 'test.svg' if name is None else name
        draw = Draw(name, width=w, height=h)
        return Canvas(x=x, y=y, w=w, h=h, draw=draw, box=box)
        
    def clip(self, x=None, y=None, w=None, h=None, pdx=None, pdy=None, box='Box'):
        pdx = 0 if pdx is None else pdx
        pdy = 0 if pdy is None else pdy
        x = (self.x if x is None else x) + pdx
        y = (self.y if y is None else y) + pdy
        w = (self.w if w is None else w)
        h = (self.h if h is None else h)
        draw=self.draw
        return Canvas(x, y, w, h, draw, box=box, 
                      parents=self.parents) if box != self.boxName else self.__init__(x, y, w, h, draw, box=box)
    
    def fill(self, c='gray'):
        self.draw.square(x=self.x, y=self.y, w=self.w, h=self.h, c=c, box=self.box)
    
    def __repr__(self):
        return str(('Canvas [{},{}] w: {} h: {}'.format(self.x, self.y, self.w, self.h)))
        
class Domain:
    def __init__(self, description, start, end, entry=None, db=None, label=None, pattern='green/'):
        self.start = start
        self.end = end
        self.db = db
        self.description = description
        self.entry = entry
        self.pattern = pattern
        self.boxH = 8   
        self.label = description.split()[0] if label is None else label
    
    def __repr__(self):
        return str((self.description, self.start, self.end))
    
class Feature:
    def __init__(self, name, start, end, gene, strand=None, contig=None, typeName='Feature', typeColor='black', boxH=6):
        self.name = name
        if start > end:
            raise Exception('Start must be greather then END ({}, {})'.format(start, end))
        self.start = start
        self.end = end
        self.size = end - start + 1
        self.gene = gene
        self.mrna = None
        self.contig = gene.contig if contig is None else contig
        self.strand = gene.strand if strand is None else strand
        self.typeName = typeName
        self.typeColor = typeColor
        self.boxH = boxH
    
    def pb2pixel(self, canvas, bp):
        u = canvas.w / self.size
        if bp < self.start or bp > self.end:
            raise Exception('Bp out of range! {} ({}, {})'.format(bp, self.start, self.end))
        diff = ((bp - self.start) if self.strand else (self.end - bp))
        x1 = canvas.x + (u * diff)
        x2 = canvas.x + (u * (diff+1))
        return x1, (x1+x2)/2, x2
    
    def pbCoordsXpixelCoords(self, canvas, a, b):
        a, b = self.pb2pixel(canvas, a), self.pb2pixel(canvas, b)
        return min([a[0], a[2], b[0], b[2]]), max([a[0], a[2], b[0], b[2]])
    
    def childCanvas(self, canvas, feature, pdy=None, h=None):
        x1, x2 = self.pbCoordsXpixelCoords(canvas, feature.start, feature.end)
        return canvas.clip(x=x1, w=x2-x1, pdy=pdy, h=h, box=feature.typeName + '_' + feature.name)
    
    def plot(self, canvas, showBox=True):
        if showBox:
            canvas.draw.square(canvas.x, canvas.half_y - (self.boxH/2), canvas.w, self.boxH, c=self.typeColor, box=canvas.box)

    def __repr__(self):
        return str((self.name, self.start, self.end, '+' if self.strand else '-', '%dbp' % self.size))

    
class Cds(Feature):
    def __init__(self, name, start, end, gene):
        super().__init__(name, start, end, gene, typeName='CDS', typeColor='purple', boxH=7)
    
    def plot(self, canvas, domains=[]):
        super().plot(canvas)
        for a, b, domain in domains:
            x1, x2 = self.pbCoordsXpixelCoords(canvas, a, b)
            canvas.draw.square(x1, canvas.half_y - (domain.boxH/2), x2-x1, 
                               domain.boxH, canvas.draw.fillPattern(domain.pattern), box=canvas.box)

            
class Intron(Feature):
    def __init__(self, name, start, end, gene):
        super().__init__(name, start, end, gene, typeName='Intron', typeColor='gray', boxH=8)
        
    def plot(self, canvas, retained=0, use=False, baixo=True, star=False):
        super().plot(canvas, showBox=False)
        white_seg = 3
        if not star:
            canvas.draw.square(canvas.x, canvas.half_y - (self.boxH / 2) - white_seg, 
                               canvas.w, self.boxH + (white_seg*2), c='white', box=canvas.box)
            if use:
                canvas.draw.seta(canvas.x, canvas.endX, canvas.half_y + (5 if baixo else - 5), box=canvas.box, baixo=baixo)
            
        if retained > 0:
            if not star:
                boxH = self.boxH * 1.5
                canvas.draw.square(canvas.x, canvas.half_y - (boxH/2), 
                                   canvas.w, boxH, c=canvas.draw.fillPattern('brown0', d=4), box=canvas.box)
            if retained > 1:
                canvas_gene = canvas.parents[self.gene.name]
                px = self.gene.pb2pixel(canvas_gene, retained)[1]
                canvas_gene.draw.star(px, canvas.half_y, fill='black', box=canvas.box)
                return True
            return False
        
        top = canvas.half_y -  self.boxH
        canvas.draw.line(x1=canvas.x, y1=canvas.half_y, 
                         x2=canvas.half_x + 0.5, y2=top, c=self.typeColor, box=canvas.box, sw=2)
        canvas.draw.line(x1=canvas.half_x - 0.5, y1=top, 
                         x2=canvas.endX, y2=canvas.half_y, c=self.typeColor, box=canvas.box, sw=2)
        return False
            
class Exon(Feature):
    def __init__(self, name, start, end, gene):
        super().__init__(name, start, end, gene, typeName='Exon', typeColor='blue', boxH=7)
        
    def getAdjIntron(self, o):
        name = self.name + '-I-' + o.name
        return Intron(name, self.end + 1, o.start -1, self.gene)
    
    def getIntrons(self):
        return [i for i in self.mrna.getIntrons() if i.end + 1 == self.start or i.start - 1 == self.end]
    
    def plot(self, canvas, skip=False, baixo=True, a3ss=0, a5ss=0):
        super().plot(canvas)
        
        if skip or ((a3ss + a5ss) > 0):
            introns = self.getIntrons()
            canvas_gene = canvas.parents[self.gene.name]
            if (skip and (len(introns) != 2)) or (((a3ss + a5ss) > 0) and (len(introns) < 1)):
                raise Exception('Number of introns diff of 2: ', introns)
            pbs = []
            for i in introns:
                pbs.extend(self.gene.pbCoordsXpixelCoords(canvas_gene, i.start, i.end))
                
            x1, x2 = min(pbs), max(pbs)
            if skip:
                canvas.draw.seta(x1, x2, canvas.half_y + (5 if baixo else - 5), box=canvas.box, baixo=baixo)
            if a3ss:
                canvas.draw.seta(x1, canvas.x, canvas.half_y - 5, box=canvas.box, baixo=False)
                _, x, _ = self.pb2pixel(canvas, a3ss)
                canvas.draw.seta(x1, x, canvas.half_y + 5, box=canvas.box, baixo=True)
            if a5ss:
                canvas.draw.seta(canvas.endX, x2, canvas.half_y - 5, box=canvas.box, baixo=False)
                _, x, _ = self.pb2pixel(canvas, a5ss)
                canvas.draw.seta(x, x2, canvas.half_y + 5, box=canvas.box, baixo=True)

    
class Mrna(Feature):
    def __init__(self, name, start, end, gene):
        super().__init__(name, start, end, gene, typeName='mRNA', typeColor='green', boxH=7)
        gene.mrnas.append(self)
        self.__exons = []
        self.__cds = []
        self.__introns = []
        
    def getExons(self):
        return self.__exons
    
    def getIntrons(self):
        return self.__introns
    
    def getCds(self):
        return sorted(self.__cds, key=lambda c: c.start if self.strand else -c.end)
    
    def cdsSize(self):
        return sum([c.size for c in self.__cds])
    
    def addCds(self, cds):
        for c in cds:
            c.mrna = self
        self.__cds = cds
    
    def addExons(self, exons):
        for e in exons:
            e.mrna = self
        self.__exons = sorted(exons, key=lambda e: e.start)
        if len(self.__exons) > 1:
            self.__introns = [e.getAdjIntron(self.__exons[self.__exons.index(e)+1]) for e in self.__exons[:-1]]
            for i in self.__introns:
                i.mrna = self
        
    def trinca2pb(self, trinca):
        cont_trinca = 0
        cont_pb = 0
        pbs = []
        for cds in self.getCds():
            for i in (range(cds.start, cds.end+1) if self.strand else range(cds.end, cds.start-1, -1)):
                if cont_pb == 0:
                    cont_pb = 1
                    cont_trinca += 1
                elif cont_pb == 3:
                    cont_pb = 0
                if cont_trinca == trinca:
                    pbs.append(i)
                elif cont_trinca > trinca:
                    return pbs
                if cont_pb > 0:
                    cont_pb += 1
        if cont_trinca < trinca:
            raise Exception('Trinca {} out of range {}'.format(trinca, cont_trinca))
        return pbs       
            
    def plot(self, canvas, domains=[], ri={}, se=[], mxe=[], a3ss={}, a5ss={}, title=False, anot=''):
        canvas_gene = canvas.parents[self.gene.name]
        if title:
            canvas.draw.text(self.name, canvas_gene.x, canvas.y, box=canvas.box, down=True)
        if anot != '':
            canvas.draw.text(anot, canvas_gene.endX, canvas.y, box=canvas.box, down=True, rg=True)
        super().plot(canvas)
        se = [se] if type(se) == str else se
        seU = []
        seD = se
        iuU = []
        iuD = []
        
        if len(mxe) == 2:
            e1, e2 = mxe[0], mxe[1]
            e1 = [e for e in self.__exons if e.name == e1][0]
            e2 = [e for e in self.__exons if e.name == e2][0]
            i1 = e1.getIntrons()
            i2 = e2.getIntrons()
            ie1 = [i for i in i1 if not i in i2][0]
            ie2 = [i for i in i2 if not i in i1][0]
            seU = [e2.name]
            iuU = [ie1.name]
            seD.append(e1.name)
            iuD = [ie2.name]
        
        for exon in self.__exons:
            c = self.childCanvas(canvas, exon)
            exon.plot(c, 
                      skip=exon.name in (seU+seD), baixo=exon.name in seD,
                      a3ss=a3ss[exon.name] if exon.name in a3ss else 0,
                      a5ss=a5ss[exon.name] if exon.name in a5ss else 0
                     )
        d_buff = {}
        for domain in domains:
            bps = self.trinca2pb(domain.start) + self.trinca2pb(domain.end)
            d_buff[(min(bps), max(bps))] = domain
        used_doms = []
        for cds in self.__cds:
            c = self.childCanvas(canvas, cds)
            doms = []
            for r, d in d_buff.items():
                if  (r[0] >= cds.start and r[0] <= cds.end and r[1] > cds.end): ##dom_comeca_cds
                    doms.append((max(r[0], cds.start), cds.end, d))
                elif  (r[1] >= cds.start and r[1] <= cds.end and r[0] < cds.start): ##dom_termina_cds
                    doms.append((cds.start, min(r[1], cds.end), d))
                elif  (r[0] >= cds.start and r[1] <= cds.end): ##dom_dentro_cds
                    doms.append((r[0], r[1], d))
                elif  (r[0] <= cds.start and r[1] >= cds.end): ##dom_cobre_cds
                    doms.append((cds.start, cds.end, d))
            cds.plot(c, doms)
            used_doms.extend(doms)
        repaint = []
        for intron in self.__introns:
            c = self.childCanvas(canvas, intron)
            if intron.plot(c, 
                        retained=ri[intron.name] if intron.name in ri else 0, 
                        use = intron.name in (iuD + iuU),
                        baixo = intron.name in iuD
                       ):
                repaint.append(intron)
        for intron in repaint:
            c = self.childCanvas(canvas, intron)
            intron.plot(c, 
                        retained=ri[intron.name] if intron.name in ri else 0, 
                        use = intron.name in (iuD + iuU),
                        baixo = intron.name in iuD,
                        star=True
                       )
            
        return used_doms
                
        
class Gene(Feature):
    def __init__(self, name, start, end, strand, contig):
        super().__init__(name, start, end, self, strand=strand, contig=contig, typeName='Gene', typeColor='gray', boxH=6)
        self.mrnas = []
        self.primary = []
    
    def getMrnas(self):
        return self.mrnas if len(self.primary) < 1 else self.primary
    
    def setPrimary(self, mrnas=None):
        mrnas = [max(self.mrnas, key=lambda e: e.size).name] if mrnas is None else mrnas
        self.primary = [m for m in self.mrnas if m.name in mrnas]
    
    def plot(self, canvas, domains={}, 
             title=True, mrna_title=True,
             scale=False, scale_strand=False, legend=True,
             marks=[], mrnas_plot=[], mrnas_anot={},
             ri={}, se={}, mxe={}, a3ss={}, a5ss={}, relative=[]):
        canvas = canvas.clip(pdx=20, w=canvas.w-40, box='PlotGene')
        if len(relative) > 0:
            maiorG = max(relative, key=lambda e: e.size).size
            u = canvas.w / maiorG
            canvas = canvas.clip(w=self.size*u, box=self.name + '_relative')
            
        canvas = canvas.clip(box=self.name)
        super().plot(canvas, showBox=False)
        
        titleH = 10 if title else 0
        scaleH = 20 if scale else 0
        legendH = 10 if legend else 0
        scale_threshold = 5000
        
        if title:
            canvas.draw.text(self.name, canvas.x, canvas.y + (titleH/2), vcenter=True, box=canvas.box, bold=True)
        
        mrnas = [m for m in self.getMrnas() if (len(mrnas_plot) < 1) or (m.name in mrnas_plot)]
        mrnaCanvas = canvas.clip(h=canvas.h-scaleH-titleH-legendH, pdy=titleH, box='Mrnas')
        mrnaH = mrnaCanvas.h / len(mrnas)
        used_doms = []
        for mrna in mrnas:
            childCanvas = self.childCanvas(mrnaCanvas, mrna, h=mrnaH, pdy=(mrnas.index(mrna) * mrnaH))
            canvas.draw.square(canvas.x, childCanvas.half_y - (self.boxH/2), 
                               canvas.w, self.boxH, self.typeColor, box=childCanvas.box)
            canvas.draw.text("5'", canvas.x-3, childCanvas.half_y, vcenter=True, rg=True, box=childCanvas.box)
            canvas.draw.text("3'", canvas.endX+3, childCanvas.half_y, vcenter=True, box=childCanvas.box)
            canvas.draw.circle(canvas.x if self.strand else canvas.endX, childCanvas.half_y, (self.boxH/2), self.typeColor, box=childCanvas.box)
            
            if self.strand:
                canvas.draw.path('M {5} {1} L {0} {1} L {2} {3} L {0} {4} L {5} {4}'.format(
                    canvas.endX, 
                    childCanvas.half_y - (self.boxH/2),
                    canvas.endX + (self.boxH/2), childCanvas.half_y,
                    childCanvas.half_y + (self.boxH/2), canvas.endX - 2), c='none', f=self.typeColor, box=childCanvas.box)
            else:
                canvas.draw.path('M {5} {1} L {0} {1} L {2} {3} L {0} {4} L {5} {4}'.format(
                    canvas.x, 
                    childCanvas.half_y - (self.boxH/2),
                    canvas.x - (self.boxH/2), childCanvas.half_y,
                    childCanvas.half_y + (self.boxH/2), canvas.x + 2), c='none', f=self.typeColor, box=childCanvas.box)
                
            used_doms.extend(mrna.plot(childCanvas, 
                      title=mrna_title, anot=mrnas_anot[mrna.name] if mrna.name in mrnas_anot else '',
                      domains=domains[mrna.name] if mrna.name in domains else [],
                      ri=ri[mrna.name] if mrna.name in ri else {},
                      se=se[mrna.name] if mrna.name in se else [],
                      mxe=mxe[mrna.name] if mrna.name in mxe else [],
                      a3ss=a3ss[mrna.name] if mrna.name in a3ss else {},
                      a5ss=a5ss[mrna.name] if mrna.name in a5ss else {}
                     ))
            for bp in marks:
                canvas.draw.seta2(self.pb2pixel(canvas, bp)[1], childCanvas.quarterBottom, box=childCanvas.box)
                
        if scale:
            s1 = 10 if self.size < scale_threshold else 100
            s2 = 100 if self.size < scale_threshold else 1000
            y = canvas.bottom-(scaleH + legendH - 2)
            cont = 1 if self.strand else self.size
            for pb in range(self.start, self.end + 1):
                pos = pb if scale_strand else cont
                if (self.size < 1000) or pos % s1 == 0:
                    it = str(pos)
                    _, px, _ = self.pb2pixel(canvas, pb)
                    py = y+(2 * ''.join(['0' if x is '0' else '1' for x in it])[::-1].index('1'))+1
                    canvas.draw.line(px, y, px, py, box=canvas.box)
                    if (self.size < 1000 and pos % s1 == 0) or pos % s2 == 0:
                        it = it if pos < 1000 else ('%.1fk' % (pos / 1000)).replace('.0k', 'k')
                        canvas.draw.text(it, px, py+8, box=canvas.box, size=8, center=True)
                cont += 1 if self.strand else -1
        
        if legend and len(used_doms) > 0:
            doms = []
            tsize = 9
            
            for _, _, dom in used_doms:
                if not dom.entry in [d.entry for d in doms]:
                    doms.append(dom)
            dom_w = canvas.w / len(doms)
            dom_y = canvas.bottom - (legendH/2)
            x = canvas.x
            pq = sum([len(d.description)*tsize for d in doms]) > canvas.w
            for dom in doms:
                xd = x + (doms.index(dom) * dom_w)
                lg = legendH/2
                canvas.draw.square(xd, dom_y - (lg/2), lg, lg, Cds('', 1, 2, self).typeColor, box=canvas.box)
                canvas.draw.square(xd, dom_y - (lg/2), lg, lg, canvas.draw.fillPattern(dom.pattern), box=canvas.box)
                canvas.draw.text(dom.label if pq else dom.description, xd + lg + 2, dom_y, vcenter=True, box=canvas.box, size=tsize)
    
    
class Gff:
    def __init__(self, gff, fasta=None, target=None):
        self.gff = gff
        self.fasta = None if fasta is None else SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
        self.data = [ x for x in list(csv.reader(open(gff), delimiter='\t')) if len(x) == 9 and not x[0].startswith('#')]
        if not target is None:
            seqs = set([x[0] for x in self.data if target in x[8]])
            self.data = [d for d in self.data if d[0] in seqs]
        self.genes = {y[0]: Gene(y[0], int(y[1][3]), int(y[1][4]), y[1][6] != '-', contig=y[1][0]) for y in [(x[8].split('ID=')[1].split(';')[0], x) 
                                                                                    for x in self.data if x[2] == 'gene']}
        self.mrnas = {y[0]: Mrna(y[0], int(y[1][3]), int(y[1][4]), self.genes[y[1][8].split('Parent=')[1].split(';')[0]]) for y in [(x[8].split('ID=')[1].split(';')[0], x) 
                                                                                    for x in self.data if x[2] == 'mRNA']}
        mes = {m: [] for m in self.mrnas}
        mcs = {m: [] for m in self.mrnas}
        for exon in [x for x in self.data if x[2] == 'exon']:
            n = exon[8].split('ID=')[1].split(';')[0] if 'ID=' in exon[8] else None
            for p in exon[8].split('Parent=')[1].split(';')[0].split(','):
                n = '{}-{}-{}-{}-{}-{}'.format(exon[0], exon[2], exon[3], exon[4], exon[6], p) if n is None else n
                mes[p].append(Exon(n, int(exon[3]), int(exon[4]), self.mrnas[p].gene))
                
        for m, exons in mes.items():
            self.mrnas[m].addExons(exons)
            
        for cds in [x for x in self.data if x[2] == 'CDS']:
            n = cds[8].split('ID=')[1].split(';')[0] if 'ID=' in exon[8] else None
            for p in cds[8].split('Parent=')[1].split(';')[0].split(','):
                n = '{}-{}-{}-{}-{}-{}'.format(exon[0], exon[2], exon[3], exon[4], exon[6], p) if n is None else n
                mcs[p].append(Cds(n, int(cds[3]), int(cds[4]), self.mrnas[p].gene))
          
        for m, cds in mcs.items():
            self.mrnas[m].addCds(cds)
    
    def plotGene(self, name, 
                 canvas=None, 
                 domains={}, 
                 title=True, mrna_title=True, mrnas_anot={},
                 scale=False, legend=True,
                 auto_primary=False,
                 scale_strand=False, 
                 marks=[], mrnas_plot=[],
                 ri={}, se={}, mxe={}, a3ss={}, a5ss={}, 
                 mrnaH=50,
                 relative=[], interpro=None, filterDb=None, storeName=None):
        legendH = 10 if legend else 0
        gene = self.genes[name]
        if auto_primary:
            gene.setPrimary()
            
        if not interpro is None:
            domains.update(self.parseInterpro(interpro, [gene], filterDb=filterDb)[gene.name])
            
        h = (mrnaH * (len(gene.getMrnas()) if len(mrnas_plot) < 1 else len(mrnas_plot))) + legendH
        canvas = Canvas.create(name=storeName, h=h) if canvas is None else canvas.clip(h=h, box='plotGene')
        
        gene.plot(canvas, domains=domains, title=title, mrna_title=mrna_title, mrnas_anot=mrnas_anot,
                  scale=scale, scale_strand=scale_strand, marks=marks, mrnas_plot=mrnas_plot,
                  ri=ri, se=se, mxe=mxe, a3ss=a3ss, a5ss=a5ss, relative=relative)
        return canvas
    
    def default(self, title=True, mrna_title=True, mrnas_anot={}, scale=False, legend=True,
                scale_strand=False, auto_primary=False, marks=[], mrnas_plot=[], ri={}, se={}, mxe={}, a3ss={}, a5ss={}):
        return {
            'title': title, 'mrna_title': mrna_title, 'mrnas_anot': mrnas_anot,
            'scale': scale, 'scale_strand': scale_strand, 'legend': legend,
            'auto_primary': auto_primary, 'marks': marks, 'mrnas_plot': mrnas_plot,
            'ri': ri, 'se': se, 'mxe': mxe, 'a3ss': a3ss, 'a5ss': a5ss
        }
    
    def defaultFromGenes(self, genes, 
                         title=True, mrna_title=True, mrnas_anot={}, legend=True,
                         scale=False, scale_strand=False, auto_primary=False, marks=[], mrnas_plot=[]):
        return {gene: self.default(title=title, mrna_title=mrna_title, mrnas_anot=mrnas_anot,
                                   scale=scale, scale_strand=scale_strand, legend=legend, mrnas_plot=mrnas_plot,
                                   auto_primary=auto_primary, marks=marks) for gene in genes}
    
    def plotGenes(self, confs, canvas=None, mrnaH=50, interpro=None, relative=True, filterDb=None, storeName=None):
        genes = [self.genes[g] for g in confs]
        domains = {g: {} for g in confs}
        if not interpro is None:
            domains = self.parseInterpro(interpro, genes, filterDb=filterDb)
                                                                
        h = sum([mrnaH * (len(self.genes[g].getMrnas()) if len(i['mrnas_plot']) == 0 else len(i['mrnas_plot'])) + (10 if i['legend'] else 0) for g, i in confs.items()])
        canvas = Canvas.create(name=storeName, h=h) if canvas is None else canvas
        y = canvas.y
        for gene, c in confs.items():
            msplot = c['mrnas_plot']
            clip = canvas.clip(y=y, h=(mrnaH * (len(self.genes[gene].getMrnas()) if len(msplot) == 0 else len(msplot))) + (10 if c['legend'] else 0), box='Plot_'+gene)
            try:
                self.plotGene(gene, 
                     canvas=clip, 
                     domains=domains[gene], 
                     title=c['title'], mrna_title=c['mrna_title'], mrnas_anot=c['mrnas_anot'],
                     scale=c['scale'], scale_strand=c['scale_strand'], legend=c['legend'],
                     auto_primary=c['auto_primary'], marks=c['marks'], mrnas_plot=msplot,
                     ri=c['ri'], se=c['se'], mxe=c['mxe'], a3ss=c['a3ss'], a5ss=c['a5ss'], 
                     mrnaH = mrnaH,
                     relative=genes if relative else [])
            except Exception as e:
                print('ERROR in plot gene ' + gene)
                print(e)
            y = clip.bottom
        return canvas
    
    def parseInterpro(self, file, genes, filterDb=None):
        interpro = list(csv.reader(open(file), delimiter='\t'))
        colors = ['white/', 'yellow\\', 'cyan-', 'black|'] + Draw('tmp').listPatterns()
        doms = {}
        for gene in genes:
            ds = {m.name: [ Domain(x[5], int(x[6]), int(x[7]), entry=x[4], db=x[3]) for x in interpro if x[0] == m.name and len(x[5].strip()) > 1 and (filterDb is None or (x[3] in filterDb))] 
                  for m in gene.mrnas}
            nms = set()
            for d in ds.values():
                nms = nms.union([x.description for x in d])
            nms = list(nms)
            for x in ds.values():
                for d in x:
                    d.pattern = colors[nms.index(d.description)]
            doms[gene.name] = ds
        return doms
    
    def showGene(self, name, interpro=None, mrnaH=40, marks=[]):
        self.plotGene(name, interpro=interpro, mrna_title=True, mrnaH=mrnaH, 
                      scale=True, scale_strand=True, marks=marks).draw.show()
        
        
        
class Rmats:
    def __init__(self, folders, gff, fasta=None, pfasta=None, fdr=0.05, interpro=None):
        
        self.gff = gff
        self.fasta = pfasta if fasta is None else SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
        self.interpro = interpro
        self.files = {}
        
        for folder in folders:
            folder = folder + ('' if folder.endswith("/") else "/") 
            self.files[folder] = {x: self.parseTable(folder + x, fdr=fdr)
                                   for x in os.listdir(folder) if '.MATS.JC' in x and x.endswith('.txt')}
            print('Do diretorio "' + folder + '":')
            for file in self.files[folder]:
                print('{} events in {} file.'.format(len(self.files[folder][file]), file.replace(".MATS.", " ")[:-4]))
    
    def parseTable(self, file, fdr=0.05):
        data = list(csv.reader(open(file).readlines(), delimiter='\t'))
        if data[0][0] != 'ID':
            raise Exception('File %s not looks like with rMATS output file' % file)
        if 'SE.' in file:
            return self.parseSE(data[1:], file_name=file, fdr=fdr)
        if 'MXE.' in file:
            return self.parseMXE(data[1:], file_name=file, fdr=fdr)
        if 'RI.' in file:
            return self.parseRI(data[1:], file=file, fdr=fdr)
        if ('A3SS' in file) or ('A5SS' in file):
            return self.parseAXSS(data[1:], file, fdr=fdr)
        return data[1:]
    
    def validateLine(self, seq, gene, strand):
        gene = self.gff.genes[gene]
        if strand != gene.strand:
            raise Exception('Strands not math! MATS: %s and GFF: %s' % ('+' if strand else '-', '+' if gene.strand else '-'))
        contig = gene.contig
        if seq[3:] == contig:
            seq = seq[3:]
        elif seq != contig:
            raise Exception('Seq diff %s and %s' % (seq, contig))
                
        return seq, gene, strand
    
    def inCds(self, exon, cds):
        for c in cds:
            exon_dentro_cds = c.start <= exon[0] and c.end >= exon[0]
            cds_dentro_exon = exon[0] <= c.start and exon[1] >= c.end
            cds_comeca_exon = c.start >= exon[0] and c.start <= exon[1]
            cds_termina_exon = c.end >= exon[0] and c.end <= exon[1]
            if exon_dentro_cds or cds_dentro_exon or cds_comeca_exon or cds_termina_exon:
                return True
        return False
    
    def parseSE(self, data, file_name, fdr=0.05, fix={'1': ['mrnaID', 'exonID']}):
        ret = {}
        for id, gene, _, seq, strand, eS, eE, upsES, upsEE, downsES, downsEE, _,_,_,_,_,_,_,_, FDR, _,_,_ in data:
            
            if float(FDR) > fdr:
                continue
                
            exon = (int(eS)+1, int(eE))
            exonU = (int(upsES)+1, int(upsEE))
            exonD = (int(downsES)+1, int(downsEE))
            seq, gene, strand = self.validateLine(seq=seq, gene=gene, strand=strand == '+')
            mrnas = []
            mrnas_com_e = []
            mrnas_com_s = []
                        
            for m in gene.mrnas:
                exons = m.getExons()
                cds = m.getCds()
                e = [ e for e in exons if e.start == exon[0] and e.end == exon[1]]
                u = [ e for e in exons if e.start == exonU[0] and e.end == exonU[1]]
                d = [ e for e in exons if e.start == exonD[0] and e.end == exonD[1]]
                if len(e) > 0:
                    mrnas_com_e.append((m.name, e[0].name, gene.name, self.inCds(exon, m.getCds())))
                if len(u) == len(d) == 1:
                    mrnas_com_s.append((m.name, None, gene.name, self.inCds(exon, m.getCds())))
                if len(e) == len(u) == len(d) == 1:
                    mrnas.append((m.name, e[0].name, gene.name, self.inCds(exon, m.getCds())))
                if (id in fix) and (m.name == fix[id][0]) and (fix[id][1] in [x.name for x in exons]):
                    mrnas.append((m.name, fix[id][1], gene.name, self.inCds(exon, m.getCds())))
                    break

            
            if len(mrnas) > 0:
                ret[id] = max([m for m in mrnas], key= lambda x: -(self.gff.mrnas[x[0]].size + self.gff.mrnas[x[0]].cdsSize()))
            else:
                if len(mrnas_com_e) > 0:
                    mv = [m for m in mrnas_com_e if len([ e for e in self.gff.mrnas[m[0]].getExons() if e.name == m[1]][0].getIntrons()) > 1]
                    if len(mv) > 0:
                        mrnas_com_e = mv
                    else:
                        print('WARN: has Exon skiping not surrouding by introns in event {}, mrnas {}'.format(id, mrnas_com_e))
                        continue
                        
                    ret[id] = max([m for m in mrnas_com_e], key= lambda x: -(self.gff.mrnas[x[0]].size + self.gff.mrnas[x[0]].cdsSize()))
                    continue
                self.gff.showGene(gene.name, marks=exon+exonU+exonD)
                print('{} Fix event: {} gene: {} mrnas: {}'.format(file_name, id, gene.name, [m.name for m in gene.mrnas]))
                raise Exception('Not anything mRNA of gene {} math whit evt [{}] upstream: {} exon: {} downstream: {}'.format(gene, id, exonU, exon,  exonD))
        return ret
    
    def parseMXE(self, data, file_name, fdr=0.05, fix={'1': ['mrnaID', 'exonAID', 'exonBID']}):
        ret = {}
        for id, gene, _,  seq, strand, aES, aEE, bES, bEE, uES, uEE, dES, dEE, _,_,_,_,_,_,_,_, FDR, _,_,_ in data:

            if float(FDR) > fdr:
                continue

            exonA = (int(aES)+1, int(aEE))
            exonB = (int(bES)+1, int(bEE))
            exonU = (int(uES)+1, int(uEE))
            exonD = (int(dES)+1, int(dEE))
            seq, gene, strand = self.validateLine(seq=seq, gene=gene, strand=strand == '+')
            mrnas = []
            mrnas_com_e = []
            mrnas_com_e_a = []
            mrnas_com_e_b = []

            for m in gene.mrnas:
                exons = m.getExons()
                cds = m.getCds()
                ea = [ e for e in exons if e.start == exonA[0] and e.end == exonA[1]]
                eb = [ e for e in exons if e.start == exonB[0] and e.end == exonB[1]]
                u = [ e for e in exons if e.start == exonU[0] and e.end == exonU[1]]
                d = [ e for e in exons if e.start == exonD[0] and e.end == exonD[1]]
                hasCds = (len(ea) > 0 and self.inCds(exonA, m.getCds())) or (len(eb) > 0 and self.inCds(exonB, m.getCds()))
                if len(ea) > 0 and len(eb) > 0:
                    mrnas_com_e.append((m.name, ea[0].name, eb[0].name, gene.name, hasCds))
                if len(ea) == len(eb) == len(u) == len(d) == 1:
                    mrnas.append((m.name, ea[0].name, eb[0].name, gene.name, hasCds))
                if len(ea) == len(u) == len(d) == 1:
                    mrnas_com_e_a.append((m, ea[0], hasCds, u[0], d[0]))
                if len(eb) == len(u) == len(d) == 1:
                    mrnas_com_e_b.append((m, eb[0], hasCds))
                if (id in fix) and (m.name == fix[id][0]) and (fix[id][1] in [x.name for x in exons]) and (fix[id][2] in [x.name for x in exons]):
                    mrnas.append((m.name, fix[id][1], fix[id][2], gene.name, hasCds))
                    break

            if len(mrnas) > 0:
                ret[id] = max([m for m in mrnas], key= lambda x: -(self.gff.mrnas[x[0]].size + self.gff.mrnas[x[0]].cdsSize()))
            else:
                if len(mrnas_com_e) > 0:
                    ret[id] = max([m for m in mrnas_com_e], key= lambda x: -(self.gff.mrnas[x[0]].size + self.gff.mrnas[x[0]].cdsSize()))
                    continue
                if len(mrnas_com_e_a) > 0  and len(mrnas_com_e_b) > 0:
                    ma, ea, ahasCds, u, d = mrnas_com_e_a[0]
                    mb, eb, bhasCds = mrnas_com_e_b[0]
                    name  = 'MXE' + ma.name + '-' + mb.name
                    if not name in self.gff.mrnas:
                        newMRNA = Mrna(name, ma.start, ma.end, gene)
                        newMRNA.addExons([Exon(x.name, x.start, x.end, gene) for x in [u, ea, eb, d]])
                        self.gff.mrnas[newMRNA.name] = newMRNA
                    ret[id] = (name, ea.name, eb.name, gene.name, ahasCds or bhasCds)
                    print('mRNA created: ' + name + ' for gene: ' + gene.name)
                    continue
                    
                self.gff.showGene(gene.name, marks=exonA+exonB+exonU+exonD)
                print('{} Fix event: {} gene: {} mrnas: {}'.format(file_name, id, gene.name, [m.name for m in gene.mrnas]))
                raise Exception('Not anything mRNA of gene {} math whit evt [{}] upstream: {} exonA: {} exonB: {} downstream: {}'.format(gene, id, exonU, exonA, exonB, exonD))
        return ret

    
    def parseAXSS(self, file, file_name, fdr=0.05):
        data, erros = {}, []
        for l in file:
            
            if float(l[19]) > fdr:
                continue
                
            id = l[0]
            contig, gene, strand = self.validateLine(seq=l[3], gene=l[1], strand=l[4] == '+')
                
            longES, longEE, shortES, shortEE, flankingES, flankingEE = int(l[5])+1, int(l[6]), int(l[7])+1, int(l[8]), int(l[9])+1, int(l[10])
            pbMudanca = (shortEE if gene.strand else shortES) if "5" in file_name else (shortES if gene.strand else shortEE)
            mrnas_flankE = set()
            mrnas_longE = set()
            mrnas_shortE = set()
            m2e = {}
            for mrna in gene.mrnas:
                exons = mrna.getExons()
                exon_flanking = [e for e in exons if e.start == flankingES and e.end == flankingEE]
                exon_long = [e for e in exons if e.start == longES and e.end == longEE]
                exon_short = [e for e in exons if e.start == shortES and e.end == shortEE]
                if len(exon_flanking) == 1:
                    mrnas_flankE.add(mrna.name)
                    m2e[mrna.name + 'F'] = exon_flanking[0].name
                if len(exon_long) == 1:
                    mrnas_longE.add(mrna.name)
                    m2e[mrna.name + 'L'] = exon_long[0].name
                    m2e[mrna.name + 'A'] = exon_long[0].name
                if len(exon_short) == 1:
                    mrnas_shortE.add(mrna.name)
                    m2e[mrna.name + 'S'] = exon_short[0].name
                    m2e[mrna.name + 'A'] = exon_short[0].name
                if len(exon_long) == 1 and len(exon_short) == 1 and exon_long[0].name != exon_short[0].name:
                    raise Exception('Exon long and short are different from same mRNA!', id, exon_long[0].name, exon_short[0].name)
                    
            mrnas = mrnas_flankE.intersection(mrnas_longE.union(mrnas_shortE))
            if len(mrnas) > 0:
                data[id] = (contig, gene.name, [(m, m2e[m + 'F'], m2e[m + 'A'], any([c.start <= pbMudanca and c.end >= pbMudanca 
                                                                                     for c in self.gff.mrnas[m].getCds()])) for m in mrnas],
                            max([self.gff.mrnas[m] for m in mrnas_longE], key=lambda e: -(e.size + e.cdsSize())).name, ## maior mrna com o long
                            list(mrnas_shortE),## tem mrna com o short
                            pbMudanca,
                            ((longEE - pbMudanca) if gene.strand else (pbMudanca - longES)) if "5" in file_name else ((pbMudanca - longES) if gene.strand else (longEE - pbMudanca))
                           )
            else:
                erros.append([(flankingES, flankingEE), (longES, longEE), (shortES, shortEE), l])
        if len(erros) > 0:
            raise Exception('Has erros in %d ids' % len(erros))
        return data
    
    def parseRI(self, data, file, fdr=0.05):
        dataR = {}
        for id, gene, _, seq, strand, riS, riE, uES, uEE, dES, dEE, _,_,_,_,_,_,_,_, FDR, _,_,_ in  data:
            
            if float(FDR) > fdr:
                continue
            
            longE = (int(riS)+1, int(riE))
            intron = (int(uEE)+1, int(dES))
            exonU = (int(uES)+1, int(uEE))
            exonD = (int(dES)+1, int(dEE))
            seq, gene, strand = self.validateLine(seq=seq, gene=gene, strand=strand == '+')

            def fromMrnas(marg=0):
                dfs = {}
                def comp(x, m):
                    y = abs(x)
                    if y <= marg:
                        if m in dfs:
                            dfs[m] = max(dfs[m], y)
                        else:
                            dfs[m] = y
                    return y
                for m in gene.mrnas:
                    n = m.name
                    exons = m.getExons()
                    introns = m.getIntrons()
                    cds = m.getCds()
                    i = [ i for i in introns if ((comp(i.start - intron[0], n) <= marg, n) and 
                         (comp(i.end - intron[1], n) <= marg)) or 
                         ((comp(i.start - intron[0], n) <= marg) and (comp(i.end - intron[1], n) <= marg))]
                    u = [ e for e in exons if (comp(e.start - exonU[0], n) <= marg) and (comp(e.end - exonU[1], n) <= marg)]
                    d = [ e for e in exons if (comp(e.start - exonD[0], n) <= marg) and (comp(e.end - exonD[1], n) <= marg)]
                    incds = ((len(u) > 0) and self.inCds((u[0].start, u[0].end), cds)) or ((len(d) > 0) and self.inCds((d[0].start, d[0].end), cds))
                    if (len(i) > 0) and (len(u) == len(d) == 1):
                        return m.name, [x.name for x in i], sum([x.size for x in i]), gene.name, incds, max(dfs.values())

            mrna = fromMrnas()
            if mrna is None:
                mrna = fromMrnas(50)
                print('WARN: %s event %s processed with %d PB difference' % (file, id, mrna[-1]))
            if mrna is None:
                raise Exception('Event %s mistake gff3 data' % id)
            else:
                dataR[id] = mrna[:-1]
        return {k: (v[0], self.findSC(v), v[2], v[3], v[4]) for k, v in dataR.items()}
    
    def findSC(self, event):
        mrna_name, intron_names, _, _, _ = event
        mrna = self.gff.mrnas[mrna_name]
        introns = {i.name: i for i in mrna.getIntrons() if i.name in intron_names}
        intronSeq = {}
        for i in introns.values():
            intronSeq[i.name] = self.fasta[mrna.contig].seq[i.start-1:i.end]
            if mrna.strand:
                intronSeq[i.name] = str(intronSeq[i.name])
            else:
                intronSeq[i.name] = str(intronSeq[i.name].reverse_complement())

        
        parts = mrna.getCds() + list(introns.values())
        partsO = mrna.getCds() + list(introns.values())
        coords = sorted([(x.start, x.end) for x in parts], key=lambda e: e[0] if mrna.strand else -e[1])
        coordsO = sorted([(x.start, x.end) for x in partsO], key=lambda e: e[0] if mrna.strand else -e[1])
        
        coords_index = [(k, coords.index((x.start, x.end))) for k, x in introns.items()]
        coords_exons = [x for x in range(len(coords)) if not x in [y[1] for y in coords_index]]
        valid_introns = [x[0] for x in coords_index if ((x[1]+1) in coords_exons) and  ((x[1]-1) in coords_exons) ]
        
        if len(valid_introns) < 1:
            return {i: 1 for i in intron_names}
        
        contig = self.fasta[mrna.contig].seq
        cds = ''.join([str(contig[p[0]-1:p[1]] if mrna.strand else contig[p[0]-1:p[1]].reverse_complement()) for p in coords])
        cdsO = ''.join([str(contig[p[0]-1:p[1]] if mrna.strand else contig[p[0]-1:p[1]].reverse_complement()) for p in coordsO])
        aa = str(Seq.Seq(cds[0: 3* (len(cds)//3)]).translate())
        aaO = str(Seq.Seq(cds[0: 3* (len(cds)//3)]).translate())
            
        if (not '*' in aa[:-1]) or (('*' in aaO) and ((aaO.index('*') < aa.index('*')))):
            return {i: 1 for i in intron_names}
        
        aaR = ''.join([x*3 for x in aa])
        cont = 0
        for a, b in coords:
            size = b - a + 1
            s = aaR[cont:cont+size]
            if '*' in s:
                return {i: a+s.index('*') for i in intron_names}
            cont += size
        raise Exception('program error')
    
    
    def plotAXSS(self, events, id, x, interpro=None, show=True, filterDb=None, storeName=None):
        interpro = self.interpro if interpro is None else interpro
        contig, gene, mrnas, mrna, mrnasWithShort, pbMudanca, asSize = events[id]
        data = [x for x in mrnas if x[0] == mrna][0]
        plot = self.gff.plotGene(gene, storeName=storeName, interpro=interpro, filterDb=filterDb, 
                                 mrna_title=True, mrnaH= 80,
                                 mrnas_anot={m.name: ('A%sSS[%s]' % (x, id)) for m in self.gff.genes[gene].mrnas},
                            mrnas_plot=[mrna]+mrnasWithShort, scale=True, scale_strand=True, marks=[pbMudanca], 
                 a3ss={mrna: {data[2]: pbMudanca}} if x == '3' else {}, a5ss={mrna: {data[2]: pbMudanca}} if x == '5' else {})
        if show:
            plot.draw.show()
            
        if not storeName is None:
            plot.draw.store()
            
        return plot
    
    def plotRI(self, events, id, interpro=None, show=True, filterDb=None, storeName=None):
        mrna, introns, _, gene, incds = events[id]
        interpro = self.interpro if interpro is None else interpro
        plot = self.gff.plotGene(gene, storeName=storeName, 
                                 interpro=interpro, filterDb=filterDb, scale=True, scale_strand=True, mrnaH= 80,
                                 mrnas_anot={mrna: ('RI[%s]' % id)}, mrnas_plot=[mrna], ri={mrna: introns})
        if show:
            plot.draw.show()
            
        if not storeName is None:
            plot.draw.store()
            
        return plot
    
    def plotSE(self, events, id, interpro=None, show=True, filterDb=None, storeName=None):
        mrna, exon, gene, _ = events[id]
        interpro = self.interpro if interpro is None else interpro
        plot = self.gff.plotGene(gene, storeName=storeName, 
                                 interpro=interpro, filterDb=filterDb, scale=True, scale_strand=True, mrnaH= 80,
                                 mrnas_anot={mrna: ('SE[%s]' % id)}, mrnas_plot=[mrna], se={mrna: [exon]})
        if show:
            plot.draw.show()
            
        if not storeName is None:
            plot.draw.store()
            
        return plot
    
    def plotMXE(self, events, id, interpro=None, show=True, filterDb=None, storeName=None):
        mrna, exonA, exonB, gene, _ = events[id]
        interpro = self.interpro if interpro is None else interpro
        exs = [exonA, exonB]
        exons = [e for e in self.gff.mrnas[mrna].getExons() if e.name in exs]
        mrnas = [mrna] if not mrna.startswith('MXE') or any([not m in self.gff.mrnas for m in mrna[3:].split('-')]) else ([mrna] + mrna[3:].split('-'))
        plot = self.gff.plotGene(gene, storeName=storeName, 
                                 interpro=interpro, filterDb=filterDb, scale=True, scale_strand=True, mrnaH= 80,
                                 marks = [int((e.start + e.end) / 2) for e in exons] if len(mrnas) > 1 else [],
                                 mrnas_anot={m: ('MXE[%s]' % id) for m in mrnas}, mrnas_plot=mrnas, mxe={mrna: exs})
        if show:
            plot.draw.show()
            
        if not storeName is None:
            plot.draw.store()
            
        return plot
    
    def plotAXSSevents(self, events, x, ids=None, interpro=None, 
                       filterDb=None, show=False, mrnaH=50, relative=False, sids=None, storeName=None, skip=False):
        sids = [] if sids is None else sids
        interpro = self.interpro if interpro is None else interpro
        data_plot = {}
        ids = list(events.keys()) if ids is None else ids
        
        g_id = {}
        
        for id in ids:
            if id in sids or not id in events:
                continue
            contig, gene, mrnas, mrna, mrnasWithShort, pbMudanca, asSize = events[id]
            
            if gene in g_id:
                if skip:
                    if VERBOSE:
                        print('INFO: gene %s skiped event %s of A%sSS' % (gene, x, id))
                    continue
                raise Exception('Unsuported plot mupliples events for same gene: {} has events {}, {}'.format(gene, g_id[gene], id))
            else:
                g_id[gene] = id
            
            data = [x for x in mrnas if x[0] == mrna]
            if len(data) < 1:
                data = mrnas[0]
                mrna = data[0]
            else:
                data = data[0]
            data_plot[gene] = self.gff.default(title=False, mrna_title=True, mrnas_anot={m.name: ('A%sSS[%s]' % (x, id)) for m in self.gff.genes[gene].mrnas},
                                          scale=False, scale_strand=False, 
                                          auto_primary=False, marks=[], mrnas_plot=[mrna]+mrnasWithShort,
                                          a3ss={mrna: {data[2]: pbMudanca}} if '3' == x else {}, a5ss={mrna: {data[2]: pbMudanca}} if '5' == x else {})
            
        plot = self.gff.plotGenes(data_plot, canvas=None, mrnaH=mrnaH, 
                                  interpro=interpro, filterDb=filterDb, relative=relative, storeName=storeName)
        if show:
            plot.draw.show()
            
        if not storeName is None:
            plot.draw.store()
            
        return plot
            
    def plotSEevents(self, events, ids=None, interpro=None, 
                     filterDb=None, show=False, mrnaH=50, relative=False, sids=None, storeName=None, skip=False):
        sids = [] if sids is None else sids
        interpro = self.interpro if interpro is None else interpro
        data_plot = {}
        ids = list(events.keys()) if ids is None else ids
        g_id = {}
        for id in ids:
            if id in sids or not id in events:
                continue
            mrna, exon, gene, _ = events[id]
            if gene in g_id:
                if skip:
                    if VERBOSE:
                        print('INFO: gene %s skiped event %s of SE' % (gene, id))
                    continue
                raise Exception('Unsuported plot mupliples events for same gene: {} has events {}, {}'.format(gene, g_id[gene], id))
            else:
                g_id[gene] = id
            
            data_plot[gene] = self.gff.default(title=False, mrna_title=True, mrnas_anot={mrna: ('SE[%s]' % id)},
                                          scale=False, scale_strand=False, mrnas_plot=[mrna], se = {mrna: [exon]})
            
        plot = self.gff.plotGenes(data_plot, canvas=None, mrnaH=mrnaH, 
                                  interpro=interpro, filterDb=filterDb, relative=relative, storeName=storeName)
        if show:
            plot.draw.show()
            
        if not storeName is None:
            plot.draw.store()
            
        return plot
    
    def plotRIevents(self, events, ids=None, interpro=None, 
                     filterDb=None, show=False, mrnaH=50, relative=False, sids=None, storeName=None, skip=False):
        sids = [] if sids is None else sids
        interpro = self.interpro if interpro is None else interpro
        data_plot = {}
        ids = list(events.keys()) if ids is None else ids
        g_id = {}
        for id in ids:
            if id in sids or not id in events:
                continue
            mrna, introns, _, gene, incds = events[id]
            if gene in g_id:
                if skip:
                    if VERBOSE:
                        print('INFO: gene %s skiped event %s of RI' % (gene, id))
                    continue
                raise Exception('Unsuported plot mupliples events for same gene: {} has events {}, {}'.format(gene, g_id[gene], id))
            else:
                g_id[gene] = id
            
            data_plot[gene] = self.gff.default(title=False, mrna_title=True, mrnas_anot={mrna: ('RI[%s]' % id)},
                                          scale=False, scale_strand=False, mrnas_plot=[mrna], ri={mrna: introns})
            
        plot = self.gff.plotGenes(data_plot, canvas=None, mrnaH=mrnaH, 
                                  interpro=interpro, filterDb=filterDb, relative=relative, storeName=storeName)
        if show:
            plot.draw.show()
            
        if not storeName is None:
            plot.draw.store()
            
        return plot 
    
    def plotMXEevents(self, events, ids=None, interpro=None, 
                      filterDb=None, show=False, mrnaH=50, relative=False, sids=None, storeName=None, skip=False):
        sids = [] if sids is None else sids
        interpro = self.interpro if interpro is None else interpro
        data_plot = {}
        ids = list(events.keys()) if ids is None else ids
        g_id = {}
        for id in ids:
            if id in sids or not id in events:
                continue
            mrna, exonA, exonB, gene, _ = events[id]
            exs = [exonA, exonB]
            exons = [e for e in self.gff.mrnas[mrna].getExons() if e.name in exs]
            mrnas = [mrna] if not mrna.startswith('MXE') or any([not m in self.gff.mrnas for m in mrna[3:].split('-')]) else ([mrna] + mrna[3:].split('-'))
  
            if gene in g_id:
                if skip:
                    if VERBOSE:
                        print('INFO: gene %s skiped event %s of MXE' % (gene, id))
                    continue
                raise Exception('Unsuported plot mupliples events for same gene: {} has events {}, {}'.format(gene, g_id[gene], id))
            else:
                g_id[gene] = id
            
            data_plot[gene] = self.gff.default(title=False, mrna_title=True, mrnas_anot={m: ('MXE[%s]' % id) for m in mrnas},
                                               marks = [int((e.start + e.end) / 2) for e in exons] if len(mrnas) > 1 else [],
                                               scale=False, scale_strand=False, mrnas_plot=mrnas, mxe={mrna: exs})
            
        plot = self.gff.plotGenes(data_plot, canvas=None, mrnaH=mrnaH, 
                                  interpro=interpro, filterDb=filterDb, relative=relative, storeName=storeName)
        if show:
            plot.draw.show()
            
        if not storeName is None:
            plot.draw.store()
            
        return plot 
    
    def plotAllEvents(self, mrnaH=50, relative=False, filterDb=None, skip=False, ids=None, sids=None):
        filterDb=None if filterDb is None else filterDb.split(',')
        for folder, files in self.files.items():
            for file, data in files.items():
                fname = folder.replace('/', '_').replace('.', '_') + '_' + file + '.svg'
                evts = list(data) if ((ids is None) or (len(ids) < 1)) else [x for x in data if x in ids]
                evts = evts if ((sids is None) or (len(sids) < 1)) else [x for x in evts if not x in sids]
                if len(evts) < 1:
                    print('WARN: %s has 0 events, skipping ...' % fname)
                    continue
                if VERBOSE:
                    print(fname)
                if file.startswith('A3SS'):
                    self.plotAXSSevents(data, '3', show=False, filterDb=filterDb, storeName=fname, relative=relative, mrnaH=mrnaH, skip=skip, ids=ids, sids=sids)
                elif file.startswith('A5SS'):
                    self.plotAXSSevents(data, '5', show=False, filterDb=filterDb, storeName=fname, relative=relative, mrnaH=mrnaH, skip=skip, ids=ids, sids=sids)
                elif file.startswith('SE'):
                    self.plotSEevents(data, show=False, filterDb=filterDb, storeName=fname, relative=relative, mrnaH=mrnaH, skip=skip, ids=ids, sids=sids)
                elif file.startswith('MXE'):
                    self.plotMXEevents(data, show=False, filterDb=filterDb, storeName=fname, relative=relative, mrnaH=mrnaH, skip=skip, ids=ids, sids=sids)
                elif file.startswith('RI'):
                    self.plotRIevents(data, show=False, filterDb=filterDb, storeName=fname, relative=relative, mrnaH=mrnaH, skip=skip, ids=ids, sids=sids)
                else:
                    print('WARN: unknown file type ' + file)

                    
if VERBOSE:
    print('events dirs: {}'.format(dirs))
    print('events ids: {}'.format(options.events))
    print('events exclude: {}'.format(options.xevents))
    print('fasta: {}'.format(options.fasta))
    print('gff: {}'.format(options.gff))
    print('interpro: {}'.format(options.interpro))
    print('interpro dbs: {}'.format(options.dbs))
    print('relative: {}'.format(options.relative))
    print('skip: {}'.format(options.skip))
    print('mrna: {}'.format(options.mrna))
    print('fdr: {}'.format(options.fdr))                    
                            
                    
print('[1/4] loading gff3 %s ...' % options.gff)
gff3 = Gff(options.gff)
print('[2/4] loading fasta %s ...' % options.fasta)
pf = SeqIO.to_dict(SeqIO.parse(options.fasta, 'fasta'))
print('[3/4] loading MATS files ...')
rmats = Rmats(dirs, gff3, interpro=options.interpro, pfasta=pf) 

print('[4/4] plotting events ...')
ids = None if options.events is None else options.events.split(',')
sids = None if options.xevents is None else options.xevents.split(',')
rmats.plotAllEvents(skip=options.skip, relative=options.relative, mrnaH=options.mrna, filterDb=options.dbs, ids=ids, sids=sids)

print('finished.')
