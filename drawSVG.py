#!/usr/bin/env python
# coding: utf-8

# In[287]:


from IPython.display import display, SVG
import svgwrite

class Draw:
    def __init__(self, file, width=600, height=200, marginX=10, marginY=10):
        
        self.file = file
        
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
        
    def show(self):
        self.drawObj.save()
        display(SVG(filename=self.file))

    def listPatterns(self):
        return ' '.join([' '.join([c+'|', c+'-', c+'\\', c+'/']) for c in self.PATTERNS_COLORS]).split()
    
    def fillPattern(self, name):
        if not name in self.listPatterns():
            raise Exception('Options are: ' + ', '.join(self.listPatterns()))
        
        name = name.replace('\\', '@')
        if not name in self.added_patterns:
            c = name[:-1]
            o = name[-1]
            s = (3, 3) if (o in ['-', '|']) else (100, 100)
            pattern = self.drawObj.pattern(size=s, id='pattern-%s' % name, patternUnits="userSpaceOnUse")
            
            if (o in ['-', '|']):
                tLine = ((0, 2), (3, 2)) if o == '-' else ((2, 3), (2, 0))
                pattern.add(self.drawObj.line(start=tLine[0], end=tLine[1], stroke=c, stroke_width=0.5))
            else:
                a = 100 if o == '/' else 0
                b = 0 if o == '/' else 100
                for i in range(-100, 200, 4):
                    pattern.add(self.drawObj.line(start=(i, a), end=(i+100, b), stroke=c, stroke_width=0.5))
            
            self.drawObj.defs.add(pattern)
            self.added_patterns.append(name)
        return "url(#pattern-%s)" % name
    
    def box(self, name=None):
        name = ('box%s' % len(self.boxes)) if name is None else name
        return self.drawObj.add(self.drawObj.g(id=name))
    
    def square(self, x, y, w, h, c, box=None):
        box = self.mainBox if box is None else box
        box.add(self.drawObj.rect(insert=(x, y), size=(w, h), fill=c))
    
    def text(self, text, x, y, box=None, italic=False, family='Times', size=12, bold=False):
        box = self.mainBox if box is None else box
        self.mainBox.add(self.drawObj.text(
            text, insert=(x, y), 
            font_style="italic" if italic else 'normal', 
            font_family=family, 
            font_size=size, 
            font_weight= 'bold' if bold else 'normal'))
        
    def star(self, x, y, box=None, fill="gray", stroke="white", r=5):
        box = self.mainBox if box is None else box
        STAR8 = "0 8 L 8 12 L 7 2 L 14 -4 L 4 -6 L 0 -15 L -4 -6 L -14 -4 L -7 2 L -8 12 L 0 8"
        STAR5 = "0 5 L 5 8 L 4 1 L 9 -3 L 2 -4 L 0 -10 L -2 -4 L -9 -3 L -4 1 L -5 8 L 0 5"
        STAR = STAR8 if r == 8 else STAR5
        d = 'M ' + ' L '.join([str(int(p.split()[0])+x) + ' ' + str(int(p.split()[1])+y) for p in STAR.split(" L ")]) + ' z'
        box.add(self.drawObj.path(d=d, fill=fill, stroke=stroke))

    def seta(self, x1, x2, y, box=None, color='black', h=20, start=False, pS=4):
        box = self.mainBox if box is None else box
        box.add(self.drawObj.path(d='M {0} {1} C {0} {2} {3} {2} {3} {1}'.format(x1, y, y+h, x2), 
                          fill="none", stroke=color, stroke_width=pS//2))
        xs = x1 if start else x2
        box.add(self.drawObj.path(
            d='M {0} {1} L {2} {1} L {3} {4} z'.format(xs-pS, y, xs+pS, xs, y-pS), stroke=color, fill=color))



# In[291]:


draw = Draw('test.svg')
draw.square(30, 30, 100, 150, draw.fillPattern('red\\'))
draw.square(30, 30, 100, 150, draw.fillPattern('green/'))
#draw.square(30, 30, 100, 150, draw.fillPattern('blue|'))
#draw.square(30, 30, 100, 150, draw.fillPattern('black-'))
draw.text('Miquéias', 30, 90)
draw.star(80, 90)
draw.seta(30, 80, 100)
draw.show()


# In[ ]:




