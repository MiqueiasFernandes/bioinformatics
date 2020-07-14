#!/usr/bin/env python
# coding: utf-8

# Copyright(c) 2020 - Miquéias Fernandes <bio@mikeias.net>

import os
import sys
import csv
from venn import venn
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#################################################
usage = "python3 plotAS.py dir1,dir2 lab1,lab2"
#################################################

try: _, dirs, labs=sys.argv
except: sys.exit("Correct usage is:\n"+usage)
    
FDR = 0.05

print('FDR: %f' % FDR)
    

def importFromFolder(folder, fdr=0.05):
    print('loading folder %s ...' % folder)
    folder = folder if folder.endswith("/") else (folder+'/')
    files = [f for f in os.listdir(folder) if '.MATS.JC' in f and f.endswith('.txt')]
    if len([y for y in files if '.JC.' in y]) > 0 and len([y for y in files if '.JCEC.' in y]) > 0:
        raise Exception('Has more than one MATS tipe in folder! choose just only one (JC or JCEC)')
    def getData(f):
        print('Load file: %s ...' % f)
        fp = list(csv.reader(open(folder + f).readlines(), delimiter='\t'))
        idx = fp[0].index('FDR') if 'FDR' in fp[0] else fp[0].index('fdr')
        res = [x[1] for x in fp[1:] if float(x[idx]) <= fdr] 
        open('genes_com_fdr.txt', 'a').writelines([f + '\t' + x + '\n' for x in set(res)])
        return res
    return {f.split('.')[0]: getData(f) for f in files}


def parseData(folders):
    folders = [x if x.endswith('/') else (x+'/') for x in folders]
    data = {f: {t: [] for t in ['A3SS', 'A5SS', 'SE', 'MXE', 'RI']} for f in folders}
    for folder in folders:
        data[folder].update(importFromFolder(folder, fdr=FDR))
    return data
        

def compareEvents(data, classes, cmap=["viridis", "plasma", "cool", list("rgb"), "Set1", 'Pastel1', 'Pastel2', 'Paired', 'Accent',
                        'Dark2', 'Set1', 'Set2', 'Set3',
                        'tab10', 'tab20', 'tab20b', 'tab20c', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'], 
                  file='vennCondition.pdf'):
    types = {}
    for t in data:
        orig = t
        for c in classes:
            t = t.replace(c, '_')
        if t in types:
            types[t].append(orig)
        else:
            types[t] = [orig]
        
    cindex = list(types)
    cmap = ["viridis"] if cmap is None else ([cmap] if type(cmap) == str else cmap)
    if len(cmap) == 1:
        cmap = [cmap[0] for i in range(len(cindex))]
    elif len(cmap) < len(cindex):
        raise Exception('Has less palets (%d) than types (%d)!' % (len(cmap), len(cindex)))
        
    _, axs = plt.subplots(ncols=len(types), nrows=len(classes), figsize=(7 * len(types), 8*len(classes)))
    if len(classes) < 2:
        axs = [axs]
    if len(types) < 2:
        axs = [axs]
    a = 0
    for classe in classes:
        b = 0
        for t, labels in types.items():
            for label in labels:
                if classe in label:
                    ax = axs[a][b]
                    ax.set_title(label)
                    b += 1
                    if len(data[label]) < 1 or sum([len(x) for x in data[label].values()]) < 1:
                        print(classe, label, ' ZEROU ')
                        continue
                    venn({k: set(v) for k, v in data[label].items()}, fontsize=8, cmap=cmap[cindex.index(t)], ax=ax)
        a += 1
            
    plt.savefig(file, dpi=150)
    
      
def compareGenes(data, classes, file='vennEvents.pdf'):
    types = {}
    evts = set()
    for t in data:
        orig = t
        evts = evts.union(data[t])
        for c in classes:
            t = t.replace(c, '')
        if t in types:
            types[t].append(orig)
        else:
            types[t] = [orig]
    
    print('Types found: "%s"' % ','.join(types)) 
    evts = list(evts)    
    cmap = ["cool", list("rgb"), "plasma", "viridis", "Set1"]
    cs = len(types) + (len(classes) if len(types) > 1 else 0)
    _, axs = plt.subplots(ncols=cs, nrows=5, figsize=(7 * cs, 20))
    
    if len(types) < 2:
        axs = [[axs[i]] for i in range(5)]
        
    col = 0
    for _, dirs in types.items():
        for evt in evts:
            fs = {d: set(data[d][evt]) for d in dirs}
            ax = axs[evts.index(evt)][col]
            ax.set_title(evt)
            if len(fs) > 1 and sum([len(x) for x in fs.values()]) > 0:
                venn(fs, fontsize=12, cmap=cmap[evts.index(evt)], legend_loc="lower right", ax=ax)
        col += 1
        
    if len(types) > 1:
        for classe in classes:
            for evt in evts:
                fs = {d: set(data[d][evt]) for d in [x for x in data if classe in x]}
                ax = axs[evts.index(evt)][col]
                ax.set_title(evt)
                if len(fs) > 1 and sum([len(x) for x in fs.values()]) > 0:
                    venn(fs, fontsize=12, cmap=cmap[evts.index(evt)], legend_loc="lower right", ax=ax)
            col += 1
            
    plt.savefig(file, dpi=150)

    
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
    

def plotHeatMap(data, file='asHeatMap.pdf'):
    evts = ['A3SS', 'A5SS', 'SE', 'MXE', 'RI']
    dirs = list(data)
    byE = np.array([[len(data[f][e]) for e  in evts] for f in dirs])
    byG = np.array([[len(set(data[f][e])) for e  in evts] for f in dirs])

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 4))

    im, cbar = heatmap(byE, dirs, evts, ax=axL, cmap="magma_r", cbarlabel="Events")
    texts = annotate_heatmap(im, valfmt="{x:d}")
    im, cbar = heatmap(byG, dirs, evts, ax=axR, cmap="magma_r", cbarlabel="Genes")
    texts = annotate_heatmap(im, valfmt="{x:d}")
    plt.savefig(file, dpi=150)
    
def compareNumEvents(data, classes, file='quantAS.pdf'):
    types = {}
    evts = set()
    for t in data:
        orig = t
        evts = evts.union(data[t])
        for c in classes:
            t = t.replace(c, '')
        if t in types:
            types[t].append(orig)
        else:
            types[t] = [orig]
    
    def autolabel(rects, ax):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')
            

    evts = list(evts) 
    _, axs = plt.subplots(ncols=len(types), nrows=5, figsize=(7 * len(types), 25))
    
    
    maxY = max([max([len(x) for x in v.values()]) for v in data.values()])
    col = 0
    width = 0.35
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    marks = ['o', 's', 'P', '*', 'X']
    for t, dirs in types.items():
        n_events = {d.replace(t, ''): [len(data[d][e]) for e in evts] for d in dirs}
        n_events2 = {d.replace(t, ''): [len(set(data[d][e])) for e in evts] for d in dirs}
        x = np.arange(len(evts))
        ax = axs[0][col] if len(types) > 1 else axs[0]
        ax2 = axs[1][col] if len(types) > 1 else axs[1]
        ax3 = axs[2][col] if len(types) > 1 else axs[2]
        ax4 = axs[3][col] if len(types) > 1 else axs[3]
        ax5 = axs[4][col] if len(types) > 1 else axs[4]
        lgT = width * len(n_events)
        c = 0
        for k, v in n_events.items():
            r1 = ax.bar(x - (lgT / 2) + (width * c) + (width / 2), v, width, label=k)
            r2 = ax2.bar(x - (lgT / 2) + (width * c) + (width / 2), n_events2[k], width, label=k)
            c += 1
            autolabel(r1, ax)
            autolabel(r2, ax2)
        ax.set_xticklabels(evts)
        ax.set_xticks(x)
        ax.set_ylabel('Events')
        ax.set_title('AS Events in ' + t)
        ax.legend()
        ax.set_ylim([0, maxY + (maxY/10)])
        
        ax2.set_xticklabels(evts)
        ax2.set_xticks(x)
        ax2.set_ylabel('Genes')
        ax2.set_title('AS Genes in ' + t)
        ax2.legend()
        ax2.set_ylim([0, maxY + (maxY/10)])
        
        z = 0
        for k in n_events:
            for i in range(5):
                ax3.scatter(x=n_events[k][i], y=n_events2[k][i], 
                            c=colors[z], label=k.replace(t, '')+evts[i], s=50, marker=marks[i], alpha=0.5)
            z += 1
            
        ax3.set_xlabel('Events')
        ax3.set_ylabel('Genes')
        ax3.set_title('AS Genes X Events in ' + t)
        ax3.legend()
        ax4.set_title('Mean AS Events in ' + ', '.join([x.replace(t, '') for x in n_events]) + ' of ' + t)
        ax4.pie([np.mean([x[e] for x in n_events.values()]) for e in range(len(evts))], labels=evts, autopct='%1.0f%%')
        ax5.set_title('Mean AS Genes in ' + ', '.join([x.replace(t, '') for x in n_events2]) + ' of ' + t)
        ax5.pie([np.mean([x[e] for x in n_events2.values()]) for e in range(len(evts))], labels=evts, autopct='%1.0f%%')
        
        col += 1
        
    plt.savefig(file, dpi=150)

dirs = dirs.split(',')
labs = labs.split(',')

data = parseData(dirs)
compareEvents(data, labs)
compareGenes(data, labs)
compareNumEvents(data, labs)
plotHeatMap(data, file='asHeatMap.pdf')

print('finished.')
