#!/usr/bin/env python
# coding: utf-8

# Copyright(c) 2020 - Miquéias Fernandes <bio@mikeias.net>

import os
import sys
import csv
import matplotlib.pyplot as plt
import numpy as np

#################################################
usage = "python3 plotASFilters.py dir1,dir2,dir3..."
#################################################

try: _, dirs=sys.argv
except: sys.exit("Correct usage is:\n"+usage)
    
FDR = 0.05

print('FDR: %f' % FDR)


def importFromFolder(folder, fdr=-1, maser=False, col='GeneID', raw=False, byIds={}):
    folder = folder if folder.endswith("/") else (folder+'/')
    files = [f for f in os.listdir(folder) if ((f.startswith('sign_events_')) if maser else ('.MATS.JC' in f and f.endswith('.txt')))]
    if len([y for y in files if '.JC.' in y]) > 0 and len([y for y in files if '.JCEC.' in y]) > 0:
        raise Exception('Has more than one MATS tipe in folder! choose just only one (JC or JCEC)')
    def getData(f):
        fp = list(csv.reader(open(folder + f).readlines(), delimiter='\t'))
        if raw:
            return [fp[0]] + [r for r in fp[1:] if r[0] in byIds[folder][f.split('.')[0]]] 
        idC, idF = fp[0].index(col), fp[0].index('FDR')
        return [x[idC] for x in fp[1:] if (fdr < 0) or (float(x[idF]) <= fdr)] 
    return {f.split('.')[0].replace('sign_events_', ''): getData(f) for f in files}


def tipoFolder(folder):
    return 'JCEC' if len([f for f in os.listdir(folder) if f.endswith('.MATS.JCEC.txt')]) > 0 else 'JC'


def parseData(folders, fdr=-1, maser=False, gene=True, raw=False, byIds={}, clear=[]):
    folders = [x if x.endswith('/') else (x+'/') for x in folders]
    data = {f: {t: [] for t in ['A3SS', 'A5SS', 'SE', 'MXE', 'RI']} for f in folders}
    for folder in [f for f in folders if not f in clear]:
        data[folder].update(importFromFolder(folder, fdr=fdr, maser=maser, col='GeneID' if gene else 'ID', raw=raw,byIds=byIds))
    return data


def plotFilt(byMats, byFdr, byCoverage, byMaser, file='maser_filters.pdf'):
    dirs = [ x for x in list(set(byMats).union(byFdr).union(byCoverage).union(byMaser)) if not x.endswith('coverage_filt/')]
    def autolabel(rects, ax):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    _, axs = plt.subplots(ncols=2, nrows=len(dirs), figsize=(16, 5 * len(dirs)))
    if len(dirs) < 2:
        axs = [axs]
        
    for i in range(len(dirs)):
        line = dirs[i]
        orig, fdr, cov, sign = byMats[line], byFdr[line], byCoverage[line + 'coverage_filt/'], byMaser[line] 
        axE = axs[i][0]
        axG = axs[i][1]
        
        width = 0.2
        lgT = width * 4
        data = {}
        evts = ['A3SS', 'A5SS', 'SE', 'MXE', 'RI']
        data['MATS_E'] = [len(orig[evt]) for evt in evts]
        data['MATS_G'] = [len(set(orig[evt])) for evt in evts]
        data['COV_E'] = [len(cov[evt]) for evt in evts]
        data['COV_G'] = [len(set(cov[evt])) for evt in evts]
        data['FDR_E'] = [len(fdr[evt]) for evt in evts]
        data['FDR_G'] = [len(set(fdr[evt])) for evt in evts]
        data['SIG_E'] = [len(sign[evt]) for evt in evts]
        data['SIG_G'] = [len(set(sign[evt])) for evt in evts]
        
        maxY = max([max(v) for v in data.values()])

        def plotBar(ax, data, order, labs, yl, t):
            c = 0
            x = np.arange(5)
            for k in order:
                r = ax.bar(x - (lgT / 2) + (width * c) + (width / 2), data[k], width, label=labs[c])
                autolabel(r, ax)
                c += 1
            ax.set_xticklabels(evts)
            ax.set_xticks(x)
            ax.set_ylabel(yl)
            ax.set_title(t[:-1])
            ax.legend()
            ax.set_ylim([0, maxY + (maxY/10)])

        plotBar(axE, data, 
                ['MATS_E', 'COV_E', 'FDR_E', 'SIG_E'], 
                ['MATS', 'COV', 'FDR', 'SIG'],
               'Events', 'Filtered Events in ' + line)
        plotBar(axG, data, 
                ['MATS_G', 'COV_G', 'FDR_G', 'SIG_G'], 
                ['MATS', 'COV', 'FDR', 'SIG'],
               'Genes', 'Filtered Genes in ' + line)
 
    plt.savefig(file, dpi=150)

    
def exportFiltered(raw): 
    for d, fs in raw.items():
        nd = d[:-1]  + '_maser/' 
        if not os.path.isdir(nd):
            os.mkdir(nd)
        for f, lns in fs.items():
            nf = f + '.MATS.' + tipoFolder(d) + '.txt'
            open(nd + nf, 'w').writelines(['\t'.join(l)+'\n' for l in lns])
        
        
dirs = dirs.split(',')
d_orig = parseData(dirs)
d_fdr = parseData(dirs, fdr=FDR)
cov_dirs = [d + ('' if d.endswith('/') else '/') + 'coverage_filt/' for d in dirs]
d_cov = parseData(cov_dirs, clear=[d for d in cov_dirs if not os.path.isdir(d)])
d_sign = parseData(dirs, maser=True)
plotFilt(d_orig, d_fdr, d_cov, d_sign)

d_export = parseData(dirs, maser=True, gene=False)
d_export_raw = parseData(dirs, raw=True, byIds=d_export)
exportFiltered(d_export_raw)

print('finished.')
