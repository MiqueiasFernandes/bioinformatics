#!/usr/bin/python

## LICENCE MIT
## REV 07/19
## Extract data from aragorn output trnas to tsv file
## Usage: aragorn2tsv.py output_file_aragorn.txt out.tsv
## www.mikeias.net
## bio@mikeias.net


import sys
import re

def process(in_file, out_file):
    arag = [l.rstrip() for l in open(in_file).readlines() ]
    inart = False
    parts = []
    tmp = []
    scf = ''
    old = ''
    for l in arag:
        if re.search("^\d+ nucleotides in sequence$", l):
            scf = old
        elif re.search("^\d+\.$", l):
            inart = True
        elif re.search("^>.+", l):
            inart = False
            tmp.append(l)
            image = ''
            for l in tmp:
                if l.startswith('  '):
                    image += l + '\\n'
                elif len(image) > 4:
                    break
            tip = tmp[-1]
            ini = tip.split('[')[1].split(',')[0]
            end = tip.split(',')[1][:-1]
            strand = '-' if tip.count('c[') > 0 else '+'
            antic = tip.split('(')[1].split(')')[0] if tip.count("?") > 0 else tip.split('-')[1].split('(')[0]
            tipo = tip.split('(')[2].split(')')[0] if tip.count("?") > 0 else tip.split('(')[1].split(')')[0]
            parts.append([scf, 'aragorn', ini, end, strand, antic, tipo, image])
            tmp = []
        if inart:
            tmp.append(l)
        old = l
    with open(out_file, 'w') as f:
        f.write('#sequence\tsource\tstart\tstop\tstrand\ttype\tanticodon\tTXTart\n') 
        f.write('\n'.join(['\t'.join(x) for x in parts]) + '\n')
    print('%d sequences parsed' % len(parts))
    print('terminado com sucesso!')
    print('by mikeias.net')


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("use: ./aragorn2tsv.py output_file_aragorn.txt out.tsv")
    else:
        process(sys.argv[1], sys.argv[2])






