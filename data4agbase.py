#!/usr/bin/python

import sys


def process(files):
    data = []
    for f in files:
        print('importando arquivo %s ...' % f)
        data.extend([l.strip().split('\t') for l in open(f).readlines()])
    cont = 0
    gs = 0
    print('processando %d registros ...' % len(data))
    with open('data4agbase.tsv', 'w') as out:
        for g in set([e[0] for e in data]):
            gs += 1
            for k in set(','.join([x[1] for x in data if x[0] == g]).split(',')):
                out.write('%s\t%s\n' % (g, k))
                cont += 1
    print('%d genes salvos' % gs)
    print('%d gos salvos' % cont)
    print('Terminado com sucesso!\nby mikeiasnet')


if __name__ == "__main__":
        if len(sys.argv) < 2:
                print("use: ./data4agbase.py tsvDATA")
        else:
                process(sys.argv[1:])