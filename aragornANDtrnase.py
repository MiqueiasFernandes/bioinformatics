#!/usr/bin/python

## LICENCE MIT
## REV 07/19
## Join results as unic for data from aragorn output and trnase of tRNAs to tsv file
## Usage: aragornANDtrnase.py aragorn.tsv out_file_trnase.tsv outFile.tsv [limit optional default 10]
## www.mikeias.net
## bio@mikeias.net


import sys

def process(in_file_aragorn, in_file_trnase, out_file, limit=10):
    arag = [l.rstrip().split('\t') for l in open(in_file_aragorn).readlines() if not l.startswith('#')]
    trnase = [l.strip().split('\t') for l in open(in_file_trnase).readlines() if not l.startswith('#')]
    if trnase[2][0].startswith('---'):
        trnase = trnase[3:]
    gfftrnas = []
    ids = []

    ## parse data from aragorn 
    for t in arag:
        dt = [
        t[0].strip(),
        'aragorn',
        t[5].strip().upper(),
        t[6].strip().upper(),
        t[2],
        t[3],
        t[4],
        t[7] if len(t) > 7 else ''
        ]
        idT = (dt[0], dt[3], int(dt[4]), int(dt[5]), dt[6])
        if len([x for x in ids if x[0] == idT[0] and  x[1] == idT[1] and x[4] == idT[4] and (abs(x[2] - idT[2]) + abs(x[3] - idT[3]) < limit) ]) < 1:
            gfftrnas.append(dt)
            ids.append(idT)
        else:
            print("redundant in aragorn => %s %s (%d - %d)[%s]" % idT)


    ## parse data from trnase 
    for t in trnase:
        a = int(t[2])
        b = int(t[3])
        dt = [
        t[0].strip(),
        'trnase',
        t[4].strip().upper(),
        t[5].strip().upper(),
        str(min(a, b)),
        str(max(a, b)),
        '+' if a < b else '-',
        t[9] if (len(t) > 9 and len(t[9]) > 2) else ''
        ]
        idT = (dt[0], dt[3], int(dt[4]), int(dt[5]), dt[6])
        if len([x for x in ids if x[0] == idT[0] and  x[1] == idT[1] and x[4] == idT[4] and (abs(x[2] - idT[2]) + abs(x[3] - idT[3]) < limit) ]) < 1:
            gfftrnas.append(dt)
            ids.append(idT)
        else:
            print("redundant in trnase => %s %s (%d - %d)[%s]" % idT)

    with open(out_file, 'w') as f:
        f.write('\n'.join(['\t'.join(x) for x in sorted(gfftrnas, key=lambda e:( int(e[0].split('.')[2]) if e[0].count('.') > 1 else e[0], int(e[4]), e[3]))]) + '\n')

    print('%d sequences parsed with limit %d\n\n' % (len(gfftrnas), limit))
    print('terminado com sucesso!')
    print('by mikeias.net')


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("use: ./aragornANDtrnase.py aragorn.tsv out_file_trnase.tsv outFile.tsv [limit optional default 10]")
    else:
        process(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 10 )






