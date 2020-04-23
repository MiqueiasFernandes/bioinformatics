##!/usr/bin/python

## LICENCE MIT
## REV 04/20
## Parse genes.read_group_tracking CUFFDIFF file
## www.mikeias.net
## bio@mikeias.net
## on UFV jupyter, before run, try
## module load python/3.7.4

import sys


#################################################
usage = "python3 parserCuffdiff.py genes.read_group_tracking"
#################################################

try: 
    file2parse=sys.argv[1]
except: sys.exit("Correct usage is:\n"+usage)


tabela = [l.strip().split("\t") 
          for l in open(file2parse).readlines()]

header = tabela[0]
lines = tabela[1:]


indexes = {h: header.index(h) for h in header}


samples = list(set(["{}_{}".format(line[indexes['condition']], line[indexes['replicate']]) 
                    for line in lines]))


parsed_lines = {}
for line in tabela[1:]:
    gene = line[indexes['tracking_id']]
    sample = "{}_{}".format(line[indexes['condition']], line[indexes['replicate']])
    if gene in parsed_lines:
        parsed_lines[gene][sample] = line[indexes['FPKM']]
    else:
        parsed_lines[gene] = {sample: line[indexes['FPKM']]}
        

transposed = []

with open('tabela_parseada.csv', 'w') as file:
    transposed = []
    header = ['gene'] + samples
    file.write(';'.join(header) + '\n')
    transposed.append(header)
    for gene in parsed_lines:
        l = [gene] + [parsed_lines[gene][h] for h in header[1:]]
        file.write(';'.join(l) + '\n')
        transposed.append(l)
        
with open('tabela_parseada_transposed.csv', 'w') as file:
    transposed.append(['sample_replicate'] + transposed[0][1:])
    transposed[0] = [x.replace('gene', 'sample').split('_')[0] for x in transposed[0]]
    for i in range(len(transposed[0])):
        file.write(';'.join([c[i] for c in transposed]) + '\n')


print('terminado com sucesso!\nby mikeias.net')

