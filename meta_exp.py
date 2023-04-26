#!/usr/bin/python3

import argparse
import pandas as pd

def extract_meta(meta, organism, controle, tratamento):
  df = pd.read_excel(meta, sheet_name='dataset')
  genoma = df.loc[df['Organismo'] == organism].loc[df['Tipo'] == 'Genoma']['Acesso'].iloc[0]
  transcriptoma = df.loc[df['Organismo'] == organism].loc[df['Tipo'] == 'Transcriptoma']['Acesso'].iloc[0]
  proteoma = df.loc[df['Organismo'] == organism].loc[df['Tipo'] == 'Proteoma']['Acesso'].iloc[0]
  gff = df.loc[df['Organismo'] == organism].loc[df['Tipo'] == 'Anotação estrutural']['Acesso'].iloc[0]
  sra_ctrl = list(df.loc[df['Organismo'] == organism]
                  .loc[df['Tipo'] == 'Amostra de RNASeq']['Acesso'].loc[df['Descrição'].apply(lambda e: controle in e and not tratamento in e)])
  sra_trat = list(df.loc[df['Organismo'] == organism]
                  .loc[df['Tipo'] == 'Amostra de RNASeq']['Acesso'].loc[df['Descrição'].apply(lambda e: tratamento in e and not controle in e)])
  open('genoma', 'w').write(genoma)
  open('transcriptoma', 'w').write(transcriptoma)
  open('proteoma', 'w').write(proteoma)
  open('genes', 'w').write(gff)
  open('controle', 'w').writelines([f'{x}\n' for x in sra_ctrl])
  open('tratamento', 'w').writelines([f'{x}\n' for x in sra_trat])
  print('finished.')

parser = argparse.ArgumentParser(description="""
Script to extract:
    1. cleaned proteome
    2. gene sequences
    3. table GENE > TRANSCRIPT > PROTEIN
from experimental design excel sheet.
""")

parser.add_argument('file')
parser.add_argument('organism')
parser.add_argument('control')
parser.add_argument('case')
args = parser.parse_args()

extract_meta(args.file, args.organism, args.control, args.case)              
