#!/bin/bash

GENOMA=$1
GFF=$2
PROTEOMA=$3
REPO=https://raw.githubusercontent.com/MiqueiasFernandes/bioinformatics/master

apt install wget curl python3 python3-pip 
pip install biopython

wget -O genoma.fna.gz $GENOMA 
wget -O proteoma.faa.gz $PROTEOMA 
wget -O genes.gff.gz $GFF 
gunzip -f *.gz

curl -s $REPO/ptnaClean.py | python3 -  proteoma.faa 
curl -s $REPO/gff.py | python3 - genes.gff new genoma.fna proteoma.faa_protein_id.txt 

[ -f new_gene2mrna2ptna.txt ] \
&& mkdir clean \
&& cp genoma.fna clean \
&& cp new_clean.gff clean/genes.gff \
&& cp new_genes.fna clean/genes.fna \
&& cp proteoma.faa_clean.faa clean/proteoma.fna \
&& cp new_mrnas.fna clean/transcripts.fna \
&& cp new_gene2mrna2ptna.txt clean/gene2mrna2ptna.txt \
&& rm new_* genes.gff genoma.fna proteoma.*

