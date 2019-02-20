#!/bin/bash

## LICENCE MIT
## REV 02/19
## Extract sequences from GFF by Fasta source to standard output
## Usage: gff2fasta.sh file.gff3 data.fasta
## www.mikeias.net
## bio@mikeias.net

if [ -z "$1" ] || [ -z "$2" ]
	then
		echo "Usage: gff2fasta.sh file.gff3 data.fasta"
		exit
fi

GFF=$1
FASTA=$2

for gene in $(grep  -P "\tgene\t" $GFF | cut -d= -f2); do
	echo \>$gene
	scaffold=$(echo $gene | cut -d. -f-3)
	echo $(for cds in $(grep -P "^[^#].+\tCDS\t.+$gene([^0-9]|$)" $GFF | tr \\t , ) ; do \
		posicao=$(echo $cds | cut -d, -f4-5 | tr , - | sort -t- -gk1); \
		grep -xm1 "^>$scaffold" $FASTA -A40000 | tail -n+2 | tr -d \\n | tr \> \\n | head -n1 | cut -b$posicao ; 
		done)| tr -d \  | sed -e "s/.\{80\}/&\n/g"
done


