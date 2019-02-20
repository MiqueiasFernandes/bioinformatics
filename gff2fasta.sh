#!/bin/bash

## LICENCE MIT
## REV 02/19
## Extract sequences from GFF by Fasta source to standard output
## Usage: gff2fasta.sh file.gff3 data.fasta
## www.mikeias.net
## bio@mikeias.net

if [ -z "$1" ] || [ -z "$2" ]
	then
		echo "Usage: gff2fasta.sh file.gff3 data.fasta feature(opcional)"
		exit
fi

GFF=$1
FASTA=$2
FILTER=${3:-.+}

for gene in $(grep  -P "\tgene\t" $GFF | cut -d= -f2); do
	for features in $(grep -P "^[^#].+\t$FILTER\t.+$gene([^0-9]|$)" $GFF | tr \\t , ) ; do
		scaffold=$(echo $features | cut -d, -f1 | head -n1)
		gene=$(echo $features | cut -d= -f2 | cut -d\; -f1 | head -n1)
		feature=$(echo $features | cut -d, -f3)
		posicao=$(echo $features | cut -d, -f4-5 | tr , - )
		id=$(echo $features | cut -d, -f9 | cut -d\; -f1 | cut -d= -f2)
		echo \>$gene\|$feature\|$posicao\|$id
		grep -xm1 "^>$scaffold" $FASTA -A40000000 | tail -n+2 | tr -d \\n | tr \> \\n | head -n1 | cut -b$posicao | sed -e "s/.\{80\}/&\n/g"
	done
done
