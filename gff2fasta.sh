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
FILTER=${3:-gene}

bedtools getfasta -fi $FASTA -bed \
	<(paste <(grep -P "\t$FILTER\t" $GFF | cut -f1,4-5) \
		<(grep -P "\t$FILTER\t" $GFF | cut -d= -f2 | cut -d\; -f1) 
	 ) \
  -name | cut -d: -f1 | sed -e "s/.\{80\}/&\n/g"
