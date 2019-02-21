#!/bin/bash

## LICENCE MIT
## REV 02/19
## Sort scaffolds of fasta and output to std with new name
## Usage: sort_scaffolds.sh file.fasta
## www.mikeias.net
## bio@mikeias.net

if [ -z "$1" ]
	then
		echo "Usage: sort_scaffolds.sh file.fasta"
		exit
fi

FASTA=$1
PREFIX=${2:-pg.scf.}
LINE_SIZE=${3:-80}
TMP_DIR=$(mktemp -d)

split -ul1 <(\
		cat $FASTA | \
		tr -d a-z0-9\\n | \
		tr \> \\n | \
		grep -P "^[ACTGN]+$" \
	) $TMP_DIR/scaffold_

paste -d\\n \
	<(  for l in $(seq 1 $(grep -c ">" $FASTA )); do echo \>$PREFIX$l ; done ) \
	<(  for l in $(wc -c $TMP_DIR/scaffold_* | sort -grk1 | grep "scaffold_" | tr -d [:blank:]0-9 | cut -d/ -f4 ) ; \
		do cat $TMP_DIR/$l ; \
	    done ) \
| sed -e "s/.\{$LINE_SIZE\}/&\n/g"


#rm $TMP_DIR -rf 1>/dev/null 2>/dev/null
