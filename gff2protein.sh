#!/bin/bash

## LICENCE MIT
## REV 02/19
## Extract aminoacids sequences of GFF to Fasta in standard output
## Usage: gff2protein.sh file.gff3
## www.mikeias.net
## bio@mikeias.net


if [ -z "$1" ] || [ -z "$2" ]
	then
		echo "Usage: gff2protein.sh file.augustus.gff3"
		exit
fi


GFF=$1

grep -P "(^# [A-Z]+(]|$))|(.+= \[\w+)|(.+start gene .+)" $GFF  | \
    tr \ \#= \\n| tr e \> | \
    grep -vP "^(pr|se|s)|^$" | \
    tr -d \\n | tr n[] \\n | \
    grep -P "^>|[A-Z]" | \
    sed -e "s/.\{80\}/&\n/g"
