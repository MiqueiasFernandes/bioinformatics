#!/bin/bash

## LICENCE MIT
## REV 04/19

## predicts location of ribosomal RNA genes in full genome sequences
## by using Hidden Markov Models based on alignments from a highly cureated dataset of 
## structurally aligned sequnces.

## Usage: rnammer.sh file.fasta
## www.mikeias.net
## bio@mikeias.net


## Deprecated 04/19 (server need long time to compute)
## Use: https://bigr.bios.cf.ac.uk/sacim/k8s-ssh/tree/master/7 e http://hmmer.org/download.html enjoy


if [ -z "$1" ]
	then
		echo "Usage: rnammer.sh file.fasta"
		exit
fi


FASTA=$1
KINGDOM={-2:"bac"}
TMP_DIR=rnatmp
SEQS_DIR=$TMP_DIR/SEQS
HEADS_DIR=$TMP_DIR/HEADS

mkdir $SEQS_DIR -p 2>/dev/null
mkdir $HEADS_DIR -p 2>/dev/null

for scf in $(grep -n ">" $FASTA | tr -s :\> ,)
        do 
                ID=$(echo $scf | cut -d, -f2)
                echo fetch $ID ...
                tail -n+$(( 1 + $(echo $scf | cut -d, -f1) )) $FASTA | tr -d \\n | cut -d\> -f1 > $SEQS_DIR/$ID
                             
                curl -X POST \
                     -F 'kingdom=bac' \
                     -F "SEQPASTE=$SEQ" \
                     -F "SEQSUB=@$SEQS_DIR/$ID" \
                     -F 'configfile=/usr/opt/www/pub/CBS/services/RNAmmer-1.2/RNAmmer.cf' \
                  http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi --dump-header $HEADS_DIR/$ID
        done



#END



