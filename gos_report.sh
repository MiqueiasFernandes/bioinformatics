#!/bin/bash

## LICENCE MIT
## REV 02/19
## Estract 3 list of GOs of list with all ocurrences 
## Usage: gos_report.sh all_gos_inline.txt gos_mapped.txt
## www.mikeias.net
## bio@mikeias.net


if [ -z "$1" ] || [ -z "$2" ]
        then
                echo "Usage: gos_report.sh all_gos_inline.txt gos_mapped.txt"
                echo "create gos_mapped.txt in https://yeastmine.yeastgenome.org/yeastmine/bag.do"
                echo "! after run this prog wait genetate UNIQ_GOS_FILE to submit on site !"
                exit
fi

ALL_GOS=$1.clear
GOS_MAPPED=$2
UNIQ_GOS=$ALL_GOS.uniq
GOS_STATS_ALL=$ALL_GOS.stats.all
GOS_STATS_BP=$ALL_GOS.stats.bp
GOS_STATS_CC=$ALL_GOS.stats.cc
GOS_STATS_MF=$ALL_GOS.stats.mf

echo "Clearing $ALL_GOS ..."

cat $1 | tr -cd GO:0-9\\n > $ALL_GOS
sort -u $ALL_GOS > $UNIQ_GOS

echo "Generate all stats in $GOS_STATS_ALL ..."

sort -grk2 <( \
                for GO in $(cat $UNIQ_GOS) ; do \
                        echo -e $GO\\t$(fgrep -xc $GO $ALL_GOS) ; \
                done) > $GOS_STATS_ALL


printf "Searching for $GOS_MAPPED ..."

while [ ! -f $GOS_MAPPED ]
        do
                printf "."
                sleep 1
done


echo -e \\n"Naming ids ..."

join -t\| \
                <(sort -k1 $GOS_STATS_ALL | tr \\t \|) \
                <(sort -k1 $GOS_MAPPED | tr \\t \|) \
        | tr \| \\t | sort -grk2 > tmp
paste <(cut -f1,3 tmp) <(cut -f2 tmp) <(cut -f4- tmp) > $GOS_STATS_ALL

echo "Spliting table $GOS_STATS_ALL ..."

grep -P "\tbiological_process\t" $GOS_STATS_ALL > $GOS_STATS_BP
grep -P "\tmolecular_function\t" $GOS_STATS_ALL > $GOS_STATS_MF
grep -P "\tcellular_component\t" $GOS_STATS_ALL > $GOS_STATS_CC

echo
echo

echo "STATS ALL Gos ($(cat $GOS_STATS_ALL | wc -l) gos)"
head -n5 $GOS_STATS_ALL | cut -f2-4
echo
echo "STATS Biological Process ($(cat $GOS_STATS_BP | wc -l) gos)"
head -n5 $GOS_STATS_BP | cut -f2,3
echo
echo "STATS Molecular Function ($(cat $GOS_STATS_MF | wc -l) gos)"
head -n5 $GOS_STATS_MF | cut -f2,3
echo
echo "STATS Cellular Component ($(cat $GOS_STATS_CC | wc -l) gos)"
head -n5 $GOS_STATS_CC | cut -f2,3

echo
echo

echo 'done'
echo 'by mikeias.net'

