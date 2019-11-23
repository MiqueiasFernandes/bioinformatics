#!/bin/bash

mkdir tmp

FILE=$1
cp $FILE tmp/main.tsv
FILE=tmp/main.tsv
OUT=$FILE.ptnas.anot

cut -f1,14 $FILE  | grep GO | cut -f1 | sort -u > $OUT.go
cut -f1,12 $FILE  | grep IPR  | cut -f1 | sort -u > $OUT.ipro
cut -d\; -f1 $FILE  | sort -u > $OUT.TOTAL
cut -f4 $FILE  | sort -u > dbs

while read DB
    do cut -f1,4 $FILE | grep $DB | cut -f1 | sort -u > $OUT.$DB
done < dbs

for PATHWAY in `cut -f1,15 $FILE  | grep : | cut -f2 | tr \| \\n | cut -d\: -f1 | sort -u`
    do cut -f1,15 $FILE  | grep $PATHWAY: | cut -f1 | sort -u > $OUT.$PATHWAY
done

echo -e DB\\tproteins\\tgenes
for f in `ls -1 $OUT.*`
    do echo `echo $f | sed 's/^.*ptnas.anot.//'`,`cat $f | wc -l`,`cut -d. -f1 $f | sort -u | wc -l`
done |  sort -t, -nk3 | tr , \\t > $1.stats

rm dbs tmp -r

cat $1.stats
