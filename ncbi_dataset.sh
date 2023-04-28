#!/bin/bash

## require wget / curl / pandas

REPO='https://raw.githubusercontent.com/MiqueiasFernandes/bioinformatics/master'
FILE="$1"
ORGANISM="$2"
CTRL="$3"
CASE="$4"
ORG=`echo $ORGANISM | tr -s \  _ `
EXP=$ORG.$CTRL.$CASE
echo $EXP > experiment.txt

curl -s $REPO/meta_exp.py | python3 - "$FILE" "$ORGANISM" "$CTRL" "$CASE"
echo wait...

wget -qO genoma.fna.gz `head -1 genoma`
wget -qO genoma.gtf.gz `head -1 genes`
wget -qO transcriptoma.fna.gz `head -1 transcriptoma`
wget -qO proteoma.faa.gz `head -1 proteoma`
gunzip -f genoma.fna.gz genoma.gtf.gz proteoma.faa.gz transcriptoma.fna.gz

curl -s $REPO/parse_ncbi.py | python3 - genoma.fna proteoma.faa genoma.gtf

mv clean.fna genes.fna
mv clean.txt gene2mrna2ptna.csv
mv clean_ptnas.faa proteoma.faa
mv controle smp_ctrl.txt
mv tratamento smp_case.txt
rm genes genoma proteoma transcriptoma
cat smp_ctrl.txt smp_case.txt > $EXP.txt
gzip -k genes.fna && cp genes.fna.gz $ORG.genes.fna.gz
gzip -k proteoma.faa && cp proteoma.faa.gz $ORG.proteoma.faa.gz

echo all done.
