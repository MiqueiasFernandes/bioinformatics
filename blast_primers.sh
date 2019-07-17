#!/bin/bash

## LICENCE MIT
## REV 07/19
## Find primers pair in fasta sequences in coverage 50% up to 100% of primers seq 
## Usage: blast_primers.sh query.tsv db.fasta COVERAGE[90]
## www.mikeias.net
## bio@mikeias.net

if [ -z "$1" ] || [ -z "$2" ]
    then
        echo "Usage: blast_primers.sh query.tsv db_blastn_name"
        echo "query format: PRIMER_NAME    SEQ_FOWARD    SEQ_REVERSE"
    exit
fi

DB=$1
QUERY=$2

rm *blast_out* primer_part_* -r 2>/dev/null

if (( `ls -1 $DB.n* 2>/dev/null | wc -l` < 3 ))
    then
    makeblastdb -in $DB -dbtype nucl
else
    echo db $DB existis ...
fi

echo database ready ...

for p in $(cut -f1 $QUERY | sort -u)
    do
    grep -P "^$p\t.+" $QUERY | cut -f1,2 | sed 's/^/>/' | sed 's/\t/_F\n/' > primer_part_$p
    grep -P "^$p\t.+" $QUERY | cut -f1,3 | sed 's/^/>/' | sed 's/\t/_R\n/' >> primer_part_$p
done

echo files ready ...


for COV in 50 60 70 80 90 100
    do
    mkdir blast_out_$COV
    echo
    echo "***** RESULTS FOR $COV% COV *****"
    for p in $(ls -1 primer_part_*)
        do
        blastn \
                -task "blastn-short" \
                -db $DB  \
                -outfmt '10 qseqid sseqid evalue bitscore length pident nident gaps gapopen qcovs qcovhsp sstart send' \
                -num_threads 24 \
                -qcov_hsp_perc $COV \
                -query $p > blast_out_$COV/$p\_blast_out
        for g in $(cut -d, -f2 blast_out_$COV/$p\_blast_out | sort -u)
            do
            if (( $( grep ",$g," blast_out_$COV/$p\_blast_out | wc -l) > 1))
                then
                if (( $(cut -f1 -d, blast_out_$COV/$p\_blast_out | grep -vP ".+_F$" | wc -l) >= 1 )) && \
                   (( $(cut -f1 -d, blast_out_$COV/$p\_blast_out | grep -vP ".+_R$" | wc -l) >= 1 )) ; 
                   then 
                       grep "$g" blast_out_$COV/$p\_blast_out
                   fi
                fi
        done > blast_out_$COV/$p\_blast_out_FILTERED
        echo $(echo $p | cut -d_ -f3-) =\> $(cut -d, -f2 blast_out_$COV/$p\_blast_out_FILTERED | sort -u | wc -l )
    done

    echo '#qseqid,sseqid,evalue,bitscore,length,pident,nident,gaps,gapopen,qcovs,qcovhsp,sstart,send' > blast_out_$COV/blast_out_FILTERED_ALL.csv
    cat blast_out_$COV/*blast_out_FILTERED >> blast_out_$COV/blast_out_FILTERED_ALL.csv
    echo =\> $( grep -v \# blast_out_$COV/blast_out_FILTERED_ALL.csv | cut -f1 -d, | sed 's/_[FR]//' | sort -u | grep -vP "^$" | wc -l) PRIMERS
    echo =\> $( grep -v \# blast_out_$COV/blast_out_FILTERED_ALL.csv | cut -f2 -d, | sort -u | grep -vP "^$" | wc -l) SEQUENCES :
    grep -v \# blast_out_$COV/blast_out_FILTERED_ALL.csv | cut -f2 -d, | sort -u | grep -vP "^$"
    echo
done
