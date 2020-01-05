#!/bin/bash

#author: Miquéias Fernandes 01/2020 bio@mikeias.net
#run filogeny and download Nwik tree from genome.jp tool

SERVER="https://www.genome.jp/tools-bin/ete"
TIMEOUT=5

FILE=$1

SEQ=${2:-"protein"}     #  nucleotide
workflow1=${3:-"mafft_default"}  #  mafft_einsi   mafft_linsi mafft_ginsi  clustalo_default  muscle_default
workflow2=${4:-"-none"} # -trimal001 -trimal01 -trimal02 -trimal05  -trimal_gappyout
workflow3=${5:-"-none"} # -prottest_default  -pmodeltest_full_ultrafast  -pmodeltest_full_fast  -pmodeltest_full_slow  -pmodeltest_soft_ultrafast  -pmodeltest_soft_fast -pmodeltest_soft_slow
workflow4=${6:-"-fasttree_default"} # -bionj_default  -fasttree_default  -fasttree_full  -phyml_default  -phyml_default_bootstrap  -raxml_default  -raxml_default_bootstrap

DATA="upload_file=@$FILE"  #DATA="sequence=`cat $FILE`"

JOB=$(curl -s \
    -F "seqtype=$SEQ" \
    -F "seqformat=unaligned" \
    -F "$DATA" \
    -F "workflow1=$workflow1" \
    -F "workflow2=$workflow2" \
    -F "workflow3=$workflow3" \
    -F "workflow4=$workflow4" \
    -F "workflow=$workflow1$workflow2$workflow3$workflow4" $SERVER | grep -m1 'ete?id=' | cut -d= -f2 | cut -d\" -f1)

if ((  `echo $JOB | awk '{print length}' ` > 10 )) ; then

    while (( `curl -s $SERVER'?id='$JOB | grep -c "Your job is still running"` > 0 )) ; do
        sleep $TIMEOUT
    done

    echo `curl -s $SERVER'?id='$JOB | grep -m1 midpoint_data | cut -d\" -f2`
fi


