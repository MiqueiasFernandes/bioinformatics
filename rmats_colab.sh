#!/bin/bash

GENOME=$1
GTF=$2
CTRL=$3
CASE=$4
RLEN=$5
ERRO=

( [ -z $GENOME ] || [ ! -f $GENOME ] ) && echo "ERROR: genome.fna is obrigatory!!!" && ERRO=1
( [ -z $GTF ] || [ ! -f $GTF ] ) && echo "ERROR: genome.gtf is obrigatory!!!" && ERRO=1
( [ -z $CTRL ] || [ ! -f $CTRL ] ) && echo "ERROR: samples control.txt is obrigatory!!!" && ERRO=1
( [ -z $CASE ] || [ ! -f $CASE ] ) && echo "ERROR: samples case.txt is obrigatory!!!" && ERRO=1
[ -z $RLEN ] && echo "ERROR: read length for rMATS is obrigatory!!!" && ERRO=1

([ ! -d hisat2 ] ||  [ ! -d rmats ]) && echo "*************************************"
([ ! -d hisat2 ] ||  [ ! -d rmats ]) && echo "************* INSTALING *************"
([ ! -d hisat2 ] ||  [ ! -d rmats ]) && echo "*************************************"

([ ! -d hisat2 ] ||  [ ! -d rmats ]) && \
   apt install wget samtools unzip curl python3 libgsl-dev -y 1> logs.install.out.txt 2> logs.install.err.txt

[ ! -d hisat2 ] && wget -qO hisat.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
  && unzip -qq -o hisat.zip && mv hisat2-2.2.1/ hisat2

## instalar rmats
[ ! -d rmats ] && mkdir rmats && cd rmats && \
wget -qO rmats "https://github.com/Xinglab/rmats-turbo/releases/download/v4.1.2/rmats_turbo_v4_1_2.tar.gz" && \
tar -xf rmats && \
cd rmats_turbo* && make 1>> logs.install.out.txt 2>> logs.install.err.txt && \
cd .. && rm rmats && ln -s $(pwd)/$(ls rmats_turbo*/rmats.py) . 

[ $ERRO ] && echo "*************************************"
[ $ERRO ] && echo "************* ABORTING *************"
[ $ERRO ] && echo "*************************************"
[ $ERRO ] && echo 
[ $ERRO ] && echo "usage on colab: bash rmats_colab.sh genome.fna genome.gtf control.txt case.txt 150" && exit -1

echo "*************************************"
echo "************* INDEXING *************"
echo "*************************************"

hisat2/hisat2-build -p 4 $GENOME idxgenoma 1> logs.idxgenoma.out.txt 2> logs.idxgenoma.err.txt

echo "*************************************"
echo "************* MAPPING *************"
echo "*************************************"


for smp in `cat $CTRL`
  do echo "starting run $smp ..." && \
    hisat2/hisat2 -x idxgenoma --sra-acc $smp -p 4 --no-unal \
      -S $smp.sam 1> logs.$smp.hisat2.out.txt 2> logs.$smp.hisat2.err.txt && \
    samtools sort -@ 4 -m 2G $smp.sam -o ctrl.$smp.sorted.bam && rm -rf $smp.sam
  done
 
for smp in `cat $CASE`
  do echo "starting run $smp ..." && \
    hisat2/hisat2 -x idxgenoma --sra-acc $smp -p 4 --no-unal \
      -S $smp.sam 1> logs.$smp.hisat2.out.txt 2> logs.$smp.hisat2.err.txt && \
    samtools sort -@ 4 -m 2G $smp.sam -o case.$smp.sorted.bam && rm -rf $smp.sam
  done

echo "*************************************"
echo "************* ANALISYiNG $RLEN ************"
echo "*************************************"

ls -1 ctrl*.bam | tr \\n , | sed 's/,$//' > control
ls -1 case*.bam | tr \\n , | sed 's/,$//' > case

python3 rmats/rmats.py \
     --b1 control --b2 case --gtf $GTF -t single \
        --od rmats_out \
        --tmp tmp_out --readLength $RLEN
        
