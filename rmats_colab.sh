#!/bin/bash

GENOME=$1
GTF=$2
CTRL=$3
CASE=$4
ERRO=

[ ! -f $GENOME ] || [ -z $GENOME ] && echo "ERROR: genome.fna is obrigatory!!!" && ERRO=1
[ ! -f $GTF ] || [ -z $GTF ] && echo "ERROR: genome.gtf is obrigatory!!!" && ERRO=1
[ ! -f $CTRL ] || [ -z $CTRL ] && echo "ERROR: samples control.txt is obrigatory!!!" && ERRO=1
[ ! -f $CASE ] || [ -z $CASE ] && echo "ERROR: samples case.txt is obrigatory!!!" && ERRO=1

[ $ERRO ] && echo "usage on colab: bash rmats_colab.sh genome.fna genome.gtf control.txt case.txt"

echo "*************************************"
echo "************* INSTALING *************"
echo "*************************************"

apt install wget samtools unzip curl python3 libgsl-dev -y 

wget -O hisat.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
  && unzip -qq -o hisat.zip && mv hisat2-2.2.1/ hisat2

## instalar rmats
mkdir rmats && cd rmats && \
wget -O rmats "https://github.com/Xinglab/rmats-turbo/releases/download/v4.1.2/rmats_turbo_v4_1_2.tar.gz" && \
tar -xvf rmats && \
cd rmats_turbo* && make && \
cd .. && rm rmats && ln -s $(pwd)/$(ls rmats_turbo*/rmats.py) . 

[ $ERRO ] && exit -1

echo "*************************************"
echo "************* INDEXING *************"
echo "*************************************"

hisat2/hisat2-build -p 4 $GENOME idxgenoma 1> logs.idxgenoma.out 2> logs.idxgenoma.err

echo "*************************************"
echo "************* MAPPING *************"
echo "*************************************"


for smp in `cat $CTRL`
  do echo "starting run $smp ..." && \
    hisat2/hisat2 -x idxgenoma --sra-acc $smp -p 4 --no-unal \
      -S $smp.sam 1> logs.$smp.sam.out.txt 2> logs.$smp.sam.err.txt && \
    samtools sort -@ 4 -m 2G $smp.sam -o ctrl.$smp.sorted.bam && rm -rf $smp.sam
  done
 
for smp in `cat $CASE`
  do echo "starting run $smp ..." && \
    hisat2/hisat2 -x idxgenoma --sra-acc $smp -p 4 --no-unal \
      -S $smp.sam 1> logs.$smp.sam.out.txt 2> logs.$smp.sam.err.txt && \
    samtools sort -@ 4 -m 2G $smp.sam -o case.$smp.sorted.bam && rm -rf $smp.sam
  done

echo "*************************************"
echo "************* ANALISYiNG ************"
echo "*************************************"

ls -1 ctrl*.bam | tr \\n , | sed 's/,$//' > control
ls -1 case*.bam | tr \\n , | sed 's/,$//' > case

python3 rmats/rmats.py \
     --b1 control --b2 case --gtf $GTF -t single \
        --od rmats_out \
        --tmp tmp_out --readLength 150
        
