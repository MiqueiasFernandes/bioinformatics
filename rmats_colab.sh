#!/bin/bash

EXP=`head -1 experiment.txt`
CORES=4 && MEM=2G && echo "$EXP on ENV: threads $CORES ram $MEM"
GENOME=$1
GTF=$2
CTRL=$3
CASE=$4
RLEN=$5
ERRO=

[ -z $EXP ] && echo "ERROR: experiment.txt is obrigatory!!!" && ERRO=1
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
cd .. && rm rmats && ln -s $(pwd)/$(ls rmats_turbo*/rmats.py) . && cd ..

( [ ! -f rmats/rmats.py ] || [ ! -f hisat2/hisat2 ] ) && echo "instalation failure ...."  && ERRO=1
[ $ERRO ] && echo "*************************************"
[ $ERRO ] && echo "************* ABORTING **************"
[ $ERRO ] && echo "*************************************"
[ $ERRO ] && echo 
[ $ERRO ] && echo "usage on colab: bash rmats_colab.sh genome.fna genome.gtf control.txt case.txt 150" && exit -1

echo "*************************************"
echo "************* INDEXING **************"
echo "*************************************"

[ ! -f indexed ] && \
hisat2/hisat2-build -p $CORES $GENOME idxgenoma 1> logs.idxgenoma.out.txt 2> logs.idxgenoma.err.txt && \
echo "indexed on `date +%d/%m\ %H:%M` ..." > indexed

echo "*************************************"
echo "************* MAPPING ***************"
echo "*************************************"


for smp in `cat $CTRL`
  do echo "starting run $smp on `date +%d/%m\ %H:%M` ..." && [ ! -f ctrl.$smp.sorted.bam ] && \
    hisat2/hisat2 -x idxgenoma --sra-acc $smp -p $CORES --no-unal \
      -S $smp.sam 1> logs.$smp.hisat2.out.txt 2> logs.$smp.hisat2.err.txt && \
    samtools sort -@ $CORES -m $MEM $smp.sam -o ctrl.$smp.sorted.bam \
    1> logs.$smp.samtools.out.txt 2> logs.$smp.samtools.err.txt && rm -rf $smp.sam
  done
 
for smp in `cat $CASE`
  do echo "starting run $smp on `date +%d/%m\ %H:%M` ..." && [ ! -f case.$smp.sorted.bam ] \
    hisat2/hisat2 -x idxgenoma --sra-acc $smp -p $CORES --no-unal \
      -S $smp.sam 1> logs.$smp.hisat2.out.txt 2> logs.$smp.hisat2.err.txt && \
    samtools sort -@ $CORES -m $MEM $smp.sam -o case.$smp.sorted.bam \
    1> logs.$smp.samtools.out.txt 2> logs.$smp.samtools.err.txt && rm -rf $smp.sam
  done

echo "********************************************"
echo "************* ANALISYiNG $RLEN *************"
echo "********************************************"

ls -1 ctrl.*.sorted.bam | tr \\n , | sed 's/,$//' > control
ls -1 case.*.sorted.bam | tr \\n , | sed 's/,$//' > case
echo "Running rMATS `date +%d/%m\ %H:%M` ...."
python3 rmats/rmats.py \
     --b1 control --b2 case --gtf $GTF -t single \
        --od rmats_out \
        --tmp tmp_out --readLength $RLEN --nthread $CORES \
        1> logs.rmats.out.txt 2> logs.rmats.err.txt
 
 tail -15 logs.rmats.out.txt
 zip -q results_$EXP.zip -r rmats_out   
 zip -q logs_$EXP.zip logs.*       
 echo "finished $EXP on `date +%d/%m\ %H:%M`."
