#!/bin/bash

base=base
filtros=filtros
var=1
START=$(date +%s.%N)

for filename in $filtros/*; do
	s1=${filename%.*}
	seqaremover=$filename
	out_dir=filtro_$((var++))_${s1#*/}
	echo processando $out_dir ...

	rm $out_dir -r 2>/tmp/a
	mkdir $out_dir

	#Input:
	illumina1=$base/illumina1.fastq
	illumina2=$base/illumina2.fastq
	pacbio=$base/pacbio.fa
	china=$base/china.fa

	bwa index -p $out_dir/refremover -a bwtsw $seqaremover

	bwa aln -t 24 $out_dir/refremover $illumina1 > $out_dir/illumina1.sai
	bwa samse -f $out_dir/illumina1.sam $out_dir/refremover $out_dir/illumina1.sai $illumina1
	samtools view --threads 24 -f 4 -bS $out_dir/illumina1.sam > $out_dir/illumina1.bam
	../rmrdna/bbmap/reformat.sh in=$out_dir/illumina1.bam out=$out_dir/illumina1.fastq ow=t

	bwa aln -t 24 $out_dir/refremover $illumina2 > $out_dir/illumina2.sai
	bwa samse -f $out_dir/illumina2.sam $out_dir/refremover $out_dir/illumina2.sai $illumina2
	samtools view --threads 24 -f 4 -bS $out_dir/illumina2.sam > $out_dir/illumina2.bam
	../rmrdna/bbmap/reformat.sh in=$out_dir/illumina2.bam out=$out_dir/illumina2.fastq ow=t

	bwa aln -t 24 $out_dir/refremover $pacbio > $out_dir/pacbio.sai
	bwa samse -f $out_dir/pacbio.sam $out_dir/refremover $out_dir/pacbio.sai $pacbio
	samtools view --threads 24 -f 4 -bS $out_dir/pacbio.sam > $out_dir/pacbio.bam
	../rmrdna/bbmap/reformat.sh in=$out_dir/pacbio.bam out=$out_dir/pacbio.fa ow=t

	bwa aln -t 24 $out_dir/refremover $china > $out_dir/china.sai
	bwa samse -f $out_dir/china.sam $out_dir/refremover $out_dir/china.sai $china
	samtools view --threads 24 -f 4 -bS $out_dir/china.sam > $out_dir/china.bam
	../rmrdna/bbmap/reformat.sh in=$out_dir/china.bam out=$out_dir/china.fa ow=t

	../rmrdna/bbmap/repair.sh in=$out_dir/illumina1.fastq in2=$out_dir/illumina2.fastq out=$out_dir/2illumina1.fastq out2=$out_dir/2illumina2.fastq
	rm $out_dir/illumina1.fastq $out_dir/illumina2.fastq
	mv $out_dir/2illumina1.fastq $out_dir/illumina1.fastq
	mv $out_dir/2illumina2.fastq $out_dir/illumina2.fastq

	rm $out_dir/*sai $out_dir/*sam $out_dir/refremover*

	base=$out_dir ##encadeado
	echo $out_dir terminado ...
done

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)


echo "CHINA:"
cat $out_dir/china.fa | grep ">" | wc -l
echo "PACBIO:"
cat $out_dir/pacbio.fa | grep ">" | wc -l
echo "ILLUMINA 1:"
cat $out_dir/illumina1.fastq | grep "@MG" | wc -l
echo "ILLUMINA 2:"
cat $out_dir/illumina2.fastq | grep "@MG" | wc -l



echo terminado com sucesso em $DIFF ...
