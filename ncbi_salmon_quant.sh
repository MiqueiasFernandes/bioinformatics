#!/bin/bash

#  Copyright (c) 2022 Miquéias Fernandes

#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:

#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.

#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.

N_ARGS=$#
if [ $N_ARGS -lt 5 ]
then
 echo "Usage:  $> bash ncbi_salmon_quant.sh http://...genoma.gz http://...gtf.gz http://...cds.gz dir_temp/ RUN,SAMPLE,FACTOR ... ... RUNn,SAMPLEn,FACTORn"
 exit 1
fi

MIN_READ_LEN=80
export TZ=America/Sao_Paulo
tid=t$(date +%s)
mkdir results$tid
TEMP_DIR=$4
if [ ! -d $TEMP_DIR ]
    then 
    echo "criando diretorio temporario: $TEMP_DIR" > results$tid/resumo.txt
    mkdir  $TEMP_DIR
fi
TEMP_DIR="../$TEMP_DIR"
cd results$tid
echo "[0    ] ARGS:  $@"
echo "[1    ] $( date +%D.%H:%M:%S) preparando o ambiente results$tid/ ..."
p=1

## sra-toolkit : https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
## trimmomatic : http://www.usadellab.org/cms/?page=trimmomatic   http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
## fastqc      : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## salmon      : https://salmon.readthedocs.io/en/latest/salmon.html
## samtools    : https://www.htslib.org/doc/samtools.html
## hisat2      : http://daehwankimlab.github.io/hisat2/manual/

for prog in sra-toolkit trimmomatic fastqc salmon samtools bamtools hisat2
    do
    if ! command -v $prog 1> /dev/null 2> /dev/null
    then
        echo "[1.$p  ] instalando o $prog ..."
        apt install $prog -y 1> _1.$p\_install.$prog.log 2> _1.$p\_install.$prog.err
        (( p=p+1 ))
    fi
done

## instalar nova versao do sra
## https://www.biostars.org/p/9527325/#9527333
## https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

mkdir sranew && cd sranew
wget -qO sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
 && tar -xzf sratoolkit.tar.gz >/dev/null \
 && ln -s ./sratoolkit.*/bin/fastq-dump
sratoolkit.*/bin/vdb-config -i --restore-defaults 1> vdb_config.log
cd ../

echo "usando o sra-toolkit NOVO ! Versão: $( sranew/fastq-dump --version )" >  _1.0_pacotes.log
echo "usando o salmon ! Versão: " >> _1.0_pacotes.log && salmon -v 2>> _1.0_pacotes.log
echo "usando o hisat2 ! Versão: $( hisat2 --version )" >>  _1.0_pacotes.log
echo "usando o bamtools ! Versão: $( bamtools --version )" >>  _1.0_pacotes.log
echo "usando o sra-toolkit ! Versão: $( fastq-dump --version )" >>  _1.0_pacotes.log
echo "usando o trimmomatic ! Versão: $( TrimmomaticPE -version )" >> _1.0_pacotes.log
echo "usando o fastqc ! Versão: $( fastqc --version )" >> _1.0_pacotes.log

## multiqc   : https://multiqc.info/docs/
## biopython : https://biopython.org/
## deeptools : https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

for pkg in multiqc biopython deeptools
    do
    if [[ ! `pip list | grep $pkg` ]]
        then 
        echo "[1.$p  ] instalando o $pkg ..."
        pip install $pkg 1> _1.$p\_install.$pkg.log 2> _1.$p\_install.$pkg.err
        (( p=p+1 ))
    fi
    if [[ `pip list | grep $pkg` ]]
        then 
        echo "usando o $pkg ! $( pip list | grep $pkg | tr -s \  \   )" >> _1.0_pacotes.log
    else
        echo ERRO: ao instalar pacote $pkg
    fi
done

if ! grep 'pysam.index(bamFile)' /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py 1>/dev/null 2>/dev/null
    then
    cp /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py .
    cp bamHandler.py /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py.old
    grep -B100000 'bam = pysam.Samfile(bamFile' bamHandler.py | grep -v 'bam = pysam.Samfile(bamFile' > xtemp
    echo '        pysam.index(bamFile)' >> xtemp
    tail bamHandler.py -n+`grep -n 'bam = pysam.Samfile(bamFile' bamHandler.py | cut -d: -f1` >> xtemp
    cp xtemp /usr/local/lib/python3.7/dist-packages/deeptools/bamHandler.py
    rm xtemp bamHandler.py
fi

echo "[2    ] $( date +%D.%H:%M:%S) importando genoma e a anotação ..."
echo '[2.1  ] baixando o genoma ...'
wget -O genoma.$tid.fa.gz $1 1> _2.1_genoma.download.log 2> _2.1_genoma.download.err
echo '[2.2  ] descompactando o genoma ...'
gunzip genoma.$tid.fa.gz 1> _2.2_genoma.unzip.log 2> _2.2_genoma.unzip.err
echo "Tamanho do genoma: $(grep -v \>  genoma.$tid.fa | tr -d '\n' | wc -c | rev | cut -c7- | rev)Mpb" >> resumo.txt
echo "Quantiade de sequencias no genoma: $(grep -c \>  genoma.$tid.fa)"  >> resumo.txt
echo '[2.3  ] indexando o genoma para validar as amostras ...'
hisat2-build genoma.$tid.fa idxgenoma.$tid 1> _2.3_genoma.index.log 2> _2.3_genoma.index.err
echo '[2.4  ] baixando o GTF ...'
wget -O gene.$tid.gtf.gz $2 1> _2.4_gtf.download.log 2> _2.4_gtf.download.err
echo '[2.5  ] descompactando o GTF ...'
gunzip gene.$tid.gtf.gz 1> _2.5_gtf.unzip.log 2> _2.5_gtf.unzip.err
echo "Quantidade de genes: $(grep -v '^#' gene.$tid.gtf | cut -f 3 | grep -c gene)"  >> resumo.txt

echo "[3    ] $( date +%D.%H:%M:%S) importando transcritos ..."
echo '[3.1  ] baixando os transcritos ...'
wget -O cds.$tid.fa.gz $3 1> _3.1_transcripts.download.log 2> _3.1_transcripts.download.err
echo '[3.2  ] descompactando os transcritos ...'
gunzip cds.$tid.fa.gz 1> _3.2_transcripts.unzip.log 2> _3.2_transcripts.unzip.err
echo "Quantidade de genes cod prot: $(grep -v '^#' gene.$tid.gtf | cut -f3,9 | grep '^CDS' | cut -f2 | tr \; '\n' | grep '^gene_id ' | uniq | wc -l)"  >> resumo.txt
echo "Quantiade de sequencias CDS: $(grep -c \>  cds.$tid.fa)"  >> resumo.txt
echo "Tamanho total da CDS: $(grep -v \>  cds.$tid.fa | tr -d '\n' | wc -c | rev | cut -c7- | rev)Mpb"  >> resumo.txt
echo '[3.3  ] indexando todos transcritos para validar a anotacao ...'
hisat2-build cds.$tid.fa idxcds.$tid 1> _3.3_cds.index.log 2> _3.3_cds.index.err

echo '[3.4  ] filtrando os transcritos ...'
echo "cds = 'cds.$tid.fa'" > script.py
cat >> script.py << EOF 
seqs = [(l.strip(), l[1:-1].split()) for l in open(cds).readlines() if l.startswith('>')]
print(len(seqs), 'sequencias de CDS')
pars = [[a, b[0], c] for a, b, c in [[x[0], [z for z in x if 'gene=' in z], y] for y, x in seqs] if len(b) == 1]
conts = {g: [0, []] for g in set([x[1] for x in pars])}
print(len(conts), 'genes')
for a, b, c in pars:
  conts[b][0]+=1
  conts[b][1].append(c)
print(len(set([k for k, v in conts.items() if v[0] < 2])), 'genes sem AS')
print(len(set([k for k, v in conts.items() if v[0] > 1])), 'genes com AS')
ok = []
for k, v in conts.items():
  if v[0] > 1:
    for seq in v[1]:
      ok.append(seq)
print(len(ok), 'CDS de genes com AS')
k=False
as_cds = []
genetrn = []
for l in open(cds).readlines():
  if l.startswith('>'):
    k = l.strip() in ok
    if k:
      genetrn.append(f"{l[1:].split()[0]},{l.split('gene=')[1].split(']')[0].strip()}\n")
  if k:
    as_cds.append(l)
open(cds, 'w').writelines(as_cds)
open('transcript_gene_mapping.csv', 'w').writelines(genetrn)
EOF
python3 script.py 1> _3.4_transcripts.filter.log 2> _3.4_transcripts.filter.err
rm script.py
echo "Genes com AS anotado: $(grep 'genes com AS' _3.4_transcripts.filter.log | head -1 | cut -d\  -f1)"  >> resumo.txt
echo "CDS de genes com AS anotado: $(grep -c \>  cds.$tid.fa)"  >> resumo.txt
echo "Tamanho total da CDS de genes com AS: $(grep -v \>  cds.$tid.fa | tr -d '\n' | wc -c | rev | cut -c7- | rev)Mpb"  >> resumo.txt

echo '[3.5  ] extraindo sequencia de genes ...'
echo "cds = 'cds.$tid.fa'" > script.py
echo "genoma = 'genoma.$tid.fa'" >> script.py
echo "gtf = 'gene.$tid.gtf'" >> script.py
cat >> script.py << EOF
from Bio import SeqIO, Seq, SeqRecord
gen_acecc = set([l.split('gene=')[1].split()[0].replace(']', '') for l in open(cds).readlines() if l.startswith('>')])
gns = [l.strip().split('\t') for l in open(gtf).readlines() if '\tgene\t' in l]
cords = [[x[0], int(x[3]), int(x[4]), x[6] == '+', 
          [z.split('"')[1] for z in [k.strip() for k in x[-1].split(";")] if z.startswith('gene_id "') or z.startswith('gene "')]
          ] for x in gns]
print(len(cords), 'genes no GTF')
cords = [x for x in cords if any([z for z in x[-1] if z in gen_acecc])]
print(len(cords), 'genes a exportar')
seqs = SeqIO.to_dict(SeqIO.parse(genoma, 'fasta'))
gseqsF = [SeqRecord.SeqRecord(seqs[s[0]].seq[s[1]-1:s[2]], id=s[-1][0], description=s[-1][1]) for s in cords if s[3]]
gseqsR = [SeqRecord.SeqRecord(seqs[s[0]].seq[s[1]:s[2]+1].reverse_complement(), id=s[-1][0], description=s[-1][1]) for s in cords if not s[3]]
SeqIO.write(gseqsF+gseqsR, 'gene_seqs.fa', 'fasta')
print('finalizado.')
EOF
python3 script.py 1> _3.5_genes.extract.log 2> _3.5_genes.extract.err
rm script.py
echo '[3.6  ] indexando sequencia de AS genes para gerar o BED ...'
hisat2-build gene_seqs.fa idxgenes 1> _3.6_genes.index.log 2> _3.6_genes.index.err

echo '[3.7  ] indexando os transcritos para quantificar ...'
salmon index -t cds.$tid.fa --index idx$tid 1> _3.7_transcripts.index.log 2> _3.7_transcripts.index.err

salvar () {
    rm logs -rf && mkdir logs
    cp *.log *.err resumo.txt logs
    zip -r $1.zip out_$1*/** 1>/dev/null 2>/dev/null
    zip -r logs.zip logs/** 1>/dev/null 2>/dev/null
    echo "Salvando $1.zip e logs.zip em  $TEMP_DIR ..." >> resumo.txt
    cp $1.zip logs.zip $TEMP_DIR
}

restaurar () {
    if [ -f $TEMP_DIR/$1.zip ]
    then 
        echo "tentando restaurar $TEMP_DIR/$1.zip" >> resumo.txt
        cp $TEMP_DIR/$1.zip .
        unzip $1.zip 1>/dev/null 2>/dev/null
        if [ -f out_$1/$1.quant.sf ] 
        then 
          return 0 
        else
          rm out_$1/ -rf
          echo "impossivel restaurar $1 ..." >> resumo.txt
          return 1
        fi
    else
        return 1
    fi
}

echo "[4    ] $( date +%D.%H:%M:%S) quantificando amostras ..."
GENE=$(grep \> gene_seqs.fa | head -1000 | tail | head -1 | tr -d \> | cut -d\   -f1)
echo "RUN,SAMPLE,FACTOR,FOLDER" > experimental_design.csv
i=1
for x in $@
    do 
        if [[ `echo $x | grep ,` ]]
        then
            RUN=`echo $x | cut -d, -f1`
            SAMPLE=`echo $x | cut -d, -f2`
            FACTOR=`echo $x | cut -d, -f3`
            
            echo "$RUN,$SAMPLE,$FACTOR,$SAMPLE" >> experimental_design.csv

            if restaurar $SAMPLE
                then 
                    echo "[4.$i  ] $SAMPLE restaurado de $TEMP_DIR/ ..."
                    continue
                else
                    echo "Rodando $RUN em $SAMPLE ..." >> resumo.txt
            fi

            echo "[4.$i.1] $( date +%D.%H:%M:%S) obtendo a amostra $SAMPLE pelo acesso $RUN no sra ..."
            fastq-dump --split-3 --minReadLen $MIN_READ_LEN $RUN 1> _4.$i.1_download.$RUN.$SAMPLE.log 2> _4.$i.1_download.$RUN.$SAMPLE.err
            
            if [[ $(ls -lh | grep -c _[12].fastq) < 1 ]]
            then
            echo "[4.$i.1] $( date +%D.%H:%M:%S) obtendo a amostra $SAMPLE pelo acesso $RUN no NOVO sra ..."
            sranew/fastq-dump --split-3 --minReadLen $MIN_READ_LEN $RUN 1> _4.$i.1_download2.$RUN.$SAMPLE.log 2> _4.$i.1_download2.$RUN.$SAMPLE.err
            fi
            
            grep 'Read' _4.$i.1_download.$RUN.$SAMPLE.log >> resumo.txt
            
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                echo "[4.$i.2] fazendo controle de qualidade da amostra $SAMPLE com o TrimmomaticPE ..."
                echo "tratando $SAMPLE como PE ..." >> resumo.txt
                TrimmomaticPE \
                    $RUN\_1.fastq $RUN\_2.fastq \
                    $SAMPLE.F.fq $SAMPLE.1.unp.fq \
                    $SAMPLE.R.fq $SAMPLE.2.unp.fq \
                    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                    1> _4.$i.2_qc.$SAMPLE.log 2> _4.$i.2_qc.$SAMPLE.err
                else
                echo "[4.$i.2] fazendo controle de qualidade da amostra $SAMPLE com o TrimmomaticSE ..."
                echo "tratando $SAMPLE como SE ..." >> resumo.txt
                TrimmomaticSE \
                    $RUN.fastq $SAMPLE.fq \
                    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
                    1> _4.$i.2_qc.$SAMPLE.log 2> _4.$i.2_qc.$SAMPLE.err
            fi
            grep 'Surviving' _4.$i.2_qc.$SAMPLE.err >> resumo.txt
            
            echo "[4.$i.3] reportando controle de qualidade da amostra $SAMPLE com fastqc ..."
            rm qc_$SAMPLE -rf && mkdir qc_$SAMPLE
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                    fastqc $SAMPLE.F.fq $SAMPLE.R.fq -o qc_$SAMPLE 1> _4.$i.3_stats.$SAMPLE.log 2> _4.$i.3_stats.$SAMPLE.err
                else
                    fastqc $SAMPLE.fq -o qc_$SAMPLE 1> _4.$i.3_stats.$SAMPLE.log 2> _4.$i.3_stats.$SAMPLE.err
            fi
            
            echo "[4.$i.4] mapeando $SAMPLE  no genoma com hisat2..."
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                    hisat2 -x idxgenoma.$tid -1 $SAMPLE.F.fq -2 $SAMPLE.R.fq --no-unal -S $SAMPLE.maped.sam 1> _4.$i.4_genomap.$SAMPLE.log 2> _4.$i.4_genomap.$SAMPLE.err
                else
                    hisat2 -x idxgenoma.$tid -U $SAMPLE.fq --no-unal -S $SAMPLE.maped.sam 1> _4.$i.4_genomap.$SAMPLE.log 2> _4.$i.4_genomap.$SAMPLE.err
            fi
            rm $SAMPLE.maped.sam
            echo "Mapeamento no genoma: $(grep 'overall' _4.$i.4_genomap.$SAMPLE.err)" >> resumo.txt
            
            echo "[4.$i.5] mapeando $SAMPLE  na CDS com hisat2..."
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                    hisat2 -x idxcds.$tid -1 $SAMPLE.F.fq -2 $SAMPLE.R.fq --no-unal -S $SAMPLE.maped.sam 1> _4.$i.5_cdsmap.$SAMPLE.log 2> _4.$i.5_cdsmap.$SAMPLE.err
                else
                    hisat2 -x idxcds.$tid -U $SAMPLE.fq --no-unal -S $SAMPLE.maped.sam 1> _4.$i.5_cdsmap.$SAMPLE.log 2> _4.$i.5_cdsmap.$SAMPLE.err
            fi
            rm $SAMPLE.maped.sam
            echo "Mapeamento na CDS: $(grep 'overall' _4.$i.5_cdsmap.$SAMPLE.err)" >> resumo.txt

            echo "[4.$i.6] quantificando com a amostra $SAMPLE com salmon ..."
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                    salmon quant -1 $SAMPLE.F.fq -2 $SAMPLE.R.fq \
            -o quant_$SAMPLE --libType IU --index idx$tid 1> _4.$i.6_quant.$SAMPLE.log 2> _4.$i.6_quant.$SAMPLE.err
                else
                    salmon quant -r $SAMPLE.fq \
            -o quant_$SAMPLE --libType IU --index idx$tid 1> _4.$i.6_quant.$SAMPLE.log 2> _4.$i.6_quant.$SAMPLE.err
            fi
            echo "CDS expressa em $SAMPLE: $(cut -f4 quant_$SAMPLE/quant.sf | tail -n+2 | grep -cv  '^0$')" >> resumo.txt

            echo "[4.$i.7] mapeando com a amostra $SAMPLE com hisat2 ..."
            if [[ $(ls -lh | grep -c _[12].fastq) > 1 ]]
                then
                    hisat2 -x idxgenes -1 $SAMPLE.F.fq -2 $SAMPLE.R.fq --no-unal -S $SAMPLE.maped.sam  1> _4.$i.7_map.$SAMPLE.log 2> _4.$i.7_map.$SAMPLE.err
                else
                    hisat2 -x idxgenes -U $SAMPLE.fq --no-unal -S $SAMPLE.maped.sam  1> _4.$i.7_map.$SAMPLE.log 2> _4.$i.7_map.$SAMPLE.err
            fi
            echo "Mapeamento nos AS genes: $(grep 'overall' _4.$i.7_map.$SAMPLE.err)" >> resumo.txt
            
            echo "[4.$i.8] gerando arquivo de cobertura para a amostra $SAMPLE com deeptools ..."
            samtools view -S -b $SAMPLE.maped.sam > $SAMPLE.maped.bam  2> _4.$i.8_bam.$SAMPLE.err
            bamtools sort -in $SAMPLE.maped.bam -out $SAMPLE.sorted.bam  1> _4.$i.8_bam.$SAMPLE.log 2>> _4.$i.8_bam.$SAMPLE.err
            bamCoverage -b $SAMPLE.sorted.bam -o $SAMPLE.bed --outFileFormat bedgraph --binSize 3 -p 2 -r $GENE 1> _4.$i.8_cov.$SAMPLE.log 2> _4.$i.8_cov.$SAMPLE.err

            echo "[4.$i.9] limpando dados de $SAMPLE ..."
            mkdir out_$SAMPLE
            cp quant_$SAMPLE/quant.sf out_$SAMPLE/$SAMPLE.quant.sf
            mv qc_$SAMPLE out_$SAMPLE
            mv quant_$SAMPLE out_$SAMPLE
            mv $SAMPLE.bed out_$SAMPLE
            mv $SAMPLE.sorted.bam out_$SAMPLE
            mv _4.$i.*.log  _4.$i.*.err out_$SAMPLE
            rm *.fastq *.fq *.bam* *.sam* -f
            (( i=i+1 )) 
            salvar $SAMPLE   
        fi
done 

echo "[5    ] $( date +%D.%H:%M:%S) executando o multiqc ..."
multiqc out_*/qc_* 1> _5_multiqc.log 2> _5_multiqc.err
rm logs -rf && mkdir logs
cp *.log *.err resumo.txt logs
zip -r logs.zip logs/** 1>/dev/null 2>/dev/null
cp logs.zip multiqc_*.html $TEMP_DIR

echo "[6    ] $( date +%D.%H:%M:%S) preparando output para o 3D-RNAseq ..."
mkdir to3d
for o in out_*
do 
    x=`echo $o|cut -c5-`
    mkdir to3d/$x
    cp $o/*.quant.sf to3d/$x/quant.sf
done
cd to3d && zip -r to3d.zip * 1>/dev/null 2>/dev/null && mv to3d.zip ../ && cd ..

echo "[7    ] $( date +%D.%H:%M:%S) compactando para RESULTS.zip ..."
zip -r RESULTS.zip out_*/** to3d.zip transcript_gene_mapping.csv experimental_design.csv multiqc_*.html *.log *.err 1>/dev/null 2>/dev/null
cp RESULTS.zip ../ && cp RESULTS.zip to3d.zip $TEMP_DIR

cd ..
echo $( date +%D.%H:%M:%S) terminado.
