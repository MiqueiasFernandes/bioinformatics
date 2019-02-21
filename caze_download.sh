#!/bin/bash

FAM=${1:-AA1_eukaryota}
PAGS=${2:-33}


echo "baixando os IDS..."

for page in $(seq 0 $PAGS)
	do 
	page=$page"00"
	wget "http://www.cazy.org/$FAM.html?debut_TAXO=$page#pagination_TAXO" -qO $page &
done

wait

cat * | grep ncbi.nlm.nih.gov/entrez | cut -d= -f5- | cut -d\  -f1 | grep -vP "^$" | sort -u > links.txt

echo $(cat links.txt | wc -l) sequencias serão baixadas ...

split -dl190 links.txt part

echo "baixando do NCBI..."

for arquivo in $(ls part* | tr \s \\n)
	do
		ids=$(cat $arquivo | tr \\n , | cut -d, -f-200)
		list=${ids: :-1}
		wget "http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id=$list" -qO fasta.$arquivo
done



cat fasta.* | grep -vP "^$" > $FAM.fa

echo "os arquivos foram salvos em $FAM.fa !"

