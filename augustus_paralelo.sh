#!/bin/bash

#Augustus PARALELO
#Copyright (C) 2018  Miquéias Fernandes

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## require BBMAP : http://sourceforge.net/projects/bbmap

## bio@mikeias.net
## http://mikeias.net

## Rodar Augustus em PARALELO
## :: Os genes vão sair com ID padrão: <NOME DO SCAFOLD><GENE ID>

GENOME=$1
SPECIE=$2
ARGS=$3  ## append args as: --AUGUSTUS_CONFIG_PATH=/my/custom/path/config in augusuts comand line
WAYS=`cat <(echo $(nproc)) <(grep -c \> $GENOME) | sort -n | head -1`

k=2
((k = $WAYS - 1))

mkdir -p parts && cd parts
echo "particionando arquivo de genoma em $WAYS parts..."
partition.sh in=$GENOME out=genome.%.fa ways=$WAYS -Xmx30g
cd ..

echo "O genoma foi particionado! ..."

mkdir -p gffs
for i in $(eval echo "{0..$k}")
        do augustus --species=$SPECIE --uniqueGeneId=true --gff3=on parts/genome.$i.fa $ARGS > gffs/augustus.$i.gff &
done;

echo "Fazendo predição com augustus com specie $SPECIE ..."
wait

echo "Predição terminada! ..."

PWD=$(pwd)
OUT=augustus.abinitio.$SPECIE.gff
touch $OUT

for i in $(eval echo "{0..$k}")
	do
		echo "### ### ### ### INICIO ARQUIVO [$i] GFF => $PWD/gffs/augusus.$i.gff ### ### ### ###" >> $OUT
		cat gffs/augustus.$i.gff >> $OUT
		echo "### ### ### ### FIM ARQUIVO $i GFF ### ### ### ###" >> $OUT
done

echo "Arquivos salvos em $PWD ..."
echo "Terminado com sucesso!"
echo "by mikeias.net"
