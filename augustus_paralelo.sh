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

## miqueias.fernandes@aluno.ufes.br
## mikeias.net

## Rodar Augustus em PARALELO
## :: Os genes vão sair com ID padrão: <NOME DO SCAFOLD><GENE ID>

GENOME=$1
SPECIE=$2
WAYS=$3

mkdir parts && cd parts
partition.sh in=../$GENOME out=genome.%.softmasked.fa ways=$WAYS
cd ..

echo "O genoma foi particionado! ..."

mkdir gffs
for i in {0..23}; 
	do augustus --species=$SPECIE --uniqueGeneId=true --gff3=on parts/genome.$i.softmasked.fa > gffs/augustus.$i.gff & done;

echo "Fazendo predição com augustus com specie $SPECIE ..."
wait

echo "Predição terminada! ..."


for i in {0..23}; 
	do gffs/augustus.$i.gff >> augustus.abinitio.gff;
done;


echo "Terminado com sucesso!"
echo "by mikeias.net"




