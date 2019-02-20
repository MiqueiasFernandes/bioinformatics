#!/bin/bash

## LICENCE MIT
## REV 02/19
## Extract sequences from GFF by Fasta source to standard output
## Usage: gff2fasta.sh file.gff3 data.fasta
## www.mikeias.net
## bio@mikeias.net

if [ -z "$1" ] || [ -z "$2" ]
	then
		echo "Usage: gff2fasta.sh file.gff3 data.fasta"
		exit
fi

GFF=$1
FASTA=$2
HAS_GENE=
START_GENE=0
END_GENE=0
END_ULT_EX=0
START_EX=0
END_EX=0
let CONT_IN
GENE=""


function setGENE {
	START_GENE=$(echo $1 | cut -d, -f4)
	END_GENE=$(echo $1 | cut -d, -f5)
	GENE=$(echo $1 | cut -d, -f9 | cut -d= -f2)
	HAS_GENE=1
	END_ULT_EX=0
	START_EX=0
	END_EX=0
	CONT_IN=1
}

function printINTRON {
	echo -e $(echo $1 | cut -d, -f-2),intron,$2,$3,$(echo $1 | cut -d, -f6-8),$(echo $1 | cut -d, -f9| cut -d\; -f1| cut -d. -f-5).intron$CONT_IN\;$(echo $1 | cut -d, -f9| cut -d\; -f2-) | tr , \\t
	CONT_IN=$(( $CONT_IN + 1 ))
}


for linha in $(grep -vP "^#" genes.gff | tr \\t ,| head -n30); do 
	feature=$(echo $linha | cut -d, -f3)
	if [ $HAS_GENE ]
		then ## Verificar gene
			if [ $(echo $linha | grep -cP ".+$GENE([^0-9]|$)") -gt 0 ]
				then
					#echo esta no gene $GENE
					if [ $feature = "exon" ]
						then
							START_EX=$(echo $linha | cut -d, -f4)
							END_EX=$(echo $linha | cut -d, -f5)
							#echo exon $START_EX $END_EX $END_ULT_EX
							if (( $END_ULT_EX < 1 ))
								then ## primeiro exon
									if (( $START_GENE < $START_EX ))
										then
											printINTRON $linha $START_GENE $(($START_EX - 1))
											#echo "intron inicio $START_GENE $(($START_EX - 1))"
									fi
							else  ### pelo menos um exon foi lido
								if (( ($START_EX - $END_ULT_EX) > 2 )) ## > 1pb entre exons
									then
										printINTRON $linha $(($END_ULT_EX + 1)) $(($START_EX - 1))
										#echo "intron meio $(($END_ULT_EX + 1)) $(($START_EX - 1))"
								fi
							fi
							END_ULT_EX=$END_EX
					fi
			elif [ $feature = "gene" ]
				then
					if (( $END_GENE > $END_ULT_EX )) ## ultimo exon do gene anterior
						then
							printINTRON $linha $(($END_ULT_EX + 1)) $FIM_GENE
							#echo "intron fim $(($END_ULT_EX + 1)) $FIM_GENE"
					fi		
					#echo terminou gene $GENE
					setGENE $linha
			fi
	else         ## Nenhim gene lido
		if [ $feature = "gene" ] 
			then 
				#echo first gene lido
				setGENE $linha
		fi
	fi
	echo $linha | tr , \\t
done




