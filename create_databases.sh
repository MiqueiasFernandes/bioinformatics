#!/bin/bash

mkdir db/sprot -p && cd db/sprot
echo baixando sprot ...
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
echo descompactando sprot ...
gunzip uniprot_sprot.fasta.gz

cd .. && mkdir trembl && cd trembl
echo baixando trembl ...
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
echo descompactando trembl ...
gunzip uniprot_trembl.fasta.gz

cd .. && mkdir nr && cd nr
echo baixando e descompactando nr ...
update_blastdb --decompress nr

cd .. && cd sprot
echo criando banco sprot ... 
makeblastdb -in uniprot_sprot.fasta -dbtype prot parse_seqids
cd .. && cd trembl
echo criando banco trembl ...
makeblastdb -in uniprot_trembl.fasta -dbtype prot parse_seqids

echo terminado com sucesso
echo by mikeias.net


