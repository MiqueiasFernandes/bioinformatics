#!/bin/bash

### Anotação Funcional manual (NR e Uniprot) à obtenção de GOs
### Author: Miquéias Fernandes
### Last Rev: jan/2018
### More: https://github.com/MiqueiasFernandes/bioinformatics
### licence MIT

### HEADS DE uniprot.hash:
### gi	Gene ontology	IDs	Gene ontology (biological process)	Gene ontology (cellular component)	Gene ontology (molecular function)

### obter uniprot.hash: 
### cat uniprot.tsv | tr \\t \# > uniprot.hash

### HEADS DE go.names:
### GO_TERM>IDENTIFIER	GO_TERM>NAME	GO_TERM>NAMESPACE	GO_TERM>DESCRIPTION	

### mapeie os GO IDs em http://www.geneontology.org/faq/how-do-i-get-term-names-my-list-go-ids
### cut -d\# -f2 uniprot.hash | tr \; \\n | perl -pe 's/^\ // and s/\ $//' | grep -P ".+" | sort -u > gos.to.map

### obter fasta da tabela inicial
### paste -d'>' <(seq 1 $(cat seq_uniq | wc -l) | tr -d [0-9]) seq_uniq| tr \\t \\n > seq_uniq.fa

### anotar sequencias
### clusterblast seq_uniq.fa "blastx -db /home/cluster/nrdb/nr -evalue 1e-5 -outfmt '10 qseqid sgi sseqid evalue bitscore stitle' -gilist /home/cluster/shared/eudicotyledons.gilist -max_target_seqs 1" blastx_temp

### obter seq2gi: 
### cut -d, -f-2 seq_uniq.fa.out | sort -u > seq2gi


rm euc/ noeuc/ all/ -rf 2>/dev/null

echo '####################### EUCALYPTUS #####################'

mkdir euc/ && cd euc/

ln ../uniprot.hash
ln ../seq2gi seq2gi
ln ../ids_euc ids

join -t\# -12 -21 <(join ids <(sort -k1 seq2gi) | tr \  \# | sort -t\# -k2) uniprot.hash | tr \# \\t > uniprot.tab
join ids <(sort -k1 seq2gi) | cut -d\  -f1 | sort -u > genicas
echo -e "# sequencias:"\\t$(cat ids | wc -l)\\n"# genicas:"\\t$(cat genicas | wc -l)\\n"# seqs com GOs:"\\t$(cut -f2 uniprot.tab | sort -u | wc -l)

cut -f3 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.all
cat uniprot.go.all | perl -pe 's/^\ // and s/\ $//' > uniprot.go.all.ids 
sort -u uniprot.go.all.ids > uniprot.go.all.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.all.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.all.ids) ; done) > uniprot.go.all.stats

cut -f4 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.bp
cut -d[ -f2 uniprot.go.bp | tr -d ] > uniprot.go.bp.ids 
sort -u uniprot.go.bp.ids > uniprot.go.bp.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.bp.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.bp.ids) ; done) > uniprot.go.bp.stats

cut -f5 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.cc
cut -d[ -f2 uniprot.go.cc | tr -d ] > uniprot.go.cc.ids 
sort -u uniprot.go.cc.ids > uniprot.go.cc.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.cc.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.cc.ids) ; done) > uniprot.go.cc.stats

cut -f6 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.mf
cut -d[ -f2 uniprot.go.mf | tr -d ] > uniprot.go.mf.ids 
sort -u uniprot.go.mf.ids > uniprot.go.mf.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.mf.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.mf.ids) ; done) > uniprot.go.mf.stats

join -t\; <(sort -k1 uniprot.go.all.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.all.tsv
join -t\; <(sort -k1 uniprot.go.bp.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.bp.tsv
join -t\; <(sort -k1 uniprot.go.cc.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.cc.tsv
join -t\; <(sort -k1 uniprot.go.mf.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.mf.tsv

head -n20 uniprot.go.*.stats

cd ..

echo '##################### PSIDIUM #####################'

mkdir noeuc/ && cd noeuc/

ln ../uniprot.hash
ln ../seq2gi
ln ../ids_noeuc ids

join -t\# -12 -21 <(join ids <(sort -k1 seq2gi) | tr \  \# | sort -t\# -k2) uniprot.hash | tr \# \\t > uniprot.tab
join ids <(sort -k1 seq2gi) | cut -d\  -f1 | sort -u > genicas
echo -e "# sequencias:"\\t$(cat ids | wc -l)\\n"# genicas:"\\t$(cat genicas | wc -l)\\n"# seqs com GOs:"\\t$(cut -f2 uniprot.tab | sort -u | wc -l)

cut -f3 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.all
cat uniprot.go.all | perl -pe 's/^\ // and s/\ $//' > uniprot.go.all.ids 
sort -u uniprot.go.all.ids > uniprot.go.all.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.all.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.all.ids) ; done) > uniprot.go.all.stats

cut -f4 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.bp
cut -d[ -f2 uniprot.go.bp | tr -d ] > uniprot.go.bp.ids 
sort -u uniprot.go.bp.ids > uniprot.go.bp.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.bp.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.bp.ids) ; done) > uniprot.go.bp.stats

cut -f5 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.cc
cut -d[ -f2 uniprot.go.cc | tr -d ] > uniprot.go.cc.ids 
sort -u uniprot.go.cc.ids > uniprot.go.cc.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.cc.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.cc.ids) ; done) > uniprot.go.cc.stats

cut -f6 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.mf
cut -d[ -f2 uniprot.go.mf | tr -d ] > uniprot.go.mf.ids 
sort -u uniprot.go.mf.ids > uniprot.go.mf.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.mf.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.mf.ids) ; done) > uniprot.go.mf.stats

join -t\; <(sort -k1 uniprot.go.all.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.all.tsv
join -t\; <(sort -k1 uniprot.go.bp.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.bp.tsv
join -t\; <(sort -k1 uniprot.go.cc.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.cc.tsv
join -t\; <(sort -k1 uniprot.go.mf.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.mf.tsv


head -n20 uniprot.go.*.stats

cd ..

echo '##################### ALL #####################'

mkdir all/ && cd all/

ln ../uniprot.hash
ln ../seq2gi
ln ../ids_all ids

join -t\# -12 -21 <(join ids <(sort -k1 seq2gi) | tr \  \# | sort -t\# -k2) uniprot.hash | tr \# \\t > uniprot.tab
join ids <(sort -k1 seq2gi) | cut -d\  -f1 | sort -u > genicas

echo -e "# sequencias:"\\t$(cat ids | wc -l)\\n"# genicas:"\\t$(cat genicas | wc -l)\\n"# seqs com GOs:"\\t$(cut -f2 uniprot.tab | sort -u | wc -l)

cut -f3 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.all
cat uniprot.go.all | perl -pe 's/^\ // and s/\ $//' > uniprot.go.all.ids 
sort -u uniprot.go.all.ids > uniprot.go.all.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.all.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.all.ids) ; done) > uniprot.go.all.stats

cut -f4 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.bp
cut -d[ -f2 uniprot.go.bp | tr -d ] > uniprot.go.bp.ids 
sort -u uniprot.go.bp.ids > uniprot.go.bp.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.bp.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.bp.ids) ; done) > uniprot.go.bp.stats

cut -f5 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.cc
cut -d[ -f2 uniprot.go.cc | tr -d ] > uniprot.go.cc.ids 
sort -u uniprot.go.cc.ids > uniprot.go.cc.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.cc.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.cc.ids) ; done) > uniprot.go.cc.stats

cut -f6 uniprot.tab | tr \; \\n | grep -P ".+" > uniprot.go.mf
cut -d[ -f2 uniprot.go.mf | tr -d ] > uniprot.go.mf.ids 
sort -u uniprot.go.mf.ids > uniprot.go.mf.ids.uniq
sort -grk2 <(for g in $(cat uniprot.go.mf.ids.uniq) ; do echo -e $g\\t$(fgrep -c $g uniprot.go.mf.ids) ; done) > uniprot.go.mf.stats


join -t\; <(sort -k1 uniprot.go.all.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.all.tsv
join -t\; <(sort -k1 uniprot.go.bp.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.bp.tsv
join -t\; <(sort -k1 uniprot.go.cc.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.cc.tsv
join -t\; <(sort -k1 uniprot.go.mf.stats | tr \\t \;) <( cut -f-2 ../go.names | tr \\t \;) |  sort -k2 -grt\; | tr \; \\t > stats.mf.tsv

head -n20 uniprot.go.*.stats

cd ..


echo 'done'
echo 'by mikeias.net'





