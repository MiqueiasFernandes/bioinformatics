#!/usr/bin/env perl

## mikeias.net 19/03/19

## 5) Criar um script perl que traduza as sequências de nucleotídeos para proteína do arquivo sequencia.txt. 
##    Dica: entrar com dois arquivos no programa, 
##    um contendo as sequências de nucleotídeos 
##    e outro contendo a relação entre códons e aminoácidos.  

## usage: perl lista3_exercicio5.pl [nucleotideos.txt] [codon_aa.tsv]

use strict;
use warnings;

my $nuc_arq = $ARGV[0];
my $cod_aa_arq = $ARGV[1];

open ARQ_NUC, "<$nuc_arq";
open ARQ_COD_AA, "<$cod_aa_arq";
open ARQ_SAIDA, ">$nuc_arq.prt";

my @linhas_nuc = <ARQ_NUC>;
my @linhas_cod_aa = <ARQ_COD_AA>;

my @cod;
my @aa;
my $imports = 0;

for ($imports = 0; $imports <= $#linhas_cod_aa ; $imports++) {
    chomp $linhas_cod_aa[$imports];
    my @cols = split("=", $linhas_cod_aa[$imports]);
    push(@cod, $cols[0]);
    push(@aa, $cols[1]);
}

print "Foram importados $imports codons\n";

for(my $i = 0 ; $i <= $#linhas_nuc; $i++) {
    my $linha = $linhas_nuc[$i];
    chomp $linha;
    if (index($linha, ">") == 0) {
        print "$linha\n";
        if ($i > 1) {
            print ARQ_SAIDA "\n";
        }
        print ARQ_SAIDA "$linha\n";
    } else {
        my $trad = $linha;
        for (my $j = 0 ; $j < $imports ; $j++) {
            $trad =~ s/$cod[$j]/$aa[$j]/g;
        }
        print "$linha => $trad\n";
        print ARQ_SAIDA "$trad";
    }
}


close or die "Impossivel fechar os arquivos corretamente!";

## END
