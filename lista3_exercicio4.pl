#!/usr/bin/env perl

## mikeias.net 18/03/19

## 4) Criar um script em perl que leia o arquivo sequencia.txt e altere timinas pela letra X. 
##    As sequências devem ser impressas no formato fasta.

## usage: perl lista3_exercicio4.pl [arquivo.tsv]

use strict;
use warnings;

my $arquivo = $ARGV[0];
open ARQ_ENTRADA, "<$arquivo";
my @linhas = <ARQ_ENTRADA>;

for (my $i = 0; $i <= $#linhas; $i++) {
    my $linha = $linhas[$i];
    $linha =~ s/T/X/g;
    $linha =~ s/\t/\n/g;
    print ">".$linha;
}

close or die "Impossivel fechar os arquivos corretamente!";

## END
