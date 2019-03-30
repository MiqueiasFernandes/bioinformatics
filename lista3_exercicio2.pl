#!/usr/bin/env perl

## mikeias.net 18/03/19

## 2) Criar o arquivo sequencia.txt com o conteúdo abaixo: 
##    dlah_C8	ACTTTATATATT
##    Criar um script em perl que imprima o arquivo acima como arquivo multifasta 
##    e salvar o resultado no arquivo sequencia.fa.

## usage: perl lista3_exercicio2.pl [arquivo.tsv]

use strict;
use warnings;

my $arquivo = $ARGV[0];

open ARQ_ORIGINAL, "<$arquivo" or die "Não foi possível abrir o arquivo $arquivo !\n";
open ARQ_NOVO, ">sequencia.fa" or die "Não foi possível abrir o arquivo sequencia.fa !\n";

my @linhas = <ARQ_ORIGINAL>;

for (my $i=0; $i <= $#linhas ; $i++ ) {
    chomp $linhas[$i];
    my @colunas = split("\t", $linhas[$i]);
    print ARQ_NOVO ">$colunas[0]\n$colunas[1]\n";
}

close || die "Impossivel fechar os arquivos corretamente!";

## END
