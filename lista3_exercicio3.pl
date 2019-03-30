#!/usr/bin/env perl

## mikeias.net 18/03/19

## 3) Criar um script perl que leia o arquivo sequencia.fa e salve o identificador em um array 
##    e a sequência em outro array com o mesmo índice dos identificadores. 
##    O programa deve imprimir primeiro a lista de identificadores e depois a lista de sequências.

## usage: perl lista3_exercicio3.pl [arquivo.tsv]

use strict;
use warnings;

my $arquivo = $ARGV[0];

open ARQ_ENTRADA, "<$arquivo" or die "Houve um erro ao abrir o arquivo $arquivo !";

my @linhas = <ARQ_ENTRADA>;
my @identificadores;
my @seqs;

for(my $i = 0; $i <= $#linhas; $i++ ) {
    chomp $linhas[$i];
    if (index($linhas[$i], ">") == 0) { ## é identificador
        push(@identificadores, $linhas[$i]."\n");
    } else { ## é sequencia
        push(@seqs, $linhas[$i]."\n");
    }
}

print @identificadores;
print @seqs;

close || die "Impossivel fechar os arquivos corretamente!";

## END
