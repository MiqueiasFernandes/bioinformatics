#!/usr/bin/env perl

## mikeias.net 19/03/19

## 6)  Desafio: Criar um arquivo sequencia2.txt com as duas sequências abaixo 
##     (a sequência deve ser mantida quebrada em várias linhas):
##     >CAJ75785.1 teste
##     MPLSYQHFRKLLLLDDGTEAGP...
##     >P03680.1
##     MKHMPRKMYSCDFETTTKV...
##     Criar um script em perl que conte quantos aminoácidos existe na sequência. 
##     Ele deve imprimir primeiro o identificador e depois o número de aminoácidos. 

## usage: perl lista3_exercicio6.pl [aminoacidos.fa]

use strict;
use warnings;

my $arquivo = $ARGV[0];
open ARQ, "<$arquivo";
my @linhas = <ARQ>;
my $qtd = -1;
for (my $i = 0; $i <= $#linhas ; $i++ ) {
    my $linha = $linhas[$i];
    chomp $linha;
    if (index($linha, ">") == 0) {
        if ($qtd > -1) {
            print "$qtd\n";
        }
        $linha =~ s/>//;
        print "$linha\t";
        $qtd = 0;
    } else {
        $qtd += length $linha;
    }
}

if ($qtd > -1) {
    print "$qtd\n";
}

close or die "Impossivel fechar os arquivos corretamente!";

## END
