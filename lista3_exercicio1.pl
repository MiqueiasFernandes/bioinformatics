#!/usr/bin/env perl

## mikeias.net 18/03/19

## 1) Criar um script em perl que calcule o fatorial de número inteiro a escolha do usuário. 
##    Lembre-se, fatorial de 4 é igual a 4 x 3 x 2 x 1 = 24.  
##    Dica: utilize a função for para calcular o fatorial de um número.

## usage: perl lista3_exercicio1.pl [Numero]

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $numero = $ARGV[0];

if (!looks_like_number($numero) || $numero < 0){
     die "Digite um número igual ou maior que zero!\n";
} else {
   if ($numero <= 1) {
        print "1\n";
        exit;
    }
    for (my $i = $numero; $i >= 2; $i-- ) {
        $numero *= $i -1;
    }
    print "$numero\n";
}

## END
