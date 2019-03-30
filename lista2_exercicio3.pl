#!/usr/bin/env perl

## mikeias.net 17/03/19

## 3) Criar um script contendo duas variáveis de string (variável 1: ATTTGC e variável 2: AGGTCC). 
##    O programa deve imprimir uma mensagem indicando se as variáveis são iguais ou não. 

use strict;
use warnings;

my $variavel1 = "ATTTGC";
my $variavel2 = "AGGTCC";


if ($variavel1 eq $variavel2) {
    print "As duas variáveis são iguais!\n";
} else {
    print "As duas variáveis são diferentes!\n";
}

## END
