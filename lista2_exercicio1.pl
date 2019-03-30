#!/usr/bin/env perl

## mikeias.net 17/03/19

## 1) Criar um script contendo quatro variáveis de números e decide qual delas é a maior.

use strict;
use warnings;

my $v1 = int($ARGV[0]);
my $v2 = int($ARGV[1]);
my $v3 = int($ARGV[2]);
my $v4 = int($ARGV[3]);

if (($v1 > $v2) && ($v1 > $v3) && ($v1 > $v4)) {         # verifica se o valor de v1 é maior entre os outros
    print "A PRIMEIRA variável é a maior: $v1\n"; # caso SIM v1 é informado
} elsif (($v2 > $v1) && ($v2 > $v3) && ($v2 > $v4)) {    # verifica se valor de v2 é maior entre os outros
    print "A SEGUNDA variável é a maior: $v2\n";  # caso SIM v2 é informado
} elsif (($v3 > $v1) && ($v3 > $v2) && ($v3 > $v4)) {    # verifica se valor de v3 é maior entre os outros
    print "A TERCEIRA variável é a maior: $v3\n"; # caso SIM v3 é informado
} elsif (($v4 > $v1) && ($v4 > $v2) && ($v4 > $v3)) {    # verifica se valor de v4 é maior entre os outros
    print "A QUARTA variável é a maior: $v4\n";   # caso SIM v4 é informado
} else {                                          # se v1 não é maior que todos, nem v2, nem v3 e nem v4 
    print "Impossível determinar o maior valor\n";# provavelmte há mais de uma variavel com mesmo valor
}

## END
