#!/usr/bin/env perl

## mikeias.net 17/03/19

## 5) Criar um programa que recebe um valor em moeda real através da linha de comando 
##    (lembre-se da variável $ARGV[]) e imprime a conversão do valor para dólar americano.

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use constant COTACAO_DOLAR => 3.81;

my $valor_real = $ARGV[0];


if (looks_like_number($valor_real)){
    my $valor_dolar = $valor_real / COTACAO_DOLAR;
    print "\$ $valor_dolar\n";
} else {
    die "Digite um valor numérico!";
}

## END
