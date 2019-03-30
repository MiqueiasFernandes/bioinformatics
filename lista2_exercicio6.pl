#!/usr/bin/env perl

## mikeias.net 17/03/19

## 6) Desafio: Criar um script com uma variável que contem a sequência de nucleotídeos  “AATCTGGTGC”. 
##    O programa deve contar quantos nucleotídeos a sequência tem e imprimir este número.

use strict;
use warnings;

my $seq = "AATCTGGTGC";
## usando a função length para determinar o tamanho da string
my $qtd_nuc = length $seq;

## pode usar split para tranformar em array e $#nome_array para obter ultimo indice do array 
## my @n = split("", $seq);
## my $qtd_nuc = $#n + 1;

## ou execute no cmd
## print `printf $seq | wc -c`;

print "$qtd_nuc\n";

## END
