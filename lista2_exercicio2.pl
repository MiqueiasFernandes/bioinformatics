#!/usr/bin/env perl

## mikeias.net 17/03/19

## 2) Criar um script em perl que possui quatro variáveis numéricas nomeadas de v1 a v4. 
##    Imprimir uma mensagem quando as variáveis pares forem maiores que as duas variáveis impares 
##    ou quando a soma das impares for maior que a soma das pares.

use strict;
use warnings;

my @PARES = ();
my @IMPARES = ();

for (my $i=0; $i<4; $i++) {
    if (($ARGV[$i] % 2)  == 0) { # testa se o "resto da divisão por 2" é igual a zero
        push(@PARES, $ARGV[$i]);
    } else {               # se o resto da divisão for diferente de zero é ímpar
        push(@IMPARES, $ARGV[$i]);
    }
}

my $quantidade_de_impares = @IMPARES;
my $quantidade_de_pares = @PARES;

if ($quantidade_de_impares != 2) { # Termina o programa caso não haja duas ímpares
    die "Não há duas variáveis ímpares";
}

## verificar se as variáveis pares forem maiores que as duas variáveis impares

my $var_par_maior_q_2_var_impar = 1;

if (($PARES[0] < $IMPARES[0]) || ($PARES[0] < $IMPARES[1])) {
    $var_par_maior_q_2_var_impar = 0; # atribui falso caso a variavel par for menor que alguma ímpar
}
if (($PARES[1] < $IMPARES[0]) || ($PARES[1] < $IMPARES[1])) {
    $var_par_maior_q_2_var_impar = 0;  # atribui falso caso a variavel par for menor que alguma ímpar
}

## Imprimir uma mensagem quando as variáveis pares forem maiores que as duas variáveis impares
if ($var_par_maior_q_2_var_impar) {
    print "As variáveis par (@PARES) são maiores que as duas variáveis impares (@IMPARES)!\n";
} else { ## verificar se soma das impares for maior que a soma das pares.
    my $soma_impares = $IMPARES[0] + $IMPARES[1];
    my $soma_pares = $PARES[0];
    if ($quantidade_de_pares > 1) {
        $soma_pares += $PARES[1];
    }
    if ($soma_impares > $soma_pares) {
        print "A soma das impares ($soma_impares) é maior que a soma das pares ($soma_pares)!\n"
    }
}

## END
