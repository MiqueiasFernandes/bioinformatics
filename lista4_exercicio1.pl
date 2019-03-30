#!/usr/bin/env perl

## mikeias.net 19/03/19

## 1)  Criar um script em perl que imprima a sequência complementar reversa do arquivo sequencia.faa.

## usage: perl lista4_exercicio1.pl [sequencia.faa]

use strict;
use warnings;

sub imprime_reverso_complemento {
    my $sequencia = $_[0];
    if (length($sequencia) > 0) {
            my @seq = split("", $sequencia);
            for (my $j = $#seq; $j >= 0 ; $j--) {
                if ($seq[$j] eq "A" || $seq[$j] eq "a") {
                    print 'T';
                } elsif ($seq[$j] eq "T" || $seq[$j] eq "t") {
                    print 'A';
                } elsif ($seq[$j] eq "C" || $seq[$j] eq "c") {
                    print 'G';
                } elsif ($seq[$j] eq "G" || $seq[$j] eq "g") {
                    print 'C';
                } else {
                    print $seq[$j];
                }
            }
            print "\n";
    }
}

my $arquivo = $ARGV[0];
open ARQ, "<$arquivo";
my @linhas = <ARQ>;
close ARQ or die "Impossivel fechar os arquivos corretamente!";

my $sequencia = "";

for (my $i =0; $i <= $#linhas; $i++) {
    my $linha = $linhas[$i];
    chomp $linha;
    if (index($linha, ">") == 0) {
        ## caso for identificador
        imprime_reverso_complemento $sequencia;
        print "$linha\n";
        $sequencia = "";
    } else {
        $sequencia = $sequencia.$linha;
    }
}

imprime_reverso_complemento $sequencia;

## END
