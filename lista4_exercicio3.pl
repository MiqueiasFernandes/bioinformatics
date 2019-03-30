#!/usr/bin/env perl

## mikeias.net 30/03/19

## 3)  Criar um script em perl que extraia as sequências em formato fasta de proteínas do arquivo.

## usage: perl lista4_exercicio3.pl [augustus.gff]

use strict;
use warnings;

my $arquivo = $ARGV[0];
open(ARQ, "<$arquivo");
my @linhas = <ARQ>;
close ARQ or die "Impossivel fechar os arquivos corretamente!";

my $em_aa = 0;

for (my $i =0; $i <= $#linhas; $i++) {
    my $linha = $linhas[$i];
    chomp $linha;
    if (index($linha, "### gene ") == 0) {
        $linha =~ s/### gene />/;
        print $linha."\n";
        $em_aa = 0;
    } elsif (index($linha, "# protein sequence = [") == 0) {
        $linha =~ s/# protein sequence = \[//;
        if (index($linha, "]") < 0) {
             $em_aa = 1;
        } else {
            $linha =~ s/\]//;
        }
        print $linha;
    } elsif ($em_aa == 1) {
        $linha =~ s/# //;
        if (index($linha, "]") > 0) {
            $linha =~ s/\]/\n/;
            $em_aa = 0;
        }
        print $linha;
    }
}

## END
