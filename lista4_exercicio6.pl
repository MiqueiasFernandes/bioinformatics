#!/usr/bin/env perl

## mikeias.net 30/03/19

## 6)  Criar um script em perl que receba um arquivo multifasta 
##     e quebre este arquivo em múltiplos arquivos contendo o 
##     número de sequências definido pelo usuário.

## usage: perl lista4_exercicio6.pl [seqs.fa] [#fragment]

use strict;
use warnings;

my $arquivoFA = $ARGV[0];
open(FASTA, "<$arquivoFA");
my @linhasFASTA = <FASTA>;

my $nArq = 1;

my $arquivoOUT = $arquivoFA.$nArq.".fa";
open(FASTA_OUT, ">$arquivoOUT");
print "salvando em $arquivoOUT ";
my $numSeqsUsuario = $ARGV[1];
my $numSeqsArquivo = 0;
my $trocar = 0;

for (my $i = 0; $i <= $#linhasFASTA; $i++) {
    my $linha = $linhasFASTA[$i];
    if (index($linha, ">") == 0) {
        if ($trocar == 1) {
            close FASTA_OUT or die "Impossivel fechar o arquivo $arquivoOUT corretamente!";
            print " => ".$numSeqsArquivo." seqs\n";
            $nArq++;
            $arquivoOUT = $arquivoFA.$nArq.".fa";
            open(FASTA_OUT, ">$arquivoOUT");
            $numSeqsArquivo = $trocar = 0;
            print "salvando em $arquivoOUT ";
        }
        $numSeqsArquivo++;
    }
    if ($numSeqsArquivo >= $numSeqsUsuario) {
        $trocar = 1;
    }
    print FASTA_OUT $linha;
}
print " => $numSeqsArquivo seqs\n";

close or die "Impossivel fechar os arquivos corretamente!";

# END

