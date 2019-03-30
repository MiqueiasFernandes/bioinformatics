#!/usr/bin/env perl

## mikeias.net 30/03/19

## 4)  Criar um script em perl que extraia as sequências em formato fasta para o gene completo.

## usage: perl lista4_exercicio4.pl [augustus.gff] [seqs.fa] gene1

use strict;
use warnings;

sub imprime_seq {
    my @linhaGFF = split("\t", $_[0]);
    my @seq = split("", $_[1]);
    if ($#linhaGFF < 8) {
        @linhaGFF[8] = $linhaGFF[0]."-".$linhaGFF[2]."-undefined";
    }
    print ">".$linhaGFF[8]."\n";
    for (my $j = $linhaGFF[3] -1; $j <= $#seq && $j < $linhaGFF[4]; $j++) {
        print $seq[$j];
    }
    print "\n";
}

my $arquivoGFF = $ARGV[0];
open(GFF, "<$arquivoGFF");
my @linhasGFF = <GFF>;

my $arquivoFA = $ARGV[1];
open(FASTA, "<$arquivoFA");
my @linhasFASTA = <FASTA>;

my $gene = $ARGV[2];
my $chr = "";
my $em_gene = 0;
my @sequence;

for (my $i = 0; $i <= $#linhasGFF; $i++) {
    my $linha = $linhasGFF[$i];
    chomp $linha;
    if ($linha =~ /.+\tAUGUSTUS\tgene\t.+$gene[^0-9]*/) {
        $em_gene = 1;
        $chr = $linha;
        $chr =~ s/\t.*//;
        ## obter a sequencia do chr
        my $em_seq = 0;
        for (my $i=0; $i <= $#linhasFASTA; $i++) {
            if ($em_seq == 0) {
                my $seq = $linhasFASTA[$i];
                chomp $seq;
                if (index($seq, ">".$chr) == 0) {
                    ## esta no chr
                    $em_seq++;
                    $chr = "";
                    ## pega as proximas linhas enquanto elas não começarem com ">"
                    while($i <= $#linhasFASTA && index($linhasFASTA[++$i], ">") < 0) {
                        chomp $linhasFASTA[$i];
                        $chr .= $linhasFASTA[$i];
                    }
                }
            }
        }
        imprime_seq $linha, $chr;
    } elsif ($em_gene == 1 && index($linha, "#") < 0) {
        imprime_seq $linha, $chr;
    } elsif ($em_gene == 1 && index($linha, "# protein sequence = [") == 0) {
        $em_gene++;
    }
}

if ($em_gene < 1) {
    print "O gene $gene nao foi encontrado!\n";
}

close or die "Impossivel fechar os arquivos corretamente!";

# END
