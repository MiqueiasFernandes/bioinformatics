#!/usr/bin/env perl

## mikeias.net 30/03/19

## 5)  Criar um script em perl que extrai as sequências em formato fasta para as cds.

## usage: perl lista4_exercicio5.pl [augustus.gff] [seqs.fa]

use strict;
use warnings;

sub imprime_seq {
    my @linhaGFF = split("\t", $_[0]);
    my @seq = split("", $_[1]);
    for (my $j = $linhaGFF[3] -1; $j <= $#seq && $j < $linhaGFF[4]; $j++) {
        print $seq[$j];
    }
}

my $arquivoGFF = $ARGV[0];
open(GFF, "<$arquivoGFF");
my @linhasGFF = <GFF>;

my $arquivoFA = $ARGV[1];
open(FASTA, "<$arquivoFA");
my @linhasFASTA = <FASTA>;

my $chr = "";
my $sequence;

for (my $i = 0; $i <= $#linhasGFF; $i++) {
    my $linha = $linhasGFF[$i];
    chomp $linha;
    if ($linha =~ /\#\#\# gene .+/) {
        $linha =~ s/\#\#\# gene //;
        print ">".$linha."\n";
    } elsif ($linha =~ /.+\tAUGUSTUS\tCDS\t.+/) {
        my $chrq = $linha;
        $chrq =~ s/\t.*//;
        if ($chrq ne $chr) {
            ## obter a sequencia do chr
            my $em_seq = 0;
            for (my $i=0; $i <= $#linhasFASTA; $i++) {
                if ($em_seq == 0) {
                    my $seq = $linhasFASTA[$i];
                    chomp $seq;
                    if (index($seq, ">".$chr) == 0) {
                        ## esta no chr
                        $em_seq++;
                        my $buff = "";
                        ## pega as proximas linhas enquanto elas não começarem com ">"
                        while($i <= $#linhasFASTA && index($linhasFASTA[++$i], ">") < 0) {
                            chomp $linhasFASTA[$i];
                            $buff .= $linhasFASTA[$i];
                        }
                        $sequence = $buff;
                        $chr = $chrq;
                    }
                }
            }
        }
        if ($chrq eq $chr) {
            imprime_seq $linha, $sequence;
        } else {
            die "sequencia não encontrada no fasta ".$chrq."\n";
        }
    } elsif ($linha =~ /\#\#\# end gene .*/) {
        print "\n";
    }
}

close or die "Impossivel fechar os arquivos corretamente!";

# END
