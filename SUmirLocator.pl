#SUmirLocator.pl –Perl script that uses SUmirPredictor.pl outputs to count the number of occurrences of
#putative miRNAs in a given genome (or any other fasta file) by conservation. Copyright © Hikmet
#Budak, 2016. Digital copies of the scripts and a usage guide are freely available from the author (e-mail:
#hikmet.budak@montana.edu)
#!/usr/bin/env perl
# mirna representation (IDs are in forms of miR156 etc.)

#################################################
use strict;
use warnings;
#################################################

my $usage = "SUmirLocator.pl genome.fasta SUmirPredictor.output.out.tbl";
my $genome;
my $filename;

if ($#ARGV<1){
    print("Correct usage is: ".$usage."\n");
    exit -1;
} else {
    $genome = $ARGV[0];
    $filename = $ARGV[1];
}


# --------------------
sub alllocations {
    my ($alist, $x) = @_;
    my @ret = ();
    while ($alist =~ /($x)/g) {
        push @ret, (pos($alist)-length $1);
    }
    return @ret;
}# ----------------------


open (my $FILE, $genome) or die "Cannot open genome file: ".$genome."\n";
my @geneID = ();
my @geneSeq = ();
my $done = 0;
my $seq = '';

while (my $line = <$FILE>) {
    chomp($line);
    if (substr($line, 0, 1) eq ">") {
        if ($seq ne '') {
            push(@geneSeq, $seq);
            $seq = '';
        }
        push(@geneID, $line);
        $done = 1;
    } else {
        if ($done == 1) {
            $seq .= $line;
        }
    }
}


push(@geneSeq, $seq);
close($FILE);


# --------------------
# my $filename = "plantmiRNAs.txt.fsa.results.tbl.hairpins.tbl.out.tbl";
open (my $FILE2, $filename);
chomp(my @alllines = <$FILE2>);
open (my $OUT1, ">".$filename.".edited");
open (my $OUT3, ">".$filename.".expression.tbl");
my %senseSeq = ();
my %antiSeq = ();
my @mirnaID = ();
my %mirnaIndex = ();
my %howmany = ();

for(my $i = 0; $i <= $#alllines; $i++){
    my $line = $alllines[$i];
    if (substr($line, 0, 3) eq "miR"){
        my @lineelems = split(/\s/, $line);
        my $mirna = $lineelems[0];
        my $mirna2 = $lineelems[-1];
        if (grep(/^$mirna$/, @mirnaID)) {
            my $a = $mirnaIndex{$mirna}+1;
            $mirnaIndex{$mirna} = $a;
        } else {
            $mirnaIndex{$mirna} = '1';
            push(@mirnaID, $mirna);
        }
        my $index = $lineelems[3];
        my $seq = $lineelems[2];
        my $newseq = $seq;
        $newseq =~ tr/UAGC/ATCG/;
        my $newseq2 = $seq;
        $newseq2 =~ tr/UAGC/TAGC/;
        while($seq =~ m/([^UAGC])/g){
            print("something is wrong with ".$seq."\n")
        }
        $newseq = reverse($newseq);
        if (grep(/^$newseq2$/, keys %senseSeq)) {
            $senseSeq{$newseq2} .= "\t".$mirna.'-'.$mirnaIndex{$mirna}.','.$index;
        } else {
            $senseSeq{$newseq2} = $mirna."-".$mirnaIndex{$mirna}.",".$index;
        }
        if (grep(/^$newseq$/, keys %antiSeq)) {
            $antiSeq{$newseq} .= "\t".$mirna.'-'.$mirnaIndex{$mirna}.','.$index;
        } else {
            $antiSeq{$newseq} = $mirna."-".$mirnaIndex{$mirna}.",".$index;
        }
        my $summary = join("\t", @lineelems[1..$#lineelems]);
        foreach my $mir2 (split(/,/, $mirna2)){
            print $OUT3 $mirna.'-'.$mirnaIndex{$mirna}."\t".$lineelems[1]."\t".$mir2."\n";
        }
        print $OUT1 $mirna."-".$mirnaIndex{$mirna}."\t".$summary."\n";

        # If you want to count mirna Isomers instead of general mirnaIDs
        # uncomment the 1st line and comment the 2nd line below
        #my $newmirna = $mirna.'-'.$mirnaIndex{$mirna};
        my $newmirna = $mirna;
        if (!(grep(/^$newmirna$/, keys %howmany))){
            $howmany{$newmirna} = [0, 0];
        }
    }
}

close($FILE2);
close($OUT1);
close($OUT3);


# --------------------------
open(my $OUT2, ">premirna-locations.csv") or die "Cannot create file: premirna-locations.csv\n";
print $OUT2 "mirnaID\tindex\treadID\tpremirna location\tpremirna length\tstrand";
eval{
    foreach my $seq2 (keys %senseSeq){
        my $geneCount = -1;
        for(my $j = 0; $j <= $#geneSeq; $j++){
            my $gene = $geneSeq[$j];
            $geneCount++;
            my @loc1 = alllocations($gene, $seq2);
            if (scalar(@loc1) > 0) {
                foreach my $loc (@loc1) {
                    my @allmirnas = split(/\t/, $senseSeq{$seq2});
                    foreach my $mirnas (@allmirnas) {
                        my @mirnaselems = split(/,/, $mirnas);
                        my @firstelems = split(/-/, $mirnaselems[0]);
                        my @geneIDelems = split(/\s/, $geneID[$geneCount]);
                        my $mirna = join("-", @firstelems[0..($#firstelems-1)]);
                        my $index = $mirnaselems[1];
                        print $OUT2 "\n".$mirnaselems[0]."\t".$index."\t".substr($geneIDelems[0], 1, length($geneIDelems[0]))."\t".$loc."\t".length($seq2)."\tSENSE";
                        if (!grep(/^$mirna$/, (keys %howmany))){
                            die;
                        }
                        $howmany{$mirna}[0]++;
                    }
                }
            }
        }
    }

    foreach my $seq3 (keys %antiSeq){
        my $geneCount = -1;
        foreach my $gene (@geneSeq){
            $geneCount++;
            my @loc1 = alllocations($gene, $seq3);
            if (scalar(@loc1) > 0) {
                foreach my $loc (@loc1) {
                    my @allmirnas = split(/\t/, $antiSeq{$seq3});
                    foreach my $mirnas (@allmirnas) {
                        my @mirnaselems = split(/,/, $mirnas);
                        my @firstelems = split(/-/, $mirnaselems[0]);
                        my @geneIDelems = split(/\s/, $geneID[$geneCount]);
                        my $mirna = join("-", @firstelems[0..($#firstelems-1)]);
                        my $index = $mirnaselems[1];
                        print $OUT2 "\n".$mirnaselems[0]."\t".$index."\t".substr($geneIDelems[0], 1, length($geneIDelems[0]))."\t".$loc."\t".length($seq3)."\tANTISENSE";
                        if (!grep(/^$mirna$/, (keys %howmany))){
                            die $mirnas;
                        }
                        $howmany{$mirna}[1]++;
                    }
                }
            }
        }
    }
}; 

if ($@) {
    print "Check if there is only one miRNA per line.\nYou should choose one from the miRNAs separated by commas.\nError in mirnaID ".$@."\n";
}

close($OUT2);
# ---------------------------------------


open(my $OUTPUT3, ">premirna-counts.csv") or die "Cannot create file:premirna-counts.csv\n";
print $OUTPUT3 "mirnaID\ton SENSE\ton ANTISENSE\ttotal";
foreach my $many (keys %howmany){
    print $OUTPUT3 "\n".$many."\t".$howmany{$many}[0]."\t".$howmany{$many}[1]."\t".($howmany{$many}[0]+$howmany{$many}[1]);
}

close($OUTPUT3);

