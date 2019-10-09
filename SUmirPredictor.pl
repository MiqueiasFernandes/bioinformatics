#Supplementary Document 1: ‘SUmirPredictor’ and ‘SUmirLocator’ Perl scripts
#SUmirPredictor.pl: A script that uses SUmirFold.pl outputs to screen structures of predicted miRNAs
#for specified folding criteria by conservation Copyright © Hikmet Budak, 2016. Digital copies of the
#scripts and a usage guide are freely available from the author (e-mail: hikmet.budak@montana.edu)
#Folding criteria:
# 1- Mismatches: max 4 for sense and max 6 for antisense strands
# 2- Dicer: mismatched sequences in the cut points of mature miRNA and miRNA* start sites
# 3- Multiloop: only one loop is acceptable between mature mirna and mirna*
# 4- Head: both mature mirna and mirna* sequences should not be involved in the head part)


use strict;
use warnings;
my $filename = "";
my $dic = "";

#################################################
my $usage = "SUmirPredictor.pl SUmirFold.hairpins.tbl.fileoutput SUmirFold.hairpins.folderoutput";
#################################################


if ($#ARGV==-1){
    print("Correct usage is: ".$usage."\n");
    exit -1;
} else {
    if ($#ARGV==1){
        $filename = $ARGV[0];$dic = $ARGV[1];
    } else {
        $filename = $ARGV[0];
        $dic = substr($filename, 0, length($filename)-4);
    }
}

if (!-d $dic) {
    print("Correct usage is: ".$usage."\n");
    exit -1;
}

my @hairpins = ();
my @start1 = ();
my @end1 = ();
my @screen = ();
my @whole_data = ();
my @blacklist = ('Unique', 'Hit', 'Hit ID');

open(my $FILE, $filename) || die "Cannot open file: ".$filename."\n";
chomp(my @lines = <$FILE>);
close $FILE;


my $x = -1;
my $lineID = 0;

for(my $i=0; $i<=$#lines; $i++){
    $lineID++;
    my $line = $lines[$i];
    my @lineelems = split(/\t/, $line);
    if ($line ne "" && substr($line, 0, 2) ne "\t" && !(grep( /^$lineelems[0]$/, @blacklist ))){
        $x++;
        my $n = scalar(@lineelems);
        my $temp = $line;
        if ($n < 20){
            $temp .= $lines[$lineID+1];
            my @tmparr = split(/\t/, $lines[$lineID+1]);
            $n += scalar(@tmparr);
            if ($n < 20){
                $temp .= $lines[$lineID+2];
            }
        }
        my @store = split(/\t/, $temp);
        push(@whole_data, $temp);
        push(@hairpins, $dic."/".$store[0].".hairpin.fsa_1.ct");
        push(@start1, $store[8]);
        push(@end1, $store[9]);push(@screen, '');
    }
}

my $er2count = 0;
my $count = -1;
my @b = ();

foreach my $i (@hairpins){
    $count++;
    open(my $FILE2, $i) || die "Cannot open file: ".$i."\n";
    chomp(my @lines2 = <$FILE2>);
    close $FILE2;
    my @a = (0);
    @b = (0);
    my @c = ('0');
    my $m = 0;
    eval {
        foreach my $line (@lines2) {
            $m++;
            if ($m > 1 && $line ne "") {
                my @store = split(/\t/, $line);push(@a, $store[0]);
                push(@b, $store[4]);
                push(@c, $store[1]);
            } elsif ($m == 1) {
                my $primirnalength = (split(/\t/, $line))[0];
            }
        }

        # 1- check for mismatches, max 4 for sense and max 6 for antisense
        my $mismatchB = 0;
        my $start2 = $b[$end1[$count] - 2];
        my $end2 = $b[$start1[$count]] + 2;
        my $mirna2 = '';
        for (my $j = $start2; $j < $end2 + 1; $j++) {
            $mirna2 .= $c[$j];
            if ($b[$j] == 0) {
                $mismatchB++;
                if ($mismatchB > 6) {
                    $screen[$count] = "MismatchB";
                }
            }
        }

        my @newdata = split(/\t/, $whole_data[$count]);
        $newdata[10] = $start2;
        $newdata[11] = $end2;
        $newdata[12] = $mirna2;
        $whole_data[$count] = join("\t", @newdata);
        if (!$screen[$count]) {
            my $mismatchA = 0;
            for (my $j = $start1[$count]; $j < $end1[$count] + 1; $j++) {
                if ($b[$j] == 0) {
                    $mismatchA++;
                    if ($mismatchA > 4) {
                        $screen[$count] = "Mismatch";
                    }
                }
            }
        }

        # 2- check for Dicer functions
        if (!$screen[$count]) {
            if ($b[$start1[$count]] == 0) {
                $screen[$count] = "Dicer";
            }
            if ($b[$end1[$count] - 2] == 0) {
                $screen[$count] = "Dicer";
            }
        }

        # determine start-end sites
        my $start = $start1[$count];
        my $end = $end2;
        if ($start2 < $start) {
            $start = $start2;
            $end = $end1[$count];
        }

        # 3- check for multiloop
        if (!$screen[$count]) {
            my $temp = $b[$start];
            for (my $k = $start + 1; $k < $end + 1; $k++) {
                my $check = $b[$k];
                if ($check != 0) {
                    if ($check < $temp) {
                        $temp = $check;
                    } elsif ($check > $temp) {
                        $screen[$count] = "Multiloop";
                    }
                }
            }
        }

        # 4- check for head
        if (!$screen[$count]) {
            my $temp = $b[$start];
            my $check = 0;
            for (my $k = $start + 1; $k < $end + 1; $k++) {
                if ($b[$k] != 0) {
                    if ($check == $b[$k]) {
                        my @u2 = ();
                        for (my $u = $k - $x; $u < $k; $u++) {
                            push(@u2, $u);
                        }
                        for (my $test = $start1[$count]; $test < $end1[$count] + 1; $test++) {
                            if (grep(/^$test$/, @u2)) {
                                $screen[$count] = "Head";
                            }
                        }
                        if (!$screen[$count]) {
                            for (my $test = $start2; $test < $end2 + 1; $test++) {
                                if (grep(/^$test$/, @u2)) {
                                    $screen[$count] = "Head";
                                }
                            }
                        }
                    }
                    $x = 0;
                } elsif ($b[$k] == 0) {
                    $x += 1;
                    $check = $a[$k - 1] if ($x == 1);
                }
            }
        }

        if (!$screen[$count]) {
            $screen[$count] = "OK";
        }
    }; 

    if($@) {
        $er2count++;
    }
}


open(my $OUTPUT, ">".$filename.".edited.tbl") || die "Cannot open file for writing: ".$filename.".edited.tbl\n";


print $OUTPUT "Screen\tUnique Hit ID\tNew miRNA ID\tNew miRNA Sequence\tNew miRNALength\tConserved miRNA ID\tConserved miRNA Sequence\tConserved miRNA Mismatch\tSequenceID\tMature Start\tMature End\tmiRNA* Start\tmiRNA* End\tmiRNA* Sequence\tHairpin Location\tPre-miRNA Length\tPre-miRNA MFE\tPre-miRNA GC%\tPre-miRNA MFEI\tPre-miRNA Start\tPre-miRNA Sequence\tmirna* Sequence";


my @mirnas = ();
my @seqs = ();
my @mismatch = ();
my @homolog = ();
my @index = ();
my @location = ();
my @mirna2 = ();

for (my $i=0; $i < scalar(@screen); $i++){
    my $data = $whole_data[$i];
    my @splitdata = split(/\t/, $data);
    print $OUTPUT "\n".$screen[$i]."\t".$data;
    my $start2 = $b[$end1[$count]-2];
    my $end2 = $b[$start1[$count]] + 2;
    if ($screen[$i] eq "OK" && (index("N", $splitdata[2]) == -1)){
        eval {
            push(@seqs, $splitdata[2]."\t".$splitdata[19]);
            my $mirID = (split(/\./, (split(/-/, $splitdata[1]))[1]))[0];
            my $mir = $mirID;
            $mir =~ s/[^\d]//g;
            push(@mismatch, $splitdata[6]);
            push(@homolog, $splitdata[4]);
            push(@index, $i+1);
            my $loc = substr($splitdata[13], 0, 1);
            push(@mirna2, $splitdata[12]);
            push(@location, $loc);
            my $homolog_loc = (split(/-/, $splitdata[4]))[-1];
            if (grep(/^$homolog_loc$/, ("3p", "5p")) && $loc ne substr($homolog_loc, 0, 1)){
                $mir = $mir;
            } else {
                $mir .= "-".$loc."p";
            }
            push(@mirnas, "miR".$mir);
        } or do {
            print("3: error in: ".$data."\n");
        };
    }
}


my %full;

for (my $i = 0; $i < scalar(@mirnas); $i++){
    if (grep(/^$seqs[$i]$/, keys(%full))){
        my $mirnas_now = $mirnas[$i];

        if ($full{$seqs[$i]}[1] > $mismatch[$i]){
            $full{$seqs[$i]} = [$mirnas_now, $mismatch[$i], $index[$i], $homolog[$i], $mirna2[$i]];
        } elsif($full{$seqs[$i]}[1] == $mismatch[$i]) {
            if (! grep(/^$mirnas[$i]$/, split(/,/, $full{$seqs[$i]}[0]))) {
                $full{$seqs[$i]}[0] .= ','.$mirnas[$i];
                $full{$seqs[$i]}[2] .= ','.$index[$i];
                if (! grep(/^$homolog[$i]$/, split(/,/, $full{$seqs[$i]}[3]))){
                    $full{$seqs[$i]}[3] .= ','.$homolog[$i];
                    $full{$seqs[$i]}[4] .= ','.$mirna2[$i];
                }
                print("similar mismatch found. look at mirnaID: ".$full{$seqs[$i]}[0]."\n");
            } elsif(grep(/^$mirnas[$i]$/, split(/,/, $full{$seqs[$i]}[0]))) {
                $full{$seqs[$i]}[2] .= ','.$index[$i];
                if (! grep(/^$homolog[$i]$/, split(/,/, $full{$seqs[$i]}[3]))){
                    $full{$seqs[$i]}[3] .= ','.$homolog[$i];
                    $full{$seqs[$i]}[4] .= ','.$mirna2[$i]
                }
            }
        }
    } else {
        $full{$seqs[$i]} = [$mirnas[$i], $mismatch[$i], $index[$i], $homolog[$i], $mirna2[$i]];
    }
}

close($OUTPUT);
open(my $OUTPUT2, ">".$filename.".out.tbl");
foreach my $x (keys(%full)){
    my @tmparr = @{$full{$x}};
    my $joinedstr = join("\t", splice(@tmparr, 2, scalar($full{$x})-2));
    print $OUTPUT2 $full{$x}[0]."\t".$x."\t".$joinedstr."\n";
}
close($OUTPUT2);

print("Number of structures with missing information: ".$er2count."\n");
print("\nFinished successfully.\n");


