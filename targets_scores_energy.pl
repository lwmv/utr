#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use autodie;
use Tie::File;
use Data::Dumper;

my @mirna = glob "./mirna/*";

mkdir "targes_scores_energy";
foreach my $mirna(@mirna){
    (my $pre = basename($mirna)) =~ s/\.\w+$//;
    my $btar = $pre . "_btar";
    my $copredict = $pre . "_copredict";
    
    my @files = map{ (my $x=basename($_))=~s/\.\w+$//; $x} glob("./$btar/*.fasta");

    my %results;
    foreach(@files){
        my @array;
        tie @array, 'Tie::File', "./$btar/$_.fasta", recsep=>"\n";
        my $name = $_;
        my $line_one = $array[0];

        my $lines = 0;
        my $utrlen = "NA";
        my $energy = 0;
        my $score = 0;
        if($line_one =~ m/^\d/) {
            $lines = scalar(@array);
            $utrlen = (split(/\s+/,$line_one))[6];
            foreach(@array) {
                my($e,$s) = (split)[7,8];
                $energy += $e;
                $score += $s;
            }
            $energy = $energy/$lines;
            $score = $score/$lines;
        }

        # compute total targets;
        my @total_targets;
        tie @total_targets,'Tie::File',"./$copredict/$_.fasta",recsep=>"\n";
        my $total_targets = scalar(@total_targets) - 10;

        $results{$name} = [$name, $utrlen, $total_targets, $lines, $energy, $score];
    }
    ## output
    open(my $fh,'>',"targes_scores_energy/$pre.txt"); 
    print{ $fh } join("\t",qw/utrname utrlen total_targets effective_targets energy score/);
    print{ $fh } "\n";
    foreach(sort(keys(%results))){
        print{$fh} join("\t",@{ $results{$_} });
        print{$fh} "\n";
    }
    close $fh;
}

