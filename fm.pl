#!/usr/bin/perl -w
# perl fm.pl file.txt
use strict;
use warnings;
use Data::Dumper;
use File::Basename;

my $dir = dirname( $ARGV[0] );
my $out = basename( $ARGV[0] );
open my $file,"<",shift or die "no input! \n";


###fold
#system "rm $dir/tmp";
mkdir "$dir/tmp";
####

my $out_file = "coe_$out";
unlink $out_file if -e $out_file;
open OUT, ">", "$dir/tmp/$out_file" or die "use sudo \n";

my  $UTR;

( $UTR ) = $out =~ /(.+)\./;

my  %end;
my  %start;
my  $total_energy;
my  $total_score;
my  $no = 0;

while (my $line = <$file>){
  $line =~ s/\s*$//;
  my @getin = split /\t/,$line;
  my $id   = $getin[1];
  my $mirna = $getin[2];
  my $start = $getin[3];
  my $end   = $getin[4];
  my $utr   = $getin[5];
#  $UTR = $utr;
  my $tag   = 0;
  my $energy= $getin[7];
  my $score = $getin[8];
  $total_energy += $energy;
  $total_score  += $score;
  $no++;
  $end{$id} = [ $id, $mirna, $start, $end, $utr, $tag]; # where 0 means no syn
  $start{$start} = [ $id, $mirna, $start, $end, $utr, $tag];
}

foreach my $end ( sort keys %end ){
  my $tag = $end{$end}->[5];
  my $new_start = $end{$end}->[2] + 17;
  my $new_end   = $end{$end}->[2] + 35;
  foreach my $start ( sort keys %start ){
    if ( $start >= $new_start && $start <= $new_end ){
      $tag++;
      push @{$end{$end}},@{$start{$start}}[1,2,3];
      next;
    }
    else {
      next;
    }
  }
  $end{$end}->[5] = $tag;
}
############################################################
# generate the co-operate mirna
############################################################

my $average_energy = $total_energy / $no;
my $average_score  = $total_score  / $no;

foreach ( keys %end ) {
  print OUT join "\t",@{$end{$_}}[1,2,3,5..( $#{$end{$_}} )],"\n";
} 

my $total_targets = 0;
my $coe1 = 1.2;
my $coe2 = 1.4;
my $coe3 = 1.6; ##
my $no_of_co = 0;
my $no_of_co_0 = 0;
my $no_of_co_1 = 0;
my $no_of_co_2 = 0;
my $no_of_co_3 = 0;

foreach my $targets ( keys %end ) {
  if ( $end{$targets}->[5] == 0 ) { $total_targets ++; $no_of_co_0++ };
  if ( $end{$targets}->[5] == 1 ) { $total_targets += 1 * $coe1; $no_of_co++ ; $no_of_co_1++ };
  if ( $end{$targets}->[5] == 2 ) { $total_targets += 1 * $coe2; $no_of_co++ ; $no_of_co_2++ };
  if ( $end{$targets}->[5] > 2 ) { $total_targets += 1*$coe3; $no_of_co++ ; $no_of_co_3++ };
}

#############################################################################################
#utrname   total_targets(coe considered)  co-effect( pairs)  energy(average) score(average)
#############################################################################################
#if ( ! -z $file ) { print "no targets!"s }
#else {
printf <<EOT;
$UTR\t$no\t$total_targets\t$no_of_co_0\t$no_of_co\t$no_of_co_1\t$no_of_co_2\t$no_of_co_3\t$average_energy\t$average_score
EOT

#}
__END__
