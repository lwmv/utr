#!/usr/bin/perl -w
use strict;
# use Tk;
use Data::Dumper;
use File::Basename;

die "illeagle usage" if $#ARGV > 2 || $#ARGV < 0;
print, die  <<'BB' if $ARGV[0] eq '-h';

  Usage:     perl copredict.pl [mirnafiles] [utrfiles]
  Author:    blinfu@gmail.com
  Attention: You need to install findtar, miranda & microTar.

BB

my $mfile = basename ( $ARGV[0] );
my $ufile = basename ( $ARGV[1] );
my $mfile_dir = dirname ( $ARGV[0] );
my $ufile_dir = dirname ( $ARGV[1] );
my ( $utr )= $ufile =~ /(.*)\./;
my $utr_length = 0;

format STDOUT =
@|||||||||||||||||||||||||||||||||||||||||||||||||||||||||
'############################################################'
@<@|||||||||||||||||||||||||||@|||||||||||||||||||||||||@>
'#','MicroRNA File',              'UTR File','#'
@<@|||||||||||||||||||||||||||@|||||||||||||||||||||||||@>
'#',$mfile,          $ufile,'#'
@<@|||||||||||||||||||||||||||@|||||||||||||||||||||||||@>
'#','Path',                        'Path','#'
@<@|||||||||||||||||||||||||||@|||||||||||||||||||||||||@>
'#',$mfile_dir,      $ufile_dir,'#'
@<@|||||||||||||||||||||||||||@|||||||||||||||||||||||||@>
'#','ALL', $utr,'#'
@|||||||||||||||||||||||||||||||||||||||||||||||||||||||||
'############################################################'
.

write ;

####################################################
#run others scripts                                #
####################################################

mkdir "tmp";
system "findtar $mfile_dir/$mfile $ufile_dir/$ufile -energy -15 -score 5 > tmp/$utr.findtar";
system "miranda $mfile_dir/$mfile $ufile_dir/$ufile > tmp/$utr.miranda";
system "microtar -q $mfile_dir/$mfile -t $ufile_dir/$ufile > tmp/$utr.microtar";


mkdir "../UTR";
open UTR,"<$ufile_dir/$ufile" or die "failed to open ! \n";


while (<UTR>){
  if ( /^>/ ){
  print;
  }
  else {
    chomp;
    s/\r\n|\n|\r//;
    $utr_length += length($_);
    print;
       }
}
close UTR;
print "\n";


####################################################
#3 files created! x.findtar, x.miranda, x.microtar #
####################################################

open FINDTAR,  "<tmp/$utr.findtar"  or die "no such file! $! \n";
open MIRANDA,  "<tmp/$utr.miranda"  or die "no such file! $! \n";
open MICROTAR, "<tmp/$utr.microtar" or die "no such file! $! \n";

my %score_findtar,  my $no_findtar  = 1;
my %score_miranda,  my $no_miranda  = 1;
my %score_microtar, my $no_microtar = 1;

#####################################################
# generate findtar status                           #
#####################################################

while (<FINDTAR>) {
  next if $.<14;
  my @score = split "\t";
  $score_findtar{$no_findtar} = {
	       'soft'  => 'findtar',
	       'mirna' => $score[0],
	       'score' => $score[1],
	       'energy'=> $score[2],
	       'start' => $score[3],
	       'end'   => $score[4],
	       'utr'   => $utr,
	       'utrlen'=> $utr_length,
	      };
  $no_findtar++;
}
close FINDTAR;

#####################################################
# generate miranda status                           #
#####################################################

while (<MIRANDA>){
  next if $.<32;
  next unless $_ =~ /^>>/;
  my ( $mirna) = $_ =~ />>(.*?)\s+/;
  my @score = split /\s+/;
  my $start = $score[-1];
  $score_miranda{$no_miranda} = {
				 'soft'  => 'miranda',
				 'mirna' => $mirna,
				 'utr'   => $utr,
				 'start' => $start,
				};
  $no_miranda++;
}
close MIRANDA;

#####################################################
# generate microTar status                          #
#####################################################

my @score = <MICROTAR>;
close MICROTAR;

foreach my $line ( 9 .. ($#score-2) ){
  chomp $score[$line];
  if (exists $score[$line+6] &&  $score[$line] =~ /.*?:\s(\w+)/ && $score[$line+6]=~/No targets found/){
    next;
  }
   elsif ( exists $score[$line+8] && $score[$line] =~ /.*?:\s(\w+)/ && $score[$line+8]=~/\s\smRNA/){
     my ( $mirna ) = $score[$line] =~ /.*?:\s(.*)/;
     $score_microtar{$no_microtar} = {
				      'soft'=> 'microtar',
				     'mirna'=> $mirna,
				      'utr' => $utr,
				     };
     $no_microtar++;
     next;
   }
}

##############################################################################
#现在已经把每个软件预测出来的miran名字，utr名字存放在3个hash里面。hash结构是编号，然后  #
#是一个存放信息的hash里面有miran名字，接下来，只需要看看，我们findtar找出来的mirna是否 #
#存在于其他2个软件之间的某一个里面就可以了。因此只需要遍历我们的findtar的含有mirna名字的 #
#hash，然后同时看这个miran是否存在于其他2个软件里面。如果有则保留，否则删除，得到我们的最 #
#后结果      %score_findtar     %score_miranda      %score_microtar           #
##############################################################################

#############################################################################
#
#   $score_findtar{$no_findtar} = {
# 	       'soft'  => 'findtar',
# 	       'mirna' => $score[0],
# 	       'score' => $score[1],
# 	       'energy'=> $score[2],
# 	       'start' => $score[3],
# 	       'end'   => $score[4],
# 	       'utr'   => $utr,
# 	      };


#   $score_miranda{$no_miranda} = {
# 				 'soft'  => 'miranda',
# 				 'mirna' => $mirna,
# 				 'utr'   => $utr,
# 				 'start' => $start,
# 				};
#   $no_miranda++;
# }

#      $score_microtar{$no_microtar} = {
# 				      'soft'=> 'microtar',
# 				     'mirna'=> $mirna,
# 				      'utr' => $utr,
# 				     };
#
###########################################################################

my %miranda;
my %microtar;

foreach my $no ( sort keys %score_miranda) {
  my $mirna = $score_miranda{$no}->{mirna};
  $miranda{$mirna} = 'yes';
}

foreach my $no ( sort keys %score_microtar ){
  my $mirna = $score_microtar{$no}->{mirna};
  $microtar{$mirna} = 'yes';
}
##################################################################
#  %miranda     %microtar                                        #
##################################################################

foreach my $no ( sort keys %score_findtar ){
  my $mirna = $score_findtar{$no}->{mirna};
  if ( exists $miranda{$mirna} ){
    $score_findtar{$no}->{'miranda'} = 'yes';
  }
  else {
    $score_findtar{$no}->{'miranda'} = 'no'; 
 };
  
  if( exists $microtar{$mirna} ){
    $score_findtar{$no}->{'microtar'} = 'yes';
  }
  else {
    $score_findtar{$no}->{'microtar'} = 'no';
  };
}

#print Dumper(%score_findtar);

# print join "\t",'Id','Microrna','Start','End','Utr','Utrlen','Energy','Score','Miranda','Microtar';
# print "\n";
# foreach my $id (sort {$a <=> $b} keys %score_findtar){
# print join "\t", $id, $score_findtar{$id}->{mirna},$score_findtar{$id}->{start},$score_findtar{$id}->{end},$score_findtar{$id}->{utr},$score_findtar{$id}->{utrlen},$score_findtar{$id}->{energy},$score_findtar{$id}->{score},$score_findtar{$id}->{miranda},$score_findtar{$id}->{microtar};
# print "\n";
# }

###############################################
#下面输出的是至少有其他一个软件预测出来有靶点的情况   #
###############################################

print '#';
print join "\t",'Id','Microrna','Start','End','Utr','Utrlen','Energy','Score','Miranda','Microtar';
print "\n";
foreach my $id (sort {$a <=> $b} keys %score_findtar){
next unless $score_findtar{$id}->{miranda} eq 'yes' || $score_findtar{$id}->{microtar} eq 'yes';
print join "\t", $id, $score_findtar{$id}->{mirna},$score_findtar{$id}->{start},$score_findtar{$id}->{end},$score_findtar{$id}->{utr},$score_findtar{$id}->{utrlen},$score_findtar{$id}->{energy},$score_findtar{$id}->{score},$score_findtar{$id}->{miranda},$score_findtar{$id}->{microtar};
print  "\n";
}


__END__
my $main = MainWindow->new;
$main->Button(-text => 'RUN', -command => \&aboutme,)->pack;
$main->Button(-text => 'About me',-command => \&aboutme,)->pack(-side => "bottom");

sub aboutme {
	my $main_about = $main -> Toplevel;
	my $aboutme = $main_about -> Text ( -width => 40, -height => 2)-> pack;
	$aboutme->insert('end', "Welcome to use this widget!\nAuthour: blinfu\@gmail.com");
	}

MainLoop;
