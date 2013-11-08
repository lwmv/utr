#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use autodie;
use Tie::File;


my @utr=glob "./utr/*";
my @mirna = glob "./mirna/*";

#### first run copredict ####
#############################
foreach my $mirna (@mirna) {
    (my $namirna = basename ($mirna)) =~ s/$/_copredict/;

    mkdir "$namirna",0755 or die "can't mkdir $mirna : $!";
    
    my $total = scalar(@utr);
    my $current = 0;
    foreach my $utr (@utr) {
        $current += 1;
        print basename($mirna),": ",$current, "/", "$total";
        my $nautr = basename ($utr);
        print `uptime`;
        system "perl copredict.pl $mirna $utr |tee ./$namirna/$nautr";
        #system "perl copredict.pl $mirna $utr > ./$namirna/$nautr";

    }
    $current = 0;
}

print `uptime`;

#### second run btar ####
#########################
foreach (@mirna){

    my $co=basename($_)."_copredict";
    my $btar=basename($_)."_btar";
    mkdir "$btar" or die "Can't mkdir $btar : $!";

    foreach (glob "./$co/*") {
        my $na=basename ($_);
        system "perl Btar.pl $_ > ./$btar/$na"
    }
}



#### third fun fm　####
#######################
mkdir "./fm" or die "Can't make ./fm : $!";
foreach(@mirna){
    my $pre = basename($_);

    my $btar="$pre"."_btar";
    open (my $fm, ">", "./fm/$pre.txt");
    #    or die "Can't open $_.fm.txt : $!";
    # my @fm = glob "./$btar/*";
    #my @result = map `perl fm.pl $_`, glob "./$btar/*";
    my @result;
    foreach(glob "./$btar/*"){
        #(my $utr_name= basename $_)=~s/\.\w+$//;
        my @array;
        tie @array, 'Tie::File', $_;
        my @temp;
        push @temp, (split /\s+/,`perl fm.pl $_`)[0,8,9];
        #print "@temp";
        
        #open(my $fh, "<", $_);
        my $length = (split /\s+/,$array[0])[6];
        #close $fh;
        splice( @temp,1,0,$length,scalar @array );

        push @result, join("\t",@temp)."\n";
    }
    
    print {$fm} "UTR\tlength\ttotalTargets\taverageEnergy\taverageScore\n";
    print {$fm} @result;
    close $fm;
}



#### fourth, make table from results of Btar.pl ####
####################################################
mkdir "./table" or die "Can't mkdir ./talbe :$!";
foreach(@mirna){
    my $pre = basename($_);
    open( my $mirna, "<", "./mirna/$pre");
    my @mirnas = sort map {chomp( my $x=$_ );$x=~s/>//; $x } grep {/^>/} <$mirna>;
    my %mirnas = mkhash( \@mirnas, glob "./${pre}_btar/*" );
    mktable( \%mirnas, \@mirnas, "${pre}_table" );
}

=pod
#读入microRNA名
open (my $rna1, "<", "./mirna/diff")
open (my $rna2, "<", "./mirna/undiff")
my @ran1 = map {chomp( my $x=$_ );$x=~s/>//; $x } grep {/^>/} <$rna1>;
my @rna2 = map {chomp( my $x=$_ );$x=~s/>//; $x } grep {/^>/} <$rna2>;
# print "$_\n" foreach @naofrna_hypo;

#构造哈希表，其结构类似 %hypoxia{utr}{microrna}=target_number
my %rna1 = mkhash( \@rna1, glob "./hypoxia_btar/*" );
my %rna2 = mkhash( \@rna2, glob "./normal_btar/*" );

#输出表
#mktable( \%hypoxia, "table_hypo");
#mktable( \%normal, "table_norm");
#存储前输出测试
#print Dumper(\%hypoxia);
###把这俩哈希表store至硬盘文件
#store [\%hypoxia, \%normal], 'table_2';

=cut

sub mkhash{
    my ($na,@file) = @_;
    my %hash;
    foreach ( @file ) {
        next if -d $_;
        my $fname = basename $_;
        my ($name) = $fname =~ /(.+)\./;
        $hash{$name}{$_} = 0 foreach @$na; ##初始化哈希表，对每个utr序列，为每个microrna建立一个键值。
        # print "$name\n";
        open (my $fh, "<", "$_") or die "Can't open $_ : $!";
        while( <$fh> ){
            chomp;
            #print "$_\n";
            my $third = (split /\s+/)[2];
            #print "$third\n";
            $hash{$name}{$third} += 1 if $third;
        }
    }
    return %hash;
}


#构造输出表
sub mktable {

    my ($h_mirna, $a_mirna, $name)=@_;

    my @utr=sort keys %$h_mirna;
    my @mirna=@$a_mirna;

    open (my $fh, ">>", "./table/${name}.txt");

    #初始化表头
    print {$fh} "$name\t";
    print {$fh} "$_\t" foreach @utr;
    print {$fh} "\n";
    
    #具体数据
    foreach my $m (@mirna){
        print {$fh} "$m\t";
        foreach my $u (@utr) {
            print {$fh} "$h_mirna->{$u}{$m}\t"
        }
        print {$fh} "\n";
    }
}



#### compute the folloing fields:
# total targets
# effective targets
# average score
# average energy
#################################


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

