#!/usr/bin/perl -w
#######################################################################################
#  only receive data format generate by copredict.pl; 
#  input is the original and generate the format like below:
#   line ->  original_line;
#   New_ID   Id     Microrna    Start  End   Utr   Utrlen Energy Score Miranda Microtar
#  '32' -> [ '44','hsa-miR-210','686','718','VEGF','1796','-22.6','20','yes','no']
# $VAR3 = '32';                                                           
# $VAR4 = [                                                               
#           '44',                                                         
#           'hsa-miR-210',                                                
#           '686',                                                        
#           '718',                                                        
#           'VEGF',                                                       
#           '1796',                                                       
#           '-22.6',                                                      
#           '20',                                                         
#           'yes',                                                        
#           'no'                                                          
#         ];    
# %outar <==== screening results
# %intar <==== original input
########################################################################################
use strict;
use warnings;
use Data::Dumper;

##选出得分最高的靶点ID
sub btar {
    my %tar = @_;
    my $max_score = 0;
    my $min_energy = 0;
    my $max_id;

    foreach my $tar_id ( keys %tar ){
        my $score = $tar{$tar_id}->[7];
        my $energy = $tar{$tar_id}->[6];

        if ( $score > $max_score ){
            #print "===$score====$max_score\n";
            $max_id = $tar_id;
            $max_score = $tar{$tar_id}->[7];
            $min_energy = $tar{$tar_id}->[6];
        }
        else{
            if ( $score < $max_score ) {
            #什么都不做
            }
            else{
                if ( $score = $max_score ) {
                    if ( $energy >= $min_energy ) {
                    #在score相等的情况下，如果energy>=最小能量值，则什么也不做
                    }else{
                        #如果energy小于最小能量值
                        $max_id = $tar_id;
                        $max_score = $tar{$tar_id}->[7];
                        $min_energy = $tar{$tar_id}->[6];
                    }
                }
            }
        }
    }
    return $max_id;
}

###############################################################################
#      Id Microrna Start End Utr Utrlen Energy Score Miranda Microtar         #
###############################################################################

open my $file,'<',shift or die $!;
my %intar;
my @targets;
my $id = 1;
##读入原始数据，构造哈希%intar，key为数字序号，值为数组引用[@mirna]
while (my $line = <$file>){ #Id Microrna Start End Utr Utrlen Energy Score Miranda Microtar
    next if $. <= 10;
    $line =~ s/\r|\n//; #＃删换行符？
    my @mirna =  split /\t/,$line;
    $intar{$id} = [ @mirna ];
    $id++;
}


my %outar;

INPUT:
my $mid   = btar(%intar);##选出得分最高的靶点
#print "mid is $mid \n";
my $left  = ( $intar{$mid}->[2] );   #left equals [2] 
$left += 1;  #2-7 seed 区开始
#print "left  is $left  \n";
my $right = ( $intar{$mid}->[3] );   #right equal [3]
$right += 6; #2-7 seed 区结束
#print "right  is $right \n";

$outar{$mid} = $intar{$mid}; #copy; ##记录得分最高的靶点
delete $intar{$mid};##删除得分最高的mirna靶点ID
##删除与得分最高靶点有竞争的靶点
foreach my $to_be_delete ( keys %intar ){
    my $start = $intar{$to_be_delete}->[2] +1;
    my $end = $intar{$to_be_delete}->[2] +6;
    if( ($start >= $left && $start <= $right) or ($end >= $left && $end <= $right) ){
        delete $intar{$to_be_delete};
    } 
#    if ( $intar{$to_be_delete}->[2] >= $left && $intar{$to_be_delete}->[2] <= $right ){
#        delete $intar{$to_be_delete}; ##如果起始区落在left 和 right 之间，则删除。left->right之间即所谓seed 区？？
#        next;
#    }
#
#    if ( $intar{$to_be_delete}->[3] >= $left && $intar{$to_be_delete}->[3] <= $right ){
#        delete $intar{$to_be_delete};
#        next;
#    }

}

my @remain = ( keys %intar );
#print join "<=>",@remain,"\n";
my $remain = $#remain + 1;
#print "remain is $remain! \n";
if ( $remain > 0) { goto INPUT } ##循环，直到%intar为空

#print "Screening finished! \n";

foreach my $screened ( keys %outar ) {
    print join "\t",($screened,@{$outar{$screened}},"\n");
}
