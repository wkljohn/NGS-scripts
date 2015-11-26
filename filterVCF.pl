#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <input vcf file> -o <out file name> -d <min depth> -q <min GQ>

INFO

my ($refName, $outname, $mindp, $mingq);
GetOptions(
"i=s"=>\$refName,
"o=s"=>\$outname,
"d=i"=>\$mindp,
"q=i"=>\$mingq,
);

die "$usage" if(!$refName || !$mindp || !$mingq); 


#my %harshALLVCFPOS;
my @refRange;/Users/wkljohn/Desktop/Other scripts/filterVCF.pl
my $refharshedcnt=0;
my $filteredcnt = 0;

open (REFFI, "<$refName");
open (OUTFI, ">$outname");
while (<REFFI>) {
    if ($_ =~ m/^\#/){
        print OUTFI $_ ;
        next;
    }
    
    my $ln = $_;
    my @splitarr = split /\t/, $_;
    chomp($splitarr[0]);
    
#$harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[0] = $splitarr[5];
#$harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[1] = $splitarr[7];
#$harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[2] = $splitarr[9]; 
#$harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[3] = 0;
    
    my @splitDP = split /;/, $splitarr[7];
    my @splitGQ = split /:/, $splitarr[9];
    my $DP = substr($splitDP[0], 3);
    my $GQ = $splitGQ[1];
    chomp($GQ);
    
    if ($mindp > $DP){
        $filteredcnt++;
    }elsif ($mingq > $GQ){
        $filteredcnt++;
    }else{
        #print $harsharr[1] . " " . $harsharr[2] . " " . "$DP $GQ";
        print OUTFI $_ ;
    }
    
    $refharshedcnt++;

}
print "filtered:$filteredcnt of $refharshedcnt\n";
close(REFFI);
close(OUTFI);




