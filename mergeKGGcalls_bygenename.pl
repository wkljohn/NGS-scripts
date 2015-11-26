#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <input KGGSeq file> -o <out file>

INFO

my ($refName, $outname);
GetOptions(
"i=s"=>\$refName,
"o=s"=>\$outname,
);

die "$usage" if(!$refName || !$outname); 


my %harshALLVCFPOS;
my @refRange;
my $refharshedcnt=0;
my $indelcnt=0;
my $qcpasscnt=0;


open (REFFI, "<$refName");
while (<REFFI>) {
    if ($_ =~ m/^\#/){
        next;
    }
    
    my $ln = $_;
    my $isindel=0;
    my $dataqcpass=0;
    my @splitarr = split /\t/, $_;
    chomp($splitarr[0]);
    chomp($_);

    if ($splitarr[2] =~ m/\+/ || $splitarr[2] =~ m/\-/){
        $indelcnt++;
        $isindel=1;
        $dataqcpass=1;
    }else{
        #missing data test for non-indel
        if ($splitarr[12] + $splitarr[13] + $splitarr[14] < 4 || $splitarr[15] + $splitarr[16] + $splitarr[17] < 6){
            $dataqcpass = 0;
        }else{
            $dataqcpass = 1;
        }
    }
    
    if ($dataqcpass == 1){
        $qcpasscnt++;
        if (!exists $harshALLVCFPOS{$splitarr[3]}){
            #12-17
            $harshALLVCFPOS{$splitarr[3]}[0] = $splitarr[13] + $splitarr[14];
            $harshALLVCFPOS{$splitarr[3]}[1] = 8-$harshALLVCFPOS{$splitarr[3]}[0];#$splitarr[12];
            $harshALLVCFPOS{$splitarr[3]}[2] = $splitarr[16] + $splitarr[17];
            $harshALLVCFPOS{$splitarr[3]}[3] = 12-$harshALLVCFPOS{$splitarr[3]}[2];#$splitarr[15];
            $harshALLVCFPOS{$splitarr[3]}[4] = $_;
        }else{
            $harshALLVCFPOS{$splitarr[3]}[0] += $splitarr[13] + $splitarr[14];
            $harshALLVCFPOS{$splitarr[3]}[1] -= $harshALLVCFPOS{$splitarr[3]}[0];#$splitarr[12];
            #$harshALLVCFPOS{$splitarr[3]}[1] += $splitarr[12];
            $harshALLVCFPOS{$splitarr[3]}[2] += $splitarr[16] + $splitarr[17];
            #$harshALLVCFPOS{$splitarr[3]}[3] += $splitarr[15];
            $harshALLVCFPOS{$splitarr[3]}[3] -= $harshALLVCFPOS{$splitarr[3]}[2];#$splitarr[15]
    	
        }
    
    }
    $refharshedcnt++;
}
print "harshed:$refharshedcnt\n";
print "filtered:$qcpasscnt\n";
close(REFFI);

my $keycnt = 0 ;
open (OUTFI, ">$outname");
print "printing res\n";
foreach my $key (keys %harshALLVCFPOS){
    if ($harshALLVCFPOS{$key}[1] < 0) { $harshALLVCFPOS{$key}[1] = 0;};
    if ($harshALLVCFPOS{$key}[3] < 0) { $harshALLVCFPOS{$key}[3] = 0;};
    
	print OUTFI $key . "\t" . $harshALLVCFPOS{$key}[0] .  "\t" . $harshALLVCFPOS{$key}[1] .  "\t" . $harshALLVCFPOS{$key}[2] .  "\t" . $harshALLVCFPOS{$key}[3] .  "\t" . $harshALLVCFPOS{$key}[4] . "\n";
	++$keycnt;
}
print "keys:$keycnt \n";
close(OUTFI);



