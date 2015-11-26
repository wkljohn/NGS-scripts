#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <input annovar with reference file> -r <reference FPKM file> -o <output> -c <cut off FPKM> -x <optional file to fix the vcf>

INFO

my ($refName, $bedfile, $output, $cutoffexp, $fixref);
GetOptions(
"i=s"=>\$bedfile,
"r=s"=>\$refName,
"o=s"=>\$output,
"c=i"=>\$cutoffexp,
"f=s"=>\$fixref,
);

die "$usage" if(!$refName || !$bedfile); 

my %genesharsh;
my %fixrefharsh;
my %harshexp;
my %harshrangeLOW;
my %harshrangeNAME;
my @refRange;
my $refindex=0;
my $harshcnt=0;

#open reference file
open (REFFI, "<$refName");
while (<REFFI>) {
    if ($_ =~ /^gene_id/){
        next;
    }
    chomp();
    my @splitarr = split /\t/, $_;
    my $genename = $splitarr[0];
    for (my $i = 3; $i <= $#splitarr; ++$i){
        push @{$harshexp{$genename}}, $splitarr[$i];
    }
    ++$harshcnt;
    
}
close (REFFI);
print "harshed $harshcnt genes\n";

if ($fixref){
    open (FIXFI, "<$fixref");
    while (<FIXFI>) {
        if ($_ =~ /^#/){
            next;
        }
        chomp();
        my @splitarr = split /\t/, $_;
    
    #print "harshed:" . $splitarr[0] . "-" . $splitarr[1] . "-" . $splitarr[3] . "\n";
        $fixrefharsh{$splitarr[0]}{$splitarr[1]}{$splitarr[3]} = $_;
        
    }
    close (FIXFI);
}

my $missingcnt = 0;
my $refcnt = 0;
my $harshfailedcnt = 0;
my $norefcnt = 0;
open (OUTFI2, ">$output");
open (MYFILE, "<$bedfile");
print "Loading annovar\n";
while (<MYFILE>) {
    chomp($_);
    my @splitarr = split /\t/, $_;
    if ($fixref){
        if (exists $fixrefharsh{$splitarr[7]}{$splitarr[8]}{$splitarr[10]}){
            my @splitfixref = split /\t/, $fixrefharsh{$splitarr[7]}{$splitarr[8]}{$splitarr[10]};
            my $cntx=0;
            for (my $i = 0; $i <= $#splitfixref; ++$i){
                $splitarr[$i + 7] = $splitfixref[$i];
                ++$cntx;
            }
            print "arr len: $cntx of " . $#splitarr . "\n";
            #@splitarr = split /\t/, $_;
        }else{
            
            for (my $i = 8; $i <= $#splitarr; ++$i){
                $splitarr[$i - 1] = $splitarr[$i];
            } 
            if (exists $fixrefharsh{$splitarr[7]}{$splitarr[8]}{$splitarr[10]}){
                my @splitfixref = split /\t/, $fixrefharsh{$splitarr[7]}{$splitarr[8]}{$splitarr[10]};
                for (my $i = 0; $i <= $#splitfixref; ++$i){
                    $splitarr[$i + 7] = $splitfixref[$i];
                }
            }else{
                print "unmapped: ".$splitarr[7] . "-" . $splitarr[8]. "-" .$splitarr[10] . $_ . "\n";
                $norefcnt++;
            }
            
            print "\n";
            for (my $i = 0; $i <= $#splitarr; ++$i){
                print $splitarr[$i] . "\t" ;
            }
            
        }
    }
    chomp($splitarr[0]);
        my $type;
        my $genename;
        my $vcfhead = "";
        my $vcfgeno = "";
        my @genotypearr;
        
        $genename = $splitarr[1];
        
        for (my $i = 16; $i <= $#splitarr; ++$i){
            chomp($splitarr[$i]);
            push @genotypearr, $splitarr[$i];
        }
        for (my $i = 7; $i <= 15; ++$i){
            $vcfhead .= $splitarr[$i] . "\t";
        }
        
        #find if isref, or ismissing 
        for (my $i = 0; $i <= $#genotypearr; ++$i){
            if ($genotypearr[$i] eq "."){
                if (!exists $harshexp{$genename}){
                    $vcfgeno .= $genotypearr[$i] . "\t";
                    ++$harshfailedcnt;
                }elsif ($harshexp{$genename}[$i] >= $cutoffexp){
                    my $score = $harshexp{$genename}[$i] * 2;
                    if ($score > 255){
                        $score = 255;
                    }
                    $vcfgeno .= "0/0:$score\t"; #$genotypearr[$i] . "\t";
                    ++$refcnt;
                }else{
                    $vcfgeno .= $genotypearr[$i] . "\t";
                    ++$missingcnt;
                }
            }else{
                $vcfgeno .= $genotypearr[$i] . "\t";
            }
        }
        
        
        print OUTFI2 $vcfhead . $vcfgeno . "\n";
        
        #if (length($type) > 0){
        #    print OUTFI2 $splitarr[0] . "\thg19_refGene\t" . $type . "\t" . $splitarr[1] . "\t"  . $splitarr[2] . "\t" . $splitarr[3] . "\t" . $splitarr[10] . "\t" . $splitarr[7] .  "\t" . $splitarr[11] . "\n";
        #}
}
close(MYFILE);
close(OUTFI2);

print "As missing: $missingcnt\n";
print "As reference: $refcnt\n";

if ($fixref){
    print "Fix ref err:$norefcnt\n";
}
