#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <input vcf file> -r <mt input> -o <out file name> -d <min depth> -q <min GQ> -v <vcf out>

INFO

my ($refName, $mtfile, $outname, $mindp, $mingq, $vcfout);
GetOptions(
"i=s"=>\$refName,
"r=s"=>\$mtfile,
"o=s"=>\$outname,
"d=i"=>\$mindp,
"q=i"=>\$mingq,
"v=s"=>\$vcfout,
);

die "$usage" if(!$refName || !$mtfile || !$mindp || !$mingq); 


my %harshALLVCFPOS;
my @refRange;
my $refharshedcnt=0;


open (REFFI, "<$refName");
while (<REFFI>) {
    if ($_ =~ m/^\#/){
        next;
    }
    
    my $ln = $_;
    my @splitarr = split /\t/, $_;
    chomp($splitarr[0]);
    
    $harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[0] = $splitarr[5];
    $harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[1] = $splitarr[7];
    $harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[2] = $splitarr[9]; 
    $harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[3] = 0;
    if ($vcfout){
        $harshALLVCFPOS{$splitarr[0]}{$splitarr[1]}[4] = $_;
    }
    
    $refharshedcnt++;
}
print "harshed:$refharshedcnt\n";
close(REFFI);

open (MYFILE, "<$mtfile");
open (OUTFI, ">$outname");
if ($vcfout){
    open (OUTFIV, ">$vcfout");
}
print "Loading mutation taster tsv\n";
my $unmappedvarientscnt = 0;
my $loadedcnt = 0;
my $filteredcnt = 0;
<MYFILE>; #ignores first line
while (<MYFILE>) {
    if (!($_ =~ m/^\#/)){
    	$loadedcnt++;
    	  chomp();
        my @splitarr = split /\t/, $_;
        my @chrposproc = split / /, $splitarr[1];
        my $chrpos = $chrposproc[0];	#strange chr pos format
        if ($splitarr[0] eq "23"){ $splitarr[0] = "X";}
        if ($splitarr[0] eq "24"){ $splitarr[0] = "Y";}
        
        my $chrname = "chr".$splitarr[0];
        my $resourcestr = "";
        my @harsharr;
        my @splitDP;
        my @splitGQ;
        my $successmap=0;
        
        if (exists $harshALLVCFPOS{$chrname}{$chrpos}){
      		@harsharr = @{$harshALLVCFPOS{$chrname}{$chrpos}};
      		$harshALLVCFPOS{$chrname}{$chrpos}[3] = 1;
      		$successmap=1;
        	#@splitDP = split /;/, $harshALLVCFPOS{$chrname}{$chrpos}[1];
        	#$resourcestr = $harshALLVCFPOS{$chrname}{$chrpos}[0] ."\t". $harshALLVCFPOS{$chrname}{$chrpos}[1] ."\t". $harshALLVCFPOS{$chrname}{$chrpos}[2];
        }elsif (exists $harshALLVCFPOS{$chrname}{$chrpos - 1}){
        	@harsharr = @{$harshALLVCFPOS{$chrname}{$chrpos - 1}};
        	$successmap=1;
        	#$resourcestr = $harshALLVCFPOS{$chrname}{$chrpos - 1}[0] ."\t". $harshALLVCFPOS{$chrname}{$chrpos - 1}[1] ."\t". $harshALLVCFPOS{$chrname}{$chrpos - 1}[2];
        }elsif (exists $harshALLVCFPOS{$chrname}{$chrpos + 1}){
        	$successmap=1;
        	@harsharr = @{$harshALLVCFPOS{$chrname}{$chrpos + 1}};
        	#$resourcestr = $harshALLVCFPOS{$chrname}{$chrpos + 1}[0] ."\t". $harshALLVCFPOS{$chrname}{$chrpos + 1}[1] ."\t". $harshALLVCFPOS{$chrname}{$chrpos + 1}[2];
        }else{
        	$unmappedvarientscnt++;
          print "$chrname\t$chrpos unmapped!\n";
        }
        
        if ($successmap > 0){
        	@splitDP = split /;/, $harsharr[1];
        	@splitGQ = split /:/, $harsharr[2];
        	my $DP = substr($splitDP[0], 3);
        	my $GQ = $splitGQ[1];
        	chomp($GQ);
        	
            if ($mindp > $DP){
                $filteredcnt++;
            }elsif ($mingq > $GQ){
                $filteredcnt++;
            }else{
                #print $harsharr[1] . " " . $harsharr[2] . " " . "$DP $GQ";
                $resourcestr = $harsharr[0] ."\t". $harsharr[1] ."\t". $harsharr[2];
                print OUTFI $_ . "\t" . $resourcestr . "\n";
                if ($vcfout){
                    print OUTFIV $harsharr[4];
                }
        	}
            
        }
    }
}
print "unmapped:$unmappedvarientscnt of $loadedcnt\n";
print "filtered:$filteredcnt \n";
close(MYFILE);
close(OUTFI);
if ($vcfout){
    close (OUTFIV);
}



