#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -f <fasta> \

INFO

my ($reffile1, $reffile2, $outfile);
GetOptions(
"f=s"=>\$reffile1,
"i=s"=>\$reffile2,
"o=s"=>\$outfile,
);

die "$usage" if(!$reffile1); 

open (INDEX, "<$reffile1");
while (<INDEX>) {
	if ($_ =~ /^>/){
		#my @spstr = split("\t", $_);
		close (CHROUT);
		chomp();
		my $chrname = substr($_, 1);
		open (CHROUT, ">$chrname.fa");
		print CHROUT $_ . "\n";
	}else{
		print CHROUT $_;
	}
    
}
close (INDEX);
