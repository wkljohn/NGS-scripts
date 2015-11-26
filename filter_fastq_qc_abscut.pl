#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <read1> -o <output file> -q <qc threshold> -p <percent over threshold 0-100>

INFO

my ($read1, $output, $qcthres, $pccut);

GetOptions(
"i=s"=>\$read1,
"o=s"=>\$output,
"q=i"=>\$qcthres,
"p=i"=>\$pccut,
);

die "$usage" if(!$read1 || !$output || !$qcthres  || !$pccut); 

print "processing: $read1 => $output\n";

if (!(-e $read1)){
	die " read file does not exist!\n";
}

my $keycnt1=0;
my $readswrote1 = 0;
my $readscnt1 = 0;
my $pcfloat = $pccut/100;

open (READFI, "<$read1");
open (READWRFI, ">$output");
while (<READFI>){
	my $ln1 = $_;
	my $ln2 = <READFI>;
	my $ln3 = <READFI>;
	my $ln4 = <READFI>;
	chomp($ln4);
	
	my @qcchars = split '', $ln4;
	my $overthrescnt = 0;
	my $redlen = $#qcchars + 1;
	for (my $i = 0; $i <= $#qcchars; ++$i){
		#$sumscore += ord(substr($ln4, $i, 1)) - 33;
		if (ord($qcchars[$i]) - 33 > $qcthres){
			++$overthrescnt;
		}
	}
	#print "\n$ln4\n";
	#print $avgscore ."\n";
	
	$readscnt1++;
	if ($overthrescnt / $redlen >= $pcfloat){
		$readswrote1++;
		print READWRFI $ln1;
		print READWRFI $ln2;
		print READWRFI $ln3;
		print READWRFI $ln4 . "\n";
	}
	
	#debug
	#if ($readscnt1 > 100){
	#	last;
	#}
}
close(READFI);
close(READWRFI);


print "Read 1 output/input/percent: $readswrote1/$readscnt1/". int($readswrote1/$readscnt1 * 100) . "\n";



