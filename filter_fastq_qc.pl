#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <read1> -o <output file> -q <qc threshold>

INFO

my ($read1, $output, $qcthres);

GetOptions(
"i=s"=>\$read1,
"o=s"=>\$output,
"q=i"=>\$qcthres,
);

die "$usage" if(!$read1 || !$output || !$qcthres); 

print "processing: $read1 => $output\n";

if (!(-e $read1)){
	die " read file does not exist!\n";
}

my $keycnt1=0;
my $readswrote1 = 0;
my $readscnt1 = 0;

open (READFI, "<$read1");
open (READWRFI, ">$output");
while (<READFI>){
	my $ln1 = $_;
	my $ln2 = <READFI>;
	my $ln3 = <READFI>;
	my $ln4 = <READFI>;
	chomp($ln4);
	
	my $sumscore = 0;
	my @qcchars = split '', $ln4;
	for (my $i = 0; $i <= $#qcchars; ++$i){
		#$sumscore += ord(substr($ln4, $i, 1)) - 33;
		$sumscore += ord($qcchars[$i]) - 33;
		#print ord(substr($ln4, $i, 1)) - 33 . ",";
	}
	my $avgscore = $sumscore / length($ln4);
	#print "\n$ln4\n";
	#print $avgscore ."\n";
	
	$readscnt1++;
	if ($avgscore > $qcthres){
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



