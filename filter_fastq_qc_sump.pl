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
"q=f"=>\$qcthres,
);

die "$usage" if(!$read1 || !$output || !$qcthres  ); 

print "stringency: probabality sum <=$qcthres\n";
print "processing: $read1 => $output\n";

if (!(-e $read1)){
	die " read file does not exist!\n";
}

my $keycnt1=0;
my $readswrote1 = 0;
my $readscnt1 = 0;
my $tmpint;

my $sumallp = 0.0;

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
	my $tmpsump = 0.0;
	for (my $i = 0; $i <= $#qcchars; ++$i){
		$tmpint = ord($qcchars[$i]) - 33;
		$tmpsump += 10 ** ($tmpint / -10);
		#print 10 ** ($tmpint / -10) . ",";
	}
	$sumallp += $tmpsump;
	#printf("sum: %.4f\n", $tmpsump);
	
	$readscnt1++;
	if ($tmpsump <= $qcthres){
		$readswrote1++;
		print READWRFI $ln1;
		print READWRFI $ln2;
		print READWRFI $ln3;
		print READWRFI $ln4 . "\n";
	}
	
	
	
	#debug
	#if ($readscnt1 % 100000 == 0){
	#	printf("average p: %.3f\n", $sumallp / $readscnt1);
	#	print "proportion wrote: ". int($readswrote1/$readscnt1 * 100) . "\n";
	#}
	#if ($readscnt1 > 100){
	#	last;
	#}
}
close(READFI);
close(READWRFI);

print "average p:" . int($sumallp / $readscnt1 ) . "\n";
print "Read 1 output/input/percent: $readswrote1/$readscnt1/". int($readswrote1/$readscnt1 * 100) . "\n";



