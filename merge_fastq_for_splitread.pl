#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -f <first file> -s <secound file> -o <out file> 

INFO

my ($reffile1, $reffile2, $outfile);
GetOptions(
"f=s"=>\$reffile1,
"s=s"=>\$reffile2,
"o=s"=>\$outfile,
);

die "$usage" if(!$reffile1 || !$reffile2); 


open (FILEONE, "<$reffile1");
open (FILESEC, "<$reffile2");
open (FOUT, ">$outfile");
my $readcnt=0;
my $combcnt=0;
my $tmpin;
while (1) {
	if ($readcnt < 4){
		$tmpin = <FILEONE>;
		++$readcnt;
	}else{
		$tmpin = <FILESEC>;
		++$readcnt;
		if ($readcnt == 8){
			$readcnt=0;
			++$combcnt;
		}
	}
	if (!$tmpin){
		last;
	}
    print FOUT $tmpin;
	#$print $combcnt;
	#if ($combcnt > 10){
	#	last;
	#}
}
close (FILEONE);
close (FILESEC);
close (FOUT);