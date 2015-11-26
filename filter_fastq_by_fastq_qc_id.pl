#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -k <key file> -i <input fastq> -o <output file>

INFO

my ($keyfile, $output, $fastqin);

GetOptions(
"k=s"=>\$keyfile,
"o=s"=>\$output,
"i=s"=>\$fastqin,
);

die "$usage" if(!$keyfile); 

print STDERR "processing: $keyfile \n";

if (!(-e $keyfile)){
	die " read file does not exist!\n";
}

my %inckeys;
my $keycnt1=0;
my $readswrote1 = 0;
my $readscnt1 = 0;


open (READFI, "<$keyfile");
print STDERR "Reading key file\n";
while (<READFI>){
	chomp();
	++$keycnt1;
	#print "K:".$_."\n";
	$inckeys{$_} = 1;
}
close (READFI);

print STDERR "$keycnt1 keys read \n";
print STDERR "Processing fastq file\n";
my $keyfdcnt = 0;
open (READWRFI, ">$output");
while (<STDIN>){
	chomp();
	my @splitarr = split(/\s/, $_);
	my $keyname = substr($splitarr[0], 1);

	#++$keyfdcnt;
	#if ($keyfdcnt % 1000 == 0){
	#	print "$keyfdcnt\n";
	#}	

	if (exists $inckeys{$keyname}){
		++$readswrote1;
		print READWRFI $_ . "\n";
		my $tmp = <STDIN>;
		print READWRFI $tmp;
		$tmp = <STDIN>;
		print READWRFI $tmp;
		$tmp = <STDIN>;
		print READWRFI $tmp;
		#print "$keyname\n";
	}else{
		<STDIN>;
		<STDIN>;
		<STDIN>;
		#print "NP:$keyname\n";
	}
}
print STDERR "done\n";
close(READWRFI);


print STDERR "keycnt: $keycnt1\n";
print STDERR "readswrote: $readswrote1\n";



