#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -k <key file> -o <output file>

INFO

my ($keyfile, $output);

GetOptions(
"k=s"=>\$keyfile,
"o=s"=>\$output,
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
	$inckeys{$_} = 1;
}
close (READFI);

print STDERR "$keycnt1 keys read \n";
print STDERR "Processing bam file\n";
open (READWRFI, ">$output");
while (<STDIN>){
	chomp();
	my @splitseq = split("\t", $_);
	
	if ($_ =~ /^@/){
		print $_ . "\n";
	}elsif (exists $inckeys{$splitseq[0]}){
		++$readswrote1;
		print $_ . "\n";
	}
}
close(READWRFI);


print STDERR "keycnt: $keycnt1\n";
print STDERR "readswrote: $readswrote1\n";



