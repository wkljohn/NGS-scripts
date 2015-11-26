#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -1 <read1> -2 <read2> -x <output ID only>

INFO

my ($read1, $read2, $idonly);
$idonly = 0;
GetOptions(
"1=s"=>\$read1,
"2=s"=>\$read2,
"x=i"=>\$idonly,
);

die "$usage" if(!$read1 || !$read2); 

print "processing: $read1, $read2\n";

if (!(-e $read1) || !(-e $read2)){
	die "either read file does not exist!\n";
}

my %harsh1;
my %harsh2;
my $keycnt1=0;
my $keycnt2=0;

open (READFI, "<$read1");
while (<READFI>){
	chomp();
	my @splitread = split(" ", $_);
	if (!(exists $harsh1{$splitread[0]})){
	    ++$keycnt1;
		$harsh1{$splitread[0]} = 1;
	}else{
		print "duplicated key 1:" . $splitread[0];
	}
	<READFI>;
	<READFI>;
	<READFI>;
}
close(READFI);

print "key1: $keycnt1 \n";
open (READFI, "<$read2");
while (<READFI>){
	chomp();
	my @splitread = split(" ", $_);
	if (!(exists $harsh2{$splitread[0]})){
	    ++$keycnt2;
		$harsh2{$splitread[0]} = 1;
	}else{
		print "duplicated key 2:" . $splitread[0];
	}
	<READFI>;
	<READFI>;
	<READFI>;
}
close(READFI);

print "key2: $keycnt2  \n";

my $readswrote1 = 0;
my $readswrote2 = 0;
my $readscnt1 = 0;
my $readscnt2 = 0; 

open (READFI, "<$read1");
if ($idonly){
	open (READWRFI, ">$read1.fixmate.idenc");
}else{
	open (READWRFI, ">$read1.fixmate.fastq");
}
while (<READFI>){
	my @splitread = split(" ", $_);
	++$readscnt1;
	
	if (exists $harsh2{$splitread[0]}){
		++$readswrote1;
		if ($idonly){
			print READWRFI substr($splitread[0],1) . "\n";
			<READFI>;
			<READFI>;
			<READFI>;
		}else{
			print READWRFI $_;
			<READFI>;
			print READWRFI $_;
			<READFI>;
			print READWRFI $_;
			<READFI>;
			print READWRFI $_;
		}
	}else{
		<READFI>;
		<READFI>;
		<READFI>;
	}
}
close(READFI);
close(READWRFI);
if ($idonly){
	print "Read 1 output/input/percent: $readswrote1/$readscnt1/". int($readswrote1/$readscnt1 * 100) . "\n";
	die "idonly finished\n";
}


open (READFI, "<$read2");
open (READWRFI, ">$read2.fixmate.fastq");
while (<READFI>){
	++$readscnt2;
	my @splitread = split(" ", $_);
	if (exists $harsh1{$splitread[0]}){
		++$readswrote2;
		print READWRFI $_;
		<READFI>;
		print READWRFI $_;
		<READFI>;
		print READWRFI $_;
		<READFI>;
		print READWRFI $_;
	}else{
		<READFI>;
		<READFI>;
		<READFI>;
	}
}
close(READFI);
close(READWRFI);

print "Read 1 output/input/percent: $readswrote1/$readscnt1/". int($readswrote1/$readscnt1 * 100) . "\n";
print "Read 2 output/input/percent: $readswrote2/$readscnt2/". int($readswrote2/$readscnt2 * 100) . "\n";



