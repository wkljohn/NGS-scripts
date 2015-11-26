#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -b  <G_B_bam> -c <G_C bam> -d <G_C2_bam> -o <outbam> 

INFO

my ($bam_gc, $bam_gb, $bam_gc2, $outfile);

GetOptions(
"c=s"=>\$bam_gc,
"b=s"=>\$bam_gb,
"d=s"=>\$bam_gc2,
"o=s"=>\$outfile,
);

die "$usage" if(!$bam_gc || !$bam_gb || !$bam_gc2 ); 

my %bamhar_gc;
my %bamhar_gc2;
my %bamhar_gb;
my %allkeys;


if (!(-e $bam_gc)){
	die " read file GC does not exist!\n";
}

if (!(-e $bam_gb)){
	die " read file GB does not exist!\n";
}

if (!(-e $bam_gc2)){
	die " read file GC2 does not exist!\n";
}

sub scorerecord{
	my $read1 = shift;
	my $read2 = shift;
	my $score = 0;
	
	my @splitr1 = split "\t", $read1;
	my @splitr2 = split "\t", $read2;
	my $flag1 = $splitr1[1];
	my $flag2 = $splitr2[1];
	
	if (!($flag1 & 0x4)){
		$score += 10;
	}
	if (!($flag2 & 0x4)){
		$score += 10;
	}
	
	return $score;
}

#my ($bam_gc, $bam_gb, $bam_gc2);
my @bamcontent_gc = `samtools view $bam_gc`;
my $vircnt = 0;
my $blkcnt = 0;
my $pair ;
foreach my $line (@bamcontent_gc){
	if ($line !~ "^@"){
		my @splitarr = split "\t", $line;
		my $flags = $splitarr[1];
		
		#if ($splitarr[2] ne "*"){ 
			if ($flags & 0x40){	#1st read
				$pair = 1;
			}else{
				$pair = 2;
			}
			
			#see if its blocked
			if ($flags & 0x100){
				$bamhar_gc{$splitarr[0] . "\\1"} = "B";
				$bamhar_gc{$splitarr[0] . "\\2"} = "B";
				++$blkcnt;
			}else{
				if (exists $bamhar_gc{$splitarr[0] . "\\$pair"} &&  $bamhar_gc{$splitarr[0] . "\\$pair"} eq "B"){
					++$blkcnt;
				}else{
					$bamhar_gc{$splitarr[0] . "\\$pair"} = $line;
					$allkeys{$splitarr[0]} = 1;
					++$vircnt;
				}
			}
			
	}
}
print "viral reads GC loaded: $vircnt, $blkcnt blocked\n";

my @bamcontent_gc2 = `samtools view $bam_gc2`;
$vircnt = 0;
$blkcnt = 0;
foreach my $line (@bamcontent_gc2){
	if ($line !~ "^@"){
		my @splitarr = split "\t", $line;
		my $flags = $splitarr[1];
		
		#if ($splitarr[2] ne "*"){ 
			if ($flags & 0x40){	#1st read
				$pair = 1;
			}else{
				$pair = 2;
			}
			
			#see if its blocked
			if ($flags & 0x100){
				$bamhar_gc2{$splitarr[0] . "\\1"} = "B";
				$bamhar_gc2{$splitarr[0] . "\\2"} = "B";
				++$blkcnt;
			}else{
				if (exists $bamhar_gc2{$splitarr[0] . "\\$pair"} &&  $bamhar_gc2{$splitarr[0] . "\\$pair"} eq "B"){
					++$blkcnt;
				}else{
					$bamhar_gc2{$splitarr[0] . "\\$pair"} = $line;
					$allkeys{$splitarr[0]} = 1;
					++$vircnt;
				}
			}
		#}
	}
}
print "viral reads GC2 loaded: $vircnt, $blkcnt blocked\n";


my @bamcontent_gb = `samtools view $bam_gb`;
$vircnt = 0;
$blkcnt = 0;
foreach my $line (@bamcontent_gb){
	if ($line !~ "^@"){
		my @splitarr = split "\t", $line;
		my $flags = $splitarr[1];
		
		#if ($splitarr[2] ne "*"){ 
			if ($flags & 0x40){	#1st read
				$pair = 1;
			}else{
				$pair = 2;
			}
			
			#see if its blocked
			if ($flags & 0x100){
				$bamhar_gb{$splitarr[0] . "\\1"} = "B";
				$bamhar_gb{$splitarr[0] . "\\2"} = "B";
				++$blkcnt;
			}else{
				if (exists $bamhar_gb{$splitarr[0] . "\\$pair"} &&  $bamhar_gb{$splitarr[0] . "\\$pair"} eq "B"){
					++$blkcnt;
				}else{
					$bamhar_gb{$splitarr[0] . "\\$pair"} = $line;
					$allkeys{$splitarr[0]} = 1;
					++$vircnt;
				}
			}
		#}
	}
}
print "viral reads GB loaded: $vircnt, $blkcnt blocked\n";

 
open(OPENFI, ">$outfile");
my $gc_cnt=0;
my $gc2_cnt=0;
my $gb_cnt=0;
foreach my $key (keys %allkeys){
	my $gc_score=0;
	my $gc2_score=0;
	my $gb_score=0;
	
	if (exists $bamhar_gc{$key . "\\1"} && $bamhar_gc{$key . "\\1"} ne "B"){
		$gc_score = scorerecord($bamhar_gc{$key . "\\1"});
	}
	if (exists $bamhar_gc2{$key . "\\1"} && $bamhar_gc2{$key . "\\1"} ne "B"){
		$gc2_score = scorerecord($bamhar_gc2{$key . "\\1"});
	}
	if (exists $bamhar_gb{$key . "\\1"} && $bamhar_gb{$key . "\\1"} ne "B"){
		$gb_score = scorerecord($bamhar_gb{$key . "\\1"});
	}
	
	#print "scores: $gc_score, $gc2_score, $gb_score\n";
	
	if ($gc_score > $gc2_score && $gc_score > $gb_score){
		print OPENFI  $bamhar_gc{$key . "\\1"};
		print OPENFI  $bamhar_gc{$key . "\\2"};
		$gc_cnt++;
		
	}elsif ($gc2_score > $gb_score){
		print OPENFI  $bamhar_gc2{$key . "\\1"};
		print OPENFI $bamhar_gc2{$key . "\\2"};
		$gc2_cnt++;
	
	}elsif ($gb_score > 0){
		print OPENFI $bamhar_gb{$key . "\\1"};
		print OPENFI $bamhar_gb{$key . "\\2"};
		$gb_cnt++;
	
	}else{
		print "SCORE NOT FOUND!\n";
	}
}
print "GC used: $gc_cnt\n";
print "GC2 used: $gc2_cnt\n";
print "GB used: $gb_cnt\n";

close(OPENFI);
