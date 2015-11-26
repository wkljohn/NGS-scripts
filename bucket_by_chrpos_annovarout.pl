#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <position  file>

INFO

my ($fin);
GetOptions(
"i=s"=>\$fin,
);

die "$usage" if(!$fin ); 

my %chrmap;
my %inmap;
my %cntmap;

open (READFI, "<$fin");
while (<READFI>){
	chomp();
	my @splitread = split("\t", $_);
	if (!(exists $inmap{$splitread[0]}{$splitread[2]})){
		print $splitread[0]  . "\t" . $splitread[2] . "\t" . $splitread[3] . "\n";
		push @{$inmap{$splitread[0]}{$splitread[2]}}, ($splitread[3]);
	}else{
		#print "ERR\n";
		push @{$inmap{$splitread[0]}{$splitread[2]}}, $splitread[3];
	}
}
close(READFI);


$chrmap{"chr1"} = 249250621;
$chrmap{"chr2"} = 243199373;
$chrmap{"chr3"} = 198022430;
$chrmap{"chr4"} = 191154276;
$chrmap{"chr5"} = 180915260;
$chrmap{"chr6"} = 171115067;
$chrmap{"chr7"} = 159138663;
$chrmap{"chr8"} = 146364022;
$chrmap{"chr9"} = 141213431;
$chrmap{"chr10"} = 135534747;
$chrmap{"chr11"} = 135006516;
$chrmap{"chr12"} = 133851895;
$chrmap{"chr13"} = 115169878;
$chrmap{"chr14"} = 107349540;
$chrmap{"chr15"} = 102531392;
$chrmap{"chr16"} = 90354753;
$chrmap{"chr17"} = 81195210;
$chrmap{"chr18"} = 78077248;
$chrmap{"chr19"} = 59128983;
$chrmap{"chr20"} = 63025520;
$chrmap{"chr21"} = 48129895;
$chrmap{"chr22"} = 51304566;
$chrmap{"chrX"} = 155270560;
$chrmap{"chrY"} = 59373566;

print "chr\tpos\t";
foreach my $samplekeymain (keys %inmap){
	print "$samplekeymain\t";
}
print "\n";

my %chrcount;
my %chrcountall;

foreach my $key (keys %chrmap){
	$chrcount{$key} = 0;
	for (my $i = 0; $i < $chrmap{$key}; $i += 1000000){
		my $rgnallcnt=0;
		print $key."\t".$i."\t";
		
		foreach my $samplekey (keys %inmap){
			$cntmap{$key}{$i}{$samplekey} = 0;
			#find the reads in region
				my $rgncnt=0;
			if (exists($inmap{$samplekey}{$key})){
				#print "proc:".$samplekey.":".$key."\n";
				my @chrarr = @{$inmap{$samplekey}{$key}};
				for (my $j = 0; $j <= $#chrarr ; ++$j){
					if ($chrarr[$j] >= $i && $chrarr[$j] < $i + 1000000){
						#only if >= 2 counts
						++$rgncnt;
						++$rgnallcnt;
						#++$cntmap{$key}{$i}{$samplekey};
						#if ($rgncnt >= 2){
						#	++$cntmap{$key}{$i};	
						#	last;	#only count once
						#}
					}
				}
			}
			print "$rgncnt\t";
			if ($rgncnt > 1){
				$chrcount{$key}++;
			}
			#if ($cntmap{$key}{$i} > 0){
			#	print "$key\t$i\t". ($i + 1000000 - 1) . "\t" . $cntmap{$key}{$i} . "\n";
			#}
		}
		if ($rgnallcnt > 1){
			$chrcountall{$key}++;
		}
		print "\n";
		#print "buc:".$key.":".$i."-".$cntmap{$key}{$i}."\n";
		#if ($cntmap{$key}{$i} > 0){
		#	print "$key\t$i\t". ($i + 1000000 - 1) . "\t" . $cntmap{$key}{$i} . "\n";
		#}
	}
}



foreach my $chrkey (sort keys %chrcount){
	print "$chrkey\t".$chrcount{$chrkey}."\n";
}
print "-----\n";
foreach my $chrkey (sort keys %chrcountall){
	print "$chrkey\t".$chrcountall{$chrkey}."\n";
}
