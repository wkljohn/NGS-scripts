#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <information rich annovar file> -o <old format=1>

INFO

my ($fin);
my $oldformat = 0;
GetOptions(
"i=s"=>\$fin,
"o=i"=>\$oldformat,
);

die "$usage" if(!$fin ); 

sub processread{
	my $cigar = shift;
	my $flags = shift;
	my $splitseqref = shift;	my @splitseq = @$splitseqref;
	
	my $pos = $splitseq[3];
	my $pos_b;
	if ($cigar =~ /^([5-9]|[1-9][0-9]+)S[0-9]+M$/ 
	   || $cigar =~ /^([5-9]|[1-9][0-9]+)S[1-9]I[0-9]+M$/  || $cigar =~ /^([5-9]|[1-9][0-9]+)S[0-9]+M[1-9]I[0-9]+M$/ 
		|| $cigar =~ /^([5-9]|[1-9][0-9]+)S[0-9]+M[1-8]S$/){
		
		my $posadd = substr($cigar, 0, index($cigar, 'S'));
		if ($flags & 0x10){	#reversed
			return "-\t$pos\t$cigar";
		}else{
			return "-\t$pos\t$cigar";
		}
		
		#$pos_b = $pos + $posadd;
		#print "PA:$posadd PB:$pos_b\n";
		#return "-\t$pos_b\t$cigar";
	}elsif($cigar =~ /^[2-9][0-9]M([5-9]|[1-9][0-9]+)S$/ || $cigar =~ /^[1-8]S[2-9][0-9]M([5-9]|[1-9][0-9]+)S$/){
		#fix the begin s
		if ($cigar =~ /^[0-5]S[2-9][0-9]M([5-9]|[1-9][0-9]+)S$/){
			$cigar = substr($cigar, 2);	
		}
		
		my $posadd = substr($cigar, 0, index($cigar, 'M'));
		#$posadd = substr($posadd, 0, length($posadd) - 1);
		
		$pos_b = $pos + $posadd;
		if ($flags & 0x10){	#reversed
			return "+\t$pos_b\t$cigar";
		}else{
			return "+\t$pos_b\t$cigar";
		}
		#$pos_b = $pos + $posadd;
		#print "PA:$posadd PB:$pos_b\n";
		
		#return "+\t$pos_b\t$cigar";
	}elsif($cigar =~ /^[0-9]+M$/ || $cigar =~ /^[0-9]+M[1-4]S$/ || $cigar =~ /^[1-4]S[0-9]+M$/  || $cigar =~ /^[1-4]S[0-9]+M[1-4]S$/ 
			|| $cigar =~ /^[0-9]+M[0-9]?[0-9]D[0-9]+M$/  || $cigar =~ /^[0-9]+M[1-9]I[0-9]+M$/){
		if ($flags & 0x40){	#1st read
			if ($flags & 0x10){	#reversed
				if ($flags & 0x20){	#reversed
					return "+>\t".$pos."\t$cigar";
				}else{				#anti strand pair
					return "-<\t".$pos."\t$cigar";
				}
			}else{				# not reversed
				if ($flags & 0x20){	#reversed		seems fixed by 308
					return "+>\t".($pos+101)."\t$cigar";
				}else{				#anti strand pair
					return "+>\t".($pos+101)."\t$cigar";
				}
			}
		}else{	#2nd read
			if ($flags & 0x10){	#reversed
				if ($flags & 0x20){	#reversed
					return "+>\t".($pos+101)."\t$cigar";
				}else{				#anti strand pair
					return "-<\t".($pos+101)."\t$cigar";
				}
			}else{				# not reversed
				if ($flags & 0x20){	#reversed		seems fixed by 308
					return "+>\t".$pos."\t$cigar";
				}else{				#anti strand pair
					return "+>\t".($pos+101)."\t$cigar";
				}
			}
		}
		
	}else{
full:
		print "Unhdl cigar:$cigar\n";
		return "";
	}
}


my %annovar_read_map;
my %redunmap;

open (READFI, "<$fin");
while (<READFI>){
	chomp();
	my @splitread = split("\t", $_);
	my $readid; 
	if ($oldformat){
		$readid = $splitread[16];
	}else{
		my @splitcolon = split(";", $splitread[12]);
		$readid = $splitcolon[7];
		print $readid . "\n";
		$readid =~ s/RID=//g;
	}
	if (!(exists $annovar_read_map{$readid})){
		$annovar_read_map{$readid} = $splitread[0] . "\t" . $splitread[1] . "\t" . $splitread[2] . "\t" .   $splitread[3];
	}else{
		print "Redundant:$readid\n";
		print $annovar_read_	map{$readid} . "\n";
		print $splitread[0] . "\t" . $splitread[1] . "\t" . $splitread[2] . "\t" .   $splitread[3] . "\n";
		#push @{$inmap{$splitread[0]}{$splitread[2]}}, $splitread[3];
	}
}
close(READFI);


while (<STDIN>){
	chomp();
	my @splitread = split("\t", $_);
	if (exists $annovar_read_map{$splitread[0]}){
		my $flags = int($splitread[1]);
		my $cigar = $splitread[5];
		#print $flags . ",";
		if (!($flags & 0x4) && !($flags & 0x100)){
			my $tmpres = processread($cigar, $flags, \@splitread);
			if ($tmpres ne ""){
				$annovar_read_map{$splitread[0]} .= "\t" . $tmpres;
			}
		}
	}
}

my $keycnt = 0;
open (OUTFI, ">$fin.jointHBV_HG");
foreach my $chrkey (sort keys %annovar_read_map){
	print OUTFI "$chrkey\t".$annovar_read_map{$chrkey}."\n";
	++$keycnt;
}
close(OUTFI);

print "Records print:$keycnt\n";
