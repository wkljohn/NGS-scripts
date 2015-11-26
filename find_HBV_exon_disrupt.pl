#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -r <ref gtf> -i <bam file> -o <output file>  -v <viral align> -l <full read length>

INFO

my ($refgtf, $bamfile, $outfile, $readlen, $viralalign);

GetOptions(
"r=s"=>\$refgtf,
"i=s"=>\$bamfile,
"o=s"=>\$outfile,
"v=s"=>\$viralalign,
"l=i"=>\$readlen,
);

die "$usage" if(!$refgtf || !$bamfile || !$outfile || !$viralalign ); 

my %exoninfostart;
my %exoninfoend;
my %exoninfogenename;
my %viraligninfo;

#only for <100bp reads
sub readlen{
	my $readlen = shift;
	my $cigar = shift;
	
	my $minclippedlen=25;
	
	if ($cigar =~ /^[2-9][0-9]S/){
		my $pos = substr($cigar, 0, 2);	#must be 20-99
		if ($pos < $minclippedlen){goto full;}
		
		return ($pos,'-');
	}elsif($cigar =~ /[2-9][0-9]S$/){
		my $cilen = length($cigar);
		#my $pos = substr($cigar, 0, 2);	#must be 20-99
		my $posx = substr($cigar, $cilen-3, 2);	#must be 20-99;
		
		return  ($posx, '+');#($readlen - $pos, '+');
	}else{
	full:
		return -1;
	}
}


sub processcigarlen{
	my $cigarstr = shift;
	my @chararr = split '', $cigarstr;
	my $oplen = 0;
	my $totallen = 0;
	
	for (my $i == 0; $i <= $#chararr; ++$i){
		#print  int($chararr[$i]);
		if ($chararr[$i] eq 'M'){
			$totallen += $oplen;
			$oplen = 0
		}
		elsif ($chararr[$i] eq 'I'){
			$oplen = 0
		}
		elsif ($chararr[$i] eq 'D') 
		{
			$totallen += $oplen;
			$oplen = 0
		}
		elsif ($chararr[$i] eq 'N')
		{
			$totallen += $oplen;
			$oplen = 0
		}
		elsif ($chararr[$i] eq 'S'){
			$oplen = 0
		}
		elsif ($chararr[$i] eq 'H'){
			$oplen = 0
		}
		elsif ($chararr[$i] eq 'P'){
			$oplen = 0
		}
		else
		{
			#print "pre: $oplen @" . int($chararr[$i]) . "\n";
			$oplen = $oplen * 10 + int($chararr[$i]);
			#print "post: $oplen\n";
		}
	}
	if ($oplen > 0){
		$totallen += $oplen;
	}
	
	#print "total: $totallen (from $cigarstr) \n";
	return $totallen;
}

sub isatsplice{
	my $chr = shift;
	my $pos = shift;
	my $datharshstart = shift;
	my $datharshend = shift;
	
	if (exists $$datharshstart{$chr}){
		for (my $i == 0; $i < $#{$$datharshstart{$chr}}; ++$i){
			if (abs($$datharshstart{$chr}[$i] - $pos) < 4){
				return abs($$datharshstart{$chr}[$i] - $pos);
			}
			if (abs($$datharshend{$chr}[$i] - $pos) < 4){
				return abs($$datharshend{$chr}[$i] - $pos);
			}
		}
		return -1;
	}else{
		print "chr: $chr not found!\n";
		return -1;
	}
	
}

print "processing: $bamfile by ref $refgtf\n";

if (!(-e $bamfile)){
	die " read file does not exist!\n";
}


open (READFI, "<$refgtf");
while (<READFI>){
	chomp();
	my @splitarr = split "\t", $_;
	if ($splitarr[2] eq "exon"){
		my @splitinfo = split ";",  $splitarr[8];
		my $genename = $splitinfo[0];
		$genename =~ s/gene_id//g;
		
		#if ($splitarr[0] eq "chr4" && int($splitarr[3]) > 120000000 && int($splitarr[3]) < 130000000){
		#	print  int($splitarr[3]) . "\t" .  int($splitarr[4]) . "\n";
		#}
		push @{$exoninfostart{$splitarr[0]}}, int($splitarr[3]);
		push @{$exoninfoend{$splitarr[0]}}, int($splitarr[4]);
		push @{$exoninfogenename{$splitarr[0]}}, $genename;
	}
}
close(READFI);


my @vircontent = `samtools view $viralalign`;
my $vircnt = 0;
foreach my $line (@vircontent){
	if ($line !~ "^@"){
		my @splitarr = split "\t", $line;
		
		if ($splitarr[2] ne "*"){ 
			push @{$viraligninfo{$splitarr[0]}}, $splitarr[5] . "\t" .  $splitarr[3] . "\t" .  length($splitarr[9]);
			++$vircnt;
		}
		
	}
}
print "viral reads loaded: $vircnt\n";

my @bamcontent = `samtools view $bamfile`;
print "HG loaded: ". $#bamcontent . "\n";
open (OUTFI, ">$outfile");
foreach my $line (@bamcontent){
	if ($line !~ "^@"){
		my @splitarr = split "\t", $line;
		my $custreadlen = $readlen;
		
		#if (exists ($viraligninfo{$splitarr[0]})){
		#	my @splitarrx = split "\t",$viraligninfo{$splitarr[0]}[0];
	#		$custreadlen = int($splitarrx[2]);
		#}
		
		if (length ($splitarr[9]) < $custreadlen - 1 && $splitarr[2] ne "*"){ 
			#print "searching " . $splitarr[2] . "-" . $splitarr[3] . "\n";
			my $cireadlen = processcigarlen($splitarr[5]);
			my $atsplice = isatsplice($splitarr[2], $splitarr[3], \%exoninfostart, \%exoninfoend);
			my $atsplice2 = isatsplice($splitarr[2], int($splitarr[3]) + $cireadlen, \%exoninfostart, \%exoninfoend);
			my $vircigar  = "N/A";
			my $vircigarbu  = "N/A";
			
			if (exists ($viraligninfo{$splitarr[0]})){
				$vircigar = "";
				foreach my $cent (@{$viraligninfo{$splitarr[0]}}){
					my @splitarrx = split "\t", $cent;
					my @cutlen = readlen($custreadlen, $splitarrx[0]);
					#print $splitarrx[0].":".$cutlen[0] . ":".$cutlen[1] . ":". length($splitarr[9]) . ",";
					if (abs($cutlen[0] - length($splitarr[9])) < 2){
						if ($cutlen[1] eq '-'){
							$vircigar .= "\t" . $splitarrx[0]."\t".$splitarrx[1];
						}else{
							$vircigar .= "\t" . $splitarrx[0]."\t".(int($splitarrx[1]) +  ($custreadlen -$cutlen[0] )) . "|" . int($splitarrx[1]);# . "|" . $readlen . "-" . int($cutlen[0]) ;
						}
					}
					$vircigarbu .= "\t$cent:".$cutlen[0].$cutlen[1];
				}
				#print  "\n";
				
				if (length($vircigar) == 0){
					$vircigar = "B\t".$vircigarbu ;
				}
			}
			
			if ($atsplice > -1 ){
				print OUTFI $splitarr[0] . "\t" . $atsplice . "\t" . $splitarr[2] . "-" . $splitarr[3] . "\t+\t" . length ($splitarr[9]) . "\t" . $splitarr[5] . "\t" . $vircigar . "\n";
				#print "Found on " . $splitarr[2] . "-" . $splitarr[3] . "\t+" . $splitarr[2] . "\t" . $vircigar . "\n";
			}elsif ($atsplice2 > -1){
				print OUTFI $splitarr[0] . "\t" . $atsplice2 . "\t" . $splitarr[2] . "-" . $splitarr[3] . "\t+\t" .  length ($splitarr[9]) . "\t" . $splitarr[5] . "\t" . $vircigar . "\n";
				#print "Found on " . $splitarr[2] . ":" . $splitarr[3] . "\t+" . $splitarr[2] . "\t" . $vircigar . "\n";
			}
		}
		
	}
}
close (OUTFI);


