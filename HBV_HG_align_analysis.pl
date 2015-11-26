#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -v <viral alignment bam> -h <human alignback bam> -o <output file> -r <max redundant filter, default=1> -n <N count thres>

INFO

my ($bamviral, $bamhg, $output, $rdf, $seqnthres);
$rdf = 1;
GetOptions(
"v=s"=>\$bamviral,
"h=s"=>\$bamhg,
"o=s"=>\$output,
"r=i"=>\$rdf,
"n=i"=>\$seqnthres,
);

die "$usage" if(!$output || !$bamhg); 

my $badread = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
my $badqc =   "#####################################################################################################";


sub countn{
	my $seqstr = shift;
	my @seqarr = split '', $seqstr;
	my $ccnt = 0;
	
	#print "ncnt of: $seqstr \n";
	for (my $i = 0; $i < $#seqarr; ++$i){
		if ($seqarr[$i] eq "N"){
			++$ccnt;
		}
	}
	
	return $ccnt;
}

sub qctop{
	my $qccharsstr = shift;
	my @qcchars = split '', $qccharsstr;
	my $tmpsump = 0.0;
	my $tmpint;
	
	for (my $i = 0; $i <= $#qcchars; ++$i){
		$tmpint = ord($qcchars[$i]) - 33;
		$tmpsump += 10 ** ($tmpint / -10);
	}
	my $tmpval = int($tmpsump * 1000);
	
	return $tmpval / 1000;
}

sub revdnacomp {
  my $dna = shift; # or   my $dna = shift @_;
  # ah, scalar context of scalar gives expected results.
  # my ($dna) = @_; # would work, too

  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

sub processread{
	my $splitseqref = shift;
	my $cigar = shift;
	my $flags = shift;
	my $pair = shift;
	my @splitseq = @$splitseqref;
	
	my $fastqstr="";
	my $minclippedlen=25;
	my $read=$splitseq[9];
	#print $cigar."\n";
	my $readqc=$splitseq[10];
	my $readorent;
	my $subtag;
	my $readid="@" . $splitseq[0];
	my $integrationorient;
	my $pos = $splitseq[3];
	if ($cigar & 0x10){
		$subtag="$pair:Y:0:GTGAAA";
		$readorent = "+";
	}else{
		$subtag="$pair:N:0:GTGAAA";
		$readorent = "+";
	}
	$readid.=" $subtag";
	
	
	if ($cigar =~ /^([5-9]|[1-9][0-9]+)S[0-9]+M$/ 
	   || $cigar =~ /^([5-9]|[1-9][0-9]+)S[1-9]I[0-9]+M$/  || $cigar =~ /^([5-9]|[1-9][0-9]+)S[0-9]+M[1-9]I[0-9]+M$/ 
		|| $cigar =~ /^([5-9]|[1-9][0-9]+)S[0-9]+M[1-8]S$/){
		
		my $posadd = substr($cigar, 0, index($cigar, 'S'));
		if ($flags & 0x10){	#reversed
			$integrationorient = "-";#\t$pos\t$cigar";
		}else{
			$integrationorient = "-";#\t$pos\t$cigar";
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
		
		my $pos_b = $pos + $posadd;
		if ($flags & 0x10){	#reversed
			$integrationorient = "+";#\t$pos_b\t$cigar";
		}else{
			$integrationorient = "+";#\t$pos_b\t$cigar";
		}
		#$pos_b = $pos + $posadd;
		#print "PA:$posadd PB:$pos_b\n";
		
		#return "+\t$pos_b\t$cigar";
	}elsif($cigar =~ /^[0-9]+M$/ || $cigar =~ /^[0-9]+M[1-4]S$/ || $cigar =~ /^[1-4]S[0-9]+M$/  || $cigar =~ /^[1-4]S[0-9]+M[1-4]S$/ 
			|| $cigar =~ /^[0-9]+M[0-9]?[0-9]D[0-9]+M$/  || $cigar =~ /^[0-9]+M[1-9]I[0-9]+M$/){
		if ($flags & 0x40){	#1st read
			if ($flags & 0x10){	#reversed
				if ($flags & 0x20){	#reversed
					$integrationorient = "<";#"+>\t".$pos."\t$cigar";
				}else{				#anti strand pair
					$integrationorient = "<";#"-<\t".$pos."\t$cigar";
				}
			}else{				# not reversed
				if ($flags & 0x20){	#reversed		seems fixed by 308
					$integrationorient = ">";#"+>\t".($pos+101)."\t$cigar";
				}else{				#anti strand pair
					$integrationorient = "?";#"?-+"; #"+>\t".($pos+101)."\t$cigar";
				}
			}
		}else{	#2nd read
			if ($flags & 0x10){	#reversed
				if ($flags & 0x20){	#reversed
					$integrationorient = ">";#"+>\t".($pos+101)."\t$cigar";
				}else{				#anti strand pair
					$integrationorient = "<";#"-<\t".($pos+101)."\t$cigar";
				}
			}else{				# not reversed
				if ($flags & 0x20){	#reversed		seems fixed by 308
					$integrationorient = ">";#"+>\t".$pos."\t$cigar";
				}else{				#anti strand pair
					$integrationorient = "?";#"?-+"; #"+>\t".($pos+101)."\t$cigar";
				}
			}
		}
		
	}else{
full:
		#print "Unhdl cigar:$cigar\n";
		$integrationorient = "NA";
	}
	print "int ori det:" . $integrationorient . "\n";
	
	
	my @returnval=();
	
	if ($cigar =~ /^[2-9][0-9]S/){
		my $pos = substr($cigar, 0, 2);	#must be 20-99
		if ($pos < $minclippedlen){goto full;}
		my $readf = substr($read, 0, $pos) ;#. substr($badread, $pos);
		my $readqcf = substr($readqc, 0, $pos);# . substr($badqc, $pos);
		
		#print $cigar."|$pos S\n";
		#print "$readid\n$readf\n$readorent\n$readqcf\n";
		$returnval[0] = $readid;
		$returnval[1] = $readf;
		$returnval[2] = $readorent;
		$returnval[3] = $readqcf;
		
	}elsif($cigar =~ /[2-9][0-9]S$/){
		my $cilen = length($cigar);
		my $rlen = length($read);
		my $posx = substr($cigar, $cilen-3, 2);	#must be 20-99
		if ($posx < $minclippedlen){goto full;}
		my $pos = $rlen - $posx;
		my $readf = substr($read, $posx);
		my $readqcf = substr($readqc, $posx);
		
		#print $cigar."|$pos E\n";
		#print "$readid\n$readf\n$readorent\n$readqcf\n";
		$returnval[0] = $readid;
		$returnval[1] = $readf;
		$returnval[2] = $readorent;
		$returnval[3] = $readqcf;
	}else{
full:
		if ($flags & 0x4){	#unmapped
			#even this is unmapped, data required for finding another end
			$returnval[0] = $readid;
			$returnval[1] = $read;
			$returnval[2] = $readorent;
			$returnval[3] = $readqc;
		}else{
			$returnval[0] = $readid;
			$returnval[1] = $read;
			$returnval[2] = $readorent;
			$returnval[3] = $readqc;
		}
	}
	$returnval[4] = $integrationorient;
	$returnval[5] = $cigar;
	
	return @returnval;
	#pair info
	#if ($flags & 0x40){ #1st in pair
	#	$firstread{$splitseq[0]} = processread(\@splitseq, $cigar, $flags, 1);
	#}else{
	#	$secondread{$splitseq[0]} = processread(\@splitseq, $cigar, $flags, 2);
	#}
}

my %firstread;
my %secondread;
my %allreads; # = ();
my @viralsam = `samtools view $bamviral`;
my %readidcnt;
my $refreadlen = -1;

print "Viral bam: " . $bamviral . " Size:" . $#viralsam . "\n";
#while (<STDIN>){
foreach (@viralsam){
	if ($_ !~ /^@/){
		chomp();
		my @splitseq = split("\t", $_);
		if ($refreadlen < 0){
			$refreadlen = length($splitseq[9]);
		}
		
		my $readid = $splitseq[0];
		my $flags = int($splitseq[1]);
		my $cigar = $splitseq[5];
		print $cigar ."\n"; 
		if (!($flags & 0x100)){	#primary alignment
			#count read id
			my $pair;
			if ($flags & 0x40){	#1st read
				$pair = 1;
			}else{
				$pair = 2;
			}
			
			#if (exists $readidcnt{$readid . "\\" . $pair}){
			#	$readidcnt{$readid . "\\" . $pair}++;
			#}else{
			#	$readidcnt{$readid . "\\" . $pair} = 1;
			#}
			
			my @tmp=processread(\@splitseq, $cigar, $flags, 1);
			if ($tmp[0] ne ""){
				#print "tmp: " . $tmp[0] . "\t" . $tmp[1] . "\n";
				print "id: " . $readid . "\\" . $pair . "\n";
				push @{$allreads{$readid . "\\" . $pair}} , @tmp;
				#print "all: " . @{$allreads{$readid . "\\" . $pair}}[0] . "\t" . @{$allreads{$readid . "\\" . $pair}}[1] . "\n";
			}else{
				print "skipped: $readid\\$pair\n";
			}
		#}else{
			#print "2nd\n";
		}
	}
}

my @hgalign = `samtools view $bamhg | cut -f 1-6,10,11`;
my $nexclcnt = 0;
print "HG bam: " . $bamhg . " Size:" . $#hgalign . "\n";

#cnt redundant
foreach (@hgalign){
	if ($_ !~ /^@/){
		chomp();
		my @splitseq = split("\t", $_);
		my $readid = $splitseq[0];
		my $revcomQ = revdnacomp($splitseq[6]);
		my $rqc = qctop($splitseq[7]);
		my $pair;
		#print "_Q:" . $splitseq[6] . "\n";
		#print "RQ:" . $revcomQ . "\n";
		if (exists $allreads{$readid . "\\1"} && 
			($splitseq[6] eq ${$allreads{$readid . "\\1"}}[1] || $revcomQ eq ${$allreads{$readid . "\\1"}}[1])){
			$pair = 1;
		}elsif (exists $allreads{$readid . "\\2"} && 
			($splitseq[6] eq ${$allreads{$readid . "\\2"}}[1] || $revcomQ eq ${$allreads{$readid . "\\2"}}[1])){
			$pair = 2;
		}
		
		#add to harsh
		if (exists $readidcnt{$readid . "\\" . $pair}){
			$readidcnt{$readid . "\\" . $pair}++;
		}else{
			$readidcnt{$readid . "\\" . $pair} = 1;
		}
	}
}
 
print "ref read len:$refreadlen\n";

open (READWRFI1, ">$output.vcf");
foreach (@hgalign){
	#print "===========================\n";
	if ($_ !~ /^@/){
		chomp();
		my @splitseq = split("\t", $_);
		my $pair;
		my $readid = $splitseq[0];
		#process this read or not?
		if ($splitseq[5] ne "*"){
			print "procID:" . $readid .  "\n";
			#print "R1:" . $allreads{$readid . "\\1"}[0] . "\t" . $allreads{$readid . "\\1"}[4] . "\t" . $allreads{$readid . "\\1"}[5]  .  "\n";
			#print "R2:" . $allreads{$readid . "\\2"}[0] . "\t" . $allreads{$readid . "\\2"}[4] . "\t" . $allreads{$readid . "\\2"}[5]  . "\n";
			#my @tmpinfo1 = @{$allreads{$readid . "\\1"}};
			#my @tmpinfo2 = @{$allreads{$readid . "\\2"}};
			my @viralinfoarr;
			my @viralinfoarralt;
			
			#det which read it's on
			#print "R1R:" . $allreads{$readid . "\\1"}[1] . "\n";
			#print "R2R:" . $allreads{$readid . "\\2"}[1] . "\n";
			
			#redundant read?
			my $flags = int($splitseq[1]);
			#my $redunt = 0;
			#if ($flags & 0x100){	#not primary alignment
			#	$redunt = 1;
			#	print "2nd\n";
			#}
			
			my $revcomQ = revdnacomp($splitseq[6]);
			my $rqc = qctop($splitseq[7]);
			my $ncontent;
			
			#print "_Q:" . $splitseq[6] . "\n";
			#print "RQ:" . $revcomQ . "\n";
			if (exists $allreads{$readid . "\\1"} && 
				($splitseq[6] eq ${$allreads{$readid . "\\1"}}[1] || $revcomQ eq ${$allreads{$readid . "\\1"}}[1])){
				$pair = 1;
				
				#some exception case will change the side of the read
				if (${$allreads{$readid . "\\1"}}[5] =~ /[2-9][0-9]M/){
					@viralinfoarralt = @{$allreads{$readid . "\\1"}};
					@viralinfoarr = @{$allreads{$readid . "\\2"}};
				}else{
					@viralinfoarr = @{$allreads{$readid . "\\1"}};
					@viralinfoarralt = @{$allreads{$readid . "\\2"}};
				}
				
				print "orient:1, " . $splitseq[2] . "\n";
			}elsif (exists $allreads{$readid . "\\2"} && 
				($splitseq[6] eq ${$allreads{$readid . "\\2"}}[1] || $revcomQ eq ${$allreads{$readid . "\\2"}}[1])){
				$pair = 2;
				
				#some exception case will change the side of the read
				if (${$allreads{$readid . "\\2"}}[5] =~ /[2-9][0-9]M/){
					#flip to itself if has substancial match
					@viralinfoarralt = @{$allreads{$readid . "\\2"}};
					@viralinfoarr = @{$allreads{$readid . "\\1"}};
				}else{
					@viralinfoarr = @{$allreads{$readid . "\\2"}};
					@viralinfoarralt = @{$allreads{$readid . "\\1"}};
				}
				
				print "orient:2, " . $splitseq[2] . "\n";
			}else{
				print "EXCEPTION! On $readid\n";
				print $splitseq[6] . "\n";
				print ${$allreads{$readid . "\\1"}}[1] . "\n";
			}
			print $viralinfoarr[0] . "\t" . $viralinfoarr[1] . "\t" . $viralinfoarr[2] . "\n";
			
			#det redundancy
			my $redunt = $readidcnt{$readid . "\\" . $pair} - 1;
			my $readlen =  length($splitseq[6]);#length($viralinfoarr[1]);
			if ($readlen < $refreadlen){
				$ncontent = countn($viralinfoarr[1]);
			}else{
				$ncontent = countn($viralinfoarralt[1]);
			}
			
			#orientation fix"
			if ($rdf >= $redunt && $ncontent <= $seqnthres){
				if ($viralinfoarralt[4] eq "+" || $viralinfoarralt[4] eq ">"){
					my $writeref = substr($splitseq[6], $readlen - 1);
					my $writealt = "T";
					if ($writeref eq $writealt){
						$writealt = "A";
					}
					my $infostr = "LN=$readlen;RD=$redunt;CI=".$splitseq[5].";VCI=".$viralinfoarralt[5].";VCI2=".$viralinfoarr[5].";EB=$rqc;DR=+;RID=" . $readid;
					print READWRFI1	$splitseq[2] . "\t" .  (int($splitseq[3]) + $readlen - 1) . "\t.\t$writeref\t$writealt\t" . $splitseq[4] . "\tPASS\t" . $infostr  . "\n";
				}elsif ($viralinfoarralt[4] eq "-" || $viralinfoarralt[4] eq "<"){
					my $writeref = substr($splitseq[6], 0, 1);
					my $writealt = "T";
					if ($writeref eq $writealt){
						$writealt = "A";
					}
					my $infostr = "LN=$readlen;RD=$redunt;CI=".$splitseq[5].";VCI=".$viralinfoarralt[5].";VCI2=".$viralinfoarr[5].";EB=$rqc;DR=-;RID=" . $readid;
					print READWRFI1	$splitseq[2] . "\t" .  int($splitseq[3]) . "\t.\t$writeref\t$writealt\t" . $splitseq[4] . "\tPASS\t" . $infostr  . "\n";
				
				}else{
					my $writeref = substr($splitseq[6], $readlen - 1);
					my $writealt = "T";
					if ($writeref eq $writealt){
						$writealt = "A";
					}
					my $infostr = "LN=$readlen;RD=$redunt;CI=".$splitseq[5].";VCI=".$viralinfoarralt[5].";VCI2=".$viralinfoarr[5].";EB=$rqc;DR=N;RID=" . $readid;
					print READWRFI1	$splitseq[2] . "\t" .  (int($splitseq[3]) + $readlen - 1) . "\t.\t$writeref\t$writealt\t" . $splitseq[4] . "\tPASS\t" . $infostr  . "\n";
				
					print "ERR integration direction:" . $viralinfoarralt[4] . "\n";
				}
			}elsif($ncontent > $seqnthres){
				$nexclcnt++;
			}
		}
	}
}
close(READWRFI1);

print "excl N: $nexclcnt\n";
#open (READWRFI1, ">$output.fastq");
#foreach my $str (@allreads){
#		print READWRFI1 $str . "\n";
#}
#close(READWRFI1);

#open (READWRFI1, ">$output.R1.fastq");
#open (READWRFI2, ">$output.R2.fastq");
#foreach my $key (keys %firstread){
#	if (!(exists $secondread{$key})){
#		die "key:$1 does not exist!\n";
#	}else{
#		print READWRFI1 $firstread{$key} . "\n";
#		print READWRFI2 $secondread{$key} . "\n";
#	}
#}
#close(READWRFI1);
#close(READWRFI2);



