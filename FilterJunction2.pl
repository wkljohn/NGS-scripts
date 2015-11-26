#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <input vcf file> -r <standard exonic regions> -d <min overlap>

INFO

my ($refName, $bedfile, $delfiltlen);
GetOptions(
"i=s"=>\$refName,
"r=s"=>\$bedfile,
"d=i"=>\$delfiltlen,
);

die "$usage" if(!$refName || !$bedfile); 


my %harshrangeUP;
my %harshrangeLOW;
my %harshrangeNAME;
my @refRange;
my $refindex=0;

#functions
sub IsInAllRange{
	my ($numcheck, $chrsearch) = @_;
	for (my $i=0; $i < $#{$harshrangeUP{$chrsearch}}; ++$i){
		#print "looping:" . $i . " of " . $#{$harshrangeUP{$chrsearch}} . "\n";
		if (@{$harshrangeUP{$chrsearch}}[$i] <= $numcheck && $numcheck <= $harshrangeLOW{$chrsearch}[$i]){
			#print "inrange!\n";
			return $i;
		}
	}
	return 0;
}


open (MYFILE, "<$bedfile");
print "Loading junction\n";
<MYFILE>; #ignores first line
while (<MYFILE>) {
    if (!($_ =~ m/^\#/)){
        my @splitarr = split /\t/, $_;
        chomp($splitarr[0]);
        if (!exists $harshrangeUP{$splitarr[0]}){
            $harshrangeUP{$splitarr[0]}[0] = $splitarr[1] + $delfiltlen;
            $harshrangeLOW{$splitarr[0]}[0] = $splitarr[2] - $delfiltlen;
            $harshrangeNAME{$splitarr[0]}[0] = $splitarr[7];
            
            print "creating harsh for '" . $splitarr[0] . "',LEN:" . length($splitarr[0]) . "\n";
        }else{
            
            #print $splitarr[2]  . " " .  $splitarr[3]  . " " .  $splitarr[7]  . "\n";
        	  push @{$harshrangeUP{$splitarr[0]}}, $splitarr[1] + $delfiltlen;
            push @{$harshrangeLOW{$splitarr[0]}},  $splitarr[2] - $delfiltlen;
            push @{$harshrangeNAME{$splitarr[0]}}, $splitarr[7];
        }
    #$refRange[$refindex][0] = $splitarr[1];
    #$refRange[$refindex][1] = $splitarr[2];
    #    print $refRange[$refindex][0] . " " . $refRange[$refindex][1] . "\n";
    }
}
close(MYFILE);

for my $key ( keys %harshrangeUP ) {
        print "$key\n";
}

open (REFFI, "<$refName");
open (OUTFI, ">$refName.junc_filter2");
open (OUTFIO, ">$refName.junc_filtout2");
while (<REFFI>) {
    if ($_ =~ m/^\#/){
        next;
        #last;
    }
    my $ln = $_;
    my @splitarr = split /\t/, $_;
    chomp($splitarr[0]);
    
    if (exists $harshrangeUP{$splitarr[0]}){
        my $startpos = int($splitarr[1]);
        my $endpos = $startpos + length($splitarr[3]);
        
        my $inrange1 = IsInAllRange($startpos, $splitarr[0]);
        my $inrange2 = IsInAllRange($endpos, $splitarr[0]);
        
        if ($inrange1){
            chomp();
            print OUTFI $_ . "\t" . $harshrangeNAME{$splitarr[0]}[$inrange1] . "\n";
        }elsif ($inrange2){
            #print "inrange\n";
            chomp();
            print OUTFI $_ . "\t" . $harshrangeNAME{$splitarr[0]}[$inrange2] . "\n";
            
        }else{
            print OUTFIO $_;
        }
    }else{
        print "ref seq:'". $splitarr[0] . "',LEN:" . length($splitarr[0]) . " doesn't exist\n";
        print OUTFIO $_;
    }
    
}
close(REFFI);
close(OUTFI);
close(OUTFIO);

