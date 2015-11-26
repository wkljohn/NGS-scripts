#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <input vcf file> -r <rum junction bed file> -d <deletion min length to filter>

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
my @refRange;
my $refindex=0;

#functions
sub IsInAllRange{
	my ($numcheck, $chrsearch) = @_;
	for (my $i=0; $i < $#{$harshrangeUP{$chrsearch}}; ++$i){
		#print "looping:" . $i . " of " . $#{$harshrangeUP{$chrsearch}} . "\n";
		if (@{$harshrangeUP{$chrsearch}}[$i] <= $numcheck && $numcheck <= $harshrangeLOW{$chrsearch}[$i]){
			print "inrange!\n";
			return 1;
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
            $harshrangeUP{$splitarr[0]}[0] = $splitarr[1];
            $harshrangeLOW{$splitarr[0]}[0] = $splitarr[2];
            print "creating harsh for '" . $splitarr[0] . "',LEN:" . length($splitarr[0]) . "\n";
        }else{
            push @{$harshrangeUP{$splitarr[0]}}, $splitarr[1];
            push @{$harshrangeLOW{$splitarr[0]}},  $splitarr[2];
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
open (OUTFI, ">$refName.junc_filter");
open (OUTFIO, ">$refName.junc_filtout");
while (<REFFI>) {
    if ($_ =~ m/^\#/){
        next;
        #last;
    }
    my $ln = $_;
    my @splitarr = split /\t/, $_;
    chomp($splitarr[0]);
    
    if (length($splitarr[3]) >= $delfiltlen && exists $harshrangeUP{$splitarr[0]}){
        my $startpos = int($splitarr[1]);
        my $endpos = $startpos + length($splitarr[3]);
        
        if (IsInAllRange($startpos, $splitarr[0]) || IsInAllRange($endpos, $splitarr[0])){
            print OUTFIO $_;
        }else{
            print OUTFI $_;
        }
    }else{
        print "ref seq:'". $splitarr[0] . "',LEN:" . length($splitarr[0]) . " doesn't exist\n";
        print OUTFI $_;
    }
    
}
close(REFFI);
close(OUTFI);
close(OUTFIO);

