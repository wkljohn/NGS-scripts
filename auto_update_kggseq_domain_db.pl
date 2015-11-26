#!/usr/bin/perl
use strict;
use POSIX;
use XML::Parser;
use XML::XPath;
use XML::XPath::XMLParser;
use LWP::Simple;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -i <uniprot tab based database file> -r <name cross-reference> -q <query list> -o <output>

INFO

my ($refxmlName, $crossref, $output, $queryfile);
GetOptions(
"i=s"=>\$refxmlName,
"r=s"=>\$crossref,
"q=s"=>\$queryfile,
"o=s"=>\$output,
);

die "$usage" if(!$refxmlName || !$crossref); 

#proc vars
my %querylist;
my $checkingID;
my $checkingID_reviewed;
my $checkingID_inquery;
my $querycnt=0;
my $unfoundcnt=0;
my $skippedcnt=0;
my $proccnt=0;
my $non_rev_cnt=0;

#feature handling
my $FTpending=0;
my $FTName;
my $FTtype;
my $FTdesc;
my $FTstart;
my $FTend;

my $addr = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";

sub getdbupdate{
	
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub checkprotinlist{
	my $refharsh = shift;
	my $protcheck = shift;
	my %harsh = %$refharsh;
	
	if (exists $harsh{$protcheck}){
		return 1;
	}else{
		return 0;
	}
}

sub ismydomaintype{
	my $refdomainnamelist = shift;
	my $domainnamequery = shift;
	my @domainaccarr = @$refdomainnamelist;
	
	for (my $i = 0; $i <= $#domainaccarr; ++$i){
		if ($domainaccarr[$i] eq $domainnamequery){
			return 1;
		}
	}
	
	return 0;
}

sub translatedomaintype{
	my $refdomainnamelist = shift;
	my $reftradomainnamelist = shift;
	my $domainnamequery = shift;
	my @domainaccarr = @$refdomainnamelist;
	my @domainacctra = @$reftradomainnamelist;
	
	for (my $i = 0; $i <= $#domainaccarr; ++$i){
		if ($domainaccarr[$i] eq $domainnamequery){
			return $domainacctra[$i];
		}
	}
	
	return 0;
}
sub writefeature{
	my $devices = shift;
	my $Status = "NA";
			
	if ($FTdesc =~ /(Potential)/){
		$Status = "potential";
		$FTdesc =~ s/\(Potential\)//;
	}
	
	$FTpending = 0;
	#print "FT: $checkingID,$FTtype,$FTstart,$FTend,$FTdesc \n";
	print $devices $querylist{$checkingID}[0] . "\t" . $querylist{$checkingID}[1] . "\t$checkingID\t$FTtype\t$FTdesc\t$Status\t$FTstart\t$FTend\n";
		
}

sub finalizeprot{
	my $devices = shift;
	$FTName = "";
	#my $reffeaturelist = shift;
	#my $refprotproplist = shift;
	if (length($checkingID) > 0){
		if ($FTpending == 1){
			writefeature($devices);
		}
	}
	
}

my %translateharsh;
my $cntduplications=0;
my $cntloaded=0;

open (HARSHFI, "<$crossref");
while (<HARSHFI>) {
    my @splitarr = split /\t/, $_;
	
	#if ($splitarr[1] eq "RefSeq" || $splitarr[1] eq "RefSeq_NT"){
	if ($splitarr[1] eq "RefSeq_NT"){
		my $tmpname = substr($splitarr[2], 0,  length($splitarr[2]) - 3);		
		if (exists $translateharsh{$tmpname}){
			print "harsh duplication:$splitarr[0]\t$tmpname\n";
			print "ORG:" . $translateharsh{$tmpname} . "\t" . $tmpname . "\n";
			++$cntduplications;
			next;
		}
		$translateharsh{$tmpname} = $splitarr[0];
		$translateharshrev{$splitarr[0]} = $tmpname;
		
		++$cntloaded;
	}
}
close(HARSHFI);

print "reference stats: dup: $cntduplications, cnt: $cntloaded \n";

#feature vars
my @domain_acc_arr = ("CA_BIND",				"COILED",				"COMPBIAS",							"DNA_BIND",				"DOMAIN","HELIX",	"INTRAMEM",				"NP_BIND",								"PEPTIDE","PROPEP",		"REGION",				"REPEAT","VARIANT",				"MOTIF",				"SITE","STRAND","TOPO_DOM",				"TRANSMEM",				"TURN","ZN_FING");
my @domain_acc_tra = ("calcium-binding region",	"coiled-coil region",	"compositionally biased region",	"DNA-binding region",	"domain","helix",	"intramembrane region",	"nucleotide phosphate-binding region",	"peptide","propeptide",	"region of interest",	"repeat","sequence variant",	"short sequence motif",	"site","strand","topological domain",	"transmembrane region",	"turn","zinc finger region");

#load the query list
open (QUERYFI, "<$queryfile");
open (NOTFOUNDFI, ">$output.notfound");
while (<QUERYFI>) {
	chomp();
	#print "line\n";
    my @splitarr = split /\t/, $_;
	if ($#splitarr >= 1){
		#print "checking:$splitarr[0]\n";
		if (exists $translateharsh{$splitarr[0]}){
			#$querylist{$translateharsh{$splitarr[0]}} = $splitarr[0];
			push @{$querylist{$translateharsh{$splitarr[0]}}}, $splitarr[0];
			push @{$querylist{$translateharsh{$splitarr[0]}}}, $splitarr[1];
			++$querycnt;
		}else{
			print NOTFOUNDFI $_ . "\n";
			++$unfoundcnt;
		}
	}
}
close(QUERYFI);

print "query stats: cnt:$querycnt, not_found: $unfoundcnt \n";

open (PROCESSEDFI, ">$output.proclist");
open (OUTFI, ">$output.results");
open (DBFI, "<$refxmlName");
#add header
print PROCESSEDFI "Gene_id\tGene_name\tUniprot_id\tType\tDescription\tStatus\tStart\tEnd\n";
#real process
while (<DBFI>) {
	chomp();
	my @splitarr;
	
	if ($_ =~ /^ID/){
		#finialize last ID here   
		#@splitarr = split("\t", $_);
		#print "stat:" . substr($_,29,3) . "\n";
		finalizeprot(*OUTFI);
		if (substr($_,29,3) eq "Rev"){
			$checkingID_reviewed = 1;
			#print "rev\n";
		}else{
			$checkingID_reviewed = 0;
		}
	}elsif ($_ =~ /^AC/){
		my $cutpropstr = substr($_, 5);
		@splitarr = split(";", $cutpropstr);
		#$splitarr[0] = substr($splitarr[0], 0, length($splitarr[0]) - 1);
		$checkingID = $splitarr[0];
		
		#print "checking:$checkingID\n";
		if (checkprotinlist(\%querylist, $checkingID) == 1){

			if ($checkingID_reviewed == 0){
				$checkingID_inquery = 0;
				$non_rev_cnt++;
				print "procnv:$checkingID\n";
			}else{
				$checkingID_inquery = 1;
				$proccnt++;
			}
			
			print "proc:$checkingID\n";
			print PROCESSEDFI $querylist{$checkingID}[0] . "\t$checkingID\n"
		}else{
			$checkingID_inquery = 0;
			$skippedcnt++;
			#print "skip:$checkingID\n";
		}
	
	}elsif ($checkingID_inquery == 1 && $_ =~ /^GN/){
		my $cutpropstr = substr($_, 5);
		if ($cutpropstr =~ /^Name=/){
			my @fieldarr = split(";", $cutpropstr);
			$FTName = substr($fieldarr[0], 5);
		}
	}elsif ($checkingID_inquery == 1 && $_ =~ /^FT/){
		my $featuretype = trim(substr($_, 5,10));
		
		if (ismydomaintype(\@domain_acc_arr, $featuretype) == 1){
			if ($FTpending == 1){
				writefeature(*OUTFI);
			}
			#my $cutpropstr = trim(substr($_, 5));
			my $start = trim(substr($_, 15,4));
			my $end = trim(substr($_, 21,6));
			my $desc;
			
			if (length($_) > 33){
				$desc = trim(substr($_, 33));
			}
			
			$FTtype=translatedomaintype(\@domain_acc_arr,\@domain_acc_tra,$featuretype);
			$FTdesc=$desc;
			$FTstart=$start;
			$FTend=$end;
			$FTpending = 1;
		}elsif(length($featuretype) == 0){
			if (length($_) > 33){
				my $desc;
				$desc = trim(substr($_, 33));
				$FTdesc.= "; $desc";
			}
		}
	}
}

print "skip: $skippedcnt, proc: $proccnt, unreli: $non_rev_cnt \n";
close(DBFI);
#output fi
close(OUTFI);
close(NOTFOUNDFI);
close(PROCESSEDFI);

