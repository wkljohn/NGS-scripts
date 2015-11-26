#!/usr/bin/perl
use strict;
use POSIX;
use XML::Parser;
use XML::XPath;
use XML::XPath::XMLParser;
use Getopt::Long;
use LWP::Simple;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
 
my $usage=<<INFO;

Usage:
perl $0 -r <name cross-reference> -q <(optional) query list> -o <output> -d <download file, 0/1>

INFO

#-i <uniprot tab based database file>

my ($refxmlName, $crossref, $output, $queryfile,$downdb);
GetOptions(
"r=s"=>\$crossref,
"q=s"=>\$queryfile,
"o=s"=>\$output,
"d=i"=>\$downdb,
);

die "$usage" if(!$crossref); 

#proc vars
my %querylist;
my %translateharsh;
my %translateharshrev;
my $checkingID;
#my $checkingID_inquery;
my @checkingIDs;
my $checkingID_reviewed;
my $checkingID_inquery;
my $isnonquery=0;
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
my $dbname = "uniprot_sprot.dat";
$refxmlName = $dbname;

sub getdbupdate{
	#getstore($addr, "uniprot_sprot.dat.gz");
	#system("gunzip uniprot_sprot.dat.gz");
	
	gunzip "uniprot_sprot.dat.gz" => "uniprot_sprot.dat"  or die "gunzip database failed: $GunzipError\n";
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
	#print $devices $translateharshrev{$checkingID} . "\t" . $FTName . "\t$checkingID\t$FTtype\t$FTdesc\t$Status\t$FTstart\t$FTend\n";
	print $devices $FTName . "\t$checkingID\t$FTtype\t$FTdesc\t$Status\t$FTstart\t$FTend\n";
		
}

sub finalizeprot{
	my $devices = shift;
	#my $reffeaturelist = shift;
	#my $refprotproplist = shift;
	if (length($checkingID) > 0){
		if ($FTpending == 1){
			writefeature($devices);
		}
	}
	
	$FTName = "";
}

my $cntduplications=0;
my $cntloaded=0;
my $cntproc=0;

if ($downdb){
	getdbupdate();
}

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
my @domain_acc_arr = ("REGION", "ACT_SITE", "CA_BIND", "CARBOHYD", "CHAIN", "COILED", "COMPBIAS", "CONFLICT", "CROSSLNK", "DISULFID", "DNA_BIND", "DOMAIN", "HELIX", "INTRAMEM", "INIT_MET", "MOD_RES", "LIPID", "METAL", "MOTIF", "MUTAGEN", "NON_CONS", "NON_STD", "NON_TER", "NP_BIND", "PEPTIDE", "PROPEP", "REPEAT", "SIGNAL", "SITE", "STRAND", "TRANSIT", "TOPO_DOM", "TRANSMEM", "TURN", "UNSURE", "VARIANT", "VAR_SEQ", "ZN_FING", "BINDING");
my @domain_acc_tra = ("region of interest", "active site", "calcium-binding region", "glycosylation site", "chain", "coiled-coil region", "compositionally biased region", "sequence conflict", "cross-link", "disulfide bond", "DNA-binding region", "domain", "helix", "intramembrane region", "initiator methionine", "modified residue", "lipid moiety-binding region", "metal ion-binding site", "short sequence motif", "mutagenesis site", "non-consecutive residues", "non-standard amino acid", "non-terminal residue", "nucleotide phosphate-binding region", "peptide", "propeptide", "repeat", "signal peptide", "site", "strand", "transit peptide", "topological domain", "transmembrane region", "turn", "unsure residue", "sequence variant", "splice variant", "zinc finger region", "binding site");

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

#see if its query based
if ($querycnt == 0){
	$isnonquery = 1;
	print "Non query mode enabled\n";
}

print "query stats: cnt:$querycnt, not_found: $unfoundcnt \n";

if (!(-s $refxmlName))
{ 
	die "reference database not found\n";
}

open (PROCESSEDFI, ">$output.proclist");
open (OUTFI, ">$output.results");
open (DBFI, "<$refxmlName");
#add header
#print PROCESSEDFI "Gene_id\tGene_name\tUniprot_id\tType\tDescription\tStatus\tStart\tEnd\n";
print PROCESSEDFI "Gene_name\tUniprot_id\tType\tDescription\tStatus\tStart\tEnd\n";
#real process
my @splitarr;
while (<DBFI>) {
	my $tag = substr($_, 0, 2);
	
	if ($tag eq "ID"){
		chomp();
		#finialize last ID here   
		#@splitarr = split("\t", $_);
		#print "stat:" . substr($_,29,3) . "\n";
		finalizeprot(*OUTFI);
		if (substr($_,29,3) eq "Rev"){
			$checkingID_reviewed = 1;
			$checkingID_inquery = 0;
			#print "rev\n";
		}else{
			$checkingID_reviewed = 0;
		}
	}elsif ($tag eq "AC"){
		chomp();
		my $cutpropstr = substr($_, 5);
		@splitarr = split(";", $cutpropstr);
		
		for (my $i = 0; $i <= $#splitarr; ++$i){
			if ($checkingID_inquery == 1){
				last;
			}
			
			#trim($splitarr[$i]);
			$checkingID = trim($splitarr[$i]);
			
			#print "checking:$checkingID\n";
			if (($isnonquery == 1 && checkprotinlist(\%translateharshrev, $checkingID) == 1) || checkprotinlist(\%querylist, $checkingID) == 1){
	
				if ($checkingID_reviewed == 0){
					$checkingID_inquery = 0;
					$non_rev_cnt++;
					#print "procnv:$checkingID\n";
				}else{
					$checkingID_inquery = 1;
					$proccnt++;
				}
				
				++$cntproc;
				#print "proc:$checkingID\n";
				if ($cntproc % 10 == 0){
					print "proc:$cntproc\n";
				}
				print PROCESSEDFI $querylist{$checkingID}[0] . "\t$checkingID\n";
				last;
			}else{
				$checkingID_inquery = 0;
				$skippedcnt++;
				#print "skip:$checkingID\n";
			}
		}
		
	}elsif ($checkingID_inquery == 1 && $tag eq "GN"){
		chomp();
		my $cutpropstr = substr($_, 5);
		#print "GN field:" . $_ . "; ID=$checkingID\n";
		if ($cutpropstr =~ /^Name=/){
			my @fieldarr = split(";", $cutpropstr);
			$FTName = substr($fieldarr[0], 5);
			#print "name=$FTName\n";
		}
	}elsif ((($isnonquery == 1 && $checkingID_inquery && length($FTName) > 0 ) || ( $isnonquery == 0 && $checkingID_inquery == 1)) && $tag eq "FT"){
		chomp();
		my $featuretype = trim(substr($_, 5,10));
		
		if (ismydomaintype(\@domain_acc_arr, $featuretype) == 1){
			if ($FTpending == 1){
				writefeature(*OUTFI);
			}
			#my $cutpropstr = trim(substr($_, 5));
			my $start = trim(substr($_, 15,5));
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
			
			if ($isnonquery == 0){
				$FTName = $querylist{$checkingID}[0];
			}
		}elsif(length($featuretype) == 0){
			if (length($_) > 33){
				my $desc;
				$desc = trim(substr($_, 33));
				$FTdesc.= "; $desc";
			}
		}
	#}else{
	#	print "Tag:$tag\n";
	}
}

print "skip: $skippedcnt, proc: $proccnt, unreli: $non_rev_cnt \n";
close(DBFI);
#output fi
close(OUTFI);
close(NOTFOUNDFI);
close(PROCESSEDFI);

