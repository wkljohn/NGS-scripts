#!/usr/bin/perl
use strict;
use POSIX;
use Getopt::Long;

my $usage=<<INFO;

Usage:
perl $0 -c <genelist path> -i <individual list split by comma> -s <file suffix> -r <Pathway Files path> -f <cut off threshold> -n <max number in one pathway>

INFO

my ($listpath, $indvlist, $filesuffix, $pathwaypath, $cutoffthres, $cutoffthresnum);
GetOptions(
"c=s"=>\$listpath,
"i=s"=>\$indvlist,
"s=s"=>\$filesuffix,
"r=s"=>\$pathwaypath,
"f=s"=>\$cutoffthres,
"n=s"=>\$cutoffthresnum,
);

die "$usage" if(!$listpath || !$indvlist); 


my %harshIndividualsGenes	;
my %harshIndividualsGenesEncoded	;
my %harshEncodings;
my %harshEncodingsRev;
my %harshPathScore;
my %harshPathList;
my %harshPathListparti;
my %AllInteractions;
my @refRange;
my @resultarray;
my $refindex=0;

my @snmarr = ("HKP1-133", "HKP1-151", "HKP1-228", "HKP1-251");
my @lnmarr = ("HKP1-74", "HKP1-315", "HKP1-234", "HKP2-107");
my @smarr = ("HKP2-308", "HKP2-84", "HKP2-383", "HKP2-254");
my @lmarr = ("HKP2-105", "HKP2-149", "HKP3-229", "HKP3-233");
my @ntarr = ("HKP3-151NT", "HKP3-233NT", "HKP3-234NT", "HKP3-308NT");
my %transharsh;
$transharsh{"HKP1-133"} = "SNM";
$transharsh{"HKP1-151"} = "SNM";
$transharsh{"HKP1-228"} = "SNM";
$transharsh{"HKP1-251"} = "SNM";

$transharsh{"HKP1-74"} = "LNM";
$transharsh{"HKP1-315"} = "LNM";
$transharsh{"HKP1-234"} = "LNM";
$transharsh{"HKP2-107"} = "LNM";

$transharsh{"HKP2-308"} = "SM";
$transharsh{"HKP2-84"} = "SM";
$transharsh{"HKP2-383"} = "SM";
$transharsh{"HKP2-254"} = "SM";

$transharsh{"HKP2-105"} = "LM";
$transharsh{"HKP2-149"} = "LM";
$transharsh{"HKP3-229"} = "LM";
$transharsh{"HKP3-233"} = "LM";

$transharsh{"HKP3-151NT"} = "NT";
$transharsh{"HKP3-233NT"} = "NT";
$transharsh{"HKP3-234NT"} = "NT";
$transharsh{"HKP3-308NT"} = "NT";


#functions
sub FindScore{
	my @listtocheck = @_;
	my $sumscore = 0;
	my $repeatscore = 0;
	my %allgenefound;
	my %indvaffected;
	my @mutposition = ();
	#my %foundFlagHarsh;
	
	foreach my $key (keys %harshIndividualsGenesEncoded){
		for (my $i = 0; $i <= $#{$harshIndividualsGenesEncoded{$key}}; ++$i){ 
			#check if the gene has been processed
			#find each gene in the list
			my $individualscore = 0;
			my %genecounted;
			#my $found = false;
			
			for (my $k = 0; $k <= $#listtocheck; ++$k){ 
				if ($harshIndividualsGenesEncoded{$key}[$i][0] eq $listtocheck[$k] && !exists $genecounted{$listtocheck[$k]}){
					#print "indv:" .  $harshIndividualsGenesEncoded{$key}[$i] . " mat " . $listtocheck[$k] . "\n";
                    push @mutposition, $harshIndividualsGenesEncoded{$key}[$i][1] . $harshIndividualsGenesEncoded{$key}[$i][2];
                    
					if ($individualscore == 0){
						$genecounted{$listtocheck[$k]} = 1;
						$allgenefound{$harshIndividualsGenes{$key}[$i]} = 1;
						$individualscore += 2;
						
					}else{
						$genecounted{$listtocheck[$k]} = 1;
						$allgenefound{$harshIndividualsGenes{$key}[$i]} = 1;
						$individualscore += 1;
					}
				}
				
				#if (exists $genecounted{$listtocheck[$k]}){
				#	print "key counted!\n";
				#}
				
			}
			
			#print "indv $key scores:" . $individualscore . "\n";
			if ($individualscore > 0){
				$indvaffected{$key} = 1;
			}
			$sumscore += $individualscore;
		}
	}
    
    #check the position of the mutation list to see if there are repeated occurances
    for (my $j = $#mutposition; $j > 0; --$j){
        for (my $k = $#mutposition - 1; $k >= 0 ; --$k){
            if ($mutposition[$j] eq $mutposition[$k]){
                ++$repeatscore;
                splice @mutposition, $k , 1;
                --$j;
            }
        }
    }
	
	return $sumscore, $repeatscore, \%allgenefound, \%indvaffected;
}


sub LoadInteractionList{
	my $interactarrref = shift;
	my $thres = shift;
	my $genecnt = 0;
	my $genecntthres = 0;
	
	<ININTERACT>; #first line uneless
	while (<ININTERACT>){
		$genecnt++;
		my @splitfields = split /\t/;
		#print $_ . "\n";
		#print $splitfields[0] . "\n";
		if (!exists $interactarrref->{$splitfields[0]}){
			$interactarrref->{$splitfields[0]} = ();
		}
		if ($splitfields[2] >= $thres){
		  $genecntthres++;
			push @{$interactarrref->{$splitfields[0]}}, [ @splitfields ];
			push @{$interactarrref->{$splitfields[1]}}, [ @splitfields ];
            
            #for (my $i = 0; $i <= $#splitfields; ++$i){ 
            #    print  $splitfields[$i] . " ";
            #}
            #print "\n";
 
		}
	}
	
    
    for my $key ( keys %{$interactarrref}){
        if ($#{$interactarrref->{$key}} > $cutoffthresnum){
            my @arr = @{$interactarrref->{$key}};
            #print "sorting $key\n";
            @{$interactarrref->{$key}} = sort { $b->[2] <=> $a->[2] } @arr;
            
            my $rmcnt = 0;
            my $orgcnt = $#arr;
            for (my $h = $cutoffthresnum; $h <= $#arr; ++$h){
                pop  @{$interactarrref->{$key}};
                #print "rm " . ${$interactarrref->{$key}}[$#{$interactarrref->{$key}} - 1][2] . " \n";
                ++$rmcnt;
            }
            #print "removed $rmcnt of \n";
        }
        push @{$interactarrref->{$key}}, [($key,$key,1)];
        
    }
    
    
	print "Loaded $genecnt genes with $genecntthres in threshold\n";
}


sub FindInInteractionList{
	my $genename = shift;
	#my $threshold = shift;
	
	my @geneinteractlist = ();
	my @geneallinteractlist = ();
	if (exists $AllInteractions{$genename}){
		print "$genename\n";
		if ($AllInteractions{$genename}){
			@geneallinteractlist = @{$AllInteractions{$genename}};
		}else{
			print "emt arr\n";
		}
	}else{
		print "inexist\n";
	}
	#for (my $i = 0; $i < $#AllInteractions; ++$i){
	#	if ($AllInteractions[$i][0] eq $genename){
		
		#for (my $i = 0; $i < $#geneallinteractlist; ++$i){
			#if ($geneallinteractlist[$i][2] > $threshold){
			#	push @geneinteractlist, $geneallinteractlist[$i][1];
			#}
		#}
		#}
	#}
	
	return @geneallinteractlist;
}

for (my $i = 0; $i <= $#resultarray; ++$i){ 
	print RPTFILE $resultarray[$i][0] . "\t" .  $resultarray[$i][1] . "\t" .  $resultarray[$i][2] . "\n";
}

#encode genes by their own scheme
open (INDENTFILE, "<$pathwaypath/identifier_mappings.txt");
<INDENTFILE>; #throw away 1st line

while (<INDENTFILE>) {
  chomp();
	my @indentfields = split /\t/, $_;
	
	$harshEncodings{$indentfields[1]} = $indentfields[0];
    if ( $indentfields[2] eq "Ensembl Gene Name"){
        $harshEncodingsRev{$indentfields[0]} = $indentfields[1];
    }
	#if (exists $harshEncodings{$indentfields[1]} && $indentfields[0] != $harshEncodings{$indentfields[1]}){
	#	print "multiple names: " . $indentfields[1] . " and  " . $indentfields[0] . "\n";
	#}
  #print "loading ident: " . $indentfields[1] . " as " . $indentfields[0] . "\n";
}
close(INDENTFILE);

#load individual list
my @individuallist = split /,/, $indvlist;

for (my $i = 0; $i <= $#individuallist; ++$i){
  print "Loading individual:" . $individuallist[$i] . " from $listpath" . $individuallist[$i] . $filesuffix . "\n";
	open (MYFILE, "<$listpath/" . $individuallist[$i] . $filesuffix);
	print "Loading individual:" . $individuallist[$i] . "\n";
	$harshIndividualsGenes{$individuallist[$i]} = ();
	$harshIndividualsGenesEncoded{$individuallist[$i]} = ();
	my $loadcnt = 0;
	while (<MYFILE>) {
        ++$loadcnt;
        #print "Loading gene\n";
        chomp();
        my @splitfield = split "\t", $_;
        
        if (index($splitfield[0],"exonic") > -1 || index($splitfield[0], "splicing") > -1){
            #chomp any comments
            if (index($splitfield[1], "(") > -1){ $splitfield[1] = substr($splitfield[1], 0, index($splitfield[1], "(")); print "c" . $splitfield[1] . "\n";};
            if (index($splitfield[1], ";") > -1){ $splitfield[1] = substr($splitfield[1], 0, index($splitfield[1], ";")); print "c" . $splitfield[1] . "\n";};
            if (index($splitfield[1], ",") > -1){ $splitfield[1] = substr($splitfield[1], 0, index($splitfield[1], ",")); print "c" . $splitfield[1] . "\n";};
            
            
            #see if it's with _HUMAN tag
	        if (!exists $harshEncodings{$splitfield[1]} && exists $harshEncodings{$splitfield[1] . "_HUMAN"}){
                $splitfield[1] .= "_HUMAN";
            }
            
	        if (exists $harshEncodings{$splitfield[1]}){
                push @{$harshIndividualsGenes{$individuallist[$i]}}, $splitfield[1] ;#[ ($splitfield[1], $splitfield[2], $splitfield[3]) ];
	        	push @{$harshIndividualsGenesEncoded{$individuallist[$i]}}, [ ($harshEncodings{$splitfield[1]}, $splitfield[2], $splitfield[3]) ];
	        }else{
	        	#push @{$harshIndividualsGenesEncoded{$individuallist[$i]}}, "NOT_FOUND";
	        	print $splitfield[0] . ", '" . $splitfield[1] . "' not found.\n";
	        }
        }
        
	}
	close(MYFILE);
    
    if ($loadcnt == 0){
        die "File: " . $individuallist[$i] . " is empty.";
    }
}

#load the whole interaction list
my @Interaction_Phy = <$pathwaypath/Physical_interactions.*>;
my @Interaction_Coexp = <$pathwaypath/Co-expression.*>;
my @Interaction_Pathway = <$pathwaypath/Pathway.*>;
my @Interaction_Domain = <$pathwaypath/Shared_protein_domains.*>;
my @Interaction_Gene = <$pathwaypath/Genetic_interactions.*>;

print "Loading list: Genetic_interactions\n";
for (my $i = 0; $i <= $#Interaction_Gene; $i++){
	open (ININTERACT, @Interaction_Gene[$i]);
	LoadInteractionList(\%AllInteractions, $cutoffthres);
	close (ININTERACT);
}
print "Loading list: Physical_interactions\n";
for (my $i = 0; $i <= $#Interaction_Phy; $i++){
	open (ININTERACT, $Interaction_Phy[$i]);
	LoadInteractionList(\%AllInteractions, $cutoffthres);
	close (ININTERACT);
}
#print "Loading list: Coexpression\n";
#for (my $i = 0; $i <= $#Interaction_Coexp; $i++){
#	open (ININTERACT, $Interaction_Coexp[$i]);
#	LoadInteractionList(\%AllInteractions, $cutoffthres);
#	close (ININTERACT);
#}
print "Loading list: Interaction_Pathway\n";
for (my $i = 0; $i <= $#Interaction_Pathway; $i++){
	open (ININTERACT, $Interaction_Pathway[$i]);
	LoadInteractionList(\%AllInteractions, $cutoffthres);
	close (ININTERACT);
}
print "Loading list: Shared_protein_domains\n";
for (my $i = 0; $i <= $#Interaction_Domain; $i++){
	open (ININTERACT, $Interaction_Domain[$i]);
	LoadInteractionList(\%AllInteractions, $cutoffthres);
	close (ININTERACT);
}


#grep from interactions list
foreach my $key (keys %harshIndividualsGenesEncoded){
    my $personaffected = $key;
	for (my $i = 0; $i <= $#{$harshIndividualsGenesEncoded{$key}}; ++$i){ 
		#check if the gene has been processed
		my $genenametocheck = ${$harshIndividualsGenesEncoded{$key}}[$i][0];
		print "checking " . $genenametocheck . "\n";
		if (!exists $harshPathScore{$genenametocheck}){
			my $tmpname =  ${$harshIndividualsGenesEncoded{$key}}[$i][0];
			#my $grepout = `grep -h $tmpname $pathwaypath/Physical_interactions.*`;
			
			#my $grepout .= `grep -h $tmpname $pathwaypath/Co-expression.*`;
			#my $grepout .= `grep -h $tmpname $pathwaypath/Shared_protein_domains.*`;
			#my $grepout .= `grep -h $tmpname $pathwaypath/Pathway.*`;
			#print "that means $tmpname\n";
			my @matchedlist = FindInInteractionList($tmpname);
            #push @matchedlist, [($tmpname,$tmpname,1)];
			#my @splitout = split /\n/, $grepout;
			my @linkedgenes = ();
			my @linkedgenesRev = ();
			my @allinpathgenesRev = ();
			my $strpath = "";
            my $strpathcmd = "[[:blank:]]";
			my $strpathparti = "";#$harshEncodingsRev{$genenametocheck} . " ";
			my $strpathparticmd = "[[:blank:]]";# . $harshEncodingsRev{$genenametocheck} . "[[:blank:]]|[[:blank:]]";
			my $pathlen = 0;
			my $pathlenparti = 0;
			my $pplparti = "";#$personaffected . " ";
			my $pplpartitrans = "";#$transharsh{$personaffected} . " " ;
			#my $rmed = 0;
			#my $rmstr = "";
            
            
			for (my $j = 0; $j <= $#matchedlist; ++$j){ 
                #see if there are redundant gene
                my $found=0;
                my $indtochk=1;
                if ($harshEncodingsRev{$genenametocheck} eq $harshEncodingsRev{$matchedlist[$j][1]}){
                    $indtochk=0;
                }
                for (my $x = 0; $x <= $#linkedgenesRev; $x++){
                    if ($linkedgenesRev[$x] eq $harshEncodingsRev{$matchedlist[$j][$indtochk]}){
                        $found=1;
                        last;
                    }
                }
					
                #if ($genenametocheck eq 'ENSG00000177600' && $found == 1){
                #    $rmed++;
                #    $rmstr .= $harshEncodingsRev{$matchedlist[$j][0]} . "/" . $harshEncodingsRev{$matchedlist[$j][1]} . " ";
                #}
                #add to list if not redundant
					if ($found == 0){
					  push @linkedgenes, $matchedlist[$j][$indtochk];
						push @linkedgenesRev, $harshEncodingsRev{$matchedlist[$j][$indtochk]};
						$strpath .= $harshEncodingsRev{$matchedlist[$j][$indtochk]} . " ";
						$strpathcmd .= $harshEncodingsRev{$matchedlist[$j][$indtochk]} . "[[:blank:]]|[[:blank:]]";
						$pathlen++;
					}else{
						#print "Gene already exists! \n";
					}
			}
            
            #if ($genenametocheck eq 'ENSG00000177600' ){
            #    die "total rm:$rmed \n $rmstr";
            #}
            
            #rint substr($strpathcmd, length($strpathcmd) - length("|[[:blank:]]")) . "\n";
            if (substr($strpathcmd, length($strpathcmd) - length("|[[:blank:]]")) eq "|[[:blank:]]"){
                $strpathcmd = substr($strpathcmd, 0, length($strpathcmd) - length("|[[:blank:]]") + 1);
            }
            
            if (substr($strpathcmd, length($strpathcmd) - 1, 1) == "|"){
                $strpathcmd = substr($strpathcmd, 0, length($strpathcmd) - 1);
            }
			
			my ($score, $repeatScore, $genesfoundharshref, $indvaffref) = FindScore(@linkedgenes);
			my %genesfoundharsh = %$genesfoundharshref;
			my %indvaff = %$indvaffref;
            my $proced_score1 = 0;
			foreach my $key (keys %genesfoundharsh){
                my $found=0;
                for (my $x = 0; $x <= $#allinpathgenesRev; $x++){
                    if ($key eq $allinpathgenesRev[$x]){
                        $found=1;
                        last;
                    }
                }
                
                if ($found == 0){
                    push @allinpathgenesRev, $key;
                    $strpathparti .= $key . " ";
                    $strpathparticmd .= $key . "[[:blank:]]|[[:blank:]]";
                    $pathlenparti++;
                }
			}
            
            if (substr($strpathparticmd, length($strpathparticmd) - length("|[[:blank:]]")) eq "|[[:blank:]]"){
                $strpathparticmd = substr($strpathparticmd, 0, length($strpathparticmd) - length("|[[:blank:]]") + 1);
            }
            
            if (substr($strpathparticmd, length($strpathparticmd) - 1, 1) == "|"){
                $strpathparticmd = substr($strpathparticmd, 0, length($strpathparticmd) - 1);
            }
            
            my $NMCnt = 0;
            my $MCnt = 0;
            my $NTCnt = 0;
            my $ratio = 0;
            #find the groupping of the pathway
			foreach my $key (keys %indvaff){
                $pplpartitrans .= $transharsh{$key} . " ";
                if ($transharsh{$key} eq "LNM" || $transharsh{$key} eq "SNM" ){
                    ++$NMCnt;
                }elsif ($transharsh{$key} eq "NT" ){
                    ++$NTCnt;
                }else{
                    ++$MCnt;
                }
				$pplparti .= $key . " ";
			}
            if ($MCnt == 0){$ratio = $NMCnt;}
            elsif ($NMCnt == 0 && $MCnt > 0) {$ratio = 1;}
            else {$ratio = $MCnt/$NMCnt;}
                
            #output the results
			print "Score: $score\n"; 
            if ($score > 0){
                $proced_score1 = $score / $pathlen;
            }
			push @resultarray, ($harshEncodingsRev{$genenametocheck},$score , $pathlen . "\t" . $strpath);
			
			$harshPathScore{$harshEncodingsRev{$genenametocheck}} = $score;
			$harshPathListparti{$harshEncodingsRev{$genenametocheck}} = $pplparti . "\t" . $pplpartitrans . "\t". $pathlenparti . "\t" . $strpathparti;
			$harshPathList{$harshEncodingsRev{$genenametocheck}} =  $pathlen . "\t" . $proced_score1  . "\t" . $strpath . "\t" .  $strpathparticmd . "\t" .  $strpathcmd . "\t`$NMCnt/$MCnt/$NTCnt\t$ratio\t$repeatScore";
			#print $grepout;
		}else{
			print ${$harshIndividualsGenesEncoded{$key}}[$i][0] . " already processed.\n";
		}
	}
}

#my @resultarray2 = sort { $a[1] <=> $b[1] } @resultarray;

open (RPTFILE, ">$listpath/pathreport.rpt");
foreach my $key (keys %harshPathScore){
	print RPTFILE "$key\t" . $harshPathScore{$key} . "\t" . $harshPathListparti{$key} . "\t" . $harshPathList{$key} . "\n";
}
close(RPTFILE);
#foreach my $list_ref ( sort { $a->[1] <=> $b->[1] } @resultarray ) {
#		my @fieldsarr = @$list_ref;
#    print $fieldsarr[0] . "\t" .  $fieldsarr[1] . "\t" .  $fieldsarr[1] . "\n";
#} 


die;

open (MYFILE, "<$pathwaypath/identifier_mappings.txt");
print "Loading junction\n";
<MYFILE>; #ignores first line
while (<MYFILE>) {
    #$refRange[$refindex][0] = $splitarr[1];
    #$refRange[$refindex][1] = $splitarr[2];
    #    print $refRange[$refindex][0] . " " . $refRange[$refindex][1] . "\n";
    
}
close(MYFILE);
