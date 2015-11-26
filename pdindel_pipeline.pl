#!/usr/bin/perl
use strict;
use POSIX;
use XML::Parser;
use XML::XPath;
use XML::XPath::XMLParser;
use LWP::Simple;
use Getopt::Long;
use threads ('yield',
			 'stack_size' => 64*4096,
			 'exit' => 'threads_only',
			 'stringify');
			 
my $usage=<<INFO;

Usage:
perl $0 -r <reference fasta> -a <alignment file> -f <0=sam, 1=bam> -p <pindel path> -d <dindel path> -s <insert size> -t <max threads> -v <ref version, eg:hg19>

example:  perl pdindel_pipeline.pl  -r /data/genome/hg19/ucsc.hg19.fasta -a /store/RNA_team/Exom_Validation_VC/A100096A/A100096A_lane1.very_sensitive.bowtie2.sam  -p ~/Programs/indel_caller/pindel -d ~/Programs/indel_caller/dindel -s 150 -t 5 -v hg19

INFO

my ($ref_fasta, $alignment, $isbam, $pipath, $dipath, $insertsz, $maxthreads, $ref_ver);
GetOptions(
"r=s"=>\$ref_fasta,
"a=s"=>\$alignment,
"f=i"=>\$isbam,
"p=s"=>\$pipath,
"d=s"=>\$dipath,
"s=i"=>\$insertsz,
"t=i"=>\$maxthreads,
"v=s"=>\$ref_ver,
);

die "$usage" if(!$ref_fasta || !$alignment || !$pipath || !$dipath || !$insertsz || !$maxthreads || !$ref_ver); 
die "Alignment file does not exists!: $alignment " if (!(-e $alignment));
die "Reference file does not exists!: $ref_fasta " if (!(-e $ref_fasta));

sub start_thread {
	my $cmdrun = shift;
	print "tCMD: $cmdrun \n";
	system($cmdrun);
}

sub wait_thread_end{
	my $threadlist = shift;
	my $thr_lim = shift;
	my $threadcnt = 999;
	
	while ($threadcnt > $thr_lim){
		$threadcnt = 0;
		for (my $i = @$threadlist - 1; $i >= 0; --$i){
			if ($$threadlist[$i]->is_running()){
				$threadcnt++;
			#}elsif ($$threadlist[$i]->is_joinable()){
				#my $thr_cnt = @$threadlist;
				#print "arrsz:" .  $thr_cnt . "\n";
				#$$threadlist[$i]->join();
				#splice(@$threadlist, $i);
			}
		}
		
		print "thr_cnt:$threadcnt\n";
		if ($threadcnt > $thr_lim){
			sleep 5;
		}
	}
}
	
sub runcmd{
	my $cmd = shift;
	print "CMD: $cmd \n";
	system($cmd);
}

sub waitthread{
	my $grepstr = shift;
	print "sleep?\n";
	my $proccnt = `ps alh | grep doDiploid | grep $grepstr | wc -l`;
	while ($proccnt > $maxthreads){
		print "waiting prc: $proccnt\n";
		$proccnt = `ps alh | grep doDiploid | grep $grepstr | wc -l`;
		sleep 5;
	}
}

sub termhdl {
	warn"Hasta la vista, baby!"
}

$SIG{'TERM'} = termhdl;
my $cmd="";
#sam/bam to pindel
my $fileprefix = substr($alignment, 0, rindex($alignment, "."));
my $workingdir = substr($alignment, 0, rindex($alignment, "/"));
my $shortfileprefix = substr($alignment,rindex($alignment, "/") + 1);
$shortfileprefix = substr($shortfileprefix, 0, rindex($shortfileprefix, "."));
my $bamfile = $fileprefix . ".bam";
my $samfile = $fileprefix . ".sam";
print "Alignment file: $alignment \n";
print "Working dir: $workingdir \n";
print "file prefixes: $fileprefix \n";
print "samfile: $samfile \n";
print "short prefix: $shortfileprefix \n";
#some checking 
if ($isbam){
	if (!(-e $bamfile)){
		die "Cannot open bam file\n";
	}
}else{
	if (!(-e $samfile)){
		die "Cannot open sam file\n";
	}
}

#goto START;
if ($isbam){
	print "Converting to sam format\n";
	$cmd = "samtools view $alignment > $samfile"
}
runcmd($cmd);

print "Converting to pindel format\n";
$cmd = "$pipath/sam2pindel/sam2pindel $samfile $fileprefix.pi $insertsz $fileprefix 0 > $fileprefix.cvt_pi.log";
runcmd($cmd);
if ($isbam){
	$cmd = "rm $samfile";
	runcmd($cmd);
}

#run pindel
$cmd = "touch emt";
runcmd($cmd);
#by chromosome
for (my $i = 1; $i < 25; ++$i){
	my $chrname;
	if ($i <= 22){
		$chrname = "chr$i";
	}elsif($i == 23){
		$chrname = "chrX";
	}elsif($i == 24){
		$chrname = "chrY";
	}
	print "Runing pindel on $chrname\n";
	$cmd = "$pipath/pindel_x86_64 $ref_fasta $fileprefix.pi $workingdir $workingdir/emt 4 $maxthreads $chrname > $workingdir/pi_$chrname.log";
	runcmd($cmd);
}

#convert pindel to vcf
#merge back the file
$cmd = "cat $fileprefix.pi_*D $fileprefix.pi_*I > $fileprefix.all.pi";
runcmd($cmd);
$cmd = "rm $fileprefix.pi_*D $fileprefix.pi_*I";
runcmd($cmd);
$cmd = "$pipath/pindel2vcf -r $ref_fasta -R $ref_ver -p $fileprefix.all.pi -d 20120426 -v $fileprefix.all.pi.vcf ";
runcmd($cmd);

#convert to dindel
$cmd = "cd $dipath";
runcmd($cmd);
$cmd = "python $dipath/convertVCFToDindel.py --inputFile $fileprefix.all.pi.vcf -o $fileprefix.all.pi.di --refFile $ref_fasta";
runcmd($cmd);
$cmd = "rm $fileprefix.pi";
runcmd($cmd);

#produce BAM for dindel
if (!$isbam){
	print "Converting to bam format\n";
	$cmd = "samtools view -bS $samfile > $bamfile";
	runcmd($cmd);
}
#check index file
if (!(-e "$bamfile.bai")){
	print "indexing bam file\n";
	$cmd = "samtools sort $bamfile $bamfile.sorted";
	runcmd($cmd);
	#$cmd = "mv $bamfile $bamfile.unsorted.bam";
	#runcmd($cmd);
	#$cmd = "mv $bamfile.sorted.bam $bamfile";	#removed afterwards
	#runcmd($cmd);
	$cmd = "samtools index $bamfile.sorted.bam";
	runcmd($cmd);
}

#import pindel df to dindel
$cmd = "$dipath/dindel --analysis getCIGARindels --bamFile $bamfile.sorted.bam --outputFile $fileprefix.candidate.di --ref $ref_fasta > $fileprefix.cigar.di.log";
runcmd($cmd);
$cmd = "$dipath/dindel --analysis realignCandidates --varFile $fileprefix.all.pi.di --outputFile $fileprefix.picandidate.di --ref $ref_fasta";
runcmd($cmd);

#dindel make windows
$cmd = "cat $fileprefix.candidate.di.variants.txt $fileprefix.picandidate.di.variants.txt > $fileprefix.pidi.candidate.variants.txt";
runcmd($cmd);
$cmd = "python $dipath/makeWindows.py --inputVarFile $fileprefix.pidi.candidate.variants.txt --windowFilePrefix $fileprefix.di.window --numWindowsPerFile 1000";
runcmd($cmd);

opendir(IMD, $workingdir) || die("Cannot open directory of dindel window");
my @dircontent = readdir(IMD);
my @threadlist;
my $proclim = 0;
print "Workdir fcnt:" . $#dircontent . "\n";
open (FLIST, ">$fileprefix.di.window.prc.list");
for (my $i = 0; $i <= $#dircontent; ++$i){
	if ($dircontent[$i] =~ /^$shortfileprefix.di.window.[0-9]{1,4}.txt$/){
		#++$proclim;
		#if ($proclim > 20){ last; }
		print $dircontent[$i] . "\n";
		#print "$fileprefix.di.window\n";
		$cmd = "$dipath/dindel --analysis indels --doDiploid --bamFile $bamfile.sorted.bam --ref $ref_fasta --varFile $workingdir/".$dircontent[$i]." --libFile $fileprefix.candidate.di.libraries.txt --outputFile $workingdir/".$dircontent[$i].".prc >> $fileprefix.windprc_log";
		my $thr = threads->create('start_thread', $cmd);
		push @threadlist, $thr;
		print FLIST $workingdir . "/" . $dircontent[$i].".prc.glf.txt\n";
		
		wait_thread_end(\@threadlist,$maxthreads);
	}
}
close FLIST;
print "Waiting to end\n";
wait_thread_end(\@threadlist,0);

$cmd = "cd $dipath";
runcmd($cmd);
$cmd = "python $dipath/mergeOutputDiploid.py --inputFiles $fileprefix.di.window.prc.list --outputFile $fileprefix.pi.di.final.vcf --ref $ref_fasta";
runcmd($cmd);

#flushing
for (my $i = 0; $i <= $#dircontent; ++$i){
	if ($dircontent[$i] =~ /^$shortfileprefix.di.window.[0-9]{1,4}.txt$/){
		unlink($workingdir . "/" . $dircontent[$i].".prc.glf.txt");
		unlink($workingdir . "/" . $dircontent[$i]);   
	}
}

if ($isbam){
	$cmd = "rm $samfile";
	runcmd($cmd);
}else{
	$cmd = "rm $bamfile";
	runcmd($cmd);
}
$cmd = "rm $fileprefix.candidate.di.variants.txt";
runcmd($cmd);
$cmd = "rm $fileprefix.candidate.di.libraries.txt";
runcmd($cmd);
$cmd = "rm $fileprefix.candidate.di.libraries.txt";
runcmd($cmd);

#$cmd = "rm $bamfile.sorted.bam";
#runcmd($cmd);