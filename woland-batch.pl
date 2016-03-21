########################################################################################################################################
## WOLAND Beta 0.1 (03-08-2016)
## woland-batch.pl
##
## WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. 
## It is implemented as a Perl and R tool using as inputs filtered unannotated or annotated SNV lists, combined with its 
## correspondent genome sequences.
## 
## For more details please read README file.
## 
## Use woland-batch to run multiple samples as in <input-table> using woland-anno.pl and build a grouped report using woland-report.pl.
##
## USAGE:
##
## woland-batch.pl <input_table> <chromosome_length_profile> <hotspot_window_length> <genome_version>
##
########################################################################################################################################

#! /usr/bin/perl
use IPC::System::Simple qw(system capture);
use Parallel::ForkManager;
use File::Copy;
use strict;
use warnings;

## variables
my ($inputTable, $inputTableLine); 
my @inputTableArray;  
my @i; my $i;
my @Group; my @sampleName;
my @time;
my ($profile, $hotspot, $genomeVersion);
my ($pm,$pid, $MAX_PROCESSES);
my @ARGS;
my @time1;

## main warning
unless (@ARGV){
	die "\nERROR : Usage: $0 <input.table> <chromosome_profile_file> <hotspot_window_length> <genome_version> \n";	
}

@time=localtime; #time characters for analysis folder name

## parsing input table
$inputTable = $ARGV[0]; #<input_table>
open (inputTable, $inputTable);
@inputTableArray=<inputTable>;

foreach $inputTableLine (@inputTableArray){ #two arrays for each category (group & sample results folder)
	@i = split (/\t/, $inputTableLine);
	chomp (@i);
	push (@Group, "$i[0]"); # array for group definition
	push (@sampleName, "$i[1]"); # array for sample folder definition
}

@i=();
$i=0;


mkdir("results-batch-$inputTable-$time[0].$time[1].$time[2].$time[3].$time[4].$time[5]", 0755) || die "Cannot create results folder - check if it already exists";
mkdir("results-batch-$inputTable-$time[0].$time[1].$time[2].$time[3].$time[4].$time[5]/samples-$inputTable", 0755) || die "Cannot create results folder- check if it already exists";

# executing batch woland_anno.pl:
$profile= $ARGV[1]; #chromosome profile <chromosome_length_profile>
$hotspot= $ARGV[2]; #natural number for hotspot window <hotspot_window_length> 
$genomeVersion = $ARGV[3]; #genome version as in genomes/genome_<genome_version>.fa and genomes/refseq_<genome_version>.txt

$pm = Parallel::ForkManager->new(30);


for my $i (0..$#sampleName){ #execution of woland-anno.pl for each sample
	@ARGS=();
	push (@ARGS, $sampleName[$i]);
	push (@ARGS, $profile);
	push (@ARGS, $hotspot);
	push (@ARGS, $genomeVersion);

	$pid=$pm->start and next;
	system ($^X, "woland-anno.pl", @ARGS);
	$pm->finish;
}

$pm->wait_all_children;

system ($^X, "woland-report.pl", $ARGV[0]);

# moving report folder and files to results-batch
move ("report-$inputTable", "results-batch-$inputTable-$time[0].$time[1].$time[2].$time[3].$time[4].$time[5]/report-$inputTable");

# moving each result sample folder and files to results-batch/results
for my $i (0..$#inputTableArray){
	move ("results-$sampleName[$i]", "results-batch-$inputTable-$time[0].$time[1].$time[2].$time[3].$time[4].$time[5]/samples-$inputTable/results-$sampleName[$i]");
}

@time1=localtime;

print "@time\n";
print "@time1\n";

exit;
