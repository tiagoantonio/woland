#######################################################################################################################################
## WOLAND BATCH BETA 0.1 (15-02-2016)
##
## WOLAND is a software package based on Perl and R for calculation of general mutation metrics, identification and
## comparison of predicted hotspots across biological samples. WOLAND uses Single Nucleotide Polymorphisms (SNPs) data
## from Next Generation Sequencing (NGS) pipelines as main input. Please read README file.
## 
## Use woland-batch to run multiple samples as in <input-table> using woland-anno.pl and build a grouped report using woland-report.pl.
##
## USAGE:
##
## woland-batch.pl <input_table> <chromosome_length_profile> <hotspot_window_length> <genome_version>
##
#######################################################################################################################################

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
my ($profile, $hotspot, $genome_version);
my ($pm,$pid, $MAX_PROCESSES);
my @ARGS;
my @time;

## main warning
unless (@ARGV){
	die "\nERROR : Usage: $0 <input.table> <chromosome_profile_file> <hotspot_window_length> <genome_version> \n";	
}

## parsing input table
$inputTable = $ARGV[0]; #<woland.input.table>
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

@time=localtime;

mkdir("results-batch-$inputTable-$time[0].$time[1].$time[2].$time[3].$time[4].$time[5]", 0755) || die "Cannot create results folder - check if it already exists";
mkdir("results-batch-$inputTable-$time[0].$time[1].$time[2].$time[3].$time[4].$time[5]/samples-$inputTable", 0755) || die "Cannot create results folder- check if it already exists";

# executing batch woland_anno.pl:
$profile= $ARGV[1];
$hotspot= $ARGV[2];
$genome_version = $ARGV[3];

for my $i (0..$#sampleName){
	push (@ARGS, $sampleName[$i]);
	push (@ARGS, $profile);
	push (@ARGS, $hotspot);
	push (@ARGS, $genome_version);
	
	$pm = Parallel::ForkManager->new($MAX_PROCESSES);
	$pid=$pm->start and next;
	system ($^X, "woland-anno.pl", @ARGS);
	$pm->finish;
	@ARGS=();
}

$pm = Parallel::ForkManager->new($MAX_PROCESSES);
$pid=$pm-> start and next;
system ($^X, "woland-report.pl", $ARGV[0]);
$pm->finish;

# moving folders and files to results-batch
move ("report-$inputTable", "results-batch-$inputTable-$time[0].$time[1].$time[2].$time[3].$time[4].$time[5]/report-$inputTable");

for my $i (0..$#inputTableArray){
	move ("results-$sampleName[$i]", "results-batch-$inputTable-$time[0].$time[1].$time[2].$time[3].$time[4].$time[5]/samples-$inputTable/results-$sampleName[$i]");
}

exit;