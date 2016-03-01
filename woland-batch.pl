
#########################################################################################################################
## WOLAND BATCH BETA 0.1 (15-02-2016)
##
## WOLAND is a software package based on Perl and R for calculation of general mutation metrics, identification and
## comparison of predicted hotspots across biological samples. WOLAND uses Single Nucleotide Polymorphisms (SNPs) data
## from Next Generation Sequencing (NGS) pipelines as main input. Please observe file requirements.
## 
## USES input.table.sample to run woland-anno.pl on samples to report.
#
# woland-batch.pl woland.input.table chromosome_profile_path hotspot_window_length
#
#########################################################################################################################


#! /usr/bin/perl


use IPC::System::Simple qw(system capture);
use Parallel::ForkManager;
use strict;
use warnings;


## variables

my $inputTable; my @inputTableArray; my $inputTableLine; 
my @i; my $i;
my @Group; my @sampleName;
my ($profile, $hotspot, $genome_version);
my $abc;
my ($pm,$pid, $MAX_PROCESSES);
my @ARGS;


## main Warning

unless (@ARGV){
	die "\nERROR : Usage: $0 <woland.input.table> <chromosome_profile_file> <hotspot_window_length> <genome_version> \n";	
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

@i=0;

#executing batch woland_annopl:

$profile= $ARGV[1];
$hotspot= $ARGV[2];
$genome_version = $ARGV[3];

for my $abc (0..$#sampleName){
	push (@ARGS, $sampleName[$abc]);
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

exit;