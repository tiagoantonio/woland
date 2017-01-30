########################################################################################################################################
## WOLAND Beta 0.2 (01-28-2016)
## woland-batch.pl
##
## WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. 
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

our $REVISION = '$Revision:  $';
our $DATE =	'$Date: 2017-01-28 00:11:04 -0800 (Sat,  28 Jan 2017) $';  
our $AUTHOR =	'$Author: Tiago A. de Souza <tiagoantonio@gmail.com> $';

## global variables
our ($inputtable, $tableline); our (@tablearray,@group,@sample); #input table and sample/group parsing
our ($profile, $hotspot, $genome); #arguments
our (@starttime,@endtime); #time tracking

## confi variables
my $MAXPROCESSES=30; #threads
my $pm; #parallel fork manager process

## warnings and checks
unless ($#ARGV==3){
	die "\nERROR : Incorrect number of arguments - Usage: $0 <input.table> <chromosome_profile_file> <hotspot_window_length> <genome_version> \n\n";	
}
for my $i (0..1) {
	unless (-r -e -f $ARGV[$i]){
    	die "\nERROR: $ARGV[$i] not exists or is not readable or not properly formatted. Please check file.\n\n";
    }
}

@starttime=localtime; #time characters for analysis folder name

## parsing input table
$inputtable = $ARGV[0]; #<input_table>
open (INPUT, $inputtable);
@tablearray=<INPUT>;

foreach $tableline (@tablearray){ #two arrays for each category (group & sample results folder)
	my @i = split (/\t/, $tableline);
	chomp (@i);
	push (@group, "$i[0]"); #array for group definition
	push (@sample, "$i[1]"); #array for sample folder definition
}

## arguments to woland_anno.pl:
$profile= $ARGV[1]; #chromosome profile <chromosome_length_profile>
$hotspot= $ARGV[2]; #natural number for hotspot window <hotspot_window_length> 
$genome = $ARGV[3]; #genome version as in genomes/genome_<genome_version>.fa and genomes/refseq_<genome_version>.txt

## creating output directories
mkdir("results-batch-$inputtable-$starttime[0].$starttime[1].$starttime[2].$starttime[3].$starttime[4].$starttime[5]", 0755) || die "Cannot create results folder - check if it already exists";
mkdir("results-batch-$inputtable-$starttime[0].$starttime[1].$starttime[2].$starttime[3].$starttime[4].$starttime[5]/samples-$inputtable", 0755) || die "Cannot create results folder- check if it already exists";

## woland-anno.pl multi-threading
$pm = Parallel::ForkManager->new($MAXPROCESSES);

for my $i (0..$#sample){ #execution of woland-anno.pl for each sample
	my @arguments;
	push (@arguments, $sample[$i]);
	push (@arguments, $profile);
	push (@arguments, $hotspot);
	push (@arguments, $genome);

	my $pid=$pm->start and next;
	system ($^X, "woland-anno.pl", @arguments);
	$pm->finish;
}
$pm->wait_all_children; #wait woland-anno.pl

## woland-report.pl
system ($^X, "woland-report.pl", $ARGV[0]);

## moving report folder and files to results-batch
move ("report-$inputtable", "results-batch-$inputtable-$starttime[0].$starttime[1].$starttime[2].$starttime[3].$starttime[4].$starttime[5]/report-$inputtable");

## moving each result sample folder and files to results-batch/results
for my $i (0..$#tablearray){
	move ("results-$sample[$i]", "results-batch-$inputtable-$starttime[0].$starttime[1].$starttime[2].$starttime[3].$starttime[4].$starttime[5]/samples-$inputtable/results-$sample[$i]");
}

@endtime=localtime;

print "Start Time: @starttime\n";
print "End Time: @endtime\n";

exit;