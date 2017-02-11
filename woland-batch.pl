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
use Getopt::ArgParse;

our $REVISION = '$Revision:  $';
our $DATE =	'$Date: 2017-01-28 00:11:04 -0800 (Sat,  28 Jan 2017) $';  
our $AUTHOR =	'$Author: Tiago A. de Souza <tiagoantonio@gmail.com> $';

## global variables
our ($inputtable, $tableline); our (@tablearray,@group,@sample); #input table and sample/group parsing
our ($profile, $hotspot, $genome); #arguments
our (@starttime,@endtime); #time tracking

## confi variables
my $pm; #parallel fork manager process

my $ap = Getopt::ArgParse->new_parser(
	prog => 'Woland',
	description => 'WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell.',
	epilog => '',
 );

$ap->add_arg('--input-table', '-i', required => 1, help => 'Help of Input table');
## It is possible to set a value as default?? default => 10
$ap->add_arg('--chromosome-length-profile', '-c', dest => 'chr_length', required => 1, help => 'Help of Chromosome length profile');
## It is possible to set a value as default?? default => 10
$ap->add_arg('--hotspot-window-length', '-w', dest => 'hotspot', required => 1, help => 'Help of Hotspot Window Length');
$ap->add_arg('--genome-version', '-g', dest => 'genome', required => 1, help => 'Help of Genome Version');
$ap->add_arg('--threads', '-t', default => 30, help => 'Help of Threads');

my $args = $ap->parse_args();

my $fastagenomeversion=sprintf("genomes/genome_%s.fa", $args->genome);
unless (-r -e -f $fastagenomeversion){
	die "\nERROR : Please check if a genome fasta file exists in genomes\/folder for <genome_version>\n"
}

unless ($args->hotspot>0){
	die "\nERROR : Please specify a natural number >0 for <hotspot_window_length>\n";
}

unless (-r -e -f $args->input_table){
	die sprintf("\nERROR: %s not exists or is not readable or not properly formatted. Please check file.\n\n",
		$args->input_table);
}

@starttime=localtime; #time characters for analysis folder name

## parsing input table
$inputtable = $args->input_table; #<input_table>
open (INPUT, $inputtable);
@tablearray=<INPUT>;

foreach $tableline (@tablearray){ #two arrays for each category (group & sample results folder)
	my @i = split (/\t/, $tableline);
	chomp (@i);
	push (@group, "$i[0]"); #array for group definition
	push (@sample, "$i[1]"); #array for sample folder definition
}

## arguments to woland_anno.pl:
$profile= $args->chr_length; #chromosome profile <chromosome_length_profile>
$hotspot= $args->hotspot; #natural number for hotspot window <hotspot_window_length> 
$genome = $args->genome; #genome version as in genomes/genome_<genome_version>.fa and genomes/refseq_<genome_version>.txt

## creating output directories
mkdir("results-batch-$inputtable-$starttime[0].$starttime[1].$starttime[2].$starttime[3].$starttime[4].$starttime[5]", 0755) || die "Cannot create results folder - check if it already exists";
mkdir("results-batch-$inputtable-$starttime[0].$starttime[1].$starttime[2].$starttime[3].$starttime[4].$starttime[5]/samples-$inputtable", 0755) || die "Cannot create results folder- check if it already exists";

## woland-anno.pl multi-threading
$pm = Parallel::ForkManager->new($args->threads);

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
system ($^X, "woland-report.pl", $args->input_table);

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