########################################################################################################################################
## WOLAND Beta 1.01 (09-30-2017)
## woland-batch.pl
##
## WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. 
## 
## For more details please read README file.
## 
########################################################################################################################################

#! /usr/bin/perl
use IPC::System::Simple qw(system capture);
use Parallel::ForkManager;
use File::Copy;
use File::Spec;
use strict;
use warnings;
use Getopt::ArgParse;


our $REVISION = '$Revision:  $';
our $DATE =	'$Date: 2017-09-30 00:11:04 -0800 (Sat,  30 Sep 2017) $';  
our $AUTHOR =	'$Author: Tiago A. de Souza <tiagoantonio@gmail.com> $';

## global variables
our ($inputtable, $tableline); our (@tablearray,@group,@sample); #input table and sample/group parsing
our ($profile, $hotspot, $genome); #arguments
our (@starttime,@endtime); #time tracking

## confi variables
my $pm; #parallel fork manager process

my $ap = Getopt::ArgParse->new_parser(
	prog => 'woland-batch.pl',
	description => 'WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing SNV data. 
	Use woland-batch to run multiple samples as in <input-table> using woland-anno.pl and build a grouped report using woland-report.pl. For more details please read README',
	epilog => 'If you used Woland in your research, we would appreciate your citation:
	de Souza TA, Defelicibus A, Menck CF',
 );

$ap->add_arg(
	'--input-table',
	'-i',
	required => 1,
	help => 'Tab-delimited file with samples in the 1st column and groups in the 2nd column');
$ap->add_arg(
	'--chromosome-length-profile',
	'-c',
	dest => 'chr_length',
	required => 1,
	help => 'Tab-delimited file with chr names in the 1st column and target sequenced length in the 2nd column');
$ap->add_arg(
	'--hotspot-window-length',
	'-w',
	dest => 'hotspot',
	required => 1,
	default => 1000,
	help => 'Natural number for hotspot window-length');
$ap->add_arg(
	'--genome-path',
	'-g',
	dest => 'genome',
	required => 1,
	help => 'String for genome path for genome and annotation files.');
$ap->add_arg(
	'--genome-name',
	'-n',
	dest => 'genome_name',
	required => 1,
	help => 'String for genome name for genome and annotation files.');
$ap->add_arg(
	'--refseq',
	'-r',
	dest => 'refseq',
	required => 1,
	help => 'String for complete path and file of refseq.');
$ap->add_arg(
	'--threads',
	'-t',
	default => 30,
	help => 'Set a number for the maximum number of threads');
$ap->add_arg(
	'--output',
	'-o',
	help => 'Output folder where all files will be created.');

my $args = $ap->parse_args();

# my $fastagenomeversion=sprintf("genomes/genome_%s.fa", $args->genome);
my $fastagenomeversion = File::Spec->catfile($args->genome, sprintf("%s.fa", $args->genome_name));
unless (-r -e -f $fastagenomeversion){
	die "\nERROR : Please check if a genome fasta file exists in $fastagenomeversion.\n"
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
my $genome_name = $args->genome_name;

## creating output directories
my $output_folder = File::Spec->catfile($args->output, "results-batch-$inputtable-$starttime[0].$starttime[1].$starttime[2].$starttime[3].$starttime[4].$starttime[5]");
mkdir($output_folder, 0755) || die "Cannot create results folder - check if it already exists";
mkdir(File::Spec->catfile($output_folder, "samples-$inputtable"), 0755) || die "Cannot create results folder- check if it already exists";

## woland-anno.pl multi-threading

$pm = Parallel::ForkManager->new($args->threads);

for my $i (0..$#sample){ #execution of woland-anno.pl for each sample
	my @arguments;
	push (@arguments, "-i");
	push (@arguments, $sample[$i]);
	push (@arguments, "-c");
	push (@arguments, $profile);
	push (@arguments, "-w");
	push (@arguments, $hotspot);
	push (@arguments, "-g");
	push (@arguments, $genome);
	push (@arguments, "-n");
	push (@arguments, $genome_name);
	push (@arguments, "-r");
	push (@arguments, $args->refseq);
	push (@arguments, "-o");
	push (@arguments, $args->output);

	if ($i==0){
		system ($^X, "woland-anno.pl", @arguments);
	}

	else{
		my $pid=$pm->start and next;
		system ($^X, "woland-anno.pl", @arguments);
		$pm->finish;
	}
}
$pm->wait_all_children; #wait woland-anno.pl

## woland-report.pl
system ($^X, "woland-report.pl", "-i", $args->input_table);

## moving report folder and files to results-batch
move ("report-$inputtable", File::Spec->catfile($output_folder, "report-$inputtable"));

## moving each result sample folder and files to results-batch/results
for my $i (0..$#tablearray){
	move ("results-$sample[$i]", File::Spec->catfile($output_folder, "samples-$inputtable", "results-$sample[$i]"));
}

@endtime=localtime;

print "Start Time: @starttime\n";
print "End Time: @endtime\n";

exit;