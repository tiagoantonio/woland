########################################################################################################################################
## WOLAND Beta 0.2 (01-28-2016)
## woland-mutageniclogo.pl
##
## WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. 
## 
## For more details please read README file.
## 
## Use woland-mutageniclogo to extract stranded context sequences FASTA files specific for C>T, G>T and G>C changes AFTER woland-report.
##
## USAGE:
##
## woland-mutageniclogo.pl <sensechange> <input-table>
##
## sense change: CT,GT or GC.
##
########################################################################################################################################

#! /usr/bin/perl
use strict;
use warnings;
use Getopt::ArgParse;

our $REVISION = '$Revision:  $';
our $DATE =	'$Date: 2017-01-28 00:11:04 -0800 (Sat,  28 Jan 2017) $';  
our $AUTHOR =	'$Author: Tiago A. de Souza <tiagoantonio@gmail.com> $';

#global variables
our ($table, $sensechange);
our ($revcomp);
our (@sensesequence);
our (@targetsequence);

sub parse_inputtable_to_sensechanges{ #parsing input-table, extract context sequence and apply reverse-complement to non-sense changes
	my $tableline;
	open (TABLE, $table);
	my @tablearray=<TABLE>;
	close(TABLE);
	my (@samplename,@group,@targetsequence,@sensesequence);

	foreach $tableline (@tablearray){ # two arrays for each category (group & sample results folder)
		my @i = split (/\t/, $tableline);
		chomp (@i);
		push (@group, "$i[0]"); # array for group definition
		push (@samplename, "$i[1]"); # array for sample folder definition
	}

	for my $i (0..$#samplename){
		open CONTEXTSEQANNO, "<results-batch-$table/samples-$table/results-$samplename[$i]/WOLAND-contextsequencesanno-$samplename[$i]";
		my @fastacontext=<CONTEXTSEQANNO>;
		close (CONTEXTSEQANNO);
		
		my @fastacontextprocessed=();

		for my $ii (0..$#fastacontext){
			if ($fastacontext[$ii] =~ /^>/){
				my @sequencesplitted = split (/\t/, "$fastacontext[$ii+1]");
				my @idsplitted = split (/\n/, "$fastacontext[$ii]");
				chomp(@sequencesplitted);
				push (@fastacontextprocessed, "$idsplitted[0]\n$sequencesplitted[0]\n");
			}
		}
		my $fastacontextprocessedline;
		@targetsequence=();
		foreach $fastacontextprocessedline (@fastacontextprocessed){
			my @iii = split (/\n/, "$fastacontextprocessedline");
			push (@targetsequence, $iii[1]);
		}

		for my $ii (0..$#targetsequence){ #C>T change
			if ($sensechange=~"CT"){
				if($targetsequence[$ii]=~"...C..."){
					push (@sensesequence, "$targetsequence[$ii]");
				}
				if($targetsequence[$ii]=~"...G..."){
					&to_reversecomplement($targetsequence[$ii]);
					push (@sensesequence, "$revcomp");
				}
			}

			if ($sensechange=~"GT"){ #G>T change
				if($targetsequence[$ii]=~"...G..."){
					push (@sensesequence, "$targetsequence[$ii]");
				}
				if($targetsequence[$ii]=~"...C..."){
					&to_reversecomplement($targetsequence[$ii]);
					push (@sensesequence, "$revcomp");
				}
			}

			if ($sensechange=~"GC"){ #G>C change
				if($targetsequence[$ii]=~"...G..."){
					push (@sensesequence, "$targetsequence[$ii]");
				}
				if($targetsequence[$ii]=~"...C..."){
					&to_reversecomplement($targetsequence[$ii]);
					push (@sensesequence, "$revcomp");
				}
			}	
		}
		open (LOGO, ">>results-batch-$table/samples-$table/results-$samplename[$i]/WOLAND-LOGO-$sensechange-$samplename[$i]");
		for my $ii (0..$#sensesequence){
			print LOGO "\>$i\n$sensesequence[$ii]\n";
		}
		close (LOGO);
		@sensesequence=();
	}
}

sub to_reversecomplement{ #reversecomplement based on sense-change
	my $dna="$_[0]";
	$dna = shift;
	$revcomp = reverse($dna); # reverse the DNA sequence
	$revcomp =~ tr/ACGTacgt/TGCAtgca/; # complement the reversed DNA sequence
	return $revcomp;
}

my $ap = Getopt::ArgParse->new_parser(
	prog => 'woland-mutageniclogo.pl',
	description => 'WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing SNV data.
	Use woland-mutageniclogo to extract stranded context sequences FASTA files specific for C>T, G>T and G>C changes AFTER woland-report. For more details please read README',
	epilog => 'If you used Woland in your research, we would appreciate your citation:
	de Souza TA, Defelicibus A, Menck CF',
 );

$ap->add_arg(
	'--sense-change',
	'-s',
	required => 1,
	choices => [ 'CT', 'GT', 'GC' ],
	help => 'Type of sense-nucleotide change to filter in extracted sequences');
$ap->add_arg(
	'--input-table',
	'-i',
	required => 1,
	help => 'Tab-delimited file with samples in the 1st column and groups in the 2nd column');

my $args = $ap->parse_args();

## main warning
# unless ($#ARGV>=1){
# 	die "\nERROR : Incorrect number of arguments/files - Usage: $0 <sensechange> <input-table> \n\n";	
# }

# unless ($ARGV[0] eq "GT" || $ARGV[0] eq "CT" || $ARGV[0] eq "GC"){
# 	die "\nERROR : Incorrect <typeofchange>. Please use CT,GT or GC\n\n";	
# }
unless (-r -e -f $args->input_table){
	die sprintf("\nERROR: %s not exists or is not readable or not properly formatted. Please check file.\n\n",
		$args->input_table);	
}

$sensechange = $args->sense_change;
$table = $args->input_table;

&parse_inputtable_to_sensechanges;

exit;