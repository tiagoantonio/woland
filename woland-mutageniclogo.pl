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

## main warning
unless ($#ARGV>=1){
	die "\nERROR : Incorrect number of arguments/files - Usage: $0 <sensechange> <input-table> \n\n";	
}

unless ($ARGV[0] eq "GT" || $ARGV[0] eq "CT" || $ARGV[0] eq "GC"){
	die "\nERROR : Incorrect <typeofchange>. Please use CT,GT or GC\n\n";	
}
unless (-r -e -f $ARGV[1]){
	die "\nERROR: $ARGV[1] not exists or is not readable or not properly formatted. Please check file.\n\n";	
}

$sensechange=$ARGV[0]; #sense change
$table = $ARGV[1]; # <input.table>

&parse_inputtable_to_sensechanges;

exit;