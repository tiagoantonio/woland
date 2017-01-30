########################################################################################################################################
## WOLAND Beta 0.2 (01-28-2016)
## woland-isectoannovar.pl
##
## WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. 
## 
## For more details please read README file.
## 
## Use woland-isectoannovar to perform multiple intersections between VCF files from **EXOME** resequencing experiments to select private variants and annotate them with ANNOVAR.
##
## USAGE:
##
## woland-isectoannovar.pl <typeofchange> <file1.vcf> <file2.vcf> <file3.vcf>...
## 
## type of change: all,CT,CG,CA,AT,AG,AC
##
########################################################################################################################################

#! /usr/bin/perl
use IPC::System::Simple qw(system capture);
use Parallel::ForkManager;
use File::Copy;
use IO::Handle;
use strict;
use warnings;

our $REVISION = '$Revision:  $';
our $DATE =	'$Date: 2017-01-28 00:11:04 -0800 (Sat,  28 Jan 2017) $';  
our $AUTHOR =	'$Author: Tiago A. de Souza <tiagoantonio@gmail.com> $';

#variables
our ($nucleotidechangeoption, $pm, $pid);
our (@ARGV,@variantfilelist,@variantfilelisttoisec);

#subroutines
sub convert_to_annovarformat {
	my (@col1,@col2,@col3,@col4,@col5);
	my $unfvcfisecoutputline;
	my $vcfisecoutput="$_[0]";
	open (VCFISECOUTPUT, $vcfisecoutput);
	my @unfvcfisecoutput=<VCFISECOUTPUT>;
	foreach $unfvcfisecoutputline (@unfvcfisecoutput){
		my @line=split(/\t/,$unfvcfisecoutputline);
		chomp(@line);
		push (@col1, "$line[0]");
		push (@col2, "$line[1]");
		push (@col3, "$line[2]");
		push (@col4, "$line[3]");
	}
	open (ANNOVARINPUT, ">>$_[0]-private-annovar.txt");
	for my $i (0 .. $#col1){
		print ANNOVARINPUT "$col1[$i]\t$col2[$i]\t$col2[$i]\t$col3[$i]\t$col4[$i]\n";
	}
	close (VCFISECOUTPUT);
	close(ANNOVARINPUT);
}

sub select_onlyexonic_or_onlysplicing {
	my (@col1,@col2,@col3,@col4,@col5,@col6,@col7);
	my $unfvariantfunctionarrayline;
	my $variantfunctionfile="$_[0]";
	open (VARIANTFUNCTIONFILE, $variantfunctionfile);
	my @unfvariantfunctionarray=<VARIANTFUNCTIONFILE>;
	foreach $unfvariantfunctionarrayline (@unfvariantfunctionarray){
		my @linevariantfunction=split(/\t/,$unfvariantfunctionarrayline);
		chomp(@linevariantfunction);
		push (@col1, "$linevariantfunction[0]");
		push (@col2, "$linevariantfunction[1]");
		push (@col3, "$linevariantfunction[2]");
		push (@col4,  $linevariantfunction[3]);
		push (@col5,  $linevariantfunction[4]);
		push (@col6, "$linevariantfunction[5]");
		push (@col7, "$linevariantfunction[6]");
	}
	if ($nucleotidechangeoption eq "all"){
		open (ALLSNVS, ">$_[0]-towoland.txt");
		for my $i (0..$#col1){
			if ($col1[$i] eq "exonic" || $col1[$i] eq "splicing"){
				print ALLSNVS "$col1[$i]\t$col2[$i]\t$col3[$i]\t$col4[$i]\t$col5[$i]\t$col6[$i]\t$col7[$i]\n";
			}
		}
		close (ALLSNVS);
	}

	if ($nucleotidechangeoption eq "CT" || $nucleotidechangeoption eq "ct" ){
		open (ONLYCT, ">$_[0]-CT_towoland.txt");
		for my $i (0..$#col1){
			if (($col1[$i] eq "exonic" || $col1[$i] eq "splicing") && (($col6[$i] eq "C" and $col7[$i] eq "T") || ($col6[$i] eq "G" and $col7[$i] eq "A"))) {
				print ONLYCT "$col1[$i]\t$col2[$i]\t$col3[$i]\t$col4[$i]\t$col5[$i]\t$col6[$i]\t$col7[$i]\n";
			}
		}
		close (ONLYCT);
	}

	if ($nucleotidechangeoption eq "CA" || $nucleotidechangeoption eq "ca" ){
		open (ONLYCA, ">$_[0]-CA_towoland.txt");
		for my $i (0..$#col1){
			if (($col1[$i] eq "exonic" || $col1[$i] eq "splicing") && (($col6[$i] eq "C" and $col7[$i] eq "A") || ($col6[$i] eq "G" and $col7[$i] eq "T"))) {
				print ONLYCA "$col1[$i]\t$col2[$i]\t$col3[$i]\t$col4[$i]\t$col5[$i]\t$col6[$i]\t$col7[$i]\n";
			}
		}
		close (ONLYCA);
	}

	if ($nucleotidechangeoption eq "CG" || $nucleotidechangeoption eq "cg" ){
		open (ONLYCG, ">$_[0]-CG_towoland.txt");
		for my $i (0..$#col1){
			if (($col1[$i] eq "exonic" || $col1[$i] eq "splicing") && (($col6[$i] eq "G" and $col7[$i] eq "C") || ($col6[$i] eq "C" and $col7[$i] eq "G"))) {
			print ONLYCG "$col1[$i]\t$col2[$i]\t$col3[$i]\t$col4[$i]\t$col5[$i]\t$col6[$i]\t$col7[$i]\n";
			}
		}
		close (ONLYCG);
	}	

	if ($nucleotidechangeoption eq "AT" || $nucleotidechangeoption eq "at" ){
		open (ONLYAT, ">$_[0]-AT_towoland.txt");
		for my $i (0..$#col1){
			if (($col1[$i] eq "exonic" || $col1[$i] eq "splicing") && (($col6[$i] eq "A" and $col7[$i] eq "T") || ($col6[$i] eq "T" and $col7[$i] eq "A"))) {
			print ONLYAT "$col1[$i]\t$col2[$i]\t$col3[$i]\t$col4[$i]\t$col5[$i]\t$col6[$i]\t$col7[$i]\n";
			}
		}
		close (ONLYAT);
	}	

	if ($nucleotidechangeoption eq "AG" || $nucleotidechangeoption eq "ag" ){
		open (ONLYAG, ">$_[0]-AG_towoland.txt");
		for my $i (0..$#col1){
			if (($col1[$i] eq "exonic" || $col1[$i] eq "splicing") && (($col6[$i] eq "A" and $col7[$i] eq "G") || ($col6[$i] eq "T" and $col7[$i] eq "C"))) {
			print ONLYAG "$col1[$i]\t$col2[$i]\t$col3[$i]\t$col4[$i]\t$col5[$i]\t$col6[$i]\t$col7[$i]\n";
			}
		}
		close (ONLYAG);
	}	

	if ($nucleotidechangeoption eq "AC" || $nucleotidechangeoption eq "ac" ){
		open (ONLYAC, ">$_[0]-AC_towoland.txt");
		for my $i (0..$#col1){
			if (($col1[$i] eq "exonic" || $col1[$i] eq "splicing") && (($col6[$i] eq "A" and $col7[$i] eq "C") || ($col6[$i] eq "T" and $col7[$i] eq "G"))) {
			print ONLYAC "$col1[$i]\t$col2[$i]\t$col3[$i]\t$col4[$i]\t$col5[$i]\t$col6[$i]\t$col7[$i]\n";
			}
		}
		close (ONLYAC);
	}	
}

## main warning
unless ($#ARGV>=2){
	die "\nERROR : Incorrect number of arguments/files - Usage: $0 <typeofchange> <file1.vcf> <file2.vcf> <file3.vcf>... \n\n";	
}

unless ($ARGV[0] eq "all" || $ARGV[0] eq "CT" || $ARGV[0] eq "CG" || $ARGV[0] eq "CA" || $ARGV[0] eq "AT" || $ARGV[0] eq "AG" || $ARGV[0] eq "AC"){
	die "\nERROR : Incorrect <typeofchange>. Please use all,CT,CG,CA,AT,AG or AC.\n\n";	
}
for my $i (1..$#ARGV) {
	unless (-r -e -f $ARGV[$i]){
    	die "\nERROR: $ARGV[$i] not exists or is not readable or not properly formatted. Please check file.\n\n";
    }
}

#parse nucleotide change option
$nucleotidechangeoption=$ARGV[0];

#processing vcf files
for my $i (1..$#ARGV){
	system ("bgzip $ARGV[$i]"); #bgzip
}

for my $i (1..$#ARGV){
	system ("tabix $ARGV[$i].gz"); #index using tabix
}

for my $i (1..$#ARGV){
	push (@variantfilelist, "$ARGV[$i].gz"); # variant file list
}

#
for my $i (0..$#variantfilelist){ #perform htscmd vcfisec using a vcf file and the other all files to intersect.
	@variantfilelisttoisec=@variantfilelist;
	splice @variantfilelisttoisec,$i,1;
	system ("~/tools/htslib-master/htscmd vcfisec -C $variantfilelist[$i] @variantfilelisttoisec > $variantfilelist[$i]-exclusive.txt");
	@variantfilelisttoisec=@variantfilelist;

}

for my $i (1..$#ARGV){ #converting to annovar format
	&convert_to_annovarformat ("$ARGV[$i].gz-exclusive.txt");
}

$pm = Parallel::ForkManager->new(30); #annotation using annovar
for my $i (1..$#ARGV){
	$pid=$pm->start and next;
	system ("perl ~/tools/annovar/annotate_variation.pl $ARGV[$i].gz-exclusive.txt-private-annovar.txt --geneanno --buildver hg19 ~/tools/annovar/humandb/");
	$pm->finish;
}
$pm->wait_all_children;

for my $i (1..$#ARGV){ #selecting only exonic or splicing site for exomes
	&select_onlyexonic_or_onlysplicing ("$ARGV[$i].gz-exclusive.txt-private-annovar.txt.variant_function");
}

exit;
