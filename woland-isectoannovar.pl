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
use Getopt::ArgParse;

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

my $ap = Getopt::ArgParse->new_parser(
	prog => 'woland-isectoannovar.pl',
	description => 'WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing SNV data.
	Use woland-isectoannovar to perform multiple intersections between VCF files from targeted resequencing experiments to select private variants and annotate them with ANNOVAR. For more details please read README',
	epilog => 'If you used Woland in your research, we would appreciate your citation:
	de Souza TA, Defelicibus A, Menck CF',
 );

$ap->add_arg(
	'--type-of-change',
	'-c',
	required => 1,
	choices => [ 'all', 'CT', 'CG', 'CA', 'AT', 'AG', 'AC' ],
	help => 'The type of nucleotide change to filter');
$ap->add_arg(
	'--vcf-files',
	'-f',
	required => 1,
	type => 'Array',
	split => ',',
	help => 'Input VCF files separated by common (,).');
$ap->add_arg(
	'--annovar-path',
	'-a',
	required => 1,
	help => 'ANNOVAR folder path');
$ap->add_arg(
	'--htslib-path',
	'-l',
	required => 1,
	help => 'HTSLib-VCFTools folder path');
$ap->add_arg(
	'--threads',
	'-t',
	default => 30,
	help => 'Set a number for the maximum number of threads');

my $args = $ap->parse_args();

## main warning
# unless ($#ARGV>=2){
# 	die "\nERROR : Incorrect number of arguments/files - Usage: $0 <typeofchange> <file1.vcf> <file2.vcf> <file3.vcf>... \n\n";	
# }

# unless ($ARGV[0] eq "all" || $ARGV[0] eq "CT" || $ARGV[0] eq "CG" || $ARGV[0] eq "CA" || $ARGV[0] eq "AT" || $ARGV[0] eq "AG" || $ARGV[0] eq "AC"){
# 	die "\nERROR : Incorrect <typeofchange>. Please use all,CT,CG,CA,AT,AG or AC.\n\n";	
# }
# print $args->vcf_files->[2];

# for my $i (1..$#ARGV) {
# 	unless (-r -e -f $ARGV[$i]){
#     	die "\nERROR: $ARGV[$i] not exists or is not readable or not properly formatted. Please check file.\n\n";
#     }
# }

my @vcf_files = $args->vcf_files;
for my $i (0..$#vcf_files) {
	unless (-r -e -f $vcf_files[$i]){
		die "\nERROR: $vcf_files[$i] not exists or is not readable or not properly formatted. Please check file.\n\n";
	}
}

unless (-r -e -f sprintf("%s/htscmd", $args->htslib_path)){
	die sprintf("\nERROR: htscmd not found at %s.\n\n",
		$args->htslib_path);
}

unless (-r -e -f sprintf("%s/annotate_variation.pl", $args->annovar_path)){
	die sprintf("\nERROR: annotate_variation.pl not found at %s.\n\n",
		$args->annovar_path);
}

#parse nucleotide change option
# $nucleotidechangeoption=$ARGV[0];
$nucleotidechangeoption=$args->type_of_change;;

#processing vcf files
# TODO: validade if the input files are bgzip
for my $i (0..$#vcf_files){
	system ("bgzip $vcf_files[$i]"); #bgzip
}

for my $i (0..$#vcf_files){
	system ("tabix $vcf_files[$i].gz"); #index using tabix
}

for my $i (0..$#vcf_files){
	push (@variantfilelist, "$vcf_files[$i].gz"); # variant file list
}

for my $i (0..$#variantfilelist){ #perform htscmd vcfisec using a vcf file and the other all files to intersect.
	@variantfilelisttoisec=@variantfilelist;
	splice @variantfilelisttoisec,$i,1;
	system (sprintf("%s/htscmd vcfisec -C $variantfilelist[$i] @variantfilelisttoisec > $variantfilelist[$i]-exclusive.txt",
		$args->htslib_path));
	@variantfilelisttoisec=@variantfilelist;
}

for my $i (0..$#vcf_files){ #converting to annovar format
	&convert_to_annovarformat ("$vcf_files[$i].gz-exclusive.txt");
}

$pm = Parallel::ForkManager->new($args->threads); #annotation using annovar
for my $i (0..$#vcf_files){
	$pid=$pm->start and next;
	system (sprintf("perl %s/annotate_variation.pl $vcf_files[$i].gz-exclusive.txt-private-annovar.txt --geneanno --buildver hg19 ~/tools/annovar/humandb/",
		$args->annovar_path));
	$pm->finish;
}
$pm->wait_all_children;

for my $i (0..$#vcf_files){ #selecting only exonic or splicing site for exomes
	&select_onlyexonic_or_onlysplicing ("$vcf_files[$i].gz-exclusive.txt-private-annovar.txt.variant_function");
}

exit;
