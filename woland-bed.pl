#####################################################################################################################################
## WOLAND Beta 0.2 (01-28-2016)
## woland-bed.pl
##
## WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. 
##
## For more details please read README file.
##
## Use woland-bed.pl to calculate nucleotide length using a .BED coordinate file as input and to build <chromosome_length_profile> 
## for other woland scripts.
##
## USAGE
##
## perl woland-bed.pl <file.bed> 
## 
## e.g.: perl woland-bed.pl mouse_exome_mm9.bed
##
## INPUT FILE REQUIREMENTS
## 
## <file.bed> : Coordinate bed file. IDs must be in the format: chr1, chr2, chr3 ...
##
######################################################################################################################################

#! /usr/bin/perl
use List::Util qw(sum); #module for sum of chromosome coordinates
use List::MoreUtils qw(uniq);
use strict;
use warnings;

our $REVISION = '$Revision:  $';
our $DATE =	'$Date: 2017-01-28 00:11:04 -0800 (Sat,  28 Jan 2017) $';  
our $AUTHOR =	'$Author: Tiago A. de Souza <tiagoantonio@gmail.com> $';

our (@chrbed,@pos1,@pos2, @sumlength,@uniquechr);
our $bedfile;

sub parse_bedfile{ # loading bed file
	$bedfile = $ARGV[0];
	open (BEDFILE, $bedfile);
	my @bedfilearray=<BEDFILE>;
	close (BEDFILE);

	my $rawbedline;

	foreach $rawbedline (@bedfilearray){
		my @i = split (/\t/, $rawbedline);
		chomp (@i);
		push (@chrbed, "$i[0]");
		push (@pos1, $i[1]);
		push (@pos2, $i[2]);
	}
	shift(@chrbed);
	@uniquechr = uniq @chrbed; #unique chromosome names
	
	for my $i (0..$#uniquechr){
		&calculate_length($uniquechr[$i]);
	}
}

sub calculate_length{ #absolute subtraction of pos values
	my @length;
	for my $i(0..$#chrbed){
		if($chrbed[$i] eq "$_[0]"){
			$length[$i]=abs($pos2[$i]-$pos1[$i]);
		}
		if($chrbed[$i] ne "$_[0]"){
			$length[$i]=0;
		}
	}
	push (@sumlength, sum(@length));
}

## main warning
unless ($#ARGV==0){
	die "\nERROR : Incorrect number of arguments - Usage: $0 <file.bed> \n\n";	
}
unless (-r -e -f $ARGV[0]){
    die "\nERROR: $ARGV[0] not exists or is not readable or not properly formatted. Please check file.\n\n";
}

&parse_bedfile;

open (PROFILE, ">>woland-bed-profile-$bedfile"); # printing profile file
for my $i (0..$#uniquechr){
	print PROFILE "$uniquechr[$i]\t";
	if ($sumlength[$i] ne 0){
		print PROFILE "$sumlength[$i]\n";
	}
	else{
		print PROFILE "\1\n";
	}
}
close(PROFILE);

exit;