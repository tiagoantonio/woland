#########################################################################################################################
## WOLAND BETA 0.1 (15-02-2016)
##
## WOLAND is a software package based on Perl and R for calculation of general mutation metrics, identification and
## comparison of predicted hotspots across biological samples. WOLAND uses Single Nucleotide Polymorphisms (SNPs) data
## from Next Generation Sequencing (NGS) pipelines as main input. Please observe file requirements.
##
##
## USAGE
##
## perl woland.pl <annnovar_variant_file> <chromosome_length_profile> <hotspot_window_length>
## 
## e.g.: perl woland.pl sample1_exome.tab mouse_chromosome.tab 1000
##
## INPUT FILE REQUIREMENTS
## 
## <annnovar_variant_file> : ANNOVAR variant file output.
##
## <chromosome_length_profile> : Tab-delimited file without header. First colum is the chromosome name; chr1. Second column is the 
## chromosome length; 45000000. WARNING: This file must contain, in the first column: chr1-chr21, chrX, chrY and chr M, in a total of
## 24 columns. If any chromosome is absent, please use 1 as chromosome length. 
##
## <hotspot_window_length> : Any natural X number correspondent to the nucleotide hotspot window length; X nucleotides distant in each
## way of SNP position in <tabular_snp_file>.
##
## WARNING: You must save a genome FASTA file as "genome.fa" for context sequence extraction. This genome file must be chromosome divided;
## >chr1, >chr2 ... as IDs.
##
######################################################################################################################### 



#! /usr/bin/perl

use Bio::DB::Fasta; # bioperl module for the extraction of sequences.
use Cwd;
use warnings;
use strict;

# variables

 
my $inputRawSNV; my $chrLengthProfile; my $hotspotWindowLength; my $contextSequenceLength;
my $datestring; 
my $rawLine; our $i; my $i2; my $i3;

our $ATTA; our $AGTC; our $ACTG; our $CGGC; our $CTGA; our $CAGT; our $refalt; my $SOMA; our $AVGATTA; our $AVGAGTC; our $AVGACTG; our $AVGCGGC; our $AVGCTGA; our $AVGCAGT;

my $transition; my $transversion;

my $chr1; my $chr2; my $chr3; my $chr4; my $chr5; my $chr6; my $chr7; my $chr8; my $chr9; my $chr10; my $chr11; 
my $chr12; my $chr13; my $chr14; my $chr15; my $chr16;my $chr17; my $chr18; my $chr19; my $chr20; my $chr21; my $chr22; 
my $chrX;my $chrY; my $chrM;

my $chrNorm1; my $chrNorm2; my $chrNorm3; my $chrNorm4; my $chrNorm5; my $chrNorm6; my $chrNorm7; my $chrNorm8; my $chrNorm9; my $chrNorm10; my $chrNorm11; 
my $chrNorm12; my $chrNorm13; my $chrNorm14; my $chrNorm15; my $chrNorm16; my $chrNorm17; my $chrNorm18; my $chrNorm19; my $chrNorm20; my $chrNorm21; my $chrNorm22;
 my $chrNormX; my $chrNormY; my $chrNormM;

my $a1; my $a2; my $a3; my $t1;

my $hotspotSNV; my $hotspotDownstreamSNV; my $hotspotUpstreamSNV; my $hsNum; 

our $config; our $count; our $notValid;

our @rawFile; our @config; our @i; our @chr; our @pos; our @alt; our @chrpos; our %posref; our %posalt; 
our @refalt; our @ref; our @chrLength;

my $hsChr1; our @hsChr1; our @hs_counts_chr1;
my $hsChr2; our @hsChr2; our @hs_counts_chr2;
my $hsChr3; our @hsChr3; our @hs_counts_chr3;
my $hsChr4; our @hsChr4; our @hs_counts_chr4;
my $hsChr5; our @hsChr5; our @hs_counts_chr5;
my $hsChr6; our @hsChr6; our @hs_counts_chr6;
my $hsChr7; our @hsChr7; our @hs_counts_chr7;
my $hsChr8; our @hsChr8; our @hs_counts_chr8;
my $hsChr9; our @hsChr9; our @hs_counts_chr9;
my $hsChr10; our @hsChr10; our @hs_counts_chr10;
my $hsChr11; our @hsChr11; our @hs_counts_chr11;
my $hsChr12; our @hsChr12; our @hs_counts_chr12;
my $hsChr13; our @hsChr13; our @hs_counts_chr13;
my $hsChr14; our @hsChr14; our @hs_counts_chr14;
my $hsChr15; our @hsChr15; our @hs_counts_chr15;
my $hsChr16; our @hsChr16; our @hs_counts_chr16;
my $hsChr17; our @hsChr17; our @hs_counts_chr17;
my $hsChr18; our @hsChr18; our @hs_counts_chr18;
my $hsChr19; our @hsChr19; our @hs_counts_chr19;
my $hsChr20; our @hsChr20; our @hs_counts_chr20;
my $hsChr21; our @hsChr21; our @hs_counts_chr21;
my $hsChr22; our @hsChr22; our @hs_counts_chr22;
my $hsChrX; our @hsChrX; our @hs_counts_chrX;
my $hsChrY; our @hsChrY; our @hs_counts_chrY;
my $hsChrM; our @hsChrM; our @hs_counts_chrM;

our @fastaContext; our $fastaContext; our @id; our @targetSequence; my $targetSequence;
our @SN1; our @DNApoln; our @oxoG; our @UVlambda; our @UVsolar; our @sixfour; our @enu; our $SNPnumber; my $mf1; my $mf2;
our @SN1counts; our @DNApolncounts; our @oxoGcounts; our @UVlambdacounts; our @UVsolarcounts; our @sixfourcounts; our @enucounts;
our $SN1; our $DNApoln; our $oxoG; our $UVlambda; our $UVsolar; our $sixfour; our $enu;
our $SN1counts; our $DNApolncounts; our $oxoGcounts; our $UVlambdacounts; our $UVsolarcounts; our $sixfourcounts; our $enucounts;
our $normSN1; our $normDNApoln; our $normoxoG; our $normUVlambda; our $normUVsolar; our $normsixfour; our $normenu;

## anno motif search variables

our (@geneClass, $geneClass, @geneName, $geneName);
our (@geneClasschr1, @geneClasschr2, @geneClasschr3, @geneClasschr4, @geneClasschr5, @geneClasschr6, @geneClasschr7, @geneClasschr8, @geneClasschr9, @geneClasschr10, @geneClasschr11, @geneClasschr12, @geneClasschr13, @geneClasschr14, @geneClasschr15, @geneClasschr16, @geneClasschr17, @geneClasschr18, @geneClasschr19, @geneClasschr20, @geneClasschr21, @geneClasschr22, @geneClasschrX, @geneClasschrY, @geneClasschrM);      
our (@geneNamechr1, @geneNamechr2, @geneNamechr3, @geneNamechr4, @geneNamechr5, @geneNamechr6, @geneNamechr7, @geneNamechr8, @geneNamechr9, @geneNamechr10, @geneNamechr11, @geneNamechr12, @geneNamechr13, @geneNamechr14, @geneNamechr15, @geneNamechr16, @geneNamechr17, @geneNamechr18, @geneNamechr19, @geneNamechr20, @geneNamechr21, @geneNamechr22, @geneNamechrX, @geneNamechrY, @geneNamechrM);
our (@chrStart, @chrEnd);
our (@fastaContextMS, $fastaContextMS, @fastaContextMS1, $fastacontextMS1, @geneClassMS, $geneClassMS, @geneClassMS1, $geneClassMS1, @geneNameMS, $geneNameMS, @geneNameMS1, $geneNameMS1);
our (@chrRaw, $chrRaw, $chrSt, $refSeqRaw);
our (@chrRefSeq, @irefSeq, );
our ($query_coord, $query_chr);
our (@txtendRefSeq, @txstartRefSeq);
our (@i2, @i3, $i4, $i5, @i7, $i7, $id);
our (@geneClass1, @geneName1, @coord, @chrSt, @refSeq, $refSeqline, @strandRefSeq);
our (@SN1plus, @SN1minus, @DNApolnplus, @DNApolnminus, @oxoGplus, @oxoGminus, @UVlambdaplus, @UVlambdaminus, @UVsolarplus, @UVsolarminus, @sixfourplus, @sixfourminus, @enuplus, @enuminus);
our ($strand_count, $strand_plus_count, $strand_transcript, $strand_score, $strand_value);
our ($db, $db1);
our $readContextSeqAnno;

# main Warning

unless (@ARGV){
	die "\nERROR : Usage: $0 <tabular_snp_file> <chromosome_length_profile> <hotspot_window_length> \n";	
}

# loading files and parameters

$inputRawSNV = $ARGV[0]; #<tabular_snp_file>
open (inputRawSNV, $inputRawSNV);
@rawFile=<inputRawSNV>;
unless (@rawFile){
	die "\nERROR : Could not load <tabular_snp_file>\n";
}

$chrLengthProfile = $ARGV[1]; #<chromosome_length_profile>
open (chrLenghtProfile, $chrLengthProfile);
@config=<chrLenghtProfile>;
unless (@config){
	die "\nERROR : Could not load <chromosome_length_profile>\n";
}

$hotspotWindowLength = $ARGV[2]; #<hotspot_window_length>
unless ($hotspotWindowLength){
	die "\nERROR : Please specify a natural number >0 for <hotspot_window_length>\n";
}

$contextSequenceLength=3; #<context_sequence_length> #default=3nt downstream & 3nt upstream

# loading of outputs
mkdir("results-$inputRawSNV", 0755) || die "Cannot mkdir results-$inputRawSNV";

open (BASECHANGE, ">>results-$inputRawSNV/WOLAND-basechange-$inputRawSNV");

open (MUTFREQ, ">>results-$inputRawSNV/WOLAND-mutfreq-$inputRawSNV");

open (CONTEXTSEQ, ">>results-$inputRawSNV/WOLAND-contextsequences-$inputRawSNV");

open (CONTEXTSEQANNO, ">>results-$inputRawSNV/WOLAND-contextsequencesanno-$inputRawSNV");

open (HOTSPOTS, ">>results-$inputRawSNV/WOLAND-hotspots-$inputRawSNV");

open (MOTIFS, ">>results-$inputRawSNV/WOLAND-motifs-$inputRawSNV");

open (NMOTIFS, ">>results-$inputRawSNV/WOLAND-norm_motifs-$inputRawSNV");

open (LOG, ">>results-$inputRawSNV/WOLAND-log-$inputRawSNV");

# Start Screen & LOG file

$datestring = localtime();

print LOG "WOLAND BETA 0.9 - 02-02-2016\n";
print LOG "Analysis started at $datestring\n";
print LOG "Tabular SNP File              : $inputRawSNV\n";
print LOG "Chromosome Length Profile File: $chrLengthProfile\n";
print LOG "Hotspot Window Length         : $hotspotWindowLength bases flanking SNP position in reference genome\n";
print LOG "Context Sequence Length       : $contextSequenceLength bases flanking SNP position in reference genome\n";

print "\n\nLoading SNP file...\n";

print LOG "\nOutputs:\n";
print LOG "Mutation Statistics:                       WOLAND-basechange-$inputRawSNV\n";
print LOG "Extracted Unannotated Context Sequences:   WOLAND-contextsequences-$inputRawSNV\n";
print LOG "Extracted Context Sequences:               WOLAND-contextsequencesanno-$inputRawSNV\n";
print LOG "Hotspots:                                  WOLAND-hotspots-$inputRawSNV\n";
print LOG "Mutational Motifs:                         WOLAND-motifs-$inputRawSNV\n";
print LOG "Mutational Motifs Normalized:              WOLAND-norm_motifs-$inputRawSNV\n";
print LOG "Strand Bias of SN1 Motifs:                 WOLAND-bias_uv-lambda-$inputRawSNV\n";
print LOG "Strand Bias of DNApoln Motifs:             WOLAND-bias_DNApoln-$inputRawSNV\n";
print LOG "Strand Bias of 8-oxoG Motifs:              WOLAND-bias_oxoG-$inputRawSNV\n";
print LOG "Strand Bias of UV-lambda Motifs:           WOLAND-bias_uv-lambda-$inputRawSNV\n";
print LOG "Strand Bias of UV Solar Motifs:             WOLAND-bias_UVsolar-$inputRawSNV\n";
print LOG "Strand Bias of 6-4 Motifs:                 WOLAND-bias_sixfour-$inputRawSNV\n";
print LOG "Strand Bias of ENU Motifs:                 WOLAND-bias_enu-$inputRawSNV\n";

# Conversion of each column in a dedicated array

foreach $rawLine (@rawFile)
	{
	@i = split (/\t/, $rawLine);
	chomp (@i);
	push (@geneClass, "$i[0]");
	push (@geneName, "$i[1]");
	push (@chrStart, "$i[2]");
	push (@chrEnd, "$i[3]");
	push (@pos, $i[4]);
	push (@ref, "$i[5]");
	push (@alt, "$i[6]");
	}
	
# Array for chr/position of each SNP - @chrpos

	for my $i3 (0 .. $#chrStart){
		push @chrpos, "$chrStart[$i3]_$pos[$i3]";
	}
	
# Hashe posref & posalt for information of each position - ALT e REF alelles

	for my $i3 (0 .. $#chrpos){
		$posref{"$chrpos[$i3]_$i3"} .= "$ref[$i3]";
	}
	for my $i3 (0 .. $#chrpos){
		$posalt{"$chrpos[$i3]_$i3"} .= "$alt[$i3]";
	}
	
#
# 

print "\nCalculating general mutation statistics and saving basechange-$inputRawSNV ...\n";

# Array for a single entry for each SNP containing a single REFALT string value. Considering  A->T # & T->A; A->G & T->C; A->C & T->G; C->G & G->C; C->T & G->A; C->A & G->T.

	foreach my $key (sort keys %posref){
  	push @refalt,"$posref{$key}$posalt{$key}";
	}


# Counting of changes.

$ATTA=0;
$AGTC=0;
$ACTG=0;
$CGGC=0;
$CTGA=0;
$CAGT=0;
	foreach $refalt (@refalt){
	if ($refalt eq "AT" || $refalt eq "TA"){
		$ATTA++;}
	if ($refalt eq "AG" || $refalt eq "TC"){
		$AGTC++;}
	if ($refalt eq "AC" || $refalt eq "TG"){
		$ACTG++;}
	if ($refalt eq "CG" || $refalt eq "GC"){
		$CGGC++;}
	if ($refalt eq "CT" || $refalt eq "GA"){
		$CTGA++;}
	if ($refalt eq "CA" || $refalt eq "GT"){
		$CAGT++;}
}

# Frequency of changes considering total amount of changes.

$SOMA=$ATTA+$AGTC+$ACTG+$CGGC+$CTGA+$CAGT;

foreach $refalt(@refalt){
$count++;}
print LOG "\nTotal SNPs valid=$count\n";
$notValid=$count-$SOMA;
print LOG "Total not valid SNPs=$notValid\n";

$AVGATTA=$ATTA/$SOMA;
$AVGAGTC=$AGTC/$SOMA;
$AVGACTG=$ACTG/$SOMA;
$AVGCGGC=$CGGC/$SOMA;
$AVGCTGA=$CTGA/$SOMA;
$AVGCAGT=$CAGT/$SOMA;


# Transitions & Transversion

$transition=($AGTC+$CTGA)/$SOMA;
$transversion=($ATTA+$ACTG+$CGGC+$CAGT)/$SOMA;



# Frequency per chromosome target

$chr1=0;$chr12=0;
$chr2=0;$chr13=0;
$chr3=0;$chr14=0;
$chr4=0;$chr15=0;
$chr5=0;$chr16=0;
$chr6=0;$chr17=0;
$chr7=0;$chr18=0;
$chr8=0;$chr19=0;
$chr9=0;$chr20=0;
$chr10=0;$chr21=0;
$chr11=0;$chr22=0;
$chrX=0;$chrY=0;
$chrM=0;
	for my $i6 (0 .. $#chrStart){
			if ($chrStart[$i6] eq "chr1"){
			++$chr1;}}
	for my $i7 (0 .. $#chrStart){
			if ($chrStart[$i7] eq "chr2"){
			++$chr2;}}
	for my $i8 (0 .. $#chrStart){
			if ($chrStart[$i8] eq "chr3"){
			++$chr3;}}
	for my $i9 (0 .. $#chrStart){
			if ($chrStart[$i9] eq "chr4"){
			++$chr4;}}
	for my $i10 (0 .. $#chrStart){
			if ($chrStart[$i10] eq "chr5"){
			++$chr5;}}
	for my $i11 (0 .. $#chrStart){
			if ($chrStart[$i11] eq "chr6"){
			++$chr6;}}
	for my $i12 (0 .. $#chrStart){
			if ($chrStart[$i12] eq "chr7"){
			++$chr7;}}
	for my $i13 (0 .. $#chrStart){
			if ($chrStart[$i13] eq "chr8"){
			++$chr8;}}
	for my $i14 (0 .. $#chrStart){
			if ($chrStart[$i14] eq "chr9"){
			++$chr9;}}
	for my $i15 (0 .. $#chrStart){
			if ($chrStart[$i15] eq "chr10"){
			++$chr10;}}
	for my $i16 (0 .. $#chrStart){
			if ($chrStart[$i16] eq "chr11"){
			++$chr11;}}
	for my $i17 (0 .. $#chrStart){
			if ($chrStart[$i17] eq "chr12"){
			++$chr12;}}
	for my $i18 (0 .. $#chrStart){
			if ($chrStart[$i18] eq "chr13"){
			++$chr13;}}
	for my $i19 (0 .. $#chrStart){
			if ($chrStart[$i19] eq "chr14"){
			++$chr14;}}
	for my $i20 (0 .. $#chrStart){
			if ($chrStart[$i20] eq "chr15"){
			++$chr15;}}
	for my $i21 (0 .. $#chrStart){
			if ($chrStart[$i21] eq "chr16"){
			++$chr16;}}
	for my $i22 (0 .. $#chrStart){
			if ($chrStart[$i22] eq "chr17"){
			++$chr17;}}
	for my $i23 (0 .. $#chrStart){
			if ($chrStart[$i23] eq "chr18"){
			++$chr18;}}
	for my $i24 (0 .. $#chrStart){
			if ($chrStart[$i24] eq "chr19"){
			++$chr19;}}
	for my $i25 (0 .. $#chrStart){
			if ($chrStart[$i25] eq "chr20"){
			++$chr20;}}
	for my $i26 (0 .. $#chrStart){
			if ($chrStart[$i26] eq "chr21"){
			++$chr21;}}
	for my $i27 (0 .. $#chrStart){
			if ($chrStart[$i27] eq "chr22"){
			++$chr22;}}
	for my $i28 (0 .. $#chrStart){
			if ($chrStart[$i28] eq "chrX"){
			++$chrX;}}
	for my $i29 (0 .. $#chrStart){
			if ($chrStart[$i29] eq "chrY"){
			++$chrY;}}
	for my $i30 (0 .. $#chrStart){
			if ($chrStart[$i30] eq "chrM"){
			++$chrM;}}


# Frequency of mutation considering chromosome length as present in chromosome_profile file.

print "\nCalculating mutation frequency and saving basechange-$inputRawSNV ...\n";

foreach $config (@config)
	{
	@i = split (/\t/, $config);
	chomp (@i);
	push (@chrLength, "$i[1]");}

$chrNorm1=$chr1/$chrLength[0];
$chrNorm2=$chr2/$chrLength[1];
$chrNorm3=$chr3/$chrLength[2];
$chrNorm4=$chr4/$chrLength[3];
$chrNorm5=$chr5/$chrLength[4];
$chrNorm6=$chr6/$chrLength[5];
$chrNorm7=$chr7/$chrLength[6];
$chrNorm8=$chr8/$chrLength[7];
$chrNorm9=$chr9/$chrLength[8];
$chrNorm10=$chr10/$chrLength[9];
$chrNorm11=$chr11/$chrLength[10];
$chrNorm12=$chr12/$chrLength[11];
$chrNorm13=$chr13/$chrLength[12];
$chrNorm14=$chr14/$chrLength[13];
$chrNorm15=$chr15/$chrLength[14];
$chrNorm16=$chr16/$chrLength[15];
$chrNorm17=$chr17/$chrLength[16];
$chrNorm18=$chr18/$chrLength[17];
$chrNorm19=$chr19/$chrLength[18];
if($chrLength[19] == 0){$chrNorm20="Cromossomo inexistente";
	}else{
	$chrNorm20=$chr20/$chrLength[19];}
if($chrLength[20] == 0){$chrNorm21="Cromossomo inexistente";
	}else{
	$chrNorm21=$chr21/$chrLength[20];}
if($chrLength[21] == 0){$chrNorm22="Cromossomo inexistente";
	}else{
	$chrNorm22=$chr22/$chrLength[21];}
$chrNormX=$chrX/$chrLength[22];
$chrNormY=$chrY/$chrLength[23];
$chrNormM=$chrM/$chrLength[24];


##### BASECHANGE Printing Format ####
print BASECHANGE "$inputRawSNV\tA>T\tA>G\tA>C\tC>G\tC>T\tC>A\n";
print BASECHANGE "Changes\t$ATTA\t$AGTC\t$ACTG\t$CGGC\t$CTGA\t$CAGT\n";
print BASECHANGE "Frequency\t$AVGATTA\t$AVGAGTC\t$AVGACTG\t$AVGCGGC\t$AVGCTGA\t$AVGCAGT\n";
print BASECHANGE "Type\tTransversion\tTransition\tTransversion\tTransversion\tTransition\tTransversion\n";
print BASECHANGE "\n";

print MUTFREQ "Chromosome\tTotalNumver\tMutPerBaseProfile\n";
print MUTFREQ "chr1\t$chr1\t$chrNorm1\n";
print MUTFREQ "chr2\t$chr2\t$chrNorm2\n";
print MUTFREQ "chr3\t$chr3\t$chrNorm3\n";
print MUTFREQ "chr4\t$chr4\t$chrNorm4\n";
print MUTFREQ "chr5\t$chr5\t$chrNorm5\n";
print MUTFREQ "chr6\t$chr6\t$chrNorm6\n";
print MUTFREQ "chr7\t$chr7\t$chrNorm7\n";
print MUTFREQ "chr8\t$chr8\t$chrNorm8\n";
print MUTFREQ "chr9\t$chr9\t$chrNorm9\n";
print MUTFREQ "chr10\t$chr10\t$chrNorm10\n";
print MUTFREQ "chr11\t$chr11\t$chrNorm11\n";
print MUTFREQ "chr12\t$chr12\t$chrNorm12\n";
print MUTFREQ "chr13\t$chr13\t$chrNorm13\n";
print MUTFREQ "chr14\t$chr14\t$chrNorm14\n";
print MUTFREQ "chr15\t$chr15\t$chrNorm15\n";
print MUTFREQ "chr16\t$chr16\t$chrNorm16\n";
print MUTFREQ "chr17\t$chr17\t$chrNorm17\n";
print MUTFREQ "chr18\t$chr18\t$chrNorm18\n";
print MUTFREQ "chr19\t$chr19\t$chrNorm19\n";
print MUTFREQ "chr20\t$chr20\t$chrNorm20\n";
print MUTFREQ "chr21\t$chr21\t$chrNorm21\n";
print MUTFREQ "chr22\t$chr22\t$chrNorm22\n";
print MUTFREQ "chrX\t$chrX\t$chrNormX\n";
print MUTFREQ "chrY\t$chrY\t$chrNormY\n";
print MUTFREQ "chrM\t$chrM\t$chrNormM\n";

############################################# HOT SPOT #####################################

print "\nCalculating hotspots and saving hotspots-$inputRawSNV.txt ...\n";

# One array of each chromosome.

for my $hs_a (0..$#chrStart){
	if ($chrStart[$hs_a] eq "chr1"){
		push @hsChr1, $pos[$hs_a];
		push @geneClasschr1, $geneClass[$hs_a];
		push @geneNamechr1, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr2"){
		push @hsChr2, $pos[$hs_a];
		push @geneClasschr2, $geneClass[$hs_a];
		push @geneNamechr2, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr3"){
		push @hsChr3, $pos[$hs_a];
		push @geneClasschr3, $geneClass[$hs_a];
		push @geneNamechr3, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr4"){
		push @hsChr4, $pos[$hs_a];
		push @geneClasschr4, $geneClass[$hs_a];
		push @geneNamechr4, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr5"){
		push @hsChr5, $pos[$hs_a];
		push @geneClasschr5, $geneClass[$hs_a];
		push @geneNamechr5, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr6"){
		push @hsChr6, $pos[$hs_a];
		push @geneClasschr6, $geneClass[$hs_a];
		push @geneNamechr6, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr7"){
		push @hsChr7, $pos[$hs_a];
		push @geneClasschr7, $geneClass[$hs_a];
		push @geneNamechr7, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr8"){
		push @hsChr8, $pos[$hs_a];
		push @geneClasschr8, $geneClass[$hs_a];
		push @geneNamechr8, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr9"){
		push @hsChr9, $pos[$hs_a];
		push @geneClasschr9, $geneClass[$hs_a];
		push @geneNamechr9, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr10"){
		push @hsChr10, $pos[$hs_a];
		push @geneClasschr10, $geneClass[$hs_a];
		push @geneNamechr10, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr11"){
		push @hsChr11, $pos[$hs_a];
		push @geneClasschr11, $geneClass[$hs_a];
		push @geneNamechr11, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr12"){
		push @hsChr12, $pos[$hs_a];
		push @geneClasschr12, $geneClass[$hs_a];
		push @geneNamechr12, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr13"){
		push @hsChr13, $pos[$hs_a];
		push @geneClasschr13, $geneClass[$hs_a];
		push @geneNamechr13, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr14"){
		push @hsChr14, $pos[$hs_a];
		push @geneClasschr14, $geneClass[$hs_a];
		push @geneNamechr14, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr15"){
		push @hsChr15, $pos[$hs_a];
		push @geneClasschr15, $geneClass[$hs_a];
		push @geneNamechr15, $geneName[$hs_a];}
		if ($chrStart[$hs_a] eq "chr16"){
		push @hsChr16, $pos[$hs_a];
		push @geneClasschr16, $geneClass[$hs_a];
		push @geneNamechr16, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr17"){
		push @hsChr17, $pos[$hs_a];
		push @geneClasschr17, $geneClass[$hs_a];
		push @geneNamechr17, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr18"){
		push @hsChr18, $pos[$hs_a];
		push @geneClasschr18, $geneClass[$hs_a];
		push @geneNamechr18, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr19"){
		push @hsChr19, $pos[$hs_a];
		push @geneClasschr19, $geneClass[$hs_a];
		push @geneNamechr19, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr20"){
		push @hsChr20, $pos[$hs_a];
		push @geneClasschr20, $geneClass[$hs_a];
		push @geneNamechr20, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr21"){
		push @hsChr21, $pos[$hs_a];
		push @geneClasschr21, $geneClass[$hs_a];
		push @geneNamechr21, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chr22"){
		push @hsChr22, $pos[$hs_a];
		push @geneClasschr22, $geneClass[$hs_a];
		push @geneNamechr22, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chrX"){
		push @hsChrX, $pos[$hs_a];
		push @geneClasschrX, $geneClass[$hs_a];
		push @geneNamechrX, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chrY"){
		push @hsChrY, $pos[$hs_a];
		push @geneClasschrY, $geneClass[$hs_a];
		push @geneNamechrY, $geneName[$hs_a];}
	if ($chrStart[$hs_a] eq "chrM"){
		push @hsChrM, $pos[$hs_a];
		push @geneClasschrM, $geneClass[$hs_a];
		push @geneNamechrM, $geneName[$hs_a];}
}

# Counting of Hotspots using one array for each chromosome.

for my $hotspotSNV (0..$#hsChr1){
	my $hotspotDownstreamSNV=$hsChr1[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr1[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr1(@hsChr1) {
		if ($hsChr1>=$hotspotDownstreamSNV && $hsChr1 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr1, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr2){
	my $hotspotDownstreamSNV=$hsChr2[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr2[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr2(@hsChr2) {
		if ($hsChr2>=$hotspotDownstreamSNV && $hsChr2 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr2, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr3){
	my $hotspotDownstreamSNV=$hsChr3[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr3[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr3(@hsChr3) {
		if ($hsChr3>=$hotspotDownstreamSNV && $hsChr3 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr3, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr4){
	my $hotspotDownstreamSNV=$hsChr4[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr4[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr4(@hsChr4) {
		if ($hsChr4>=$hotspotDownstreamSNV && $hsChr4 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr4, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr5){
	my $hotspotDownstreamSNV=$hsChr5[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr5[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr5(@hsChr5) {
		if ($hsChr5>=$hotspotDownstreamSNV && $hsChr5 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr5, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr6){
	my $hotspotDownstreamSNV=$hsChr6[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr6[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr6(@hsChr6) {
		if ($hsChr6>=$hotspotDownstreamSNV && $hsChr6 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr6, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr7){
	my $hotspotDownstreamSNV=$hsChr7[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr7[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr7(@hsChr7) {
		if ($hsChr7>=$hotspotDownstreamSNV && $hsChr7 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr7, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr8){
	my $hotspotDownstreamSNV=$hsChr8[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr8[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr8(@hsChr8) {
		if ($hsChr8>=$hotspotDownstreamSNV && $hsChr8 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr8, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr9){
	my $hotspotDownstreamSNV=$hsChr9[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr9[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr9(@hsChr9) {
		if ($hsChr9>=$hotspotDownstreamSNV && $hsChr9 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr9, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr10){
	my $hotspotDownstreamSNV=$hsChr10[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr10[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr10(@hsChr10) {
		if ($hsChr10>=$hotspotDownstreamSNV && $hsChr10 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr10, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr11){
	my $hotspotDownstreamSNV=$hsChr11[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr11[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr11(@hsChr11) {
		if ($hsChr11>=$hotspotDownstreamSNV && $hsChr11 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr11, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr12){
	my $hotspotDownstreamSNV=$hsChr12[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr12[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr12(@hsChr12) {
		if ($hsChr12>=$hotspotDownstreamSNV && $hsChr12 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr12, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr13){
	my $hotspotDownstreamSNV=$hsChr13[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr13[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr13(@hsChr13) {
		if ($hsChr13>=$hotspotDownstreamSNV && $hsChr13 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr13, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr14){
	my $hotspotDownstreamSNV=$hsChr14[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr14[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr14(@hsChr14) {
		if ($hsChr14>=$hotspotDownstreamSNV && $hsChr14 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr14, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr15){
	my $hotspotDownstreamSNV=$hsChr15[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr15[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr15(@hsChr15) {
		if ($hsChr15>=$hotspotDownstreamSNV && $hsChr15 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr15, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr16){
	my $hotspotDownstreamSNV=$hsChr16[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr16[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr16(@hsChr16) {
		if ($hsChr16>=$hotspotDownstreamSNV && $hsChr16 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr16, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr17){
	my $hotspotDownstreamSNV=$hsChr17[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr17[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr17(@hsChr17) {
		if ($hsChr17>=$hotspotDownstreamSNV && $hsChr17 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr17, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr18){
	my $hotspotDownstreamSNV=$hsChr18[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr18[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr18(@hsChr18) {
		if ($hsChr18>=$hotspotDownstreamSNV && $hsChr18 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr18, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr19){
	my $hotspotDownstreamSNV=$hsChr19[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr19[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr19(@hsChr19) {
		if ($hsChr19>=$hotspotDownstreamSNV && $hsChr19 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr19, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr20){
	my $hotspotDownstreamSNV=$hsChr20[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr20[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr20(@hsChr20) {
		if ($hsChr20>=$hotspotDownstreamSNV && $hsChr20 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr20, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChr21){
	my $hotspotDownstreamSNV=$hsChr21[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr21[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr21(@hsChr21) {
		if ($hsChr21>=$hotspotDownstreamSNV && $hsChr21 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr21, $hsNum;
	$hsNum=0;

}

for my $hotspotSNV (0..$#hsChr22){
	my $hotspotDownstreamSNV=$hsChr22[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChr22[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChr22(@hsChr22) {
		if ($hsChr22>=$hotspotDownstreamSNV && $hsChr22 <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chr22, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChrX){
	my $hotspotDownstreamSNV=$hsChrX[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChrX[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChrX(@hsChrX) {
		if ($hsChrX>=$hotspotDownstreamSNV && $hsChrX <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chrX, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChrY){
	my $hotspotDownstreamSNV=$hsChrY[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChrY[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChrY(@hsChrY) {
		if ($hsChrY>=$hotspotDownstreamSNV && $hsChrY <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chrY, $hsNum;
	$hsNum=0;
}

for my $hotspotSNV (0..$#hsChrM){
	my $hotspotDownstreamSNV=$hsChrM[$hotspotSNV]-$hotspotWindowLength;
	my $hotspotUpstreamSNV=$hsChrM[$hotspotSNV]+$hotspotWindowLength;
	foreach $hsChrM(@hsChrM) {
		if ($hsChrM>=$hotspotDownstreamSNV && $hsChrM <=$hotspotUpstreamSNV){++$hsNum;}
	}
	push @hs_counts_chrM, $hsNum;
	$hsNum=0;
}

##### HOTSPOTS Printing Format ####
print HOTSPOTS "geneClass\tGeneName\tCHR\tBP\tHotspotCount\n";
for my $t1 (0 .. $#hsChr1){
print HOTSPOTS "$geneClasschr1[$t1]\t$geneNamechr1[$t1]\t1\t$hsChr1[$t1]\t$hs_counts_chr1[$t1]\n";}
for my $t1 (0 .. $#hsChr2){
print HOTSPOTS "$geneClasschr2[$t1]\t$geneNamechr2[$t1]\t2\t$hsChr2[$t1]\t$hs_counts_chr2[$t1]\n";}
for my $t1 (0 .. $#hsChr3){
print HOTSPOTS "$geneClasschr3[$t1]\t$geneNamechr3[$t1]\t3\t$hsChr3[$t1]\t$hs_counts_chr3[$t1]\n";}
for my $t1 (0 .. $#hsChr4){
print HOTSPOTS "$geneClasschr4[$t1]\t$geneNamechr4[$t1]\t4\t$hsChr4[$t1]\t$hs_counts_chr4[$t1]\n";}
for my $t1 (0 .. $#hsChr5){
print HOTSPOTS "$geneClasschr5[$t1]\t$geneNamechr5[$t1]\t5\t$hsChr5[$t1]\t$hs_counts_chr5[$t1]\n";}
for my $t1 (0 .. $#hsChr6){
print HOTSPOTS "$geneClasschr6[$t1]\t$geneNamechr6[$t1]\t6\t$hsChr6[$t1]\t$hs_counts_chr6[$t1]\n";}
for my $t1 (0 .. $#hsChr7){
print HOTSPOTS "$geneClasschr7[$t1]\t$geneNamechr7[$t1]\t7\t$hsChr7[$t1]\t$hs_counts_chr7[$t1]\n";}
for my $t1 (0 .. $#hsChr8){
print HOTSPOTS "$geneClasschr8[$t1]\t$geneNamechr8[$t1]\t8\t$hsChr8[$t1]\t$hs_counts_chr8[$t1]\n";}
for my $t1 (0 .. $#hsChr9){
print HOTSPOTS "$geneClasschr9[$t1]\t$geneNamechr9[$t1]\t9\t$hsChr9[$t1]\t$hs_counts_chr9[$t1]\n";}
for my $t1 (0 .. $#hsChr10){
print HOTSPOTS "$geneClasschr10[$t1]\t$geneNamechr10[$t1]\t10\t$hsChr10[$t1]\t$hs_counts_chr10[$t1]\n";}
for my $t1 (0 .. $#hsChr11){
print HOTSPOTS "$geneClasschr11[$t1]\t$geneNamechr11[$t1]\t11\t$hsChr11[$t1]\t$hs_counts_chr11[$t1]\n";}
for my $t1 (0 .. $#hsChr12){
print HOTSPOTS "$geneClasschr12[$t1]\t$geneNamechr12[$t1]\t12\t$hsChr12[$t1]\t$hs_counts_chr12[$t1]\n";}
for my $t1 (0 .. $#hsChr13){
print HOTSPOTS "$geneClasschr13[$t1]\t$geneNamechr13[$t1]\t13\t$hsChr13[$t1]\t$hs_counts_chr13[$t1]\n";}
for my $t1 (0 .. $#hsChr14){
print HOTSPOTS "$geneClasschr14[$t1]\t$geneNamechr14[$t1]\t14\t$hsChr14[$t1]\t$hs_counts_chr14[$t1]\n";}
for my $t1 (0 .. $#hsChr15){
print HOTSPOTS "$geneClasschr15[$t1]\t$geneNamechr15[$t1]\t15\t$hsChr15[$t1]\t$hs_counts_chr15[$t1]\n";}
for my $t1 (0 .. $#hsChr16){
print HOTSPOTS "$geneClasschr16[$t1]\t$geneNamechr16[$t1]\t16\t$hsChr16[$t1]\t$hs_counts_chr16[$t1]\n";}
for my $t1 (0 .. $#hsChr17){
print HOTSPOTS "$geneClasschr17[$t1]\t$geneNamechr17[$t1]\t17\t$hsChr17[$t1]\t$hs_counts_chr17[$t1]\n";}
for my $t1 (0 .. $#hsChr18){
print HOTSPOTS "$geneClasschr18[$t1]\t$geneNamechr18[$t1]\t18\t$hsChr18[$t1]\t$hs_counts_chr18[$t1]\n";}
for my $t1 (0 .. $#hsChr19){
print HOTSPOTS "$geneClasschr19[$t1]\t$geneNamechr19[$t1]\t19\t$hsChr19[$t1]\t$hs_counts_chr19[$t1]\n";}
for my $t1 (0 .. $#hsChr20){
print HOTSPOTS "$geneClasschr20[$t1]\t$geneNamechr20[$t1]\t20\t$hsChr20[$t1]\t$hs_counts_chr20[$t1]\n";}
for my $t1 (0 .. $#hsChr21){
print HOTSPOTS "$geneClasschr21[$t1]\t$geneNamechr21[$t1]\t21\t$hsChr21[$t1]\t$hs_counts_chr21[$t1]\n";}
for my $t1 (0 .. $#hsChr22){
print HOTSPOTS "$geneClasschr22[$t1]\t$geneNamechr22[$t1]\t22\t$hsChr22[$t1]\t$hs_counts_chr22[$t1]\n";}
for my $t1 (0 .. $#hsChrX){
print HOTSPOTS "$geneClasschrX[$t1]\t$geneNamechrX[$t1]\t23\t$hsChrX[$t1]\t$hs_counts_chrX[$t1]\n";}
for my $t1 (0 .. $#hsChrY){
print HOTSPOTS "$geneClasschrY[$t1]\t$geneNamechrY[$t1]\t24\t$hsChrY[$t1]\t$hs_counts_chrY[$t1]\n";}
for my $t1 (0 .. $#hsChrM){
print HOTSPOTS "$geneClasschrM[$t1]\t$geneNamechrM[$t1]\t25\t$hsChrM[$t1]\t$hs_counts_chrM[$t1]\n";}


### Bioperl module for extraction of context sequences in reference genomes of each SNP.

print "\nExtracting context sequences and saving weblogo-$inputRawSNV ..\n";

$db = Bio::DB::Fasta->new('genome.fa');

for my $a1 (0..$#hsChr1){
	my $a2=$hsChr1[$a1]-$contextSequenceLength;
	my $a3=$hsChr1[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr1', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr1_$hsChr1[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr2){
	my $a2=$hsChr2[$a1]-$contextSequenceLength;
	my $a3=$hsChr2[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr2', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr2_$hsChr2[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr3){
	my $a2=$hsChr3[$a1]-$contextSequenceLength;
	my $a3=$hsChr3[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr3', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr3_$hsChr3[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr4){
	my $a2=$hsChr4[$a1]-$contextSequenceLength;
	my $a3=$hsChr4[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr4', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr4_$hsChr4[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr5){
	my $a2=$hsChr5[$a1]-$contextSequenceLength;
	my $a3=$hsChr5[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr5', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr5_$hsChr5[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr6){
	my $a2=$hsChr6[$a1]-$contextSequenceLength;
	my $a3=$hsChr6[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr6', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr6_$hsChr6[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr7){
	my $a2=$hsChr7[$a1]-$contextSequenceLength;
	my $a3=$hsChr7[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr7', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr7_$hsChr7[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr8){
	my $a2=$hsChr8[$a1]-$contextSequenceLength;
	my $a3=$hsChr8[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr8', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr8_$hsChr8[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr9){
	my $a2=$hsChr9[$a1]-$contextSequenceLength;
	my $a3=$hsChr9[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr9', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr9_$hsChr9[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr10){
	my $a2=$hsChr10[$a1]-$contextSequenceLength;
	my $a3=$hsChr10[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr10', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr10_$hsChr10[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr11){
	my $a2=$hsChr11[$a1]-$contextSequenceLength;
	my $a3=$hsChr11[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr11', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr11_$hsChr11[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr12){
	my $a2=$hsChr12[$a1]-$contextSequenceLength;
	my $a3=$hsChr12[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr12', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr12_$hsChr1[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr13){
	my $a2=$hsChr13[$a1]-$contextSequenceLength;
	my $a3=$hsChr13[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr13', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr13_$hsChr13[$a1]\n$seq\n";
}

for my $a1 (0..$#hsChr14){
	my $a2=$hsChr14[$a1]-$contextSequenceLength;
	my $a3=$hsChr14[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr14', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr14_$hsChr14[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChr15){
	my $a2=$hsChr15[$a1]-$contextSequenceLength;
	my $a3=$hsChr15[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr15', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr15_$hsChr15[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChr16){
	my $a2=$hsChr16[$a1]-$contextSequenceLength;
	my $a3=$hsChr16[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr16', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr16_$hsChr16[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChr17){
	my $a2=$hsChr17[$a1]-$contextSequenceLength;
	my $a3=$hsChr17[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr17', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr17_$hsChr17[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChr18){
	my $a2=$hsChr18[$a1]-$contextSequenceLength;
	my $a3=$hsChr18[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr18', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr18_$hsChr18[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChr19){
	my $a2=$hsChr19[$a1]-$contextSequenceLength;
	my $a3=$hsChr19[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr19', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr19_$hsChr19[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChr20){
	my $a2=$hsChr20[$a1]-$contextSequenceLength;
	my $a3=$hsChr20[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr20', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr20_$hsChr20[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChr21){
	my $a2=$hsChr21[$a1]-$contextSequenceLength;
	my $a3=$hsChr21[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr21', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr21_$hsChr21[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChr22){
	my $a2=$hsChr22[$a1]-$contextSequenceLength;
	my $a3=$hsChr22[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chr22', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chr22_$hsChr22[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChrX){
	my $a2=$hsChrX[$a1]-$contextSequenceLength;
	my $a3=$hsChrX[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chrX', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chrX_$hsChrX[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChrY){
	my $a2=$hsChrY[$a1]-$contextSequenceLength;
	my $a3=$hsChrY[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chrY', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chrY_$hsChrY[$a1]\n$seq\n";
}
for my $a1 (0..$#hsChrM){
	my $a2=$hsChrM[$a1]-$contextSequenceLength;
	my $a3=$hsChrM[$a1]+$contextSequenceLength;
	my $seq = $db->seq('chrM', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQ ">chrM_$hsChrM[$a1]\n$seq\n";
}

print "\nExtracting context sequences and saving extracted_sequences-$inputRawSNV ..\n";

$db1 = Bio::DB::Fasta->new('genome.fa');

for my $a1 (0..$#hsChr1){
	my $a2=$hsChr1[$a1]-$contextSequenceLength;
	my $a3=$hsChr1[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr1', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr1_$hsChr1[$a1]\n$seq\t$geneClasschr1[$a1]\t$geneNamechr1[$a1]\n";
}

for my $a1 (0..$#hsChr2){
	my $a2=$hsChr2[$a1]-$contextSequenceLength;
	my $a3=$hsChr2[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr2', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr2_$hsChr2[$a1]\n$seq\t$geneClasschr2[$a1]\t$geneNamechr2[$a1]\n";
}

for my $a1 (0..$#hsChr3){
	my $a2=$hsChr3[$a1]-$contextSequenceLength;
	my $a3=$hsChr3[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr3', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr3_$hsChr3[$a1]\n$seq\t$geneClasschr3[$a1]\t$geneNamechr3[$a1]\n";
}

for my $a1 (0..$#hsChr4){
	my $a2=$hsChr4[$a1]-$contextSequenceLength;
	my $a3=$hsChr4[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr4', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr4_$hsChr4[$a1]\n$seq\t$geneClasschr4[$a1]\t$geneNamechr4[$a1]\n";
}

for my $a1 (0..$#hsChr5){
	my $a2=$hsChr5[$a1]-$contextSequenceLength;
	my $a3=$hsChr5[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr5', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr5_$hsChr5[$a1]\n$seq\t$geneClasschr5[$a1]\t$geneNamechr5[$a1]\n";
}

for my $a1 (0..$#hsChr6){
	my $a2=$hsChr6[$a1]-$contextSequenceLength;
	my $a3=$hsChr6[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr6', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr6_$hsChr6[$a1]\n$seq\t$geneClasschr6[$a1]\t$geneNamechr6[$a1]\n";
}

for my $a1 (0..$#hsChr7){
	my $a2=$hsChr7[$a1]-$contextSequenceLength;
	my $a3=$hsChr7[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr7', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr7_$hsChr7[$a1]\n$seq\t$geneClasschr7[$a1]\t$geneNamechr7[$a1]\n";
}

for my $a1 (0..$#hsChr8){
	my $a2=$hsChr8[$a1]-$contextSequenceLength;
	my $a3=$hsChr8[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr8', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr8_$hsChr8[$a1]\n$seq\t$geneClasschr8[$a1]\t$geneNamechr8[$a1]\n";
}

for my $a1 (0..$#hsChr9){
	my $a2=$hsChr9[$a1]-$contextSequenceLength;
	my $a3=$hsChr9[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr9', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr9_$hsChr9[$a1]\n$seq\t$geneClasschr9[$a1]\t$geneNamechr9[$a1]\n";
}

for my $a1 (0..$#hsChr10){
	my $a2=$hsChr10[$a1]-$contextSequenceLength;
	my $a3=$hsChr10[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr10', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr10_$hsChr10[$a1]\n$seq\t$geneClasschr10[$a1]\t$geneNamechr10[$a1]\n";
}

for my $a1 (0..$#hsChr11){
	my $a2=$hsChr11[$a1]-$contextSequenceLength;
	my $a3=$hsChr11[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr11', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr11_$hsChr11[$a1]\n$seq\t$geneClasschr11[$a1]\t$geneNamechr11[$a1]\n";
}

for my $a1 (0..$#hsChr12){
	my $a2=$hsChr12[$a1]-$contextSequenceLength;
	my $a3=$hsChr12[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr12', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr12_$hsChr1[$a1]\n$seq\t$geneClasschr12[$a1]\t$geneNamechr12[$a1]\n";
}

for my $a1 (0..$#hsChr13){
	my $a2=$hsChr13[$a1]-$contextSequenceLength;
	my $a3=$hsChr13[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr13', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr13_$hsChr13[$a1]\n$seq\t$geneClasschr13[$a1]\t$geneNamechr13[$a1]\n";
}

for my $a1 (0..$#hsChr14){
	my $a2=$hsChr14[$a1]-$contextSequenceLength;
	my $a3=$hsChr14[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr14', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr14_$hsChr14[$a1]\n$seq\t$geneClasschr14[$a1]\t$geneNamechr14[$a1]\n";
}
for my $a1 (0..$#hsChr15){
	my $a2=$hsChr15[$a1]-$contextSequenceLength;
	my $a3=$hsChr15[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr15', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr15_$hsChr15[$a1]\n$seq\t$geneClasschr15[$a1]\t$geneNamechr15[$a1]\n";
}
for my $a1 (0..$#hsChr16){
	my $a2=$hsChr16[$a1]-$contextSequenceLength;
	my $a3=$hsChr16[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr16', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr16_$hsChr16[$a1]\n$seq\t$geneClasschr16[$a1]\t$geneNamechr16[$a1]\n";
}
for my $a1 (0..$#hsChr17){
	my $a2=$hsChr17[$a1]-$contextSequenceLength;
	my $a3=$hsChr17[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr17', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr17_$hsChr17[$a1]\n$seq\t$geneClasschr17[$a1]\t$geneNamechr17[$a1]\n";
}
for my $a1 (0..$#hsChr18){
	my $a2=$hsChr18[$a1]-$contextSequenceLength;
	my $a3=$hsChr18[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr18', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr18_$hsChr18[$a1]\n$seq\t$geneClasschr18[$a1]\t$geneNamechr18[$a1]\n";
}
for my $a1 (0..$#hsChr19){
	my $a2=$hsChr19[$a1]-$contextSequenceLength;
	my $a3=$hsChr19[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr19', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr19_$hsChr19[$a1]\n$seq\t$geneClasschr19[$a1]\t$geneNamechr19[$a1]\n";
}
for my $a1 (0..$#hsChr20){
	my $a2=$hsChr20[$a1]-$contextSequenceLength;
	my $a3=$hsChr20[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr20', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr20_$hsChr20[$a1]\n$seq\t$geneClasschr20[$a1]\t$geneNamechr20[$a1]\n";
}
for my $a1 (0..$#hsChr21){
	my $a2=$hsChr21[$a1]-$contextSequenceLength;
	my $a3=$hsChr21[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr21', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr21_$hsChr21[$a1]\n$seq\t$geneClasschr21[$a1]\t$geneNamechr21[$a1]\n";
}
for my $a1 (0..$#hsChr22){
	my $a2=$hsChr22[$a1]-$contextSequenceLength;
	my $a3=$hsChr22[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chr22', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chr22_$hsChr22[$a1]\n$seq\t$geneClasschr22[$a1]\t$geneNamechr22[$a1]\n";
}
for my $a1 (0..$#hsChrX){
	my $a2=$hsChrX[$a1]-$contextSequenceLength;
	my $a3=$hsChrX[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chrX', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chrX_$hsChrX[$a1]\n$seq\t$geneClasschrX[$a1]\t$geneNamechrX[$a1]\n";
}
for my $a1 (0..$#hsChrY){
	my $a2=$hsChrY[$a1]-$contextSequenceLength;
	my $a3=$hsChrY[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chrY', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chrY_$hsChrY[$a1]\n$seq\t$geneClasschrY[$a1]\t$geneNamechrY[$a1]\n";
}
for my $a1 (0..$#hsChrM){
	my $a2=$hsChrM[$a1]-$contextSequenceLength;
	my $a3=$hsChrM[$a1]+$contextSequenceLength;
	my $seq = $db1->seq('chrM', $a2 => $a3);
	$a2=0;
	$a3=0;
	print CONTEXTSEQANNO ">chrM_$hsChrM[$a1]\n$seq\t$geneClasschrM[$a1]\t$geneNamechrM[$a1]\n";
}


## closing of outputs

close (">>results-$inputRawSNV/WOLAND-basechange-$inputRawSNV");

close (">>results-$inputRawSNV/WOLAND-mutfreq-$inputRawSNV");

close (">>results-$inputRawSNV/WOLAND-contextsequences-$inputRawSNV");

close (">>results-$inputRawSNV/WOLAND-hotspots-$inputRawSNV");

close (">>results-$inputRawSNV/WOLAND-log-$inputRawSNV");

close (">>results-$inputRawSNV/WOLAND-contextsequencesanno-$inputRawSNV");


####################### WOLAND MOTIF SEARCHER ###################################

#outputs:

open (OUTPUTbiasSN1, ">>results-$inputRawSNV/WOLAND-bias_SN1-$inputRawSNV");
open (OUTPUTbiasDNApoln, ">>results-$inputRawSNV/WOLAND-bias_DNApoln-$inputRawSNV");
open (OUTPUTbiasoxoG, ">>results-$inputRawSNV/WOLAND-bias_oxoG-$inputRawSNV");
open (OUTPUTbiasUV, ">>results-$inputRawSNV/WOLAND-bias_UV-lambda-$inputRawSNV");
open (OUTPUTbiasUVsolar, ">>results-$inputRawSNV/WOLAND-bias_UVsolar-$inputRawSNV");
open (OUTPUTbiassixfour, ">>results-$inputRawSNV/WOLAND-bias_sixfour-$inputRawSNV");
open (OUTPUTbiasenu, ">>results-$inputRawSNV/WOLAND-bias_enu-$inputRawSNV");

# Open Context Sequences generated

open CONTEXTSEQANNO, "<results-$inputRawSNV/WOLAND-contextsequencesanno-$inputRawSNV";
@fastaContext=<CONTEXTSEQANNO>;

# Array processing

foreach $fastaContext (@fastaContext)
	{
	@i7=0;
	@i7 = split (/\t/, "$fastaContext");
	chomp (@i7);
	push (@fastaContextMS, "$i7[0]");
	push (@geneClassMS, "$i7[1]");
	push (@geneNameMS, "$i7[2]");
	}

foreach $geneClassMS (@geneClassMS){
	if ($geneClassMS =~ /[\S]/){
	push (@geneClass1, $geneClassMS);}
}

foreach $geneNameMS (@geneNameMS){
	if ($geneNameMS =~ /[\S]/){
	push (@geneName1, $geneNameMS);}
}
	
foreach $fastaContextMS (@fastaContextMS){
	if ($fastaContextMS=~ /^>+/){
	chomp (@fastaContextMS);
	push (@id, $fastaContextMS);}
	else{
	push (@targetSequence, $fastaContextMS);}
	}

foreach $id (@id)
	{
	@i2 = split (/_/, $id);
	chomp (@i2);
	push (@chrRaw, "$i2[0]");
	push (@coord, "$i2[1]");
	}

foreach $chrRaw (@chrRaw)
	{
	@i3 = split (/>/, $chrRaw);
	chomp (@i3);
	push (@chrSt, "$i3[1]");
	}

############################# STRAND BIAS REFSEQ #############################################
$refSeqRaw = 'refseq.txt';
open (REFSEQRAW, $refSeqRaw);
@refSeq=<REFSEQRAW>;
unless (@refSeq){
	die "\nERROR : Could not load <refseq.txt>\n";
}

foreach $refSeqline (@refSeq)
	{
	@irefSeq = split (/\t/, $refSeqline);
	chomp (@irefSeq);
	push (@chrRefSeq, "$irefSeq[2]");
	push (@strandRefSeq, "$irefSeq[3]");
	push (@txstartRefSeq, $irefSeq[4]);
	push (@txtendRefSeq, "$irefSeq[5]");
	}


# Search for Mutable Motifs

print "\nSearching for mutable motifs\n";

## SN1 MOTIF

foreach $targetSequence (@targetSequence){
	if ($targetSequence =~ "..AG..." || $targetSequence =~ "..GG..."|| $targetSequence =~ "...CT.."|| $targetSequence =~ "...CC.."){
		push (@SN1, "1");}
	else {
		push (@SN1, "0");
	}
	if ($targetSequence =~ "..AG..." || $targetSequence =~ "..GG..."){
		push (@SN1plus, "1");}
	else {
		push (@SN1plus, "0")
	}
	if ($targetSequence =~ "...CT.."|| $targetSequence =~ "...CC.."){
		push (@SN1minus, "1");}
	else {
		push (@SN1minus, "0")
	}
}

#### Strand Concordance Score SN1

for my $iSN1 (0..$#SN1){
	if ($SN1[$iSN1] eq "1"){
		$strand_value=$SN1plus[$iSN1];
	$query_chr=$chrSt[$iSN1];
	$query_coord=$coord[$iSN1];
	$strand_plus_count=0;
    $strand_count=0;
	for my $i4 (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_count;}
			if ($strandRefSeq[$i4] eq "+" && $chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_plus_count;} 

	}
	
	$i4=0;
	if ($strand_count >= 1){
	 $strand_transcript=$strand_plus_count/$strand_count;
     $strand_score=$strand_transcript - ($strand_value);
     print OUTPUTbiasSN1 "$chrSt[$iSN1]\t$coord[$iSN1]\t$strand_score\n";


	$strand_transcript=0;
    $strand_plus_count=0;
    $strand_count=0;
    $strand_score=0;
	}
   
    
   
                     }
 
	}

## DNA POL ETA MOTIF

foreach $targetSequence (@targetSequence){
	if ($targetSequence =~ "..AA..." || $targetSequence =~ "..TA..."|| $targetSequence =~ "...TT.."|| $targetSequence =~ "...TA.."){
		push (@DNApoln, "1");}
	else {
		push (@DNApoln, "0");
	}
		if ($targetSequence =~ "..AA..." || $targetSequence =~ "..TA..."){
		push (@DNApolnplus, "1");}
	else {
		push (@DNApolnplus, "0")
	}
	if ($targetSequence =~ "...TT.."|| $targetSequence =~ "...TA.."){
		push (@DNApolnminus, "1");}
	else {
		push (@DNApolnminus, "0")
	}
}

#### Strand Concordance Score DNAPOLeta

for my $iDNApoln (0..$#DNApoln){
	if ($DNApoln[$iDNApoln] eq "1"){
		$strand_value=$DNApolnplus[$iDNApoln];
	$query_chr=$chrSt[$iDNApoln];
	$query_coord=$coord[$iDNApoln];
	$strand_plus_count=0;
    $strand_count=0;
	for my $i4 (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_count;}
			if ($strandRefSeq[$i4] eq "+" && $chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_plus_count;} 

	}
	
	$i4=0;
	if ($strand_count >= 1){
	 $strand_transcript=$strand_plus_count/$strand_count;
     $strand_score=$strand_transcript - ($strand_value);
     print OUTPUTbiasDNApoln "$chrSt[$iDNApoln]\t$coord[$iDNApoln]\t$strand_score\n";


	$strand_transcript=0;
    $strand_plus_count=0;
    $strand_count=0;
    $strand_score=0;
	}
   
    
   
                     }
 
	}



## OXO G MOTIF

foreach $targetSequence (@targetSequence){
	if ($targetSequence =~ "..AGA.." || $targetSequence =~ "..GGG.." || $targetSequence =~ "..AGG.." || $targetSequence =~ "..GGA.."|| $targetSequence =~ "..TCT.." || $targetSequence =~ "..CCC.." || $targetSequence =~ "..CCT.."|| $targetSequence =~ "..TCC.."){
		push (@oxoG, "1");}
	else {
		push (@oxoG, "0");
	}
		if ($targetSequence =~ "..AGA.." || $targetSequence =~ "..GGG.." || $targetSequence =~ "..AGG.." || $targetSequence =~ "..GGA.."){
		push (@oxoGplus, "1");}
	else {
		push (@oxoGplus, "0")
	}
	if ($targetSequence =~ "..TCT.." || $targetSequence =~ "..CCC.." || $targetSequence =~ "..CCT.."|| $targetSequence =~ "..TCC.."){
		push (@oxoGminus, "1");}
	else {
		push (@oxoGminus, "0")
	}
}

#### Strand Concordance Score 8-oxoG

for my $ioxoG (0..$#oxoG){
	if ($oxoG[$ioxoG] eq "1"){
		$strand_value=$oxoGplus[$ioxoG];
	$query_chr=$chrSt[$ioxoG];
	$query_coord=$coord[$ioxoG];
	$strand_plus_count=0;
    $strand_count=0;
	for my $i4 (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_count;}
			if ($strandRefSeq[$i4] eq "+" && $chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_plus_count;} 

	}
	$i4=0;
	if ($strand_count >= 1){
	 $strand_transcript=$strand_plus_count/$strand_count;
     $strand_score=$strand_transcript - ($strand_value);
     print OUTPUTbiasoxoG "$chrSt[$ioxoG]\t$coord[$ioxoG]\t$strand_score\n";
	
	$strand_transcript=0;
    $strand_plus_count=0;
    $strand_count=0;
    $strand_score=0;
	}
   
    
   
                     }
 
	}




## UV LAMBDA MOTIF

foreach $targetSequence (@targetSequence){
	if ($targetSequence =~ "...TC.." || $targetSequence =~ "..TC..." || $targetSequence =~ "...CT.." || $targetSequence =~ "..CT..." ||
	    $targetSequence =~ "..TT..." || $targetSequence =~ "...TT.." || $targetSequence =~ "..CC..." || $targetSequence =~ "...CC.."|| 
	    $targetSequence =~ "...GA.." || $targetSequence =~ "..GA..." || $targetSequence =~ "...AA.." || $targetSequence =~ "..AA..." ||
	    $targetSequence =~ "...GG.." || $targetSequence =~ "..GG..." || $targetSequence =~ "..AG..." || $targetSequence =~ "...AG.."){
		push (@UVlambda, "1");}
	else {
		push (@UVlambda, "0");
	}
		if ($targetSequence =~ "...TC.." || $targetSequence =~ "..TC..." || $targetSequence =~ "...CT.." || $targetSequence =~ "..CT..." ||
		    $targetSequence =~ "..TT..." || $targetSequence =~ "...TT.." || $targetSequence =~ "..CC..." || $targetSequence =~ "...CC.."){
		push (@UVlambdaplus, "1");}
	else {
		push (@UVlambdaplus, "0")
	}
		if ($targetSequence =~ "...GA.." || $targetSequence =~ "..GA..." || $targetSequence =~ "...AA.." || $targetSequence =~ "..AA..." ||
		    $targetSequence =~ "...GG.." || $targetSequence =~ "..GG..." || $targetSequence =~ "..AG..." || $targetSequence =~ "...AG.."){
		push (@UVlambdaminus, "1");}
	else {
		push (@UVlambdaminus, "0")
	}
}

#### Strand Concordance Score UV-lambda


for my $iuv (0..$#UVlambda){
	if ($UVlambda[$iuv] eq "1"){
		$strand_value=$UVlambdaplus[$iuv];
	$query_chr=$chrSt[$iuv];
	$query_coord=$coord[$iuv];
	$strand_plus_count=0;
   $strand_count=0;
	for my $i5 (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i5] eq $query_chr and (($query_coord <= $txstartRefSeq[$i5] and $query_coord >= $txtendRefSeq[$i5]) or ($query_coord >= $txstartRefSeq[$i5] and $query_coord <= $txtendRefSeq[$i5]))){
				++$strand_count;}
			if ($strandRefSeq[$i5] eq "+" && $chrRefSeq[$i5] eq $query_chr and (($query_coord <= $txstartRefSeq[$i5] and $query_coord >= $txtendRefSeq[$i5]) or ($query_coord >= $txstartRefSeq[$i5] and $query_coord <= $txtendRefSeq[$i5]))){
				++$strand_plus_count;} 

	}
	$i5=0;
	if ($strand_count >= 1){
	 $strand_transcript=$strand_plus_count/$strand_count;
     $strand_score=$strand_transcript - ($strand_value);
     print OUTPUTbiasUV "$chrSt[$iuv]\t$coord[$iuv]\t$strand_score\n";
	
	$strand_transcript=0;
    $strand_plus_count=0;
    $strand_count=0;
    $strand_score=0;
	}
   
    
   
                     }
 
	}

## UV SOLAR UVA-1


foreach $targetSequence (@targetSequence){
	if ($targetSequence =~ "..TCG.."  || $targetSequence =~ "...TCG."|| $targetSequence =~ ".TCG..." || $targetSequence =~ "..CCG.." ||$targetSequence =~ "...CCG." || $targetSequence =~ ".CCG..." || 
	    $targetSequence =~ "..CGA.." || $targetSequence =~ "...CGA."|| $targetSequence =~ ".CGA..." || $targetSequence =~ "..CGG.." || $targetSequence =~ "...CGG."|| $targetSequence =~ ".CGG..."){
		push (@UVsolar, "1");}
	else {
		push (@UVsolar, "0");
	}
	
	if ($targetSequence =~ "..TCG.."  || $targetSequence =~ "...TCG."|| $targetSequence =~ ".TCG..." || 
		    $targetSequence =~ "..CCG.." ||$targetSequence =~ "...CCG." || $targetSequence =~ ".CCG..."){
		push (@UVsolarplus, "1");}
	else {
		push (@UVsolarplus, "0")
	}
		if ($targetSequence =~ "..CGA.." || $targetSequence =~ "...CGA."|| $targetSequence =~ ".CGA..." ||
		    $targetSequence =~ "..CGG.." || $targetSequence =~ "...CGG."|| $targetSequence =~ ".CGG...") {
		push (@UVsolarminus, "1");}
	else {
		push (@UVsolarminus, "0")
	}
}

#### Strand Concordance Score UVsolar

for my $iUVsolar (0..$#UVsolar){
	if ($UVsolar[$iUVsolar] eq "1"){
		$strand_value=$UVsolarplus[$iUVsolar];
	$query_chr=$chrSt[$iUVsolar];
	$query_coord=$coord[$iUVsolar];
	$strand_plus_count=0;
    $strand_count=0;
	for my $i4 (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_count;}
			if ($strandRefSeq[$i4] eq "+" && $chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_plus_count;} 

	}
	$i4=0;
	if ($strand_count >= 1){
	 $strand_transcript=$strand_plus_count/$strand_count;
     $strand_score=$strand_transcript - ($strand_value);
     print OUTPUTbiasUVsolar "$chrSt[$iUVsolar]\t$coord[$iUVsolar]\t$strand_score\n";
	
	$strand_transcript=0;
    $strand_plus_count=0;
    $strand_count=0;
    $strand_score=0;
	}
   
    
   
                     }
 
	}


## SIX FOUR MOTIF

foreach $targetSequence (@targetSequence){
	if ($targetSequence =~ "TTCA" || $targetSequence =~ "CTCA" || $targetSequence =~ "TGAA" || $targetSequence =~ "TGAG"){
		push (@sixfour, "1");}
	else {
		push (@sixfour, "0");
	}
		if ($targetSequence =~ "TTCA" || $targetSequence =~ "CTCA"){
		push (@sixfourplus, "1");}
	else {
		push (@sixfourplus, "0")
	}
		if ($targetSequence =~"TGAA" || $targetSequence =~ "TGAG") {
		push (@sixfourminus, "1");}
	else {
		push (@sixfourminus, "0")
	}
}

#### Strand Concordance Score 6-4

for my $isixfour (0..$#sixfour){
	if ($sixfour[$isixfour] eq "1"){
		$strand_value=$sixfourplus[$isixfour];
	$query_chr=$chrSt[$isixfour];
	$query_coord=$coord[$isixfour];
	$strand_plus_count=0;
    $strand_count=0;
	for my $i4 (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_count;}
			if ($strandRefSeq[$i4] eq "+" && $chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_plus_count;} 

	}
	$i4=0;
	if ($strand_count >= 1){
	 $strand_transcript=$strand_plus_count/$strand_count;
     $strand_score=$strand_transcript - ($strand_value);
     print OUTPUTbiassixfour "$chrSt[$isixfour]\t$coord[$isixfour]\t$strand_score\n";
	
	$strand_transcript=0;
    $strand_plus_count=0;
    $strand_count=0;
    $strand_score=0;
	}
   
    
   
                     }
 
	}

## ENU MOTIF

foreach $targetSequence (@targetSequence){
	if ($targetSequence =~ "..CAG.." || $targetSequence =~ "..GAC.."|| $targetSequence =~ "..GAG.."|| $targetSequence =~ "..CAC.."|| $targetSequence =~ "..CTG.." || $targetSequence =~ "..GTC.."|| $targetSequence =~ "..GTG.."|| $targetSequence =~ "..CTC.."){
		push (@enu, "1");}
	else {
		push (@enu, "0");
	}
			if ($targetSequence =~ "..CAG.." || $targetSequence =~ "..GAC.."|| $targetSequence =~ "..GAG.."|| $targetSequence =~ "..CAC.."){
		push (@enuplus, "1");}
	else {
		push (@enuplus, "0")
	}
		if ($targetSequence =~ "..CTG.." || $targetSequence =~ "..GTC.."|| $targetSequence =~ "..GTG.."|| $targetSequence =~ "..CTC..") {
		push (@enuminus, "1");}
	else {
		push (@enuminus, "0")
	}
	
}

#### Strand Concordance Score ENU

for my $ienu (0..$#enu){
	if ($enu[$ienu] eq "1"){
		$strand_value=$enuplus[$ienu];
	$query_chr=$chrSt[$ienu];
	$query_coord=$coord[$ienu];
	$strand_plus_count=0;
    $strand_count=0;
	for my $i4 (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_count;}
			if ($strandRefSeq[$i4] eq "+" && $chrRefSeq[$i4] eq $query_chr and (($query_coord <= $txstartRefSeq[$i4] and $query_coord >= $txtendRefSeq[$i4]) or ($query_coord >= $txstartRefSeq[$i4] and $query_coord <= $txtendRefSeq[$i4]))){
				++$strand_plus_count;} 

	}
	$i4=0;
	if ($strand_count >= 1){
	 $strand_transcript=$strand_plus_count/$strand_count;
     $strand_score=$strand_transcript - ($strand_value);
     print OUTPUTbiasenu "$chrSt[$ienu]\t$coord[$ienu]\t$strand_score\n";
	
	$strand_transcript=0;
    $strand_plus_count=0;
    $strand_count=0;
    $strand_score=0;
	}
   
    
   
                     }
 
	}

## CALCULATION OF NORMALIZED NUMBER OF MOTIFS FOUND

for my $a1 (0 .. $#id){
		++$SNPnumber;
	}

foreach $SN1 (@SN1){
	if ($SN1 eq "1"){
		++$SN1counts;}
}

foreach $DNApoln (@DNApoln){
	if ($DNApoln eq "1"){
		++$DNApolncounts;}
}

foreach $oxoG (@oxoG){
	if ($oxoG eq "1"){
		++$oxoGcounts;}
}
foreach $UVlambda (@UVlambda){
	if ($UVlambda eq "1"){
		++$UVlambdacounts;}
}

foreach $sixfour (@sixfour){
	if ($sixfour eq "1"){
		++$sixfourcounts;}
}

foreach $enu (@enu){
	if ($enu eq "1"){
		++$enucounts;}
}

foreach $UVsolar (@UVsolar){
	if ($UVsolar eq "1"){
		++$UVsolarcounts;}
}


$normSN1=$SN1counts/$SNPnumber;
$normDNApoln=$DNApolncounts/$SNPnumber;
$normoxoG=$oxoGcounts/$SNPnumber;
$normUVlambda=$UVlambdacounts/$SNPnumber;
$normsixfour=$sixfourcounts/$SNPnumber;
$normenu=$enucounts/$SNPnumber;
$normUVsolar=$UVsolarcounts/$SNPnumber;




#### PRINTING OF RAW MUTABLE MOTIFS FILE #####

print "\nSaving mutational motifs file\n";
	print MOTIFS "Chr\tPos\ttargetSequence\tClass\tGene\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUVA-solar\n";
for $a2 (0..$#targetSequence){
	print MOTIFS "$chrSt[$a2]\t$coord[$a2]\t$targetSequence[$a2]\t$geneClass1[$a2]\t$geneName1[$a2]\t$SN1[$a2]\t$DNApoln[$a2]\t$oxoG[$a2]\t$UVlambda[$a2]\t$sixfour[$a2]\t$enu[$a2]\t$UVsolar[$a2]\n";
}

#### PRINTING OF NORMALIZED MUTABLE MOTIFS FILE #####
print "\nSaving normalized mutational motifs";
	print NMOTIFS "$inputRawSNV\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUVA-solar\n";
	print NMOTIFS "Number of Total SNPs\t$SNPnumber\t$SNPnumber\t$SNPnumber\t$SNPnumber\t$SNPnumber\t$SNPnumber\t$SNPnumber\n";
	print NMOTIFS "Total Raw Number of Motifs Found\t$SN1counts\t$DNApolncounts\t$oxoGcounts\t$UVlambdacounts\t$sixfourcounts\t$enucounts\t$UVsolarcounts\n";
	print NMOTIFS "Normalized Number of Motifs Found\t$normSN1\t$normDNApoln\t$normoxoG\t$normUVlambda\t$normsixfour\t$normenu\t$normUVsolar";
	
	exit;
exit;

