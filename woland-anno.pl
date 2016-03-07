#########################################################################################################################
## WOLAND BETA 0.1 (15-02-2016)
##
## WOLAND is a software package based on Perl and R for calculation of general mutation metrics, identification and
## comparison of predicted hotspots across biological samples. WOLAND uses Single Nucleotide Polymorphisms (SNPs) data
## from Next Generation Sequencing (NGS) pipelines as main input. Please read README file.
##
## USAGE
##
## perl woland.pl <annnovar_variant_file> <chromosome_length_profile> <hotspot_window_length> <genome_version>
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
## <genome_version> : Genome version as found in genomes/genome_<genome_version>.fa and genomes/refseq_<genome_version>.txt
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
my ($inputRawSNV, $chrLengthProfile, $hotspotWindowLength, $contextSequenceLength, $genome_version);
my $datestring;
my $rawLine; 
my $i; my $i2; our $i3; my @i; my @i2; our @i3;
my ($ATTA, $AGTC, $ACTG, $CGGC, $CTGA, $CAGT, $refalt, $SOMA, $AVGATTA, $AVGAGTC, $AVGACTG, $AVGCGGC, $AVGCTGA, $AVGCAGT);
my ($transition, $transversion);
my ($chr1,$chr2,$chr3,$chr4,$chr5,$chr6,$chr7,$chr8,$chr9,$chr10,$chr11,
	$chr12,$chr13,$chr14,$chr15,$chr16,$chr17,$chr18,$chr19,$chr20,$chr21,
	$chr22,$chrX,$chrY,$chrM);
my @chrList;
my $countChr;
my @chrCountFreq;
my @chrNormList;
my $chrNorm;
my ($hotspotSNV, $hotspotDownstreamSNV, $hotspotUpstreamSNV, $hsNum); 
my $config; my $count; my $notValid;
my @rawFile; my @config; my @chr; my @pos; my @alt; my @chrpos; my %posref; my %posalt; 
my @refalt; my @ref; my @chrLength;
my $hsChr1; my @hsChr1; my @hs_counts_chr1;
my $hsChr2; my @hsChr2; my @hs_counts_chr2;
my $hsChr3; my @hsChr3; my @hs_counts_chr3;
my $hsChr4; my @hsChr4; my @hs_counts_chr4;
my $hsChr5; my @hsChr5; my @hs_counts_chr5;
my $hsChr6; my @hsChr6; my @hs_counts_chr6;
my $hsChr7; my @hsChr7; my @hs_counts_chr7;
my $hsChr8; my @hsChr8; my @hs_counts_chr8;
my $hsChr9; my @hsChr9; my @hs_counts_chr9;
my $hsChr10; my @hsChr10; my @hs_counts_chr10;
my $hsChr11; my @hsChr11; my @hs_counts_chr11;
my $hsChr12; my @hsChr12; my @hs_counts_chr12;
my $hsChr13; my @hsChr13; my @hs_counts_chr13;
my $hsChr14; my @hsChr14; my @hs_counts_chr14;
my $hsChr15; my @hsChr15; my @hs_counts_chr15;
my $hsChr16; my @hsChr16; my @hs_counts_chr16;
my $hsChr17; my @hsChr17; my @hs_counts_chr17;
my $hsChr18; my @hsChr18; my @hs_counts_chr18;
my $hsChr19; my @hsChr19; my @hs_counts_chr19;
my $hsChr20; my @hsChr20; my @hs_counts_chr20;
my $hsChr21; my @hsChr21; my @hs_counts_chr21;
my $hsChr22; my @hsChr22; my @hs_counts_chr22;
my $hsChrX; my @hsChrX; my @hs_counts_chrX;
my $hsChrY; my @hsChrY; my @hs_counts_chrY;
my $hsChrM; my @hsChrM; my @hs_counts_chrM;
my @fastaContext; my $fastaContext; my @id; my @targetSequence; my $targetSequence;
my @SN1; my @DNApoln; my @oxoG; my @UVlambda; my @UVsolar; my @sixfour; my @enu; my $SNPnumber; my $mf1; my $mf2;
my @SN1counts; my @DNApolncounts; my @oxoGcounts; my @UVlambdacounts; my @UVsolarcounts; my @sixfourcounts; my @enucounts;
my $SN1; my $DNApoln; my $oxoG; my $UVlambda; my $UVsolar; my $sixfour; my $enu;
my $SN1counts; my $DNApolncounts; my $oxoGcounts; my $UVlambdacounts; my $UVsolarcounts; my $sixfourcounts; my $enucounts;
my $normSN1; my $normDNApoln; my $normoxoG; my $normUVlambda; my $normUVsolar; my $normsixfour; my $normenu;
## anno motif search variables
my (@geneClass, $geneClass, @geneName, $geneName);
my (@geneClasschr1, @geneClasschr2, @geneClasschr3, @geneClasschr4, @geneClasschr5, @geneClasschr6, @geneClasschr7, @geneClasschr8, @geneClasschr9, @geneClasschr10, @geneClasschr11, @geneClasschr12, @geneClasschr13, @geneClasschr14, @geneClasschr15, @geneClasschr16, @geneClasschr17, @geneClasschr18, @geneClasschr19, @geneClasschr20, @geneClasschr21, @geneClasschr22, @geneClasschrX, @geneClasschrY, @geneClasschrM);      
my (@geneNamechr1, @geneNamechr2, @geneNamechr3, @geneNamechr4, @geneNamechr5, @geneNamechr6, @geneNamechr7, @geneNamechr8, @geneNamechr9, @geneNamechr10, @geneNamechr11, @geneNamechr12, @geneNamechr13, @geneNamechr14, @geneNamechr15, @geneNamechr16, @geneNamechr17, @geneNamechr18, @geneNamechr19, @geneNamechr20, @geneNamechr21, @geneNamechr22, @geneNamechrX, @geneNamechrY, @geneNamechrM);
my (@chrStart, @chrEnd);
my (@fastaContextMS, $fastaContextMS, @fastaContextMS1, $fastacontextMS1, @geneClassMS, $geneClassMS, @geneClassMS1, $geneClassMS1, @geneNameMS, $geneNameMS, @geneNameMS1, $geneNameMS1);
my (@chrRaw, $chrRaw, $chrSt, $refSeqRaw);
my (@chrRefSeq, @irefSeq, );
my ($query_coord, $query_chr);
my (@txtendRefSeq, @txstartRefSeq);
my $id;
my (@geneClass1, @geneName1, @coord, @chrSt, @refSeq, $refSeqline, @strandRefSeq);
my (@SN1plus, @SN1minus, @DNApolnplus, @DNApolnminus, @oxoGplus, @oxoGminus, @UVlambdaplus, @UVlambdaminus, @UVsolarplus, @UVsolarminus, @sixfourplus, @sixfourminus, @enuplus, @enuminus);
my ($strand_count, $strand_plus_count, $strand_transcript, $strand_score, $strand_value);
my ($db, $db1);
my $readContextSeqAnno;

# main warning
unless (@ARGV){
	die "\nERROR : Usage: $0 <tabular_snp_file> <chromosome_length_profile> <hotspot_window_length> <genome_version> \n";	
}

# subroutines
sub chrListCount {
	for my $i3 (0 .. $#chrStart){
		if ($chrStart[$i3] eq "$_[0]"){
			++$countChr;
		}
	}
	push @chrCountFreq, $countChr;
	$countChr=0;
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

$genome_version = $ARGV[3]; #<genome_version>
unless ($genome_version){
	die "\nERROR : Please specify a genome version in genomes\/folder for <genome_version>\n";
}

$contextSequenceLength=3; #<context_sequence_length> #default=3nt downstream & 3nt upstream

# loading of outputs
mkdir("results-$inputRawSNV", 0755) || die "Cannot mkdir results-$inputRawSNV - folder already exists, please delete it or change samplename";

open (BASECHANGE, ">>results-$inputRawSNV/WOLAND-basechange-$inputRawSNV");
open (MUTFREQ, ">>results-$inputRawSNV/WOLAND-mutfreq-$inputRawSNV");
open (CONTEXTSEQ, ">>results-$inputRawSNV/WOLAND-contextsequences-$inputRawSNV");
open (CONTEXTSEQANNO, ">>results-$inputRawSNV/WOLAND-contextsequencesanno-$inputRawSNV");
open (HOTSPOTS, ">>results-$inputRawSNV/WOLAND-hotspots-$inputRawSNV");
open (MOTIFS, ">>results-$inputRawSNV/WOLAND-motifs-$inputRawSNV");
open (NMOTIFS, ">>results-$inputRawSNV/WOLAND-norm_motifs-$inputRawSNV");

# start Screen & log file
$datestring = localtime();

open (LOG, ">>results-$inputRawSNV/WOLAND-log-$inputRawSNV");

print LOG "WOLAND BETA 0.9 - 02-02-2016\n";
print LOG "Analysis started at $datestring\n";
print LOG "Tabular SNP File              : $inputRawSNV\n";
print LOG "Chromosome Length Profile File: $chrLengthProfile\n";
print LOG "Hotspot Window Length         : $hotspotWindowLength bases flanking SNP position in reference genome\n";
print LOG "Context Sequence Length       : $contextSequenceLength bases flanking SNP position in reference genome\n";
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
print LOG "Strand Bias of UV Solar Motifs:            WOLAND-bias_UVsolar-$inputRawSNV\n";
print LOG "Strand Bias of 6-4 Motifs:                 WOLAND-bias_sixfour-$inputRawSNV\n";
print LOG "Strand Bias of ENU Motifs:                 WOLAND-bias_enu-$inputRawSNV\n";


# conversion of each column in a dedicated array
print "\nLoading SNP file...\n";

foreach $rawLine (@rawFile){
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
for my $i (0 .. $#chrStart){
	push @chrpos, "$chrStart[$i]_$pos[$i]";
}
	
# Hashe posref & posalt for information of each position - ALT e REF alelles
for my $i (0 .. $#chrpos){
	$posref{"$chrpos[$i]_$i"} .= "$ref[$i]";
}
for my $i (0 .. $#chrpos){
	$posalt{"$chrpos[$i]_$i"} .= "$alt[$i]";
}

print "\nCalculating general mutation statistics and saving basechange-$inputRawSNV ...\n";

# Array for a single entry for each SNP containing a single REFALT string value. Considering  A->T # & T->A; A->G & T->C; A->C & T->G; C->G & G->C; C->T & G->A; C->A & G->T.
foreach my $key (sort keys %posref){
	push @refalt,"$posref{$key}$posalt{$key}";
}


# counting of nucleotide changes.
$ATTA=0;
$AGTC=0;
$ACTG=0;
$CGGC=0;
$CTGA=0;
$CAGT=0;

foreach $refalt (@refalt){
	if ($refalt eq "AT" || $refalt eq "TA"){
		$ATTA++;
	}
	if ($refalt eq "AG" || $refalt eq "TC"){
		$AGTC++;
	}
	if ($refalt eq "AC" || $refalt eq "TG"){
		$ACTG++;
	}
	if ($refalt eq "CG" || $refalt eq "GC"){
		$CGGC++;
	}
	if ($refalt eq "CT" || $refalt eq "GA"){
		$CTGA++;
	}
	if ($refalt eq "CA" || $refalt eq "GT"){
		$CAGT++;
	}
}

# frequency of changes considering total amount of changes.
$SOMA=$ATTA+$AGTC+$ACTG+$CGGC+$CTGA+$CAGT;
if ($SOMA== 0){
	die "Please review input file format";
}

foreach $refalt(@refalt){
	$count++;
}

print LOG "\nTotal SNPs valid=$count\n";
$notValid=$count-$SOMA;
print LOG "Total not valid SNPs=$notValid\n";

$AVGATTA=$ATTA/$SOMA;
$AVGAGTC=$AGTC/$SOMA;
$AVGACTG=$ACTG/$SOMA;
$AVGCGGC=$CGGC/$SOMA;
$AVGCTGA=$CTGA/$SOMA;
$AVGCAGT=$CAGT/$SOMA;

# transitions & transversions
$transition=($AGTC+$CTGA)/$SOMA;
$transversion=($ATTA+$ACTG+$CGGC+$CAGT)/$SOMA;

# frequency per chromosome target
$countChr=0;
$chr1=0;$chr2=0;$chr3=0;$chr4=0;$chr5=0;$chr6=0;$chr7=0;$chr8=0;$chr9=0;$chr10=0;
$chr11=0;$chr12=0;$chr13=0;$chr14=0;$chr15=0;$chr16=0;$chr17=0;$chr18=0;$chr19=0;$chr20=0;$chr21=0;$chr22=0;$chrX=0;$chrY=0;$chrM=0;

@chrList = ("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
	"chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM");

for my $i (0..$#chrList){
	&chrListCount("$chrList[$i]");
}

# frequency of mutation considering chromosome length as present in chromosome_profile file.
print "\nCalculating mutation frequency and saving basechange-$inputRawSNV ...\n";

foreach $config (@config){
	@i = split (/\t/, $config);
	chomp (@i);
	push (@chrLength, "$i[1]");
}

for my $i (0..$#chrCountFreq){
	if ($chrLength[$i]==0){
		$chrNorm=0;
	}
	else {
		$chrNorm=$chrCountFreq[$i]/$chrLength[$i];
	}
	push @chrNormList, $chrNorm;
}

##### nucleotide type changes and frequency
print BASECHANGE "$inputRawSNV\tA>T\tA>G\tA>C\tC>G\tC>T\tC>A\n";
print BASECHANGE "Changes\t$ATTA\t$AGTC\t$ACTG\t$CGGC\t$CTGA\t$CAGT\n";
print BASECHANGE "Frequency\t$AVGATTA\t$AVGAGTC\t$AVGACTG\t$AVGCGGC\t$AVGCTGA\t$AVGCAGT\n";
print BASECHANGE "Type\tTransversion\tTransition\tTransversion\tTransversion\tTransition\tTransversion\n";
print BASECHANGE "\n";

print MUTFREQ "Chromosome\tTotalNumver\tMutPerBaseProfile\n";
print MUTFREQ "chr1\t$chrCountFreq[0]\t$chrNormList[0]\n";
print MUTFREQ "chr2\t$chrCountFreq[1]\t$chrNormList[1]\n";
print MUTFREQ "chr3\t$chrCountFreq[2]\t$chrNormList[2]\n";
print MUTFREQ "chr4\t$chrCountFreq[3]\t$chrNormList[3]\n";
print MUTFREQ "chr5\t$chrCountFreq[4]\t$chrNormList[4]\n";
print MUTFREQ "chr6\t$chrCountFreq[5]\t$chrNormList[5]\n";
print MUTFREQ "chr7\t$chrCountFreq[6]\t$chrNormList[6]\n";
print MUTFREQ "chr8\t$chrCountFreq[7]\t$chrNormList[7]\n";
print MUTFREQ "chr9\t$chrCountFreq[8]\t$chrNormList[8]\n";
print MUTFREQ "chr10\t$chrCountFreq[9]\t$chrNormList[9]\n";
print MUTFREQ "chr11\t$chrCountFreq[10]\t$chrNormList[10]\n";
print MUTFREQ "chr12\t$chrCountFreq[11]\t$chrNormList[11]\n";
print MUTFREQ "chr13\t$chrCountFreq[12]\t$chrNormList[12]\n";
print MUTFREQ "chr14\t$chrCountFreq[13]\t$chrNormList[13]\n";
print MUTFREQ "chr15\t$chrCountFreq[14]\t$chrNormList[14]\n";
print MUTFREQ "chr16\t$chrCountFreq[15]\t$chrNormList[15]\n";
print MUTFREQ "chr17\t$chrCountFreq[16]\t$chrNormList[16]\n";
print MUTFREQ "chr18\t$chrCountFreq[17]\t$chrNormList[17]\n";
print MUTFREQ "chr19\t$chrCountFreq[18]\t$chrNormList[18]\n";
print MUTFREQ "chr20\t$chrCountFreq[19]\t$chrNormList[19]\n";
print MUTFREQ "chr21\t$chrCountFreq[20]\t$chrNormList[20]\n";
print MUTFREQ "chr22\t$chrCountFreq[21]\t$chrNormList[21]\n";
print MUTFREQ "chrX\t$chrCountFreq[22]\t$chrNormList[22]\n";
print MUTFREQ "chrY\t$chrCountFreq[23]\t$chrNormList[23]\n";
print MUTFREQ "chrM\t$chrCountFreq[24]\t$chrNormList[24]\n";

##### hotspots 
print "\nCalculating hotspots and saving hotspots-$inputRawSNV.txt ...\n";

# one array of each chromosome.
for my $i (0..$#chrStart){

	if ($chrStart[$i] eq "chr1"){
		push @hsChr1, $pos[$i];
		push @geneClasschr1, $geneClass[$i];
		push @geneNamechr1, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr2"){
		push @hsChr2, $pos[$i];
		push @geneClasschr2, $geneClass[$i];
		push @geneNamechr2, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr3"){
		push @hsChr3, $pos[$i];
		push @geneClasschr3, $geneClass[$i];
		push @geneNamechr3, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr4"){
		push @hsChr4, $pos[$i];
		push @geneClasschr4, $geneClass[$i];
		push @geneNamechr4, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr5"){
		push @hsChr5, $pos[$i];
		push @geneClasschr5, $geneClass[$i];
		push @geneNamechr5, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr6"){
		push @hsChr6, $pos[$i];
		push @geneClasschr6, $geneClass[$i];
		push @geneNamechr6, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr7"){
		push @hsChr7, $pos[$i];
		push @geneClasschr7, $geneClass[$i];
		push @geneNamechr7, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr8"){
		push @hsChr8, $pos[$i];
		push @geneClasschr8, $geneClass[$i];
		push @geneNamechr8, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr9"){
		push @hsChr9, $pos[$i];
		push @geneClasschr9, $geneClass[$i];
		push @geneNamechr9, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr10"){
		push @hsChr10, $pos[$i];
		push @geneClasschr10, $geneClass[$i];
		push @geneNamechr10, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr11"){
		push @hsChr11, $pos[$i];
		push @geneClasschr11, $geneClass[$i];
		push @geneNamechr11, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr12"){
		push @hsChr12, $pos[$i];
		push @geneClasschr12, $geneClass[$i];
		push @geneNamechr12, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr13"){
		push @hsChr13, $pos[$i];
		push @geneClasschr13, $geneClass[$i];
		push @geneNamechr13, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr14"){
		push @hsChr14, $pos[$i];
		push @geneClasschr14, $geneClass[$i];
		push @geneNamechr14, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr15"){
		push @hsChr15, $pos[$i];
		push @geneClasschr15, $geneClass[$i];
		push @geneNamechr15, $geneName[$i];
	}
		if ($chrStart[$i] eq "chr16"){
		push @hsChr16, $pos[$i];
		push @geneClasschr16, $geneClass[$i];
		push @geneNamechr16, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr17"){
		push @hsChr17, $pos[$i];
		push @geneClasschr17, $geneClass[$i];
		push @geneNamechr17, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr18"){
		push @hsChr18, $pos[$i];
		push @geneClasschr18, $geneClass[$i];
		push @geneNamechr18, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr19"){
		push @hsChr19, $pos[$i];
		push @geneClasschr19, $geneClass[$i];
		push @geneNamechr19, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr20"){
		push @hsChr20, $pos[$i];
		push @geneClasschr20, $geneClass[$i];
		push @geneNamechr20, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr21"){
		push @hsChr21, $pos[$i];
		push @geneClasschr21, $geneClass[$i];
		push @geneNamechr21, $geneName[$i];
	}
	if ($chrStart[$i] eq "chr22"){
		push @hsChr22, $pos[$i];
		push @geneClasschr22, $geneClass[$i];
		push @geneNamechr22, $geneName[$i];
	}
	if ($chrStart[$i] eq "chrX"){
		push @hsChrX, $pos[$i];
		push @geneClasschrX, $geneClass[$i];
		push @geneNamechrX, $geneName[$i];
	}
	if ($chrStart[$i] eq "chrY"){
		push @hsChrY, $pos[$i];
		push @geneClasschrY, $geneClass[$i];
		push @geneNamechrY, $geneName[$i];
	}
	if ($chrStart[$i] eq "chrM"){
		push @hsChrM, $pos[$i];
		push @geneClasschrM, $geneClass[$i];
		push @geneNamechrM, $geneName[$i];
	}
}

# counting of hotspots using one array for each chromosome.
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

##### Saving hotspots output file
print HOTSPOTS "geneClass\tGeneName\tCHR\tBP\tHotspotCount\n";
for my $i (0 .. $#hsChr1){
	print HOTSPOTS "$geneClasschr1[$i]\t$geneNamechr1[$i]\t1\t$hsChr1[$i]\t$hs_counts_chr1[$i]\n";
}
for my $i (0 .. $#hsChr2){
	print HOTSPOTS "$geneClasschr2[$i]\t$geneNamechr2[$i]\t2\t$hsChr2[$i]\t$hs_counts_chr2[$i]\n";
}
for my $i (0 .. $#hsChr3){
	print HOTSPOTS "$geneClasschr3[$i]\t$geneNamechr3[$i]\t3\t$hsChr3[$i]\t$hs_counts_chr3[$i]\n";
}
for my $i (0 .. $#hsChr4){
	print HOTSPOTS "$geneClasschr4[$i]\t$geneNamechr4[$i]\t4\t$hsChr4[$i]\t$hs_counts_chr4[$i]\n";
}
for my $i (0 .. $#hsChr5){
	print HOTSPOTS "$geneClasschr5[$i]\t$geneNamechr5[$i]\t5\t$hsChr5[$i]\t$hs_counts_chr5[$i]\n";
}
for my $i (0 .. $#hsChr6){
	print HOTSPOTS "$geneClasschr6[$i]\t$geneNamechr6[$i]\t6\t$hsChr6[$i]\t$hs_counts_chr6[$i]\n";
}
for my $i (0 .. $#hsChr7){
	print HOTSPOTS "$geneClasschr7[$i]\t$geneNamechr7[$i]\t7\t$hsChr7[$i]\t$hs_counts_chr7[$i]\n";
}
for my $i (0 .. $#hsChr8){
	print HOTSPOTS "$geneClasschr8[$i]\t$geneNamechr8[$i]\t8\t$hsChr8[$i]\t$hs_counts_chr8[$i]\n";
}
for my $i (0 .. $#hsChr9){
	print HOTSPOTS "$geneClasschr9[$i]\t$geneNamechr9[$i]\t9\t$hsChr9[$i]\t$hs_counts_chr9[$i]\n";
}
for my $i (0 .. $#hsChr10){
	print HOTSPOTS "$geneClasschr10[$i]\t$geneNamechr10[$i]\t10\t$hsChr10[$i]\t$hs_counts_chr10[$i]\n";
}
for my $i (0 .. $#hsChr11){
	print HOTSPOTS "$geneClasschr11[$i]\t$geneNamechr11[$i]\t11\t$hsChr11[$i]\t$hs_counts_chr11[$i]\n";
}
for my $i (0 .. $#hsChr12){
	print HOTSPOTS "$geneClasschr12[$i]\t$geneNamechr12[$i]\t12\t$hsChr12[$i]\t$hs_counts_chr12[$i]\n";
}
for my $i (0 .. $#hsChr13){
	print HOTSPOTS "$geneClasschr13[$i]\t$geneNamechr13[$i]\t13\t$hsChr13[$i]\t$hs_counts_chr13[$i]\n";
}
for my $i (0 .. $#hsChr14){
	print HOTSPOTS "$geneClasschr14[$i]\t$geneNamechr14[$i]\t14\t$hsChr14[$i]\t$hs_counts_chr14[$i]\n";
}
for my $i (0 .. $#hsChr15){
	print HOTSPOTS "$geneClasschr15[$i]\t$geneNamechr15[$i]\t15\t$hsChr15[$i]\t$hs_counts_chr15[$i]\n";
}
for my $i (0 .. $#hsChr16){
	print HOTSPOTS "$geneClasschr16[$i]\t$geneNamechr16[$i]\t16\t$hsChr16[$i]\t$hs_counts_chr16[$i]\n";
}
for my $i (0 .. $#hsChr17){
	print HOTSPOTS "$geneClasschr17[$i]\t$geneNamechr17[$i]\t17\t$hsChr17[$i]\t$hs_counts_chr17[$i]\n";
}
for my $i (0 .. $#hsChr18){
	print HOTSPOTS "$geneClasschr18[$i]\t$geneNamechr18[$i]\t18\t$hsChr18[$i]\t$hs_counts_chr18[$i]\n";
}
for my $i (0 .. $#hsChr19){
	print HOTSPOTS "$geneClasschr19[$i]\t$geneNamechr19[$i]\t19\t$hsChr19[$i]\t$hs_counts_chr19[$i]\n";
}
for my $i (0 .. $#hsChr20){
	print HOTSPOTS "$geneClasschr20[$i]\t$geneNamechr20[$i]\t20\t$hsChr20[$i]\t$hs_counts_chr20[$i]\n";
}
for my $i (0 .. $#hsChr21){
	print HOTSPOTS "$geneClasschr21[$i]\t$geneNamechr21[$i]\t21\t$hsChr21[$i]\t$hs_counts_chr21[$i]\n";
}
for my $i (0 .. $#hsChr22){
	print HOTSPOTS "$geneClasschr22[$i]\t$geneNamechr22[$i]\t22\t$hsChr22[$i]\t$hs_counts_chr22[$i]\n";
}
for my $i (0 .. $#hsChrX){
	print HOTSPOTS "$geneClasschrX[$i]\t$geneNamechrX[$i]\t23\t$hsChrX[$i]\t$hs_counts_chrX[$i]\n";
}
for my $i (0 .. $#hsChrY){
	print HOTSPOTS "$geneClasschrY[$i]\t$geneNamechrY[$i]\t24\t$hsChrY[$i]\t$hs_counts_chrY[$i]\n";
}
for my $i (0 .. $#hsChrM){
	print HOTSPOTS "$geneClasschrM[$i]\t$geneNamechrM[$i]\t25\t$hsChrM[$i]\t$hs_counts_chrM[$i]\n";
}

### bioperl module for extraction of context sequences in reference genomes of each SNP.
print "\nExtracting context sequences and saving weblogo-$inputRawSNV ..\n";

$db = Bio::DB::Fasta->new("genomes/genome_$genome_version.fa");

for my $i (0..$#hsChr1){
	my $i2=$hsChr1[$i]-$contextSequenceLength;
	my $i3=$hsChr1[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr1', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr1_$hsChr1[$i]\n$seq\n";
}

for my $i (0..$#hsChr2){
	my $i2=$hsChr2[$i]-$contextSequenceLength;
	my $i3=$hsChr2[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr2', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr2_$hsChr2[$i]\n$seq\n";
}

for my $i (0..$#hsChr3){
	my $i2=$hsChr3[$i]-$contextSequenceLength;
	my $i3=$hsChr3[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr3', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr3_$hsChr3[$i]\n$seq\n";
}

for my $i (0..$#hsChr4){
	my $i2=$hsChr4[$i]-$contextSequenceLength;
	my $i3=$hsChr4[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr4', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr4_$hsChr4[$i]\n$seq\n";
}

for my $i (0..$#hsChr5){
	my $i2=$hsChr5[$i]-$contextSequenceLength;
	my $i3=$hsChr5[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr5', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr5_$hsChr5[$i]\n$seq\n";
}

for my $i (0..$#hsChr6){
	my $i2=$hsChr6[$i]-$contextSequenceLength;
	my $i3=$hsChr6[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr6', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr6_$hsChr6[$i]\n$seq\n";
}

for my $i (0..$#hsChr7){
	my $i2=$hsChr7[$i]-$contextSequenceLength;
	my $i3=$hsChr7[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr7', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr7_$hsChr7[$i]\n$seq\n";
}

for my $i (0..$#hsChr8){
	my $i2=$hsChr8[$i]-$contextSequenceLength;
	my $i3=$hsChr8[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr8', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr8_$hsChr8[$i]\n$seq\n";
}

for my $i (0..$#hsChr9){
	my $i2=$hsChr9[$i]-$contextSequenceLength;
	my $i3=$hsChr9[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr9', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr9_$hsChr9[$i]\n$seq\n";
}

for my $i (0..$#hsChr10){
	my $i2=$hsChr10[$i]-$contextSequenceLength;
	my $i3=$hsChr10[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr10', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr10_$hsChr10[$i]\n$seq\n";
}

for my $i (0..$#hsChr11){
	my $i2=$hsChr11[$i]-$contextSequenceLength;
	my $i3=$hsChr11[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr11', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr11_$hsChr11[$i]\n$seq\n";
}

for my $i (0..$#hsChr12){
	my $i2=$hsChr12[$i]-$contextSequenceLength;
	my $i3=$hsChr12[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr12', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr12_$hsChr1[$i]\n$seq\n";
}

for my $i (0..$#hsChr13){
	my $i2=$hsChr13[$i]-$contextSequenceLength;
	my $i3=$hsChr13[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr13', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr13_$hsChr13[$i]\n$seq\n";
}

for my $i (0..$#hsChr14){
	my $i2=$hsChr14[$i]-$contextSequenceLength;
	my $i3=$hsChr14[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr14', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr14_$hsChr14[$i]\n$seq\n";
}
for my $i (0..$#hsChr15){
	my $i2=$hsChr15[$i]-$contextSequenceLength;
	my $i3=$hsChr15[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr15', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr15_$hsChr15[$i]\n$seq\n";
}
for my $i (0..$#hsChr16){
	my $i2=$hsChr16[$i]-$contextSequenceLength;
	my $i3=$hsChr16[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr16', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr16_$hsChr16[$i]\n$seq\n";
}
for my $i (0..$#hsChr17){
	my $i2=$hsChr17[$i]-$contextSequenceLength;
	my $i3=$hsChr17[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr17', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr17_$hsChr17[$i]\n$seq\n";
}
for my $i (0..$#hsChr18){
	my $i2=$hsChr18[$i]-$contextSequenceLength;
	my $i3=$hsChr18[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr18', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr18_$hsChr18[$i]\n$seq\n";
}
for my $i (0..$#hsChr19){
	my $i2=$hsChr19[$i]-$contextSequenceLength;
	my $i3=$hsChr19[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr19', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr19_$hsChr19[$i]\n$seq\n";
}
for my $i (0..$#hsChr20){
	my $i2=$hsChr20[$i]-$contextSequenceLength;
	my $i3=$hsChr20[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr20', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr20_$hsChr20[$i]\n$seq\n";
}
for my $i (0..$#hsChr21){
	my $i2=$hsChr21[$i]-$contextSequenceLength;
	my $i3=$hsChr21[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr21', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr21_$hsChr21[$i]\n$seq\n";
}
for my $i (0..$#hsChr22){
	my $i2=$hsChr22[$i]-$contextSequenceLength;
	my $i3=$hsChr22[$i]+$contextSequenceLength;
	my $seq = $db->seq('chr22', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chr22_$hsChr22[$i]\n$seq\n";
}
for my $i (0..$#hsChrX){
	my $i2=$hsChrX[$i]-$contextSequenceLength;
	my $i3=$hsChrX[$i]+$contextSequenceLength;
	my $seq = $db->seq('chrX', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chrX_$hsChrX[$i]\n$seq\n";
}
for my $i (0..$#hsChrY){
	my $i2=$hsChrY[$i]-$contextSequenceLength;
	my $i3=$hsChrY[$i]+$contextSequenceLength;
	my $seq = $db->seq('chrY', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chrY_$hsChrY[$i]\n$seq\n";
}
for my $i (0..$#hsChrM){
	my $i2=$hsChrM[$i]-$contextSequenceLength;
	my $i3=$hsChrM[$i]+$contextSequenceLength;
	my $seq = $db->seq('chrM', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQ ">chrM_$hsChrM[$i]\n$seq\n";
}

print "\nExtracting context sequences and saving extracted_sequences-$inputRawSNV ..\n";

$db1 = Bio::DB::Fasta->new("genomes/genome_$genome_version.fa");

for my $i (0..$#hsChr1){
	my $i2=$hsChr1[$i]-$contextSequenceLength;
	my $i3=$hsChr1[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr1', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr1_$hsChr1[$i]\n$seq\t$geneClasschr1[$i]\t$geneNamechr1[$i]\n";
}

for my $i (0..$#hsChr2){
	my $i2=$hsChr2[$i]-$contextSequenceLength;
	my $i3=$hsChr2[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr2', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr2_$hsChr2[$i]\n$seq\t$geneClasschr2[$i]\t$geneNamechr2[$i]\n";
}

for my $i (0..$#hsChr3){
	my $i2=$hsChr3[$i]-$contextSequenceLength;
	my $i3=$hsChr3[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr3', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr3_$hsChr3[$i]\n$seq\t$geneClasschr3[$i]\t$geneNamechr3[$i]\n";
}

for my $i (0..$#hsChr4){
	my $i2=$hsChr4[$i]-$contextSequenceLength;
	my $i3=$hsChr4[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr4', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr4_$hsChr4[$i]\n$seq\t$geneClasschr4[$i]\t$geneNamechr4[$i]\n";
}

for my $i (0..$#hsChr5){
	my $i2=$hsChr5[$i]-$contextSequenceLength;
	my $i3=$hsChr5[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr5', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr5_$hsChr5[$i]\n$seq\t$geneClasschr5[$i]\t$geneNamechr5[$i]\n";
}

for my $i (0..$#hsChr6){
	my $i2=$hsChr6[$i]-$contextSequenceLength;
	my $i3=$hsChr6[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr6', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr6_$hsChr6[$i]\n$seq\t$geneClasschr6[$i]\t$geneNamechr6[$i]\n";
}

for my $i (0..$#hsChr7){
	my $i2=$hsChr7[$i]-$contextSequenceLength;
	my $i3=$hsChr7[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr7', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr7_$hsChr7[$i]\n$seq\t$geneClasschr7[$i]\t$geneNamechr7[$i]\n";
}

for my $i (0..$#hsChr8){
	my $i2=$hsChr8[$i]-$contextSequenceLength;
	my $i3=$hsChr8[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr8', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr8_$hsChr8[$i]\n$seq\t$geneClasschr8[$i]\t$geneNamechr8[$i]\n";
}

for my $i (0..$#hsChr9){
	my $i2=$hsChr9[$i]-$contextSequenceLength;
	my $i3=$hsChr9[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr9', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr9_$hsChr9[$i]\n$seq\t$geneClasschr9[$i]\t$geneNamechr9[$i]\n";
}

for my $i (0..$#hsChr10){
	my $i2=$hsChr10[$i]-$contextSequenceLength;
	my $i3=$hsChr10[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr10', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr10_$hsChr10[$i]\n$seq\t$geneClasschr10[$i]\t$geneNamechr10[$i]\n";
}

for my $i (0..$#hsChr11){
	my $i2=$hsChr11[$i]-$contextSequenceLength;
	my $i3=$hsChr11[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr11', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr11_$hsChr11[$i]\n$seq\t$geneClasschr11[$i]\t$geneNamechr11[$i]\n";
}

for my $i (0..$#hsChr12){
	my $i2=$hsChr12[$i]-$contextSequenceLength;
	my $i3=$hsChr12[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr12', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr12_$hsChr1[$i]\n$seq\t$geneClasschr12[$i]\t$geneNamechr12[$i]\n";
}

for my $i (0..$#hsChr13){
	my $i2=$hsChr13[$i]-$contextSequenceLength;
	my $i3=$hsChr13[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr13', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr13_$hsChr13[$i]\n$seq\t$geneClasschr13[$i]\t$geneNamechr13[$i]\n";
}

for my $i (0..$#hsChr14){
	my $i2=$hsChr14[$i]-$contextSequenceLength;
	my $i3=$hsChr14[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr14', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr14_$hsChr14[$i]\n$seq\t$geneClasschr14[$i]\t$geneNamechr14[$i]\n";
}
for my $i (0..$#hsChr15){
	my $i2=$hsChr15[$i]-$contextSequenceLength;
	my $i3=$hsChr15[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr15', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr15_$hsChr15[$i]\n$seq\t$geneClasschr15[$i]\t$geneNamechr15[$i]\n";
}
for my $i (0..$#hsChr16){
	my $i2=$hsChr16[$i]-$contextSequenceLength;
	my $i3=$hsChr16[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr16', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr16_$hsChr16[$i]\n$seq\t$geneClasschr16[$i]\t$geneNamechr16[$i]\n";
}
for my $i (0..$#hsChr17){
	my $i2=$hsChr17[$i]-$contextSequenceLength;
	my $i3=$hsChr17[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr17', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr17_$hsChr17[$i]\n$seq\t$geneClasschr17[$i]\t$geneNamechr17[$i]\n";
}
for my $i (0..$#hsChr18){
	my $i2=$hsChr18[$i]-$contextSequenceLength;
	my $i3=$hsChr18[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr18', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr18_$hsChr18[$i]\n$seq\t$geneClasschr18[$i]\t$geneNamechr18[$i]\n";
}
for my $i (0..$#hsChr19){
	my $i2=$hsChr19[$i]-$contextSequenceLength;
	my $i3=$hsChr19[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr19', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr19_$hsChr19[$i]\n$seq\t$geneClasschr19[$i]\t$geneNamechr19[$i]\n";
}
for my $i (0..$#hsChr20){
	my $i2=$hsChr20[$i]-$contextSequenceLength;
	my $i3=$hsChr20[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr20', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr20_$hsChr20[$i]\n$seq\t$geneClasschr20[$i]\t$geneNamechr20[$i]\n";
}
for my $i (0..$#hsChr21){
	my $i2=$hsChr21[$i]-$contextSequenceLength;
	my $i3=$hsChr21[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr21', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr21_$hsChr21[$i]\n$seq\t$geneClasschr21[$i]\t$geneNamechr21[$i]\n";
}
for my $i (0..$#hsChr22){
	my $i2=$hsChr22[$i]-$contextSequenceLength;
	my $i3=$hsChr22[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chr22', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chr22_$hsChr22[$i]\n$seq\t$geneClasschr22[$i]\t$geneNamechr22[$i]\n";
}
for my $i (0..$#hsChrX){
	my $i2=$hsChrX[$i]-$contextSequenceLength;
	my $i3=$hsChrX[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chrX', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chrX_$hsChrX[$i]\n$seq\t$geneClasschrX[$i]\t$geneNamechrX[$i]\n";
}
for my $i (0..$#hsChrY){
	my $i2=$hsChrY[$i]-$contextSequenceLength;
	my $i3=$hsChrY[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chrY', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chrY_$hsChrY[$i]\n$seq\t$geneClasschrY[$i]\t$geneNamechrY[$i]\n";
}
for my $i (0..$#hsChrM){
	my $i2=$hsChrM[$i]-$contextSequenceLength;
	my $i3=$hsChrM[$i]+$contextSequenceLength;
	my $seq = $db1->seq('chrM', $i2 => $i3);
	$i2=0;
	$i3=0;
	print CONTEXTSEQANNO ">chrM_$hsChrM[$i]\n$seq\t$geneClasschrM[$i]\t$geneNamechrM[$i]\n";
}

## closing of outputs
close (">>results-$inputRawSNV/WOLAND-basechange-$inputRawSNV");
close (">>results-$inputRawSNV/WOLAND-mutfreq-$inputRawSNV");
close (">>results-$inputRawSNV/WOLAND-contextsequences-$inputRawSNV");
close (">>results-$inputRawSNV/WOLAND-hotspots-$inputRawSNV");
close (">>results-$inputRawSNV/WOLAND-log-$inputRawSNV");
close (">>results-$inputRawSNV/WOLAND-contextsequencesanno-$inputRawSNV");


### motif search
# outputs
open (OUTPUTbiasSN1, ">>results-$inputRawSNV/WOLAND-bias_SN1-$inputRawSNV");
open (OUTPUTbiasDNApoln, ">>results-$inputRawSNV/WOLAND-bias_DNApoln-$inputRawSNV");
open (OUTPUTbiasoxoG, ">>results-$inputRawSNV/WOLAND-bias_oxoG-$inputRawSNV");
open (OUTPUTbiasUV, ">>results-$inputRawSNV/WOLAND-bias_UV-lambda-$inputRawSNV");
open (OUTPUTbiasUVsolar, ">>results-$inputRawSNV/WOLAND-bias_UVsolar-$inputRawSNV");
open (OUTPUTbiassixfour, ">>results-$inputRawSNV/WOLAND-bias_sixfour-$inputRawSNV");

# opening context sequences generated
open CONTEXTSEQANNO, "<results-$inputRawSNV/WOLAND-contextsequencesanno-$inputRawSNV";
@fastaContext=<CONTEXTSEQANNO>;

# array processing
for my $i (0..$#fastaContext){
	if ($fastaContext[$i] =~ /^>/){
		@i2 = split (/\t/, "$fastaContext[$i+1]");
		@i = split (/\n/, "$fastaContext[$i]");
		chomp(@i2);
		push (@fastaContextMS, "$i[0]\n$i2[0]\n");
		push (@geneClassMS, "$i2[1]");
		push (@geneNameMS, "$i2[2]");
		@i2=();
		@i=();
	}
}

foreach $geneClassMS (@geneClassMS){
	if ($geneClassMS =~ /[\S]/){
		push (@geneClass1, $geneClassMS);
	}
}

foreach $geneNameMS (@geneNameMS){
	if ($geneNameMS =~ /[\S]/){
		push (@geneName1, $geneNameMS);
	}
}

foreach $fastaContextMS (@fastaContextMS){
	@i = split (/\n/, "$fastaContextMS");
	push (@id, $i[0]);
	push (@targetSequence, $i[1]);
}

foreach $id (@id){
	@i = split (/_/, $id);
	chomp (@i);
	push (@chrRaw, "$i[0]");
	push (@coord, "$i[1]");
}

foreach $chrRaw (@chrRaw){
	@i = split (/>/, $chrRaw);
	chomp (@i);
	push (@chrSt, "$i[1]");
}

### SC score (strand bias)
$refSeqRaw = "genomes/refseq_$genome_version.txt";
open (REFSEQRAW, $refSeqRaw);
@refSeq=<REFSEQRAW>;
unless (@refSeq){
	die "\nERROR : Could not load <refseq_$genome_version.txt>\n";
}

foreach $refSeqline (@refSeq){
	@irefSeq = split (/\t/, $refSeqline);
	chomp (@irefSeq);
	push (@chrRefSeq, "$irefSeq[2]");
	push (@strandRefSeq, "$irefSeq[3]");
	push (@txstartRefSeq, $irefSeq[4]);
	push (@txtendRefSeq, "$irefSeq[5]");
}


# search for motifs
print "\nSearching for mutable motifs\n";

## SN1 MOTIF
foreach $targetSequence (@targetSequence){
	if ($targetSequence =~ "..AG..." || $targetSequence =~ "..GG..."|| $targetSequence =~ "...CT.."|| $targetSequence =~ "...CC.."){
		push (@SN1, "1");
	}
	else {
		push (@SN1, "0");
	}
	if ($targetSequence =~ "..AG..." || $targetSequence =~ "..GG..."){
		push (@SN1plus, "1");
	}
	else {
		push (@SN1plus, "0")
	}
	if ($targetSequence =~ "...CT.."|| $targetSequence =~ "...CC.."){
		push (@SN1minus, "1");
	}
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
		for my $i (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_count;
			}
			if ($strandRefSeq[$i] eq "+" && $chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_plus_count;
			} 
		}
		$i=0;
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
		push (@DNApoln, "1");
	}
	else {
		push (@DNApoln, "0");
	}
	if ($targetSequence =~ "..AA..." || $targetSequence =~ "..TA..."){
		push (@DNApolnplus, "1");
	}
	else {
		push (@DNApolnplus, "0")
	}
	if ($targetSequence =~ "...TT.."|| $targetSequence =~ "...TA.."){
		push (@DNApolnminus, "1");
	}
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
	for my $i (0 .. $#chrRefSeq){
		if ($chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
			++$strand_count;
		}
		if ($strandRefSeq[$i] eq "+" && $chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
			++$strand_plus_count;
		} 
	}
	$i=0;
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
		for my $i (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_count;
			}
			if ($strandRefSeq[$i] eq "+" && $chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_plus_count;
			} 
		}
		$i=0;
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
	if($targetSequence =~ "...TC.." || $targetSequence =~ "..TC..." || $targetSequence =~ "...CT.." || $targetSequence =~ "..CT..." ||
		$targetSequence =~ "..TT..." || $targetSequence =~ "...TT.." || $targetSequence =~ "..CC..." || $targetSequence =~ "...CC.."|| 
		$targetSequence =~ "...GA.." || $targetSequence =~ "..GA..." || $targetSequence =~ "...AA.." || $targetSequence =~ "..AA..." ||
		$targetSequence =~ "...GG.." || $targetSequence =~ "..GG..." || $targetSequence =~ "..AG..." || $targetSequence =~ "...AG.."){
		push (@UVlambda, "1");
	}
	else{
		push (@UVlambda, "0");
	}
	if($targetSequence =~ "...TC.." || $targetSequence =~ "..TC..." || $targetSequence =~ "...CT.." || $targetSequence =~ "..CT..." ||
		$targetSequence =~ "..TT..." || $targetSequence =~ "...TT.." || $targetSequence =~ "..CC..." || $targetSequence =~ "...CC.."){
		push (@UVlambdaplus, "1");
	}
	else{
		push (@UVlambdaplus, "0")
	}
	if ($targetSequence =~ "...GA.." || $targetSequence =~ "..GA..." || $targetSequence =~ "...AA.." || $targetSequence =~ "..AA..." ||
		$targetSequence =~ "...GG.." || $targetSequence =~ "..GG..." || $targetSequence =~ "..AG..." || $targetSequence =~ "...AG.."){
		push (@UVlambdaminus, "1");
	}
	else{
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
		for my $i (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_count;
			}
			if ($strandRefSeq[$i] eq "+" && $chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_plus_count;
			} 
		}
		$i=0;
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
		push (@UVsolar, "1");
	}
	else{
		push (@UVsolar, "0");
	}
	if ($targetSequence =~ "..TCG.."  || $targetSequence =~ "...TCG."|| $targetSequence =~ ".TCG..." || 
		    $targetSequence =~ "..CCG.." ||$targetSequence =~ "...CCG." || $targetSequence =~ ".CCG..."){
		push (@UVsolarplus, "1");
	}
	else {
		push (@UVsolarplus, "0")
	}
		if ($targetSequence =~ "..CGA.." || $targetSequence =~ "...CGA."|| $targetSequence =~ ".CGA..." ||
			$targetSequence =~ "..CGG.." || $targetSequence =~ "...CGG."|| $targetSequence =~ ".CGG...") {
		push (@UVsolarminus, "1");
	}
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
		for my $i (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_count;
			}
			if ($strandRefSeq[$i] eq "+" && $chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_plus_count;
			} 
		}
		$i=0;
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
		push (@sixfour, "1");
	}
	else {
		push (@sixfour, "0");
	}
		if ($targetSequence =~ "TTCA" || $targetSequence =~ "CTCA"){
		push (@sixfourplus, "1");
	}
	else {
		push (@sixfourplus, "0")
	}
		if ($targetSequence =~"TGAA" || $targetSequence =~ "TGAG") {
		push (@sixfourminus, "1");
	}
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
		for my $i (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_count;
			}
			if ($strandRefSeq[$i] eq "+" && $chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_plus_count;
			} 
		}
		$i=0;
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
		push (@enu, "1");
	}
	else {
		push (@enu, "0");
	}
			if ($targetSequence =~ "..CAG.." || $targetSequence =~ "..GAC.."|| $targetSequence =~ "..GAG.."|| $targetSequence =~ "..CAC.."){
		push (@enuplus, "1");
	}
	else {
		push (@enuplus, "0")
	}
		if ($targetSequence =~ "..CTG.." || $targetSequence =~ "..GTC.."|| $targetSequence =~ "..GTG.."|| $targetSequence =~ "..CTC..") {
		push (@enuminus, "1");
	}
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
		for my $i (0 .. $#chrRefSeq){
			if ($chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_count;
			}
			if ($strandRefSeq[$i] eq "+" && $chrRefSeq[$i] eq $query_chr and (($query_coord <= $txstartRefSeq[$i] and $query_coord >= $txtendRefSeq[$i]) or ($query_coord >= $txstartRefSeq[$i] and $query_coord <= $txtendRefSeq[$i]))){
				++$strand_plus_count;
			} 
		}
		$i=0;
		if ($strand_count >= 1){
			$strand_transcript=$strand_plus_count/$strand_count;
			$strand_score=$strand_transcript - ($strand_value);
			#print OUTPUTbiasenu "$chrSt[$ienu]\t$coord[$ienu]\t$strand_score\n";
			$strand_transcript=0;
			$strand_plus_count=0;
   			$strand_count=0;
			$strand_score=0;
		}
	}
}

## calculation of normalized motifs
for my $i (0 .. $#id){
	++$SNPnumber;
}

$SN1counts=0;
foreach $SN1 (@SN1){
	if ($SN1 eq "1"){
		++$SN1counts;
	}
}

$DNApolncounts=0;
foreach $DNApoln (@DNApoln){
	if ($DNApoln eq "1"){
		++$DNApolncounts;
	}
}

$oxoGcounts=0;
foreach $oxoG (@oxoG){
	if ($oxoG eq "1"){
		++$oxoGcounts;
	}
}

$UVlambdacounts=0;
foreach $UVlambda (@UVlambda){
	if ($UVlambda eq "1"){
		++$UVlambdacounts;
	}
}

$sixfourcounts=0;
foreach $sixfour (@sixfour){
	if ($sixfour eq "1"){
		++$sixfourcounts;
	}
}

$enucounts=0;
foreach $enu (@enu){
	if ($enu eq "1"){
		++$enucounts;
	}
}

$UVsolarcounts=0;
foreach $UVsolar (@UVsolar){
	if ($UVsolar eq "1"){
		++$UVsolarcounts;
	}
}

$normSN1=$SN1counts/$SNPnumber;
$normDNApoln=$DNApolncounts/$SNPnumber;
$normoxoG=$oxoGcounts/$SNPnumber;
$normUVlambda=$UVlambdacounts/$SNPnumber;
$normsixfour=$sixfourcounts/$SNPnumber;
$normenu=$enucounts/$SNPnumber;
$normUVsolar=$UVsolarcounts/$SNPnumber;

#### output for motif number
print "\nSaving mutational motifs file\n";

print MOTIFS "Chr\tPos\ttargetSequence\tClass\tGene\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUVA-solar\n";
for $i (0..$#targetSequence){
	print MOTIFS "$chrSt[$i]\t$coord[$i]\t$targetSequence[$i]\t$geneClass1[$i]\t$geneName1[$i]\t$SN1[$i]\t$DNApoln[$i]\t$oxoG[$i]\t$UVlambda[$i]\t$sixfour[$i]\t$enu[$i]\t$UVsolar[$i]\n";
}

#### output for motif number normalized
print "\nSaving normalized mutational motifs\n";
print NMOTIFS "$inputRawSNV\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUVA-solar\n";
print NMOTIFS "Number of Total SNPs\t$SNPnumber\t$SNPnumber\t$SNPnumber\t$SNPnumber\t$SNPnumber\t$SNPnumber\t$SNPnumber\n";
print NMOTIFS "Total Raw Number of Motifs Found\t$SN1counts\t$DNApolncounts\t$oxoGcounts\t$UVlambdacounts\t$sixfourcounts\t$enucounts\t$UVsolarcounts\n";
print NMOTIFS "Normalized Number of Motifs Found\t$normSN1\t$normDNApoln\t$normoxoG\t$normUVlambda\t$normsixfour\t$normenu\t$normUVsolar";

exit;