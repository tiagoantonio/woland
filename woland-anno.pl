########################################################################################################################################
## WOLAND Beta 0.2 (01-28-2017)
## woland-anno.pl
##
## WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. .
## 
## For more details please read README file.
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
###########################################################################################################################################

#! /usr/bin/perl
use Bio::DB::Fasta; # bioperl module for the extraction of sequences.
use Cwd;
use warnings;
use strict;
use List::MoreUtils qw(uniq);

our $REVISION = '$Revision:  $';
our $DATE =	'$Date: 2017-01-28 00:11:04 -0800 (Sat,  28 Jan 2017) $';  
our $AUTHOR =	'$Author: Tiago A. de Souza <tiagoantonio@gmail.com> $';

# global variables
our ($variantfile,$chrlengthprofile,$hotspotwindowlength,$genomeversion,$datestring,$contextsequencelength, $genomesequenceinstance);
our (@geneclass,@genename,@chrstart,@chrend,@pos,@ref,@alt,@refalt,@chrnames,@chrpos,@chrlength,@chrnormalizedlist, @chrcountfreq,@hotspotsdetected,@chrnameformatted);
our (@extractedsequenceformotifsearch,@chrrefseq, @strandrefseq,@txstartrefseq,@txtendrefseq,@uniquechrinsamples);

# subroutines
sub parse_inputs { #parsing .variant.function file
	$variantfile=$_[0];
	open (VARIANTFILE, $variantfile);
	my @rawvariantarray=<VARIANTFILE>;
	close (VARIANTFILE);
	foreach my $rawvariantarrayline (@rawvariantarray){ #of each column in a dedicated array
		my @i = split (/\t/, $rawvariantarrayline);
		chomp (@i);
		push (@geneclass, "$i[0]"); push (@genename, "$i[1]");#gene clss and gene names
		push (@chrstart, "$i[2]"); push (@chrend, "$i[3]"); #chromosome @chrstart
		push (@pos, $i[4]); #position 
		push (@ref, "$i[5]"); push (@alt, "$i[6]"); #ref and alt alleles
	}
	for my $i (0 .. $#chrstart){ # Array for chr/position of each SNP - @chrpos
		push @chrpos, "$chrstart[$i]_$pos[$i]";
	}
	my %posref;
	for my $i (0 .. $#chrpos){ # Hashe posref & posalt for information of each position - ALT e REF alelles
		$posref{"$chrpos[$i]_$i"} .= "$ref[$i]";
	}
	my %posalt;
	for my $i (0 .. $#chrpos){
		$posalt{"$chrpos[$i]_$i"} .= "$alt[$i]";
	}
	foreach my $key (sort keys %posref){ # Array for a single entry for each SNP containing a single REFALT string value. Considering  A->T # & T->A; A->G & T->C; A->C & T->G; C->G & G->C; C->T & G->A; C->A & G->T.
		push (@refalt,"$posref{$key}$posalt{$key}");
	}
	$chrlengthprofile=$_[1]; #parsing chromosome profile file
	open (CHRLENGTHPROFILEFILE, $chrlengthprofile);
	my @chrprofilearray=<CHRLENGTHPROFILEFILE>;
	close (CHRLENGTHPROFILEFILE);
	my $chrprofilearrayline;
	foreach $chrprofilearrayline (@chrprofilearray){
		my @i = split (/\t/, $chrprofilearrayline);
		chomp (@i);
		push (@chrnames, "$i[0]");
		push (@chrlength, "$i[1]");
	}

	$hotspotwindowlength=$_[2]; #parsing hotspot window length
	unless ($hotspotwindowlength>0){
		die "\nERROR : Please specify a natural number >0 for <hotspot_window_length>\n";
	}
	
	$genomeversion = $_[3]; #parsing genome version
	my $fastagenomeversion="genomes/genome_$genomeversion.fa";
	unless ($genomeversion){
		die "\nERROR : Please specify a genome version in genomes\/folder for <genome_version>\n";
	}
	unless (-r -e -f $fastagenomeversion){
		die "\nERROR : Please check if a genome fasta file exists in genomes\/folder for <genome_version>\n"
	}
}

sub process_refseq{ #parsing refseq info for transcriptional strand concordance
	my $refseqdata = "genomes/refseq_$genomeversion.txt";
	my @refseqprocessedarrayline;
	my $refseqarrayline;
	unless (-r -e -f $refseqdata){
		die "\nERROR : Could not load <refseq_$genomeversion.txt>. Please review RefSeq annotation.\n";
	}
	open (REFSEQDATA, $refseqdata);
	my @refseqrawarray=<REFSEQDATA>;
	close (REFSEQDATA);
	
	foreach $refseqarrayline (@refseqrawarray){
		my @refseqprocessedarrayline = split (/\t/, $refseqarrayline);
		chomp (@refseqprocessedarrayline);
		push (@chrrefseq, "$refseqprocessedarrayline[2]");
		push (@strandrefseq, "$refseqprocessedarrayline[3]");
		push (@txstartrefseq, $refseqprocessedarrayline[4]);
		push (@txtendrefseq, "$refseqprocessedarrayline[5]");
	}
}

sub count_nucleotidechanges{ #nucleotide-type changes, frequency,
	my $ATTA=0;my $AGTC=0;my $ACTG=0; my $CGGC=0; my $CTGA=0; my $CAGT=0;
	my $refaltline;
	foreach $refaltline (@refalt){
		if ($refaltline eq "AT" || $refaltline eq "TA"){
			$ATTA++;
		}
		if ($refaltline eq "AG" || $refaltline eq "TC"){
			$AGTC++;
		}
		if ($refaltline eq "AC" || $refaltline eq "TG"){
			$ACTG++;
		}
		if ($refaltline eq "CG" || $refaltline eq "GC"){
			$CGGC++;
		}
		if ($refaltline eq "CT" || $refaltline eq "GA"){
			$CTGA++;
		}
		if ($refaltline eq "CA" || $refaltline eq "GT"){
			$CAGT++;
		}
	}
	my $validtotalsnvs=$ATTA+$AGTC+$ACTG+$CGGC+$CTGA+$CAGT; # frequency of changes considering total amount of changes.
	if ($validtotalsnvs== 0){
		die "Please review input file format";
	}
	my $allsnvs;
	foreach $refaltline(@refalt){
		$allsnvs++;
	}

	print "\nTotal SNPs valid=$validtotalsnvs\n";
	my $notvalidsnvs=$allsnvs-$validtotalsnvs;
	print "Total not valid SNPs=$notvalidsnvs\n";

	my $averageATTA=$ATTA/$validtotalsnvs; my $averageAGTC=$AGTC/$validtotalsnvs;
	my $averageACTG=$ACTG/$validtotalsnvs; my $averageCGGC=$CGGC/$validtotalsnvs;
	my $averageCTGA=$CTGA/$validtotalsnvs; my $averageCAGT=$CAGT/$validtotalsnvs;

	my $transitions=($AGTC+$CTGA)/$validtotalsnvs; # transitionss & transversions
	my $transversions=($ATTA+$ACTG+$CGGC+$CAGT)/$validtotalsnvs; # transitionss & transversions

	open (BASECHANGE, ">>results-$variantfile/WOLAND-basechange-$variantfile");
	print BASECHANGE "$variantfile\tA>T\tA>G\tA>C\tC>G\tC>T\tC>A\n";
	print BASECHANGE "Changes\t$ATTA\t$AGTC\t$ACTG\t$CGGC\t$CTGA\t$CAGT\n";
	print BASECHANGE "Frequency\t$averageATTA\t$averageAGTC\t$averageACTG\t$averageCGGC\t$averageCTGA\t$averageCAGT\n";
	print BASECHANGE "Type\tTransversion\tTransition\tTransversion\tTransversion\tTransition\tTransversion\n";
	print BASECHANGE "\n";
}

sub count_chrhits{ # number of mutations per chromosome target
	my $countchr=0;
	for my $ii (0 .. $#chrstart){
		if ($chrstart[$ii] eq "$_[0]"){
			++$countchr;
		}
	}
	push (@chrcountfreq, $countchr);
	$countchr=0;
}

sub count_chrfreq{ # frequency of mutations across chromosomes
	my $chrnormalized;
	for my $i (0..$#chrcountfreq){
		$chrnormalized=0;
		if ($chrcountfreq[$i]==0){
			$chrnormalized=0;
		}
		if ($chrcountfreq[$i]!=0) {
			$chrnormalized=$chrcountfreq[$i]/$chrlength[$i];
		}
		push (@chrnormalizedlist, $chrnormalized);
	}
	open (MUTFREQ, ">>results-$variantfile/WOLAND-mutfreq-$variantfile");
	print MUTFREQ "Chromosome\tTotalNumver\tMutPerBaseProfile\n";
	for my $i(0..$#chrnames){
		print MUTFREQ "$chrnames[$i]\t$chrcountfreq[$i]\t$chrnormalizedlist[$i]\n";
	}
}

sub search_for_hotspots{ #hotspot search
	for my $ii (0..$#chrstart){
		if ($chrstart[$ii] eq $_[0]){
			my $hotspotwindowdownstream=$pos[$ii]-$hotspotwindowlength;
			my $hotspotwindowupstream=$pos[$ii]+$hotspotwindowlength;
			my $counthotspot=0;
			for my $iii (0..$#chrstart){
				if ($chrstart[$iii] eq $_[0]){
					if ($pos[$iii]>=$hotspotwindowdownstream && $pos[$iii]<=$hotspotwindowupstream){
						++$counthotspot;
					}
				}
			}
			my $chrwithoutletters = substr($_[0],3);
			if ($chrwithoutletters eq "X"){
				$chrwithoutletters=23;
			}
			if ($chrwithoutletters eq "Y"){
				$chrwithoutletters=24;
			}
			if ($chrwithoutletters eq "M" || $chrwithoutletters eq "MT"){
				$chrwithoutletters=25;
			}
			if ($chrwithoutletters =~ /^[0-9]+$/){
			}
			else{
				$chrwithoutletters="30";
			}
			push (@hotspotsdetected, $counthotspot);
			push (@chrnameformatted, $chrwithoutletters);
		}
		else{

		}
	}
}
sub extract_contextsequences{ #context sequences n=3
	my $extractedsequence;
	if (grep(/^$_[0]$/, @chrnames)){
		my $downstreamcoordinate=$_[1]-$contextsequencelength;
		my $upstreamcoordinate=$_[1]+$contextsequencelength;
		$extractedsequence=$_[2]->seq($_[0], $downstreamcoordinate=>$upstreamcoordinate);
	}else{
		$extractedsequence="NNNNNNN";
	}
	push (@extractedsequenceformotifsearch, $extractedsequence);
	print CONTEXTSEQ ">$_[0]_$_[1]\n$extractedsequence\n";
	print CONTEXTSEQANNO ">$_[0]_$_[1]\n$extractedsequence\t$_[3]\t$_[4]\n";
}

sub motif_search{ #motif search counting, normalized, strand concordance
	my $extractedsequence;
	my (@sn1,@sn1plus,@sn1minus,@dnapoln,@dnapolnplus,@dnapolnminus, @oxog,@oxogplus,@oxogminus, @uvlambda,@uvlambdaplus,@uvlambdaminus,@sixfour,@sixfourplus,@sixfourminus,@enu,@enuplus,@enuminus,@uvsolar,@uvsolarplus,@uvsolarminus);
	#sn1 motif
	foreach $extractedsequence (@extractedsequenceformotifsearch){
		if ($extractedsequence =~ "..AG..." || $extractedsequence =~ "..GG..."|| $extractedsequence =~ "...CT.."|| $extractedsequence =~ "...CC.."){
			push (@sn1, "1");
		}else{
			push (@sn1, "0");
		}
		if ($extractedsequence =~ "..AG..." || $extractedsequence =~ "..GG..."){
			push (@sn1plus, "1");
		}
		else{
			push (@sn1plus, "0")
		}
		if ($extractedsequence =~ "...CT.."|| $extractedsequence =~ "...CC.."){
			push (@sn1minus, "1");
		}
		else{
			push (@sn1minus, "0")
		}

		#dna poln motif
		if ($extractedsequence =~ "..AA..." || $extractedsequence =~ "..TA..."|| $extractedsequence =~ "...TT.."|| $extractedsequence =~ "...TA.."){
			push (@dnapoln, "1");
		}
		else {
			push (@dnapoln, "0");
		}
		if ($extractedsequence =~ "..AA..." || $extractedsequence =~ "..TA..."){
			push (@dnapolnplus, "1");
		}
		else {
			push (@dnapolnplus, "0")
		}
		if ($extractedsequence =~ "...TT.."|| $extractedsequence =~ "...TA.."){
			push (@dnapolnminus, "1");
		}
		else {
			push (@dnapolnminus, "0")
		}

		#8-oxoG motif
		if ($extractedsequence =~ "..AGA.." || $extractedsequence =~ "..GGG.." || $extractedsequence =~ "..AGG.." || $extractedsequence =~ "..GGA.."|| 
			$extractedsequence =~ "..TCT.." || $extractedsequence =~ "..CCC.." || $extractedsequence =~ "..CCT.."|| $extractedsequence =~ "..TCC.."){
			push (@oxog, "1");}
		else {
			push (@oxog, "0");
		}
			if ($extractedsequence =~ "..AGA.." || $extractedsequence =~ "..GGG.." || $extractedsequence =~ "..AGG.." || $extractedsequence =~ "..GGA.."){
			push (@oxogplus, "1");}
		else {
			push (@oxogplus, "0")
		}
		if ($extractedsequence =~ "..TCT.." || $extractedsequence =~ "..CCC.." || $extractedsequence =~ "..CCT.."|| $extractedsequence =~ "..TCC.."){
			push (@oxogminus, "1");}
		else {
			push (@oxogminus, "0")
		}

		# uv-lambda motif
		if($extractedsequence =~ "...TC.." || $extractedsequence =~ "..TC..." || $extractedsequence =~ "...CT.." || $extractedsequence =~ "..CT..." ||
			$extractedsequence =~ "..TT..." || $extractedsequence =~ "...TT.." || $extractedsequence =~ "..CC..." || $extractedsequence =~ "...CC.."|| 
			$extractedsequence =~ "...GA.." || $extractedsequence =~ "..GA..." || $extractedsequence =~ "...AA.." || $extractedsequence =~ "..AA..." ||
			$extractedsequence =~ "...GG.." || $extractedsequence =~ "..GG..." || $extractedsequence =~ "..AG..." || $extractedsequence =~ "...AG.."){
			push (@uvlambda, "1");
		}
		else{
			push (@uvlambda, "0");
		}
		if($extractedsequence =~ "...TC.." || $extractedsequence =~ "..TC..." || $extractedsequence =~ "...CT.." || $extractedsequence =~ "..CT..." ||
			$extractedsequence =~ "..TT..." || $extractedsequence =~ "...TT.." || $extractedsequence =~ "..CC..." || $extractedsequence =~ "...CC.."){
			push (@uvlambdaplus, "1");
		}
		else{
			push (@uvlambdaplus, "0")
		}
		if ($extractedsequence =~ "...GA.." || $extractedsequence =~ "..GA..." || $extractedsequence =~ "...AA.." || $extractedsequence =~ "..AA..." ||
			$extractedsequence =~ "...GG.." || $extractedsequence =~ "..GG..." || $extractedsequence =~ "..AG..." || $extractedsequence =~ "...AG.."){
			push (@uvlambdaminus, "1");
		}
		else{
			push (@uvlambdaminus, "0")
		}

		#uv-solar motif
		if ($extractedsequence =~ "..TCG.."  || $extractedsequence =~ "...TCG."|| $extractedsequence =~ ".TCG..." || $extractedsequence =~ "..CCG.." ||$extractedsequence =~ "...CCG." || $extractedsequence =~ ".CCG..." || 
			$extractedsequence =~ "..CGA.." || $extractedsequence =~ "...CGA."|| $extractedsequence =~ ".CGA..." || $extractedsequence =~ "..CGG.." || $extractedsequence =~ "...CGG."|| $extractedsequence =~ ".CGG..."){
			push (@uvsolar, "1");
		}
		else{
			push (@uvsolar, "0");
		}
		if ($extractedsequence =~ "..TCG.."  || $extractedsequence =~ "...TCG."|| $extractedsequence =~ ".TCG..." || 
			    $extractedsequence =~ "..CCG.." ||$extractedsequence =~ "...CCG." || $extractedsequence =~ ".CCG..."){
			push (@uvsolarplus, "1");
		}
		else {
			push (@uvsolarplus, "0")
		}
			if ($extractedsequence =~ "..CGA.." || $extractedsequence =~ "...CGA."|| $extractedsequence =~ ".CGA..." ||
				$extractedsequence =~ "..CGG.." || $extractedsequence =~ "...CGG."|| $extractedsequence =~ ".CGG...") {
			push (@uvsolarminus, "1");
		}
		else {
			push (@uvsolarminus, "0")
		}

		#6-4 motif
		if ($extractedsequence =~ "TTCA" || $extractedsequence =~ "CTCA" || $extractedsequence =~ "TGAA" || $extractedsequence =~ "TGAG"){
			push (@sixfour, "1");
		}
		else {
			push (@sixfour, "0");
		}
			if ($extractedsequence =~ "TTCA" || $extractedsequence =~ "CTCA"){
			push (@sixfourplus, "1");
		}
		else {
			push (@sixfourplus, "0")
		}
			if ($extractedsequence =~"TGAA" || $extractedsequence =~ "TGAG") {
			push (@sixfourminus, "1");
		}
		else {
			push (@sixfourminus, "0")
		}

		#enu motif
		if ($extractedsequence =~ "..CAG.." || $extractedsequence =~ "..GAC.."|| $extractedsequence =~ "..GAG.."|| $extractedsequence =~ "..CAC.."|| $extractedsequence =~ "..CTG.." || $extractedsequence =~ "..GTC.."|| $extractedsequence =~ "..GTG.."|| $extractedsequence =~ "..CTC.."){
			push (@enu, "1");
		}
		else {
			push (@enu, "0");
		}
				if ($extractedsequence =~ "..CAG.." || $extractedsequence =~ "..GAC.."|| $extractedsequence =~ "..GAG.."|| $extractedsequence =~ "..CAC.."){
			push (@enuplus, "1");
		}
		else {
			push (@enuplus, "0")
		}
			if ($extractedsequence =~ "..CTG.." || $extractedsequence =~ "..GTC.."|| $extractedsequence =~ "..GTG.."|| $extractedsequence =~ "..CTC..") {
			push (@enuminus, "1");
		}
		else {
			push (@enuminus, "0")
		}
	}
	## calculation of normalized motifs
	my $totalsnvs=0;
	for my $i (0 .. $#chrstart){
		++$totalsnvs;
	}
	my $sn1counts=0;
	for my $i (0..$#sn1){
		if ($sn1[$i] eq "1"){
			++$sn1counts;
		}
	}
	my $dnapolncounts=0;
	for my $i (0..$#dnapoln){
		if ($dnapoln[$i] eq "1"){
			++$dnapolncounts;
		}
	}
	my $oxogcounts=0;
	for my $i (0..$#oxog){
		if ($oxog[$i] eq "1"){
			++$oxogcounts;
		}
	}
	my $uvlambdacounts=0;
	for my $i (0..$#uvlambda){
		if ($uvlambda[$i] eq "1"){
			++$uvlambdacounts;
		}
	}
	my $sixfourcounts=0;
	for my $i (0..$#sixfour){
		if ($sixfour[$i] eq "1"){
			++$sixfourcounts;
		}
	}
	my $enucounts=0;
	for my $i (0..$#enu){
		if ($enu[$i] eq "1"){
			++$enucounts;
		}
	}
	my $uvsolarcounts=0;
	for my $i (0..$#uvsolar){
		if ($uvsolar[$i] eq "1"){
			++$uvsolarcounts;
		}
	}

	my $normsn1=$sn1counts/$totalsnvs;
	my $normdnapoln=$dnapolncounts/$totalsnvs;
	my $normoxog=$oxogcounts/$totalsnvs;
	my $normuvlambda=$uvlambdacounts/$totalsnvs;
	my $normsixfour=$sixfourcounts/$totalsnvs;
	my $normenu=$enucounts/$totalsnvs;
	my $normuvsolar=$uvsolarcounts/$totalsnvs;

	# strand-concordance calculation
	open (OUTPUTBIASSN1, ">>results-$variantfile/WOLAND-bias_SN1-$variantfile");
	for my $i (0..$#sn1){
		if ($sn1[$i] eq "1"){
			my $strandvalue=$sn1plus[$i];
			my $querychr=$chrstart[$i];
			my $querycoord=$pos[$i];
			my $strandpluscount=0;
    		my $strandcount=0;
			for my $ii (0 .. $#chrrefseq){
				if ($chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
					++$strandcount;
				}
				if ($strandrefseq[$ii] eq "+" && $chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
				++$strandpluscount;
				} 
			}
			if ($strandcount >= 1){
				my $strandtranscript=$strandpluscount/$strandcount;
    			my $strandscore=$strandtranscript - ($strandvalue);
    			print OUTPUTBIASSN1 "$chrstart[$i]\t$pos[$i]\t$strandscore\n";
    		}
    	}
	} 
	close (OUTPUTBIASSN1);

	open (OUTPUTBIASDNAPOLN, ">>results-$variantfile/WOLAND-bias_DNApoln-$variantfile");
	for my $i (0..$#dnapoln){
		if ($dnapoln[$i] eq "1"){
			my $strandvalue=$dnapolnplus[$i];
			my $querychr=$chrstart[$i];
			my $querycoord=$pos[$i];
			my $strandpluscount=0;
    		my $strandcount=0;
			for my $ii (0 .. $#chrrefseq){
				if ($chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
					++$strandcount;
				}
				if ($strandrefseq[$ii] eq "+" && $chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
				++$strandpluscount;
				} 
			}
			if ($strandcount >= 1){
				my $strandtranscript=$strandpluscount/$strandcount;
    			my $strandscore=$strandtranscript - ($strandvalue);
    			print OUTPUTBIASDNAPOLN "$chrstart[$i]\t$pos[$i]\t$strandscore\n";
    		}
    	}
	} 
	close (OUTPUTBIASDNAPOLN);


	open (OUTPUTBIASOXOG, ">>results-$variantfile/WOLAND-bias_oxoG-$variantfile");
	for my $i (0..$#oxog){
		if ($oxog[$i] eq "1"){
			my $strandvalue=$oxogplus[$i];
			my $querychr=$chrstart[$i];
			my $querycoord=$pos[$i];
			my $strandpluscount=0;
    		my $strandcount=0;
			for my $ii (0 .. $#chrrefseq){
				if ($chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
					++$strandcount;
				}
				if ($strandrefseq[$ii] eq "+" && $chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
				++$strandpluscount;
				} 
			}
			if ($strandcount >= 1){
				my $strandtranscript=$strandpluscount/$strandcount;
    			my $strandscore=$strandtranscript - ($strandvalue);
    			print OUTPUTBIASOXOG "$chrstart[$i]\t$pos[$i]\t$strandscore\n";
    		}
    	}
	} 
	close (OUTPUTBIASOXOG);

	open (OUTPUTBIASUVLAMBDA, ">>results-$variantfile/WOLAND-bias_UV-lambda-$variantfile");
	for my $i (0..$#uvlambda){
		if ($uvlambda[$i] eq "1"){
			my $strandvalue=$uvlambdaplus[$i];
			my $querychr=$chrstart[$i];
			my $querycoord=$pos[$i];
			my $strandpluscount=0;
    		my $strandcount=0;
			for my $ii (0 .. $#chrrefseq){
				if ($chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
					++$strandcount;
				}
				if ($strandrefseq[$ii] eq "+" && $chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
				++$strandpluscount;
				} 
			}
			if ($strandcount >= 1){
				my $strandtranscript=$strandpluscount/$strandcount;
    			my $strandscore=$strandtranscript - ($strandvalue);
    			print OUTPUTBIASUVLAMBDA "$chrstart[$i]\t$pos[$i]\t$strandscore\n";
    		}
    	}
	} 
	close (OUTPUTBIASUVLAMBDA);

	open (OUTPUTBIASUVSOLAR, ">>results-$variantfile/WOLAND-bias_UVsolar-$variantfile");
	for my $i (0..$#uvsolar){
		if ($uvsolar[$i] eq "1"){
			my $strandvalue=$uvsolarplus[$i];
			my $querychr=$chrstart[$i];
			my $querycoord=$pos[$i];
			my $strandpluscount=0;
    		my $strandcount=0;
			for my $ii (0 .. $#chrrefseq){
				if ($chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
					++$strandcount;
				}
				if ($strandrefseq[$ii] eq "+" && $chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
				++$strandpluscount;
				} 
			}
			if ($strandcount >= 1){
				my $strandtranscript=$strandpluscount/$strandcount;
    			my $strandscore=$strandtranscript - ($strandvalue);
    			print OUTPUTBIASUVSOLAR "$chrstart[$i]\t$pos[$i]\t$strandscore\n";
    		}
    	}
	} 
	close (OUTPUTBIASUVSOLAR);

	open (OUTPUTBIASSIXFOUR, ">>results-$variantfile/WOLAND-bias_sixfour-$variantfile");
	for my $i (0..$#sixfour){
		if ($sixfour[$i] eq "1"){
			my $strandvalue=$sixfourplus[$i];
			my $querychr=$chrstart[$i];
			my $querycoord=$pos[$i];
			my $strandpluscount=0;
    		my $strandcount=0;
			for my $ii (0 .. $#chrrefseq){
				if ($chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
					++$strandcount;
				}
				if ($strandrefseq[$ii] eq "+" && $chrrefseq[$ii] eq $querychr and (($querycoord <= $txstartrefseq[$ii] and $querycoord >= $txtendrefseq[$ii]) or ($querycoord >= $txstartrefseq[$ii] and $querycoord <= $txtendrefseq[$ii]))){
				++$strandpluscount;
				} 
			}
			if ($strandcount >= 1){
				my $strandtranscript=$strandpluscount/$strandcount;
    			my $strandscore=$strandtranscript - ($strandvalue);
    			print OUTPUTBIASSIXFOUR "$chrstart[$i]\t$pos[$i]\t$strandscore\n";
    		}
    	}
	} 
	close (OUTPUTBIASSIXFOUR);

	#output for motif number
	open (MOTIFS, ">>results-$variantfile/WOLAND-motifs-$variantfile");
	print MOTIFS "Chr\tPos\ttargetSequence\tClass\tGene\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUVA-solar\n";
	for my $i (0..$#extractedsequenceformotifsearch){
		print MOTIFS "$chrstart[$i]\t$pos[$i]\t$extractedsequenceformotifsearch[$i]\t$geneclass[$i]\t$genename[$i]\t$sn1[$i]\t$dnapoln[$i]\t$oxog[$i]\t$uvlambda[$i]\t$sixfour[$i]\t$enu[$i]\t$uvsolar[$i]\n";
	}
	close (MOTIFS);

	#output for motif number normalized
	open (NMOTIFS, ">>results-$variantfile/WOLAND-norm_motifs-$variantfile");
	print NMOTIFS "$variantfile\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUVA-solar\n";
	print NMOTIFS "Number of Total SNPs\t$totalsnvs\t$totalsnvs\t$totalsnvs\t$totalsnvs\t$totalsnvs\t$totalsnvs\t$totalsnvs\n";
	print NMOTIFS "Total Raw Number of Motifs Found\t$sn1counts\t$dnapolncounts\t$oxogcounts\t$uvlambdacounts\t$sixfourcounts\t$enucounts\t$uvsolarcounts\n";
	print NMOTIFS "Normalized Number of Motifs Found\t$normsn1\t$normdnapoln\t$normoxog\t$normuvlambda\t$normsixfour\t$normenu\t$normuvsolar";
	close (NMOTIFS);

}

# start Screen & log file
$datestring = localtime();

#main warning
unless ($#ARGV==3){
	die "\nERROR : Incorrect number of arguments - Usage: $0 <tabular_snp_file> <chromosome_length_profile> <hotspot_window_length> <genome_version>  \n\n";	
}
for my $i (0,1) {
	unless (-r -e -f $ARGV[$i]){
    	die "\nERROR: $ARGV[$i] not exists or is not readable or not properly formatted. Please check file.\n\n";
    }
}

#loading files and parameters
print "\nLoading SNP file...\n";
$contextsequencelength=3; #<context_sequence_length> #default=3nt downstream & 3nt upstream
&parse_inputs (@ARGV);
&process_refseq;

#creating results folder
mkdir("results-$variantfile", 0755) || die "Cannot mkdir results-$variantfile - folder already exists, please delete it or change samplename";

#counting of nucleotide changes.
print "\nCalculating general mutation statistics and saving basechange-$variantfile ...\n";
&count_nucleotidechanges;

#counting of nucleotide frequency across chromosomes
for my $i (0..$#chrnames){
	&count_chrhits("$chrnames[$i]");
}
&count_chrfreq;

#hotspots 
print "\nCalculating hotspots and saving hotspots-$variantfile.txt ...\n";
@uniquechrinsamples = uniq @chrstart;

for my $i(0..$#uniquechrinsamples){
	&search_for_hotspots($uniquechrinsamples[$i]);
}

open (HOTSPOTS, ">>results-$variantfile/WOLAND-hotspots-$variantfile");
print HOTSPOTS "geneclass\tgenename\tCHR\tBP\tHotspotCount\n";
for my $i (0 .. $#chrstart){
	print HOTSPOTS "$geneclass[$i]\t$genename[$i]\t$chrnameformatted[$i]\t$pos[$i]\t$hotspotsdetected[$i]\n";
}

#module for extraction of context sequences in reference genomes of each SNP.
$genomesequenceinstance = Bio::DB::Fasta->new("genomes/genome_$genomeversion.fa");
open (CONTEXTSEQ, ">>results-$variantfile/WOLAND-contextsequences-$variantfile");
open (CONTEXTSEQANNO, ">>results-$variantfile/WOLAND-contextsequencesanno-$variantfile");
print "\nExtracting context sequences and saving extracted_sequences-$variantfile ..\n";
for my $i(0..$#chrstart){
	&extract_contextsequences($chrstart[$i], $pos[$i], $genomesequenceinstance, $geneclass[$i], $genename[$i]);
}
close(CONTEXTSEQ);
close(CONTEXTSEQANNO);

print "\nPerforming motif search & transcriptional strand concordance score calculation on $variantfile\n";

#motif search & analysis
&motif_search;

#printing log file info
open (LOG, ">>results-$variantfile/WOLAND-log-$variantfile");
print LOG "WOLAND BETA 0.2 - 26-01-2017\n";
print LOG "Analysis started at $datestring\n";
print LOG "Tabular SNP File              : $variantfile\n";
print LOG "Chromosome Length Profile File: $chrlengthprofile\n";
print LOG "Hotspot Window Length         : $hotspotwindowlength bases flanking SNP position in reference genome\n";
print LOG "Context Sequence Length       : $contextsequencelength bases flanking SNP position in reference genome\n";
print LOG "\nOutputs:\n";
print LOG "Mutation Statistics:                       WOLAND-basechange-$variantfile\n";
print LOG "Extracted Unannotated Context Sequences:   WOLAND-contextsequences-$variantfile\n";
print LOG "Extracted Context Sequences:               WOLAND-contextsequencesanno-$variantfile\n";
print LOG "Hotspots:                                  WOLAND-hotspots-$variantfile\n";
print LOG "Mutational Motifs:                         WOLAND-motifs-$variantfile\n";
print LOG "Mutational Motifs Normalized:              WOLAND-norm_motifs-$variantfile\n";
print LOG "Strand Bias of SN1 Motifs:                 WOLAND-bias_uv-lambda-$variantfile\n";
print LOG "Strand Bias of DNApoln Motifs:             WOLAND-bias_DNApoln-$variantfile\n";
print LOG "Strand Bias of 8-oxoG Motifs:              WOLAND-bias_oxoG-$variantfile\n";
print LOG "Strand Bias of UV-lambda Motifs:           WOLAND-bias_uv-lambda-$variantfile\n";
print LOG "Strand Bias of UV Solar Motifs:            WOLAND-bias_UVsolar-$variantfile\n";
print LOG "Strand Bias of 6-4 Motifs:                 WOLAND-bias_sixfour-$variantfile\n";

print "\nDONE\n";

exit;