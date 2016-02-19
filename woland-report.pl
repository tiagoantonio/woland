########################################################################################################################
## WOLAND REPORT BETA 0.1 (15-02-2016)
##
## WOLAND is a software package based on Perl and R for calculation of general mutation metrics, identification and
## comparison of predicted hotspots across biological samples. WOLAND uses Single Nucleotide Polymorphisms (SNPs) data
## from Next Generation Sequencing (NGS) pipelines as main input. Please observe file requirements.
##
##
######################################################################################################################## 



#! /usr/bin/perl

use Statistics::R; # module for the bridge between Perl and R.
use strict;
use warnings;

## variables

my ($inputTable, $inputTableLine, @inputTableArray);

my ($i, $i2, $i3, $i4); my (@i, @i2, @i3); #variables for indices

my (@Group, @sampleName, $sampleNameLine);

my ($inputChromosomeProfile, $inputChromosomeProfileArrayLine, @inputChromosomeProfileArray);

my ($inputSampleBaseChange, @CHR_P, @MUT_P, @inputSampleBaseChangeArray, @AT,@AG,@AC,@CG,@CT,@CA,@ATF,@AGF,@ACF,@CGF,@CTF,@CAF);

my ($TransversionF, $TransitionF, @TransversionF, @TransitionF);

my ($inputSampleMotif, $inputSampleMotifArray, @inputSampleMotifArray);

my ($inputSampleHotSpot, $inputSampleHotSpotArrayLine, $item, @SN1,@DNApol,@oxoG,@UV,@SixFour,@ENU,@UVAsolar,@SN1F,@DNApolF,@oxoGF,@UVF,@SixFourF,@ENUF,@UVAsolarF);

my ($uniqueGroup, $input, $HofA, @inputSampleHotSpotArray, @geneClass, @GeneName, @CHR, @BP, @HotspotCount, @uniqgroup);

my ($inputStrandScoreArrayLine, $inputStrandScore, $SC_0, $SC_1,$SC_m1,$SC_SC,$totalSC,$nConcordance,$nDiscordance,$ratioConcordanceDiscordance,$n, @inputStrandScoreArray, @CHR_SC, @BP_SC, @SC_SC, @colorGauss);

my ($R, $colorGausspick, $uniqGroupLine, $pick, $legendltyV, $legendlwdV, $lastgroup, @legendname, @legendlty, @legendlwd, @legendcol);


## main Warning

unless (@ARGV){
	die "\nERROR : Usage: $0 woland.input.table ... \n";	
}


## parsing input table

$inputTable = $ARGV[0]; # <input.table>
open (inputTable, $inputTable);
@inputTableArray=<inputTable>;

foreach $inputTableLine (@inputTableArray){ # two arrays for each category (group & sample results folder)
	@i = split (/\t/, $inputTableLine);
	chomp (@i);
	push (@Group, "$i[0]"); # array for group definition
	push (@sampleName, "$i[1]"); # array for sample folder definition
}

@i=0;


##Create Report Folder

mkdir("report-$inputTable", 0755);


### Frequency Histogram across chromosomes 


for my $i3 (0..$#sampleName){
	$inputChromosomeProfile = "results-$sampleName[$i3]/WOLAND-mutfreq-$sampleName[$i3]";
	open (inputChromosomeProfile, $inputChromosomeProfile);
	@inputChromosomeProfileArray=<inputChromosomeProfile>;

	foreach $inputChromosomeProfileArrayLine (@inputChromosomeProfileArray){
		@i2 = split (/\t/, $inputChromosomeProfileArrayLine);
		chomp (@i2);
		push (@CHR_P, "$i2[0]");
		push (@MUT_P, "$i2[2]");
	}

	open (MUTFREQ, ">>report-$inputTable/mutfreq-$inputTable.tmp");

	for my $i (1..$#CHR_P){
		print MUTFREQ "$CHR_P[$i]\t$MUT_P[$i]\t$Group[$i3]\n";
	}

	$i=0;

	close (MUTFREQ);
	$inputChromosomeProfile=0;
	$inputChromosomeProfileArrayLine=0;
	@inputChromosomeProfileArray=();
	@i = 0;
	@CHR_P =();
	@MUT_P=();
}

$i3=0;

#### Nucleotide-type change - BOXPLOT & and transition/transvertion ratio PIE CHART####

## building grouped sample input arrays for R using parsed input table

foreach $sampleNameLine (@sampleName){
	$inputSampleBaseChange = "results-$sampleNameLine/WOLAND-basechange-$sampleNameLine";
	open (inputSampleBaseChange, $inputSampleBaseChange);
	@inputSampleBaseChangeArray=<inputSampleBaseChange>;
	@i = split (/\t/, $inputSampleBaseChangeArray[1]);
	chomp (@i);
	push (@AT, "$i[1]");
	push (@AG, "$i[2]");
	push (@AC, "$i[3]");
	push (@CG, "$i[4]");
	push (@CT, "$i[5]");
	push (@CA, "$i[6]");
	@i = split (/\t/, $inputSampleBaseChangeArray[2]);
	chomp (@i);
	push (@ATF, "$i[1]");
	push (@AGF, "$i[2]");
	push (@ACF, "$i[3]");
	push (@CGF, "$i[4]");
	push (@CTF, "$i[5]");
	push (@CAF, "$i[6]");
}
$inputSampleBaseChange=0;
@inputSampleBaseChangeArray=0;
@i=0;

## temporary output file for R analysis of Nucleotide-type changes

open (BASECHANGE, ">>report-$inputTable/nucleotide_type_change-$inputTable.tmp");

print BASECHANGE "X\tA.T\tA.G\tA.C\tC.G\tC.T\tC.A\n";
for my $i (0 .. $#Group){
	print BASECHANGE "$Group[$i]\t$AT[$i]\t$AG[$i]\t$AC[$i]\t$CG[$i]\t$CT[$i]\t$CA[$i]\n";
}

@i=0;
close (BASECHANGE);

## temporary output file for R analysis of Nucleotide-type changes frequency

open (BASECHANGEF, ">>report-$inputTable/nucleotide_type_changeF-$inputTable.tmp");
print BASECHANGEF "X\tA.T\tA.G\tA.C\tC.G\tC.T\tC.A\n";

for my $i (0 .. $#Group){
	print BASECHANGEF "$Group[$i]\t$ATF[$i]\t$AGF[$i]\t$ACF[$i]\t$CGF[$i]\t$CTF[$i]\t$CAF[$i]\n";
}
@i=0;
close (BASECHANGEF);

## temporary output file for Transvertion Transition Ratio

for my $i (0..$#sampleName){
	$TransversionF=$ATF[$i]+$ACF[$i]+$CGF[$i]+$CAF[$i];
	$TransitionF=$AGF[$i]+$CTF[$i];
	push (@TransversionF, "$TransversionF");
	push (@TransitionF, "$TransitionF");
}

open (TRANSITIONTRANSVERSION, ">>report-$inputTable/transitiontransversionF-$inputTable.tmp");
for my $i (0 .. $#sampleName){
	print TRANSITIONTRANSVERSION "$TransversionF[$i]\tTransversion\t$sampleName[$i]\t$Group[$i]\n$TransitionF[$i]\tTransition\t$sampleName[$i]\t$Group[$i]\n";
}

@i=0;
close (TRANSITIONTRANSVERSION);

#### Motif-Search Number BOXPLOT #####

## building grouped sample input arrays for R using parsed input table

foreach $sampleNameLine (@sampleName){
	$inputSampleMotif = "results-$sampleNameLine/WOLAND-norm_motifs-$sampleNameLine";
	open (inputSampleMotif, $inputSampleMotif);
	@inputSampleMotifArray=<inputSampleMotif>;
	@i = split (/\t/, $inputSampleMotifArray[2]);
	chomp (@i);
	push (@SN1, "$i[1]");
	push (@DNApol, "$i[2]");
	push (@oxoG, "$i[3]");
	push (@UV, "$i[4]");
	push (@SixFour, "$i[5]");
	push (@ENU, "$i[6]");
	push (@UVAsolar, "$i[7]");
	@i = split (/\t/, $inputSampleMotifArray[3]);
	chomp (@i);
	push (@SN1F, "$i[1]");
	push (@DNApolF, "$i[2]");
	push (@oxoGF, "$i[3]");
	push (@UVF, "$i[4]");
	push (@SixFourF, "$i[5]");
	push (@ENUF, "$i[6]");
	push (@UVAsolarF, "$i[7]");
}

$inputSampleMotif=0;
@inputSampleMotifArray=0;
@i=0;


## temporary output file for R analysis of MOTIF number

open (MOTIFNUMBER, ">>report-$inputTable/motif_number-$inputTable.tmp");
print MOTIFNUMBER "X\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUV-solar\n";

for my $i (0 .. $#Group){
	print MOTIFNUMBER "$Group[$i]\t$SN1[$i]\t$DNApol[$i]\t$oxoG[$i]\t$UV[$i]\t$SixFour[$i]\t$ENU[$i]\t$UVAsolar[$i]\n";
}
@i=0;
close (MOTIFNUMBER);

## temporary output file for R analysis of MOTIF normalized
open (MOTIFNUMBERNORM, ">>report-$inputTable/motif_numberNorm-$inputTable.tmp");
print MOTIFNUMBERNORM "X\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUV-solar\n";

for my $i (0 .. $#Group){
	print MOTIFNUMBERNORM "$Group[$i]\t$SN1F[$i]\t$DNApolF[$i]\t$oxoGF[$i]\t$UVF[$i]\t$SixFourF[$i]\t$ENUF[$i]\t$UVAsolarF[$i]\n";
}
@i=0;
close (MOTIFNUMBERNORM);

#### Hotspots Manhattan Plot #####

##Grouping Samples using parsed input.table

for my $i3 (0..$#sampleName){
	$inputSampleHotSpot = "results-$sampleName[$i3]/WOLAND-hotspots-$sampleName[$i3]";
	open (inputSampleHotSpot, $inputSampleHotSpot);
	@inputSampleHotSpotArray=<inputSampleHotSpot>;

	foreach $inputSampleHotSpotArrayLine (@inputSampleHotSpotArray){
		@i = split (/\t/, $inputSampleHotSpotArrayLine);
		chomp (@i);
		push (@geneClass, "$i[0]");
		push (@GeneName, "$i[1]");
		push (@CHR, "$i[2]");
		push (@BP, "$i[3]");
		push (@HotspotCount, "$i[4]");
	}

	open (HOTSPOT, ">>report-$inputTable/hotspot-$Group[$i3].tmp");

	for my $i2 (1..$#geneClass){
		print HOTSPOT "$geneClass[$i2]\t$GeneName[$i2]\t$CHR[$i2]\t$BP[$i2]\t$HotspotCount[$i2]\n";
	}

	$i2=0;
	close (HOTSPOT);
	$inputSampleHotSpot=0;
	$inputSampleHotSpotArrayLine=0;
	@inputSampleHotSpotArray=();
	@i = 0;
	@geneClass =();
	@GeneName =();
	@CHR =();
	@BP =();
	@HotspotCount =();
}

$i3=0;

## Parsing unique group names for Hotspots

my %uniqueGroup = ();

foreach $item (@Group) {
	push(@uniqgroup, $item) unless $uniqueGroup{$item}++;
}


#### Box Plot for Concordance/Discordace SC ratio #####
##Parsing Grouped Files

open (CONCORDANCEDISCORDANCERATIO, ">>report-$inputTable/SC_concordance_ratio-$inputTable.tmp");
print CONCORDANCEDISCORDANCERATIO "X\tSN1\tDNApol\t8-oxoG\tUVlambda\tSixFour\tUV-solar\n";

$input = 6;
$HofA = ();

my %HofA=();
for (0..$input) {
   $HofA{$_} = [];
}

$n=0;

sub SCRatioGroup{

	for my $i3 (0..$#sampleName){
		$inputStrandScore = "results-$sampleName[$i3]/WOLAND-bias_$_[0]-$sampleName[$i3]";
		open (inputStrandScore, $inputStrandScore);
		@inputStrandScoreArray=<inputStrandScore>;

		foreach $inputStrandScoreArrayLine (@inputStrandScoreArray){
			@i = split (/\t/, $inputStrandScoreArrayLine);
			chomp (@i);
			push (@CHR_SC, "$i[0]");
			push (@BP_SC, "$i[1]");
			push (@SC_SC, "$i[2]");
		}

		$SC_0=0;
		$SC_1=0;
		$SC_m1=0;

		for my $i4 (0..$#SC_SC){
			if ($SC_SC[$i4] == 0){ ++$SC_0;}
			if ($SC_SC[$i4] == 1){ ++$SC_1;}
			if ($SC_SC[$i4] eq "-1"){ ++$SC_m1;}
		}

		$totalSC=$SC_0+$SC_1+$SC_m1;
			if ($totalSC==0){
				$totalSC=1;
			}
		$nConcordance=$SC_0/$totalSC;
		$nDiscordance=($SC_1+$SC_m1)/$totalSC;
			if ($nDiscordance==0){
				$nDiscordance=1;
			}
		$ratioConcordanceDiscordance=$nConcordance/$nDiscordance;

		push @{ $HofA{$n} } , $ratioConcordanceDiscordance;

		$SC_0=0;
		$SC_1=0;
		$SC_m1=0;

		@CHR_SC =();
		@BP_SC =();
		@SC_SC =();
	}
	++$n;
}

&SCRatioGroup ("SN1");
&SCRatioGroup ("DNApoln");
&SCRatioGroup ("oxoG");
&SCRatioGroup ("UV-lambda");
&SCRatioGroup ("sixfour");
#&SCRatioGroup ("enu");
&SCRatioGroup ("UVsolar");

for my $i (0 .. $#Group){
	print CONCORDANCEDISCORDANCERATIO "$Group[$i]\t$HofA{0}[$i]\t$HofA{1}[$i]\t$HofA{2}[$i]\t$HofA{3}[$i]\t$HofA{4}[$i]\t$HofA{5}[$i]\n";
}

@i=0;
close (CONCORDANCEDISCORDANCERATIO);
$i3=0;


### R bridge

$R = Statistics::R->new() ;
$R->start_sharedR ;
$R->send('library (reshape2)');
$R->send('library (ggplot2)');
$R->send('library(qqman)');
$R->send('library (RColorBrewer)');
$R->send('library (plyr)');
$R->send(qq'setwd(dir = "./report-$inputTable/")');

#Box Plots of Grouped Samples

sub BoxPlotGroup{	

	$R->send(qq'WOLAND.$_[0].boxplot<-read.table("$_[0]-$inputTable.tmp", sep = "\t", header=TRUE)');
	$R->send(qq'WOLAND.$_[0].boxplot.m <- melt (WOLAND.$_[0].boxplot)');
	$R->send(qq'pdf(file="$_[0]_number_boxplot_$inputTable.pdf", width=11.692, height=8.267)');
	$R->send(qq'ggplot (WOLAND.$_[0].boxplot.m, aes(x=variable, y=value, fill = X))+
	ggtitle("$_[0]_$inputTable")+
	theme(plot.title = element_text(size=16, vjust=1.1))+
	theme(axis.title.y = element_text(size=16))+
	scale_x_discrete(name="")+
	geom_boxplot(lwd=0.25)+
	scale_fill_brewer(name="Group", palette="Spectral")');
	$R->send('dev.off()');
}

&BoxPlotGroup("nucleotide_type_changeF");
&BoxPlotGroup("nucleotide_type_change");
&BoxPlotGroup("motif_number");
&BoxPlotGroup("motif_numberNorm");


#Manhattan Plot Hotspots

foreach $uniqGroupLine (@uniqgroup){
	$R->send(qq'WOLAND.hotspot.manhattan.$uniqGroupLine <- read.delim ("hotspot-$uniqGroupLine.tmp", comment.char="#", header=FALSE)');
	$R->send(qq'pdf("manhattan.hotspot.$uniqGroupLine-$inputTable.pdf")');
	$R->send(qq'manhattan.$uniqGroupLine <- manhattan(x = WOLAND.hotspot.manhattan.$uniqGroupLine, chr="V3", bp="V4", p = "V5", logp = FALSE, ylab= "Mutations in $ARGV[2] bp window", genomewideline = FALSE, suggestiveline = FALSE, main = "Sample:$uniqGroupLine:hotspotwindow:$ARGV[2]", ylim = c(0,10), col = c("blue4", "orange3"))');
	$R->send(q'dev.off()');
}

#Gaussian Density Strand Concordance

@colorGauss=("\"#9E0142\"", "\"#5E4FA2\"","\"#D53E4F\"","\"#3288BD\"", "\"#F46D43\"","\"#66C2A5\"","\"#FDAE61\"","\"#ABDDA4\"","\"#FEE08B\"",
			 "\"#E6F598\"" );

			
sub GaussGraphPlot{

	$R->send(qq'pdf(file="gauss_$_[0]-$inputTable.pdf", width=11.692, height=8.267)');
	$R->send(qq'plot(5,
			5,
			main="Kernel Density Estimation for SC score of $_[0]-$inputTable",
			xlab="SC score",
			ylab="SC score for bw=0.05",
			xlim=c(-1.25, 1.25),
			ylim = c(0,5))'
			);
	$pick=0;
	foreach $uniqGroupLine(@uniqgroup){

		$colorGausspick=$colorGauss[$pick];

		for my $i3 (0..$#sampleName){
			if($Group[$i3] eq $uniqGroupLine){
				$R->send(qq'FirstPlot <- read.delim("../results-$sampleName[$i3]/WOLAND-bias_$_[0]-$sampleName[$i3]", header=FALSE)');
				$R->send(q'fp<-density (FirstPlot$V3, bw=0.05)');
				$R->send(qq'lines (fp, col=$colorGausspick, lwd=2, ylim=c(0,5))');
			}
		}

		$pick++;
	}

	$legendltyV="1\,";
	$legendlwdV="2.5\,";
	$lastgroup=$#uniqgroup;

	for my $i3(0..$#uniqgroup){

		push (@legendname, "\"$uniqgroup[$i3]\"\,");
		push (@legendlty, "$legendltyV");
		push (@legendlwd, "$legendlwdV");
		push (@legendcol, "$colorGauss[$i3]\,");

		if ($i3==$lastgroup){
			chop ($legendname[$i3]);
			chop ($legendlty[$i3]);
			chop ($legendlwd[$i3]);
			chop ($legendcol[$i3]);
		}
	}
	
	$R->send(qq'legend("topright",c(@legendname),lty=c(@legendlty), lwd=c(@legendlwd), col=c(@legendcol))');
	$R->send(q'dev.off()');
	$pick=0;
	@legendname=();
	@legendlty=();
	@legendlwd=();
	@legendcol=();
}

&GaussGraphPlot("UV-lambda");
&GaussGraphPlot("sixfour");
&GaussGraphPlot("DNApoln");
&GaussGraphPlot("oxoG");
#&GaussGraphPlot("enu");
&GaussGraphPlot("UVsolar");
&GaussGraphPlot("SN1");

## BoxPlot for Concordance/Discordance Strand Score ratio

&BoxPlotGroup("SC_concordance_ratio");

## Chromosome Profile Mean Mutation Rate

$R->send(qq'mutfreq <- read.delim ("mutfreq-$inputTable.tmp", comment.char="#", header=FALSE)');
$R->send(qq'svg("mutfreq-$inputTable.svg", width=11.692, height=8.267)');


$R->send(qq'melted <- melt(mutfreq,id.vars=c("V1","V2","V3"))');
$R->send(qq'means = ddply(melted,c("V1","V3"),summarise, mean=mean(V2))');
$R->send(qq'means.sem = ddply(melted,c("V1","V3"),summarise, mean=mean(V2), sem=sd(V2)/sqrt(length(V2)))');
$R->send(qq'means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)');
$R->send(qq'means.group<-ddply(means,"V3",summarise, mean=mean(mean))');
$R->send(qq'colnames(means.group)<-c("Average","mean")');

$R->send('ggplot(data=means, aes(x = V1, y = mean, fill=V3))+
	geom_bar(stat="identity", position ="dodge")+
	xlab("Chromosome") + ylab("Mutations per base sequenced") +
	theme(plot.title = element_text(size=16, vjust=1.1))+
	theme(axis.title.y = element_text(size=16))+
	geom_hline(aes(yintercept=mean), data=means.group)+
	scale_fill_brewer(name="Group", palette="Spectral")+ 
	geom_errorbar(aes(ymax=upper,ymin=lower), size=.3, width=.2, position=position_dodge(0.9),data=means.sem)+ 
	geom_hline(aes(yintercept=mean, col=Average), data=means.group)');
$R->send('dev.off()');

## Transition e Transversion Barplot Pie
$R->send(qq'pdf(file="barplot_pie_TransTransv-$inputTable.pdf", width=11.692, height=8.267)');
$R->send(qq'transitiontransversion <- read.delim("transitiontransversionF-$inputTable.tmp", header=FALSE)');
$R->send(qq'ggplot(transitiontransversion, aes(x=V3, y=V1, fill=V2))+
		ggtitle("Transversion & Transition frequency-$inputTable")+
		theme(plot.title = element_text(size=16, vjust=1.1))+
		xlab("Sample") + ylab("Frequency") +
		theme(axis.text.x = element_text(size=8))+
		geom_bar(position="fill", stat = "identity")+
		scale_fill_brewer(name="Type", palette="Spectral")');
$R->send(q'dev.off()');

$R->stop;

exit;
