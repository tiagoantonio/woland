########################################################################################################################
## WOLAND Beta 1.01 (09-30-2017)
## woland-report.pl 
##
## WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. 
## 
## For more details please read README file.
##
######################################################################################################################## 

#! /usr/bin/perl
use Statistics::R; # module for the bridge between Perl and R.
use strict;
use warnings;
use Getopt::ArgParse;

our $REVISION = '$Revision:  $';
our $DATE =	'$Date: 2017-09-30 00:11:04 -0800 (Sat,  30 Sep 2017) $';  
our $AUTHOR =	'$Author: Tiago A. de Souza <tiagoantonio@gmail.com> $';

## variables
our @tablearray; our ($table, $tableline); 
our %uniquegroup; our (@samplename, @group, @uniqgroup); our ($item, $uniqgroupline); 
our (@AT,@AG,@AC,@CG,@CT,@CA,@ATfrequency,@AGfrequency,@ACfrequency,@CGfrequency,@CTfrequency,@CAfrequency);
our (@transversionfrequency,@transitionfrequency);
our (@SN1,@DNApol,@oxoG,@UV,@SixFour,@ENU,@UVAsolar,@SN1F,@DNApolF,@oxoGF,@UVF,@SixFourF,@ENUF,@UVAsolarF);
our %motifconcordancediscordanceratio;  our @i; our ($input, $i, $pick, $Rbridge);

## subs

sub calculate_mutationfrequency{ #mutational frequency grouping across chromosomes
	my $chromosomeprofile = "results-$_[0]/WOLAND-mutfreq-$_[0]";
	open (CHROMOSOMEPROFILE, $chromosomeprofile);
	my @chromosomeprofilearray=<CHROMOSOMEPROFILE>;
	close (CHROMOSOMEPROFILE);
	my @chrprofilevalue;
	my @mutprofilevalue;

	foreach my $chromosomeprofilearrayline (@chromosomeprofilearray){
		my @ii = split (/\t/, $chromosomeprofilearrayline);
		chomp (@ii);
		push (@chrprofilevalue, "$ii[0]");
		push (@mutprofilevalue, "$ii[2]");
	}

	open (MUTATIONFREQUENCY, ">>report-$table/mutfreq-$table.txt");

	for my $ii (1..$#chrprofilevalue){
		print MUTATIONFREQUENCY "$chrprofilevalue[$ii]\t$mutprofilevalue[$ii]\t$group[$_[1]]\n";
	}
	close (MUTATIONFREQUENCY);
}

sub calculate_nucleotidebasechanges{ #nucleotide change grouping
	my $basechanges = "results-$_[0]/WOLAND-basechange-$_[0]";
	open (BASECHANGES, $basechanges);
	my @basechangesarray=<BASECHANGES>;
	close (BASECHANGES);
	my @ii = split (/\t/, $basechangesarray[1]);
	chomp (@ii);
	push (@AT, "$ii[1]"); push (@AG, "$ii[2]"); push (@AC, "$ii[3]");
	push (@CG, "$ii[4]"); push (@CT, "$ii[5]"); push (@CA, "$ii[6]");
	
	@ii = split (/\t/, $basechangesarray[2]);
	chomp (@ii);
	push (@ATfrequency, "$ii[1]"); push (@AGfrequency, "$ii[2]");
	push (@ACfrequency, "$ii[3]"); push (@CGfrequency, "$ii[4]");
	push (@CTfrequency, "$ii[5]"); push (@CAfrequency, "$ii[6]");
}

sub calculate_transtransvratio{ #transition transversion ratio
	my $transversionfrequency=$ATfrequency[$_[0]]+$ACfrequency[$_[0]]+$CGfrequency[$_[0]]+$CAfrequency[$_[0]];
	my $transitionfrequency=$AGfrequency[$_[0]]+$CTfrequency[$_[0]];
	push (@transversionfrequency, "$transversionfrequency");
	push (@transitionfrequency, "$transitionfrequency");
}

sub find_motifs{ #group motifs
	my $motif = "results-$_[0]/WOLAND-norm_motifs-$_[0]";
	open (MOTIF, $motif);
	my @motifarray=<MOTIF>;
	
	my @ii = split (/\t/, $motifarray[2]);
	chomp (@ii);
	push (@SN1, "$ii[1]"); push (@DNApol, "$ii[2]");
	push (@oxoG, "$ii[3]"); push (@UV, "$ii[4]");
	push (@SixFour, "$ii[5]"); push (@ENU, "$ii[6]");
	push (@UVAsolar, "$ii[7]");

	@ii = split (/\t/, $motifarray[3]);
	chomp (@ii);
	push (@SN1F, "$ii[1]"); push (@DNApolF, "$ii[2]");
	push (@oxoGF, "$ii[3]"); push (@UVF, "$ii[4]");
	push (@SixFourF, "$ii[5]"); push (@ENUF, "$ii[6]");
	push (@UVAsolarF, "$ii[7]");

	close(MOTIF);	
}

sub to_merged_hotspots{ #merging hotspots

	my $samplehotspot = "results-$_[0]/WOLAND-hotspots-$_[0]";
	open (SAMPLEHOTSPOT, $samplehotspot);
	my @samplehotspotarray=<SAMPLEHOTSPOT>;
	my $samplehotspotarrayline;
	my (@geneClass, @GeneName, @CHR, @BP, @HotspotCount);

	foreach $samplehotspotarrayline (@samplehotspotarray){
		my @ii = split (/\t/, $samplehotspotarrayline);
		chomp (@ii);
		push (@geneClass, "$ii[0]");
		push (@GeneName, "$ii[1]");
		push (@CHR, "$ii[2]");
		push (@BP, "$ii[3]");
		push (@HotspotCount, "$ii[4]");
	}

	open (HOTSPOT, ">>report-$table/hotspot-$group[$_[1]].txt");

	for my $ii (1..$#geneClass){
		if ($HotspotCount[$ii] ne ""){
			print HOTSPOT "$geneClass[$ii]\t$GeneName[$ii]\t$CHR[$ii]\t$BP[$ii]\t$HotspotCount[$ii]\n";
		}
		else{
		}
	}
	close (HOTSPOT);
}

sub calculate_strandscore{ #strand score calculation
	for my $ii (0..$#samplename){
		my $scorevalue = "results-$samplename[$ii]/WOLAND-bias_$_[0]-$samplename[$ii]";
		open (SCOREVALUE, $scorevalue);
		my @scorevaluearray=<SCOREVALUE>;
		close(SCOREVALUE);
		my @scorechromosome =();
		my @scorecoordinate =();
		my @scoretranscript =();

		foreach my $scorevaluearrayline (@scorevaluearray){
			my @strandscoreline = split (/\t/, $scorevaluearrayline);
			chomp (@strandscoreline);
			push (@scorechromosome, "$strandscoreline[0]");
			push (@scorecoordinate, "$strandscoreline[1]");
			push (@scoretranscript, "$strandscoreline[2]");
		}

		my $scorezero=0;
		my $scoreone=0;
		my $scoreminusone=0;

		for my $iii (0..$#scoretranscript){
			if ($scoretranscript[$iii] == 0){ ++$scorezero;}
			if ($scoretranscript[$iii] == 1){ ++$scoreone;}
			if ($scoretranscript[$iii] eq "-1"){ ++$scoreminusone;}
		}

		my $totalscore=$scorezero+$scoreone+$scoreminusone;

		if ($totalscore==0){
			$totalscore=1;
		}

		my $concordants=$scorezero/$totalscore;
		my $discordants=($scoreone+$scoreminusone)/$totalscore;

		if ($discordants==0){
			$discordants=1;
		}

		my $concordancediscordanceratio=$concordants/$discordants;
		push @{ $motifconcordancediscordanceratio{$pick} } , $concordancediscordanceratio;
	}
	++$pick;

}

sub plot_boxplot_pvalue{ #box plot graph for nucleotide changes and motifs
	$Rbridge->send(qq'WOLAND.$_[0].boxplot<-read.table("$_[0]-$table.txt", sep = "\t", header=TRUE)');
	$Rbridge->send(qq'WOLAND.$_[0].boxplot.m <- melt (WOLAND.$_[0].boxplot)');
	$Rbridge->send(qq'svg("$_[0]_number_boxplot_$table.svg", width=11.692, height=8.267)');
	$Rbridge->send(qq'ggplot (WOLAND.$_[0].boxplot.m, aes(x=variable, y=value, fill = X))+
		ggtitle("$_[0]_$table")+
		theme(plot.title = element_text(size=16, vjust=1.1))+
		theme(axis.title.y = element_text(size=16))+
		scale_x_discrete(name="")+
		geom_boxplot(lwd=0.25)+
		scale_fill_brewer(name="Group", palette="Spectral")');
	$Rbridge->send('dev.off()');
	
	# #p-value nucleotide number absolute
	$Rbridge->send(qq'WOLAND.$_[0].boxplot=read.table("$_[0]-$table.txt", sep = "\t", header=TRUE)');
	$Rbridge->send(qq'melted.boxplot.m <- melt (WOLAND.$_[0].boxplot)');
	$Rbridge->send(q'melted.boxplot.m$X=as.factor(melted.boxplot.m$X)');
	$Rbridge->send(q'melted.boxplot.m$variable=as.factor(melted.boxplot.m$variable)');

	$Rbridge->send(q'groups=permutations(n=length(levels(melted.boxplot.m$X)),r=2,v=levels(melted.boxplot.m$X))');
	$Rbridge->send(q'changes=permutations(n=length(levels(melted.boxplot.m$variable)),r=2,v=levels(melted.boxplot.m$variable), repeats.allowed=TRUE)');

	$Rbridge->send(qq'pvaluetable=matrix(NA,nrow = nrow(changes),ncol = nrow(groups))');

	$Rbridge->send(qq'groupsasnames<-rep(NA,nrow(groups))');
	$Rbridge->send(qq'for (i in 1:nrow(groups)){
		groupsasnames[i]=paste(groups[i,1],groups[i,2],sep = "::")
	}');

	$Rbridge->send(qq'changesasnames<-rep(NA,nrow(changes))');
	$Rbridge->send(qq'for (i in 1:nrow(changes)){
  		changesasnames[i]=paste(changes[i,1],changes[i,2], sep="::")}');

	$Rbridge->send(qq'colnames(pvaluetable)=groupsasnames');
	$Rbridge->send(qq'rownames(pvaluetable)=changesasnames');

	$Rbridge->send(q'for (i in 1:nrow(changes)){
		for (ii in 1:nrow(groups)){
			condition1=(melted.boxplot.m$value[melted.boxplot.m$X==groups[ii,1] & melted.boxplot.m$variable==changes[i,1]])
    		condition2=(melted.boxplot.m$value[melted.boxplot.m$X==groups[ii,2] & melted.boxplot.m$variable==changes[i,2]])
    		result<-rep(NA,5000)
    		result[1]<-diff(append(mean(condition2),mean(condition1)))
    		for(iii in 2:5000){
      			dif.dados=diff(append(mean(sample(append(condition2,condition1),length(condition2))),mean(sample(append(condition2,condition1),length(condition1)))))
      			result[iii]<-dif.dados
    		}
    	unicaudal=sum(result>=result[1])
    	p.uni=unicaudal/length(result)
    	pvaluetable[i,ii]=p.uni
  		}
  	}');

  	if ("$_[0]" eq "motif_number" || "$_[0]" eq "motif_numberNorm"){
  		$Rbridge->send(qq'include_list = c("DNApol::DNApol", "ENU::ENU", "SixFour::SixFour", "SN1::SN1","UV.lambda::UV.lambda","X8.oxoG::X8.oxoG","UV.solar::UV.solar")');
  	}
  	if ("$_[0]" eq "SC_concordance_ratio"){
  		$Rbridge->send(qq'include_list = c("DNApol::DNApol", "SixFour::SixFour", "SN1::SN1","UVlambda::UVlambda","X8.oxoG::X8.oxoG","UV.solar::UV.solar")');
  	}
  	if ("$_[0]" eq "nucleotide_type_change" || "$_[0]" eq "nucleotide_type_changeF"){
		$Rbridge->send(qq'include_list = c("C.A::C.A", "C.T::C.T", "C.G::C.G", "A.C::A.C","A.T::A.T","A.G::A.G")');
	}

	$Rbridge->send(qq'changes_matrix=pvaluetable[include_list, ]');

	$Rbridge->send(qq'my_palette <- colorRampPalette(rev(brewer.pal(3, "RdYlGn")), space="Lab")(n = 3)');
	$Rbridge->send(qq'col_breaks = c(seq(0,0.01,length=1),seq(0.011,0.05,length=1),seq(0.051,0.1,length=1),seq(0.101,1,length=1))');        

	$Rbridge->send(qq'svg("$_[0]_pvalue_samechanges_$table.svg", width=11.692, height=8.267)');
	$Rbridge->send(qq'par(cex.main=0.9)');
	$Rbridge->send(qq'heatmap.2(changes_matrix,
          cellnote = changes_matrix,
          main = "p-value of $_[0] comparisons",
          cexRow=0.9,
          cexCol = 0.9,
          notecol="black",
          density.info="none",
          key.xlab="one sided p-value",
          trace="none",
          margins =c(15,15),
          col=my_palette,
          breaks=col_breaks)');            #turn off column clustering
	$Rbridge->send(qq'dev.off()');
}

sub plot_gaussian{ #kernel density estimation ofr strand scores

	my @color=("\"#9E0142\"", "\"#5E4FA2\"","\"#D53E4F\"","\"#3288BD\"", "\"#F46D43\"","\"#66C2A5\"","\"#FDAE61\"","\"#ABDDA4\"","\"#FEE08B\"",
			 "\"#E6F598\"" );
	my $legendltyvalue="1\,";
	my $legendlwdvalue="2.5\,";
	my $pickgroup=0;

	$Rbridge->send(qq'jpeg("gauss_$_[0]-$table.jpg",width=1024, height=768,units = "px", pointsize = 12, quality = 100)');
	$Rbridge->send(qq'plot(5,5,
			main="Kernel Density Estimation for SC score of $_[0]-$table",
			xlab="SC score",
			ylab="SC score for bw=0.05",
			xlim=c(-1.25, 1.25),
			ylim = c(0,15))');
	
	foreach $uniqgroupline(@uniqgroup){
		my $colorpick=$color[$pickgroup];
		for my $ii (0..$#samplename){
			unless (-z "results-$samplename[$ii]/WOLAND-bias_$_[0]-$samplename[$ii]"){
				if($group[$ii] eq $uniqgroupline){
					$Rbridge->send(qq'FirstPlot <- read.table("../results-$samplename[$ii]/WOLAND-bias_$_[0]-$samplename[$ii]", header=FALSE)');
					$Rbridge->send(q'fp<-density (FirstPlot$V3, bw=0.05)');
					$Rbridge->send(qq'lines (fp, col=$colorpick, lwd=2, ylim=c(0,5))');
				}
			}
		}

		$pickgroup++;
	}

	my $lastgroup=$#uniqgroup;
	my (@legendname, @legendlty,@legendlwd,@legendcol);

	for my $ii(0..$#uniqgroup){
		push (@legendname, "\"$uniqgroup[$ii]\"\,");
		push (@legendlty, "$legendltyvalue");
		push (@legendlwd, "$legendlwdvalue");
		push (@legendcol, "$color[$ii]\,");

		if ($ii==$lastgroup){
			chop ($legendname[$ii]);
			chop ($legendlty[$ii]);
			chop ($legendlwd[$ii]);
			chop ($legendcol[$ii]);
		}
	}
	
	$Rbridge->send(qq'legend("topright",c(@legendname),lty=c(@legendlty), lwd=c(@legendlwd), col=c(@legendcol))');
	$Rbridge->send(q'dev.off()');
}

my $ap = Getopt::ArgParse->new_parser(
	prog => 'woland-report.pl',
	description => 'WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing SNV data.
	Use woland-report to build a grouped report using results-folder of each woland-anno.pl analyzed sample. For more details please read README',
	epilog => 'If you used Woland in your research, we would appreciate your citation:
	de Souza TA, Defelicibus A, Menck CF',
 );

$ap->add_arg(
	'--input-table',
	'-i',
	required => 1,
	help => 'Tab-delimited file with samples in the 1st column and groups in the 2nd column');

my $args = $ap->parse_args();

## main warning
# unless ($#ARGV==0){
# 	die "\nERROR : Incorrect number of arguments - Usage: $0 <input.table> \n\n";	
# }
# unless (-r -e -f $ARGV[0]){
#     die "\nERROR: $ARGV[0] not exists or is not readable or not properly formatted. Please check file.\n\n";
# }
unless (-r -e -f $args->input_table){
	die sprintf("\nERROR: %s not exists or is not readable or not properly formatted. Please check file.\n\n",
		$args->input_table);
}

## parsing input table
$table = $args->input_table; # <input.table>
open (TABLE, $table);
@tablearray=<TABLE>;

foreach $tableline (@tablearray){ # two arrays for each category (group & sample results folder)
	my @i = split (/\t/, $tableline);
	chomp (@i);
	push (@group, "$i[0]"); # array for group definition
	push (@samplename, "$i[1]"); # array for sample folder definition
}

## check if output folder exists & creating report folder

if (-d "report-$table"){
	die "\nERROR: report-$table folder already exists. Check if analysis was performed or remove this folder to repeat analysis.\n\n";
}
else{
	mkdir("report-$table", 0755);
}

print "\nGrouping samples and building reports\n";

###frequency histogram of mutations across chromosomes
for my $i (0..$#samplename){
	&calculate_mutationfrequency ($samplename[$i], $i)
}

#### nucleotide-type change - boxplot & and transition/transvertion ratio pie chart
## building grouped sample input arrays for R using parsed input table
for my $i (0..$#samplename){
	&calculate_nucleotidebasechanges ($samplename[$i]);
}

# ## grouped output file for transvertion transition Ratio
for my $i (0..$#samplename){
	&calculate_transtransvratio($i);
}

## building grouped sample input arrays for R using parsed input table
for my $i (0 .. $#samplename){
	& find_motifs ($samplename[$i]);
}

## grouping hotspot data of samples using parsed input.table
for my $i (0..$#samplename){
	&to_merged_hotspots($samplename[$i], $i);
}

## parsing unique group names for hotspots
%uniquegroup = ();
foreach $item (@group) {
	push(@uniqgroup, $item) unless $uniquegroup{$item}++;
}

#### concordance/discordance SC ratio calc

$input = 6;
for (0..$input) {
   $motifconcordancediscordanceratio{$_} = [];
}

$pick=0;

&calculate_strandscore ("SN1");
&calculate_strandscore ("DNApoln");
&calculate_strandscore ("oxoG");
&calculate_strandscore ("UV-lambda");
&calculate_strandscore ("sixfour");
&calculate_strandscore ("UVsolar");


# ### text outputs

##grouped output file for R analysis of nucleotide-type changes
open (BASECHANGE, ">>report-$table/nucleotide_type_change-$table.txt");
print BASECHANGE "X\tA.T\tA.G\tA.C\tC.G\tC.T\tC.A\n";
for my $i (0 .. $#group){
	print BASECHANGE "$group[$i]\t$AT[$i]\t$AG[$i]\t$AC[$i]\t$CG[$i]\t$CT[$i]\t$CA[$i]\n";
}
close (BASECHANGE);

##grouped output file for R analysis of nucleotide-type changes frequency
open (BASECHANGEFREQUENCY, ">>report-$table/nucleotide_type_changeF-$table.txt");
print BASECHANGEFREQUENCY "X\tA.T\tA.G\tA.C\tC.G\tC.T\tC.A\n";
for my $i (0 .. $#group){
	print BASECHANGEFREQUENCY "$group[$i]\t$ATfrequency[$i]\t$AGfrequency[$i]\t$ACfrequency[$i]\t$CGfrequency[$i]\t$CTfrequency[$i]\t$CAfrequency[$i]\n";
}
close (BASECHANGEFREQUENCY);

##transition-transversion rate calc
open (TRANSITIONTRANSVERSION, ">>report-$table/transitiontransversionF-$table.txt");
for my $i (0 .. $#samplename){
	print TRANSITIONTRANSVERSION "$transversionfrequency[$i]\tTransversion\t$samplename[$i]\t$group[$i]\n$transitionfrequency[$i]\tTransition\t$samplename[$i]\t$group[$i]\n";
}
close (TRANSITIONTRANSVERSION);

## grouped output file for R analysis of number of motifs
open (MOTIFNUMBER, ">>report-$table/motif_number-$table.txt");
print MOTIFNUMBER "X\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUV-solar\n";

for my $i (0 .. $#group){
	print MOTIFNUMBER "$group[$i]\t$SN1[$i]\t$DNApol[$i]\t$oxoG[$i]\t$UV[$i]\t$SixFour[$i]\t$ENU[$i]\t$UVAsolar[$i]\n";
}
close (MOTIFNUMBER);

## grouped output file for R analysis of number of motifs normalized
open (MOTIFNUMBERNORM, ">>report-$table/motif_numberNorm-$table.txt");
print MOTIFNUMBERNORM "X\tSN1\tDNApol\t8-oxoG\tUV-lambda\tSixFour\tENU\tUV-solar\n";

for my $i (0 .. $#group){
	print MOTIFNUMBERNORM "$group[$i]\t$SN1F[$i]\t$DNApolF[$i]\t$oxoGF[$i]\t$UVF[$i]\t$SixFourF[$i]\t$ENUF[$i]\t$UVAsolarF[$i]\n";
}
close (MOTIFNUMBERNORM);

# concordance discordance ratio using strand scores 
open (CONCORDANCEDISCORDANCERATIO, ">>report-$table/SC_concordance_ratio-$table.txt");
print CONCORDANCEDISCORDANCERATIO "X\tSN1\tDNApol\t8-oxoG\tUVlambda\tSixFour\tUV-solar\n";

for my $i (0 .. $#group){
	print CONCORDANCEDISCORDANCERATIO "$group[$i]\t$motifconcordancediscordanceratio{0}[$i]\t$motifconcordancediscordanceratio{1}[$i]\t$motifconcordancediscordanceratio{2}[$i]\t$motifconcordancediscordanceratio{3}[$i]\t$motifconcordancediscordanceratio{4}[$i]\t$motifconcordancediscordanceratio{5}[$i]\n";
}
close (CONCORDANCEDISCORDANCERATIO);

print "\nBuilding graphs...\n";

### R-bridge to graphical outputs
$Rbridge = Statistics::R->new() ;
$Rbridge->start_sharedR ;
$Rbridge->send('library (reshape2)');
$Rbridge->send('library (ggplot2)');
$Rbridge->send('library(qqman)');
$Rbridge->send('library (RColorBrewer)');
$Rbridge->send('library (plyr)');
$Rbridge->send('library (gtools)');
$Rbridge->send('library (gplots)');
$Rbridge->send(qq'setwd(dir = "./report-$table/")');

# box plot of grouped samples
&plot_boxplot_pvalue("nucleotide_type_changeF");
&plot_boxplot_pvalue("nucleotide_type_change");
&plot_boxplot_pvalue("motif_number");
&plot_boxplot_pvalue("motif_numberNorm");

# manhattan plot hotspots
foreach $uniqgroupline(@uniqgroup){

	$Rbridge->send(qq'WOLAND.hotspot.manhattan.$uniqgroupline<- read.delim ("hotspot-$uniqgroupline.txt", comment.char="#", header=FALSE)');
	$Rbridge->send(qq'hotspotdata=WOLAND.hotspot.manhattan.$uniqgroupline');
	$Rbridge->send(q'topvalues=round(mean(hotspotdata$V5)*10)');
	$Rbridge->send(q'df=hotspotdata[order(hotspotdata$V5, decreasing = TRUE),]');
	$Rbridge->send(q'df=df[ave(df$V5, FUN = seq_along) <= topvalues, ]');
	$Rbridge->send(q'topgenes=df$V2');
	$Rbridge->send(qq'topgenes=as.character(topgenes)');
	$Rbridge->send(qq'topgenes=as.factor(topgenes)');
	$Rbridge->send(qq'uniquetopgenes=levels(topgenes)');
	$Rbridge->send(qq'uniquetopgenes=as.character(uniquetopgenes)');

	$Rbridge->send(qq'jpeg("manhattan.hotspot.$uniqgroupline-$table.jpg",width=1024, height=768,units = "px", pointsize = 12, quality = 75)');
	$Rbridge->send(qq'manhattan.$uniqgroupline<- manhattan(x = WOLAND.hotspot.manhattan.$uniqgroupline, chr="V3", bp="V4", p = "V5", logp = FALSE, ylab= "Mutations per bp in the hotspot window", genomewideline = FALSE, suggestiveline = FALSE, main = "Sample:$uniqgroupline", ylim = c(0,50), col = c("blue4", "orange3"))');
	
	$Rbridge->send(qq'mtext("Unique genes above hotspot threshold cutoff",side=3,line=0,cex=0.9)');
	$Rbridge->send(qq'for (i in 1:length(uniquetopgenes)){
		mtext(uniquetopgenes[i],side=3,line=-(i),cex=0.8)
	}');

	$Rbridge->send(q'dev.off()');
}

# gaussian kernel density of SC scores
&plot_gaussian("UV-lambda");
&plot_gaussian("sixfour");
&plot_gaussian("DNApoln");
&plot_gaussian("oxoG");
&plot_gaussian("UVsolar");
&plot_gaussian("SN1");

## boxplot of concordance/discordance SC ratio
&plot_boxplot_pvalue("SC_concordance_ratio");

## chromosome profile mean mutation rate
$Rbridge->send(qq'mutfreq <- read.delim ("mutfreq-$table.txt", comment.char="#", header=FALSE)');
$Rbridge->send(qq'svg("mutfreq-$table.svg", width=11.692, height=8.267)');
$Rbridge->send(qq'melted <- melt(mutfreq,id.vars=c("V1","V2","V3"))');
$Rbridge->send(qq'means = ddply(melted,c("V1","V3"),summarise, mean=mean(V2))');
$Rbridge->send(qq'means.sem = ddply(melted,c("V1","V3"),summarise, mean=mean(V2), sem=sd(V2)/sqrt(length(V2)))');
$Rbridge->send(qq'means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)');
$Rbridge->send(qq'means.group<-ddply(means,"V3",summarise, mean=mean(mean))');
$Rbridge->send(qq'colnames(means.group)<-c("Average","mean")');
$Rbridge->send('ggplot(data=means, aes(x = V1, y = mean, fill=V3))+
	geom_bar(stat="identity", position ="dodge")+
	xlab("Chromosome") + ylab("Mutations per base sequenced") +
	theme(plot.title = element_text(size=16, vjust=1.1))+
	theme(axis.title.y = element_text(size=16))+
	geom_hline(aes(yintercept=mean), data=means.group)+
	scale_fill_brewer(name="group", palette="Spectral")+ 
	geom_errorbar(aes(ymax=upper,ymin=lower), size=.3, width=.2, position=position_dodge(0.9),data=means.sem)+ 
	geom_hline(aes(yintercept=mean, col=Average), data=means.group)');
$Rbridge->send('dev.off()');

## transition and traversion barplot pie
$Rbridge->send(qq'svg("barplot_pie_TransTransv-$table.svg", width=11.692, height=8.267)');
$Rbridge->send(qq'transitiontransversion <- read.delim("transitiontransversionF-$table.txt", header=FALSE)');
$Rbridge->send(qq'ggplot(transitiontransversion, aes(x=V3, y=V1, fill=V2))+
		ggtitle("Transversion & Transition frequency-$table")+
		theme(plot.title = element_text(size=16, vjust=1.1))+
		xlab("Samples") + ylab("Frequency") +
		geom_bar(position="fill", stat = "identity")+
		theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
		facet_grid(.~ V4, scales = "free_x")+
		scale_fill_brewer(name="Type", palette="Spectral")');
$Rbridge->send(q'dev.off()');

$Rbridge->stop;

print "\nDONE\n";

exit;
