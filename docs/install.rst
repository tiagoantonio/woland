Install
========

WOLAND is a multiplatform tool based on Perl and R. Please observe prerequisites and modules and libraries need. 
They must be installed before WOLAND installation. 

Prerequisites
-------------

The following software must be present before installing Woland:

--- Perl (Minimal recommended Perl version: 5.17)

To check Perl version type:

$ perl -v 

The following Perl modules should be present:

- Bio::DB::Fasta
- Cwd
- List::Util
- IPC::System::Simple
- IPC::Run
- Parallel::ForkManager
- Regexp::Common
- Text::Balanced(>=1.97)
- Text::Wrap
- Statistics::R

--- R (Minimal recommended R version: 3.1)

The following R packages must be installed:

- reshape2
- ggplot2
- qqman
- RColorBrewer
- plyr


Installation
-------------

--- Installing Perl modules

You can use CPAN to get any missing module. First install cpanm:
$ sudo cpan App::cpanminus

Them you can install each module:
$ sudo cpanm Module::Name

--- Installing R packages

You can use this following command in R to install each missing library:

$R install.packages("packagename")

--- Installing source WOLAND files

Download last Woland source release in <link>. Woland is provided as a tar.gz file which could be extracted using, for example:

$ tar vxzf <file>

Woland will be installed inside woland folder.

--- Copying genome sequences and annotation

Woland needs genome reference sequence and its gene annotation for each organism. User should download two files in order to perform its analysis:

- Woland-1.0/genomes/genome_<genome_version>.fa - Can be downloaded using UCSC genome sequences database.
- Woland-1.0/refGene/<refseq_<genome_version>.txt - Can be downloaded using UCSC annotation database.

You must rename genome.fa and refGene.txt according to <genome_version>. For example:

- hg19

$install_dir/genomes/genome_hg19.fa

$install_dir/genomes/refseq_hg19.fa