# ----------------------------------------------------------
Woland       : Version [1.01]

Author       : Tiago A. de Souza tiagoantonio@gmail.com

Last Modified: 09-30-2017

Credits      : de Souza TA, Defelicibus A, Menck CF

# ----------------------------------------------------------

# LICENSE
Woland is released under GNU Lesser General Public License version 3.0 (LGPLv3) for academic and research use only. Commercial licenses are available to legal entities, including companies and organizations (both for-profit and non-profit), requiring the software for general commercial use. To obtain a commercial license
please, contact author at the address below.

Tiago Antonio de Souza
e-mail: <tiagoantonio@gmail.com>
Institute of Biomedical Sciences, University of Sao Paulo (Brazil)
Av Prof Lineu Prestes 1374
05508 Sao Paulo SP Brazil
Phone :+55 11 30917499

# DISCLAIMER
This software is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.

# DOCUMENTATION AND GETTING HELP
The latest version of the documentation and examples of usage can be found online: [readthedocs](http://woland.readthedocs.io/en/latest).

Please subscribe the woland google forum for questions and help: [google-forum](https://groups.google.com/forum/#!forum/woland-forum).

# WOLAND PREREQUISITES
The following software must be present before installing Woland.

### PERL
Perl is required to run Woland. Minimal recommended Perl version is 5.17.
To check Perl version type:
```sh
$ perl -v 
```
Set your perl environment variable PERL5LIB for modules:

```sh
$ export PERL5LIB=/home/user/perl5/lib
```
Edit ~/.bashrc to include PERL5LIB variable.

### R
You also need R in order to generate reports.
Minimal recommended R version is 3.1.

# INSTALLATION

## Step 1: Download and install source files

### Option 1: Fork woland repository from GitHub
Using Git, clone woland repo:
```sh
$ git clone https://github.com/tiagoantonio/woland
```
### Option 2: Download last Woland source release [here](https://github.com/tiagoantonio/woland/blob/master/woland-source.1.01.tar.gz). Woland is provided as a tar.gz file which could be extracted using:
```sh
$ tar vxzf woland-source.1.01.tar.gz
```
## Step 2: Installing required Perl modules and R libraries

### Option 1: Run INSTALL script 
```sh
./INSTALL [woland install folder complete path]
```
### Option 2:  Manual installation of modules/packages:
Perl modules:
```sh
$ cpan App::cpanminus
$ cpanm Bio::DB::Fasta
$ cpanm Cwd
$ cpanm List::Util
$ cpanm IPC::System::Simple
$ cpanm IPC::Run
$ cpanm Parallel::ForkManager
$ cpanm Regexp::Common
$ cpanm Text::Balanced
(Text::Balanced version should be higher than 1.97)
$ cpanm Text::Wrap
$ cpanm Statistics::R
$ cpanm Moo
$ cpanm Getopt::ArgParse
$ cpanm List::MoreUtils
```
 R packages:
```sh
$ R
$R install.packages("reshape2")
$R install.packages("ggplot2")
$R install.packages("qqman")
$R install.packages("RColorBrewer")
$R install.packages("plyr")
$R install.packages("gtools")
$R q();
```

### Copying genome sequences and annotation
Woland needs genome reference sequence and its gene annotation for each organism. User should download two files in order to perform its analysis.

Download genome reference sequence, move to genomes/ folder and rename it (if necessary):
```sh
$ mv hg19.fa woland/genomes/hg19.fa
```
Download annotation, move to genomes/ folder and rename it:
```sh
$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
$ gunzip refGene.txt.gz
$ mv refGene.txt woland/genomes/refseq_hg19.txt
```

# WOLAND USAGE

### Recommended usage: woland-batch.pl

Woland.batch enables batch submission of multiple samples as provided by `input.table` file. This script runs `woland-anno.pl` for each sample followed by `woland-report.pl` generating one result folder for each sample file provided and a grouped report folder for whole analysis as provided by `input.table`:

```sh
$ woland-batch.pl -i <input.table file> -c <chromosome.profile file> -g <genomes.folder> -n <genome.version> -r <refseq.file> -w <hotspot.window length> -t <number.of.threads> -o <target.output folder>
```
### Parameters explained:

* -i 
input.table file: Path for a regular tabular file without header. First column is group name. Second column is file sample name of annovar annotated.variant file. Samples files MUST be located in the woland folder . 

* -c
chromosome.profile file: Path for a regular tabular file without header . First column is chromosome name in chr format (e.g. chr13). Second column is chromosome length sequenced. Useful for targeted-enrichment sequencing projects. User can build this file with target .BED file using `woland-bed.pl`.

* -g
genomes.folder: Path for genome folder where FASTA sequences are located.

* -n
genome.version: Genome version of genomes.folder/[genome.version].fa and genomes/refseq_[genome.version].txt files.

* -r 
refseq.file: Complete path and file of refseq annotation.

* -w
hotspot.window: A natural number N (N>1), for hotspot window length. Hotspot window corresponds to N nucleotides flanking each SNV.

* -t 
number.of.threads: Default is 30.

### Example:
```sh
$ perl woland-batch.pl -i input.table.tgca.csv -c profiles/chromosome.profile.hg19.bed.exons.txt -w 1000 -g genomes/ -n hg19 -r genomes/refseq_hg19.txt -o .
```
### Auxiliary perl scripts usage:

* woland-bed.pl

This script generates a <chromosome_profile> file using a bed file from a targeted-sequencing experiment. The <chromosome_profile> file is used in other Woland scripts.
```sh
$ perl woland-bed.pl -b <bed-file>
```
* woland-isectoannovar.pl

Use woland-isectoannovar to perform multiple intersections between VCF files from ** HUMAN EXOME** resequencing experiments to select private variants and annotate them using ANNOVAR.

```sh
$ perl woland-isectoannovar.pl -c <type of change> -f <vcf files> -a <annovar folder path> -l <HTSLib/VCFTools folder path> -t <number of threads>
```

* woland-anno.pl

woland-anno.pl uses a single `annovar.variant_function file` to calculate all patterns in a single result folder for `annovar.variant_function file` provided.

```sh
$ perl woland-batch.pl -i <annovar.variant_function> -c <chromosome.profile file> -g <genomes.folder> -n <genome.version> -r <refseq.file> -w <hotspot.window length> -o <target.output folder>
```

* woland-report.pl

This script uses a group of samples to perform `woland-anno.pl` script for each sample provided to build an unique grouped report output folder for the analysis. Sample names and groups are provided by <input.table> and one result folder is also generated for each sample provided. 
```sh
$ perl woland-report.pl -i <input.table file>
```

* woland-mutageniclogo.pl

Use `woland-mutageniclogo.pl` to extract stranded context sequences FASTA files specific for C>T, G>T and G>C changes AFTER woland-report. For more details please read README,

```sh
$ perl woland-mutageniclogo.pl -s <type of 5'3' change> -i <input.table file>
```

# WOLAND OUTPUTS

* results-samplename folder : All results from a single sample file.
* report-samplename folder  : Grouped report files [.txt] and graphic reports from a single [input.table].

# CITATION

If you used Woland in your research, we would appreciate your citation:
de Souza TA, Defelicibus A, Menck CF