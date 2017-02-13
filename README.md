# ----------------------------------------------------------
Woland       : Version [version]

Author       : Tiago A. de Souza tiagoantonio@gmail.com

Last Modified: [data]

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
The latest version of the documentation and examples of usage can be found online: <link>.
Please subscribe our mailing list for questions and help: <link>

# WOLAND PREREQUISITES
The following software must be present before installing Woland.

### PERL
Perl is required to run Woland. Minimal recommended Perl version is 5.17.
To check Perl version type:
```sh
$ perl -v 
```

If you do not have them installed you can use CPAN to get them.
The following Perl modules should be installed:
```sh
$ cpan App::cpanminus
$ cpanm Bio::DB::Fasta
$ cpanm Cwd
$ cpanm List::Util
$ cpanm IPC::System::Simple
$ cpanm IPC::Run
$ cpanm Parallel::ForkManager
$ cpanm Regexp::Common
# Text::Balanced version should be higher than 1.97
$ cpanm Text::Balanced
$ cpanm Text::Wrap
$ cpanm Statistics::R
```

### R
You also need R in order to generate reports.
Minimal recommended R version is 3.1.

The following R packages must be installed:
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

# INSTALLATION

# Installing source Woland files
Download last Woland source release in [link]. Woland is provided as a tar.gz file which could be extracted using, for example:
```sh
$ tar vxzf <file>
```

Woland will be installed inside woland/ folder.

### Copying genome sequences and annotation
Woland needs genome reference sequence and its gene annotation for each organism. User should download two files in order to perform its analysis:

```sh
# Can be downloaded using UCSC genome sequences database.
$ woland/genomes/genome_<genome_version>.fa
# Can be downloaded using UCSC annotation database.
$ woland/refGene/<refseq_<genome_version>.txt
```

# WOLAND USAGE

### WOLAND scripts
1) woland-anno
```sh
# woland-anno.pl uses a single <annovar.variant_function> to calculate all patterns in a single result folder for <annovar.variant_function> provided.
$ woland-anno.pl <annovar.variant_function> <chromosome_profile> <hotspot_window> <genome_version>
```

2) woland-report
```sh
# This script uses a group of samples to perform Woland-anno.pl script for each sample provided to build an unique grouped report output folder for the analysis. Sample names and groups are provided by <input.table> and one result folder is also generated for each sample provided. 
$ woland-report.pl <input.table> <chromosome_profile> <hotspot_window> <genome_version>
```

3) woland-batch
```sh
# Woland.batch enables batch submission of multiple samples as provided by <input.table> file. This script runs Woland-anno.pl for each sample followed by Woland-report.pl generating one result folder for each sample file provided and a grouped report folder for whole analysis as provided by <input.table>
$ woland-batch.pl <input.table> <chromosome_profile> <hotspot_window> <genome_version>
```

4) woland-bed
```sh
# This script generates a <chromosome_profile> file using a .bed file from a targeted-sequencing experiment, for example. The <chromosome_profile> file could be used in other Woland scripts.
$ woland-bed.pl <coordinates.bed>
```

# WOLAND INPUTS

* annovar.variant_function: Annotated .variant_function file which can be obtained using annotate-variation script from ANNOVAR.

* chromosome_profile: Regular tabular file without header . First column is chromosome name in chr format (e.g. chr13). Second column is chromosome length sequenced. User can build this file with target .BED file using woland-bed-to-profile.

* hotspot_window: A natural number N (N>1), for hotspot window length. Hotspot window corresponds to N nucleotides flanking each SNV.

* genome_version: Genome version of genomes/genome_[genome_version] and genomes/refseq_[genome_version] files.

* input.table: Regular tabular file without header. First column is group name. Second column is file sample name of annovar annotated.variant file. Samples files MUST be located in the Woland root folder. 

* coordinates.bed: Coordinates of target regions used in sequencing experiment in BED format.

# WOLAND OUTPUTS

* results-samplename folder : All results from a single sample file.
* report-samplename folder  : Grouped report files [.tmp] and graphic reports from a single [input.table].

# CITATION

If you used Woland in your research, we would appreciate your citation:
de Souza TA, Defelicibus A, Menck CF
