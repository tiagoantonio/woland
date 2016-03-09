Tutorial
========

The easiest way to perform a WOLAND analysis is through a single batch submission using ``woland-batch.pl``. 

First Step - Preparing VCF files
--------------------------------

Filtering
~~~~~~~~~

In most cases, a raw .vcf file containing SNVs from a resequencing pipeline is not suitable for a point mutation analysis. First, you have to filter polymorphisms and false-positives from each sample using, for example, ``vcftools`` (<http://vcftools.sourceforge.net/>) and/or ``ANNOVAR`` (<http://annovar.openbioinformatics.org/en/>).

Annotating
~~~~~~~~~~

Several tools are available to annotate .vcf files. However, WOLAND accepts only ``ANNOVAR`` (<http://annovar.openbioinformatics.org/en/>) gene annotation. It is easy to use ANNOVAR and you can find information about downloading, installing and using it at its website. Here is an example how to use annovar to annotate a .vcf file using ''annotate_variation.pl``from ANNOVAR::

	$ perl annotate_variation.pl -geneanno -buildver hg19 example/ex1.avinput humandb/

.. warning:: WOLAND accepts ONLY ``.variant_function``files from ANNOVAR. It is not possible to use ``exonic_variant_function``output.

At this time you have a ``.variant_function`` for each sample to be analyzed. Now you have to build a tabular ``input-table``file to assign samples into a group name - a "Control" or a "Treated" group, for example.

Second Step - Grouping samples
------------------------------

At this step you must create a simple tabular file (``input-table``). Each line must corresponds to each file sample name in the first column and its group in the second column. Let's see an example:

 


The following software must be present before installing Woland:

Perl (Minimal recommended Perl version: 5.17)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To check Perl version type::

	$ perl -v 

The following Perl modules should be present:

	* Bio::DB::Fasta
	* Cwd
	* List::Util
	* IPC::System::Simple
	* IPC::Run
	* Parallel::ForkManager
	* Regexp::Common
	* Text::Balanced(>=1.97)
	* Text::Wrap
	* Statistics::R

R (Minimal recommended R version: 3.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To check R version type (in R)::

	$ R.Version()

The following R packages must be installed:

	* reshape2
	* ggplot2
	* qqman
	* RColorBrewer
	* plyr


Installation
-------------

Installing Perl modules
~~~~~~~~~~~~~~~~~~~~~~~

You can use CPAN to get any missing module. First install cpanm::

	$ sudo cpan App::cpanminus

Them you can install each module::

	$ sudo cpanm Module::Name

Installing R packages
~~~~~~~~~~~~~~~~~~~~~

You can use this following command in R to install each missing library:

$R install.packages(``packagename``)

Installing source WOLAND files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download last Woland source release in <link>. Woland is provided as a tar.gz file which could be extracted using, for example::

	$ tar -xvzf woland-<version>-install.tar.gz

Woland will be installed inside woland installation folder ``install_dir``.

Copying genome sequences and annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Woland needs genome reference sequence and its gene annotation for each organism in the ``$install_dir/genomes`` folder. User should download two files in order to perform its analysis:

.. note:: These files can be downloaded using UCSC database (http://hgdownload.cse.ucsc.edu/downloads.html). See below how to rename them::

	$install_dir/genomes/genome_<genome_version>.fa
	$install_dir/genomes/refseq_<genome_version>.txt

.. warning:: You must rename ``<genome_version>.fa`` and refGene.txt according to ``<genome_version>``. For example:

- hg19::

	$install_dir/genomes/genome_hg19.fa

	$install_dir/genomes/refseq_hg19.fa

- mm10::

	$install_dir/genomes/genome_mm10.fa

	$install_dir/genomes/refseq_mm10.fa
