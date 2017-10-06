Tutorial
========

The easiest way to perform a WOLAND analysis is through a single batch submission using ``woland-batch.pl``. This will envolve initial 4-step preparation but in the next-time that you will use WOLAND with other samples (and we believe that you do!) you will use only the first step. It is easy no? Each step will prepare the inputs for this script::

$ perl woland-batch.pl -i <input.table file> -c <chromosome.profile file> -g <genomes.folder> -n <genome.version> -r <refseq.file> -w <hotspot.window length> -t <number.of.threads> -o <target.output folder>

First Step - Preparing input-table
----------------------------------

Filtering
~~~~~~~~~

In most cases, a raw .vcf file containing SNVs from a resequencing pipeline is not suitable for a point mutation analysis. First, you have to filter polymorphisms and false-positives from each sample using, for example, ``vcftools`` (<http://vcftools.sourceforge.net/>) and/or ``ANNOVAR`` (<http://annovar.openbioinformatics.org/en/>).

Annotating
~~~~~~~~~~

Several tools are available to annotate .vcf files. However, WOLAND accepts only ``ANNOVAR`` (<http://annovar.openbioinformatics.org/en/>) gene annotation. It is easy to use ANNOVAR and you can find information about downloading, installing and using it at its website. Here is an example how to use annovar to annotate a .vcf file using ``annotate_variation.pl`` from ANNOVAR::

	$ perl annotate_variation.pl -geneanno -buildver hg19 example/ex1.avinput humandb/

.. warning:: WOLAND accepts ONLY ``.variant_function`` files from ANNOVAR. It is not possible to use ``exonic_variant_function`` output.

At this time you have a ``.variant_function`` for each sample to be analyzed. You can manually annotate a file (think twice before) or force annotation when gene information are not avaialble (or not necessary). Let's take a look at a ``variant-function`` file from annovar:

+------------+--------+------+---------+---------+---+---+
|    exonic  |  Lrp1b | chr2 | 3432131 | 3432131 | A | G |
+------------+--------+------+---------+---------+---+---+
| intergenic |  Rbpj  | chr5 |  25465  |  25465  | T | A |
+------------+--------+------+---------+---------+---+---+
|  intronic  | Cmklr1 | chr5 | 4234231 | 4234231 | C | T |
+------------+--------+------+---------+---------+---+---+
|  intronic  |  Setd8 | chr5 | 8423415 | 8423415 | G | C |
+------------+--------+------+---------+---------+---+---+
|     ...    |   ...  |  ... |   ...   |   ...   | . | . |
+------------+--------+------+---------+---------+---+---+


Now you have to build a tabular ``input-table`` file to assign samples into a group name - a "Control" or a "Treated" group, for example.

Grouping samples
~~~~~~~~~~~~~~~~

At this step you must create a simple tabular file (``input-table``). Each line must corresponds to each file sample name in the first column and its group in the second column. Let's see an example:

+----------+------------------------------+
| Control  | Sample1.txt.variant_function |
+----------+------------------------------+
| Control  | Sample2.txt.variant_function |
+----------+------------------------------+
| Treated  | Sample3.txt.variant_function |
+----------+------------------------------+
| Treated  | Sample4.txt.variant_function |
+----------+------------------------------+
| Treated  | Sample5.txt.variant_function |
+----------+------------------------------+

This file ``input-table`` must be saved as a tabular text file and it will be used as the first argument in ``woland-batch.pl`` script.

.. note:: You can provide a path for each file in ``input-table`` if it does not rely on WOLAND ``$install_dir``.

Second Step - Chromosome profile
--------------------------------

At this step you must check you chromosome profile file. This file contains the length of each chromosome of your genome and it is used to calculate frequency of mutational changes. You can manually create your own chromosome profile or use the ``woland-bed.pl`` script if you have a .BED file from you targeted resequencing experiment (exome, for example). Let's see an example of a ``chr_profile`` file:

+------+-----------+
| chr1 | 195471971 |
+------+-----------+
| chr2 | 182113224 |
+------+-----------+
| chr3 | 160039680 |
+------+-----------+
| chr4 | 156508116 |
+------+-----------+
| ...  |    ...    |
+------+-----------+

.. note:: If you have a .BED file from you experiment you can use ``woland-bed.pl``. For example::

	$ perl woland-bed.pl hg19-exome-enrichment.bed

This will create a ``WOLAND-BED-PROFILE-hg19-exome-enrichment.bed`` file which can be used as ``chr_profile`` argument.

Third Step - Genome information
-------------------------------

WOLAND uses genome sequences in FASTA format to extract context sequences and RefSeq annotation to obtain gene and transcriptional information. So you must provide two files for each genome. It is easy to obtain them and you MUST rename the according to <genome_version> parameter and move them to ``$install_dir/genomes/`` folder.


A lot of genome sequences are avaialble nowadays. We advise you to use UCSC genome database to obtain your genome sequence and your RefSeq annotation file. Please check http://hgdownload.cse.ucsc.edu/downloads.html.

Genome sequence in FASTA format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The genome sequence must contain all chromosomes in ``chr`` format. For example::

	>chr1
	AGCATCGATCGGCATGCATGCTAGCTAGCTACGATGCTAGCAT (...)
	>chr2
	GCATGCATCGTACGTACGATCGATCGATCGATCGATCGATCGA (...)
	(...)

Please rename the FASTA file to ``genome_<genome_version>.fa`` and move it to ``$install_dir/genomes/``. For example::

	$ mv hg19.fa $install_dir/genomes/hg19.fa

RefSeq annotation
~~~~~~~~~~~~~~~~~

The RefSeq annotation can be obtained through http://hgdownload.cse.ucsc.edu/downloads.html . 

.. note:: You MUST download the RefGene file - usually provided as ``refGene.txt``.

Please rename the RefGene file to ``refseq_<genome_version>.txt` and move it to ``$install_dir/genomes/``. For example::

	$ mv RefGene $install_dir/genomes/refseq_hg19.txt

Fourth Step - Choosing hotspot window length and running!
---------------------------------------------------------

Now you can choose a natural number >1 for the hotspot window length ``<hotspot_window>``, for example: 1000. Now, voilà, you can run ``woland-batch.pl``!::

	$ perl woland-batch.pl -i input.table.tgca.csv -c profiles/chromosome.profile.hg19.bed.exons.txt -w 1000 -g genomes/ -n hg19 -r genomes/refseq_hg19.txt -o .

