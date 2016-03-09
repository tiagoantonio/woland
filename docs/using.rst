Using
========

WOLAND has three main scripts: ``woland-anno.pl``, ``woland-report.pl`` and ``woland-batch.pl``) plus a acessory script ``woland-bed.pl``. Please observe descriptions of each script, input file file requirements, usage and outputs.


Scripts
--------

woland-anno.pl
~~~~~~~~~~~~~~

Script used to calculate mutational patterns of a single annotated variant file. Uses a single ``<variant_function>`` to calculate all patterns in a single result folder for ``<results-variant_function>`` provided::

	$ perl woland-anno.pl <variant_function> <chr_profile> <hs_window> <genome_version>

woland-report.pl
~~~~~~~~~~~~~~~~

Script used to build a report of multiple annotated variant files assigned as groups.

woland-batch.pl
~~~~~~~~~~~~~~~

Script which automatically runs multiple instances of woland-anno.pl and build a single report using woland-report.pl  


Inputs
--------

- ``<annovar.variant_function>``: Annotated .variant_function file which can be obtained using annotate-variation.pl script from ANNOVAR.

- ``<chromosome_profile>``: Regular tabular file without header . First column is chromosome name in chr format (e.g. chr13). Second column is chromosome length sequenced. User can build this file with target .BED file using ``woland-bed-to-profile.pl``.

- ``<hotspot_window>``: A natural number N (N>1), for hotspot window length. Hotspot window corresponds to N nucleotides flanking each SNV.

- ``<genome_version>``:Genome version of genomes/genome_<genome_version> and genomes/refseq_<genome_version> files.

- ``<input.table>``: Regular tabular file without header. First column is group name. Second column is file sample name of annovar annotated.variant file. Samples files MUST be located in the Woland root folder. 

- ``<coordinates.bed>``:Coordinates of target regions used in sequencing experiment in BED format.


Usage
--------

woland-anno.pl
~~~~~~~~~~~~~~

Uses a single <annovar.variant_function> to calculate all patterns in a single result folder for <annovar.variant_function> provided::

	$ woland-anno.pl <annovar.variant_function> <chromosome_profile> <hotspot_window> <genome_version>

woland-report.pl
~~~~~~~~~~~~~~~~

This script uses a group of samples to perform Woland-anno.pl script for each sample provided to build an unique grouped report output folder for the analysis. Sample names and groups are provided by <input.table> and one result folder is also generated for each sample provided::

	$ woland-report.pl <input.table> <chromosome_profile> <hotspot_window> <genome_version>

woland-batch.pl
~~~~~~~~~~~~~~~

woland.batch enables batch submission of multiple samples as provided by <input.table> file. This script runs Woland-anno.pl for each sample followed by Woland-report.pl generating one result folder for each sample file provided and a grouped report folder for whole analysis as provided by <input.table>::

	$ woland-batch.pl <input.table> <chromosome_profile> <hotspot_window> <genome_version>

woland-bed.pl
~~~~~~~~~~~~~ 

This script generates a <chromosome_profile> file using a .bed file from a targeted-sequencing experiment, for example. The <chromosome_profile> file could be used in other Woland scripts::

$ woland-bed.pl <coordinates.bed>




Outputs
--------

<results-samplename> folder : All results from a single sample file.

<report-samplename> folder  : Grouped report files <.tmp> and graphic reports from a single <input.table> 