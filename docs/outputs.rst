Outputs
=======

WOLAND provide different set of outputs depending on the script used - and how far you want to go in your mutational pattern analysis. ``woland-batch.pl`` runs the most complete analysis but you must have at least two groups of samples (a control and a treated group, for example). 

.. note:: ``woland-batch.pl`` provides the faster and straighforward way to analyze samples. It runs ``woland-anno.pl`` for each sample then runs ``woland-report.pl``.

woland-anno.pl outputs
----------------------

This script will analyze each ANNOVAR ``variant_function`` files and will provide a total of 13 tabular text files + 1 log file. We consider these output files as the most raw type of WOLAND analysis. Let's take a look at each class and its output files explanations:

.. note:: All ``woland-anno.pl`` outputs were saved in a ``results-$sample_name/`` folder.

Nucleotide type-changes and frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``WOLAND-basechange-$sample_name``: 

``WOLAND-mutfreq-$sample_name``: 

Extracted context sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``WOLAND-contextsequences-$sample_name``: 

``WOLAND-contextsequencesanno-$sample_name``: 

Hotspots
~~~~~~~~

``WOLAND-hotspots-$sample_name``: 

Mutational motifs
~~~~~~~~~~~~~~~~~

``WOLAND-motifs-$sample_name``: 

``WOLAND-norm_motifs-$sample_name``: 

Transcriptional strand bias
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``WOLAND-bias_motif-$sample_name``: 

woland-report.pl outputs
------------------------

``woland-report.pl``uses outputs from ``woland-anno.pl`` to build a comprehensive grouped analysis with some graphical data. Those outputs are saved as tabular text files (``.tmp``) and graphical data as .SVG files.

.. note:: All ``woland-report.pl`` outputs were saved in a ``report-$input_table/`` folder.

Nucleotide type-changes and frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nucleotide_type_change-$input_table.tmp``:

``nucleotide_type_changeF-$input_table.tmp``:

.. image:: /img/nucleotide_freq.svg

``mutfreq-$input_table.tmp``:

.. image:: /img/mutfreq.svg

``transitiontransversionF-$input_table.tmp``:

.. image:: /img/barplot_trans_transv.svg

Extracted context sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hotspots
~~~~~~~~

``hotspots-$group.tmp``:

Mutational motifs
~~~~~~~~~~~~~~~~~

``motif_number-$input_table.tmp`` :

``motif_numberNorm-$input_table.tmp`` :

.. image:: /img/motif_norm.svg

Transcriptional strand bias
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SC_concordance_ratio-$input_table.tmp``:

woland-batch.pl outputs
-----------------------