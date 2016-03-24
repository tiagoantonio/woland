Galaxy
======

WOLAND is implemented in a web-friendly Galaxy environment at <link>. This implementation covers only batch submissions (more than two sample groups) and a limited number of genomes (hg19, hg18, mm10, mm9). If you are interested to use WOLAND in a custom-user way please consider to install it locally.

Using WOLAND @ Galaxy
--------------------- 

You must access our Galaxy server at <link>. You can enter WOLAND page clicking in the left-side panel. Let's take a look at WOLAND main page:

<image>

Adding ANNOVAR .variant_function files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you already have ANNOVAR .variant_function files you must submit them to Galaxy. Thus, you must assign each file submitted to a group name.

Adding chromosome profile
~~~~~~~~~~~~~~~~~~~~~~~~~

You must submit a chromosome profile file containing the total length of each chromosome to be considered in frequency calculations.

Selecting genome
~~~~~~~~~~~~~~~~

Please select the genome among organisms available. If you did not find a genome for your organism of interest please consider to install a local version of WOLAND.

Choosing hotspot window length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now you must enter a natural number > 1 which will correspond to the radius of the hotspot window length.

Starting WOLAND
---------------

Click EXECUTE and wait analysis to finish

WOLAND web-report
-----------------

After finishing you may able to access a tar.gz result file containing all woland outputs. In addition you can view and save a small web-report containing some report graphs.
