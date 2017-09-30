.. woland documentation master file, created by
   sphinx-quickstart on Tue Mar  8 10:18:24 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to WOLAND's documentation!
==================================

WOLAND is a tool to analyze point mutation patterns using resequencing data from any organism or cell. There are many suitable applications for WOLAND such as the study of the genome-wide impact of both endogenous and exogenous mutagens on organisms and cells and analysis of mutagenic profiles of DNA repair-associated diseases. Also, mechanisms of molecular mutational processes involving specific target proteins and identification of potential hazardous mutagens in environmental samples can be profiled using WOLAND.

WOLAND retrieves a comprehensive visual report with R-based graphics to enable fast and reliable interpretation by the user.

Given one or more list of SNVs in ANNOVAR variant_function format, WOLAND can perform:

- **Nucleotide type changes identification and context-sequence extration**: identify frequency of nucleotide type changes and extract context-sequences around each point mutation.
- **Search for mutagen-associated motifs**: Retrieve the number, frequency and localization of mutagen-associated motif sequences such as UV-light, 8-oxoguanine, 6-4 photoproducts, ENU, among others.
- **Search for mutational hotspots**: Identify mutational hotspots using a user-defined sliding window which considers each SNV as the window center.
- **Transcriptional strand bias**: Retrieve scores associated with the strand of each mutational motif found using RefSeq annotation data.

You can use Galaxy-WOLAND directly through your web navigator or install a local version. Galaxy-WOLAND is a friendly and easy interface while installing a local version will allow you to customize all analyses steps. Please see docs for further information!

The documentation of WOLAND is organized into:

.. toctree::
   :maxdepth: 2

   README
   install
   using
   tutorial
   outputs
   report
   faq
   galaxy
   contact

References
==========


