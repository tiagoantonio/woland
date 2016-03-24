README
========

WOLAND is a multiplatform tool to analyze point mutation patterns using resequencing data from any organism or cell. 
It is implemented as a Perl and R tool using as inputs filtered unannotated or annotated SNV lists, combined with its 
correspondent genome sequences.

What do you need to use woland:

    annotated SNV list(s)
    genome sequence and annotation
    length of each chromosome or a bed file from a targeted-resequencing experiment (exome, for example)

Features
--------

WOLAND can provide:

- the number and frequency of nucleotide type changes
- detection of regions enriched in mutations alongside the genome (hotspots)
- extraction of sequence-context sequences of each SNV. 
- count established mutational motifs associated with environmental mutagens and DNA-repair mechanisms
- calculation of transcriptional strand bias of mutations linked to the mutational motifs found.

Contribute
----------

- Issue Tracker: github.com/woland
- Source Code: github.com/woland

Support
-------

If you are having issues, please let us know.
We have a mailing list located at: woland@google-groups.com

License
-------

Woland is released under GNU Lesser General Public License version 3.0 (LGPLv3) for academic and research use only. Commercial licenses are available to legal entities, including companies and organizations (both for-profit and non-profit), requiring the software for general commercial use. To obtain a commercial license
please, contact author **tiagoantonio[at]gmail.com**
