![diamond](http://www.diamondsearch.org/diamond_white_95px.png)

Introduction
============

DIAMOND is a sequence aligner for protein and translated DNA searches,
designed for high performance analysis of big sequence data. The key
features are:

-   Pairwise alignment of proteins and translated DNA at 100x-10,000x
    speed of BLAST.
-   [Protein clustering of up to tens of billions of proteins](https://github.com/bbuchfink/diamond/wiki/Clustering)
-   Frameshift alignments for long read analysis.
-   Low resource requirements and suitable for running on standard
    desktops or laptops.
-   Various output formats, including BLAST pairwise, tabular and XML,
    as well as taxonomic classification.

[![Build](https://github.com/bbuchfink/diamond/actions/workflows/cmake.yml/badge.svg)](https://github.com/bbuchfink/diamond/actions/workflows/cmake.yml)
[![image](https://img.shields.io/cirrus/github/bbuchfink/diamond)](https://cirrus-ci.com/github/bbuchfink/diamond/master)
[![image](https://img.shields.io/github/downloads/bbuchfink/diamond/total)](https://github.com/bbuchfink/diamond/releases)
[![image](https://anaconda.org/bioconda/diamond/badges/version.svg)](https://anaconda.org/bioconda/diamond)
[![image](https://anaconda.org/bioconda/diamond/badges/downloads.svg)](https://anaconda.org/bioconda/diamond)
[![image](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fraw.githubusercontent.com%2Fbbuchfink%2Fdiamond-data%2Frefs%2Fheads%2Fmain%2Fcitations.json&query=%24.citations&style=flat&label=Citations&color=%23a020f0
)](https://scholar.google.com/citations?user=kjPIF1cAAAAJ)
    
Documentation
=============
The online documentation is located at the [GitHub Wiki](https://github.com/bbuchfink/diamond/wiki).

Support
=======
DIAMOND is actively supported and developed software. Please use the [issue tracker](https://github.com/bbuchfink/diamond/issues) for malfunctions and the [GitHub discussions](https://github.com/bbuchfink/diamond/discussions) for questions, comments, feature requests, etc.

About
=====
DIAMOND is developed by Benjamin Buchfink at the Max Planck Institute for Biology
Tübingen in collaboration with the Drost lab at the University of Dundee. From 2019-2024,
it was developed by Benjamin Buchfink at the Drost lab, Max Planck Institute for Biology
Tübingen. From 2018-2019, its development was supported by the German Federal Ministry
for Economic Affairs and Energy through an EXIST grant. From 2016-2018, it was developed
by Benjamin Buchfink as an independent researcher. From 2013-2015, the initial version
was developed by Benjamin Buchfink at the Huson lab, University of Tübingen, Germany.

\[[:email:Email](mailto:buchfink@gmail.com)\]
\[[X](https://x.com/bbuchfink)\]
\[[Bluesky](https://bsky.app/profile/bbuchfink.bsky.social)\]
\[[LinkedIn](https://www.linkedin.com/in/benjamin-buchfink-875692105/)\]
\[[Google Scholar](https://scholar.google.de/citations?user=kjPIF1cAAAAJ)\]
\[[Drost lab](https://drostlab.com/)\]
\[[MPI-BIO](https://www.bio.mpg.de/)\]

**When using the tool in published research, please cite:**

-   Buchfink B, Reuter K, Drost HG, \"Sensitive protein alignments at tree-of-life
    scale using DIAMOND\", *Nature Methods* **18**, 366–368 (2021).
    [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)

For sequence clustering:

-   Buchfink B, Ashkenazy H, Reuter K, Kennedy JA, Drost HG, \"Sensitive clustering
    of protein sequences at tree-of-life scale using DIAMOND DeepClust\", *bioRxiv*
    2023.01.24.525373; doi: https://doi.org/10.1101/2023.01.24.525373 

Original publication to cite DIAMOND until v0.9.25:

-   Buchfink B, Xie C, Huson DH, \"Fast and sensitive protein alignment
    using DIAMOND\", *Nature Methods* **12**, 59-60 (2015).
    [doi:10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)
