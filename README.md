<img src="/../../../../bbuchfink/diamond-data/blob/main/diamond_coat_of_arms2.png" height=280>

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
[![Conda](https://img.shields.io/conda/dn/bioconda/diamond.svg)](https://anaconda.org/bioconda/diamond/files)
[![image](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fraw.githubusercontent.com%2Fbbuchfink%2Fdiamond-data%2Frefs%2Fheads%2Fmain%2Fcitations.json&query=%24.citations&style=flat&label=Citations&color=%23a020f0
)](https://scholar.google.com/citations?user=kjPIF1cAAAAJ)
[![image](https://img.shields.io/badge/Live%20chat-Discord-red)](https://discord.gg/ptJnz3GSCy)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=bg_diamond)


Documentation
=============
The online documentation is located at the [GitHub Wiki](https://github.com/bbuchfink/diamond/wiki).

Support
=======
DIAMOND is actively supported and developed software. Please use the [issue tracker](https://github.com/bbuchfink/diamond/issues) for malfunctions and the [GitHub discussions](https://github.com/bbuchfink/diamond/discussions) for questions, comments, feature requests, etc. I also provide live support on [Discord](https://discord.gg/ptJnz3GSCy). Don't be too shy to ask!

About
=====
DIAMOND is developed by Dr. Benjamin J. Buchfink, independent scientist, Tübingen, Germany,
supported by the Max Planck Society for the Advancement of Science, in collaboration with
the Drost lab at the University of Dundee. From 2019-2024,
it was developed by Benjamin Buchfink at the Drost lab, Max Planck Institute for Biology
Tübingen. From 2018-2019, its development was supported by the German Federal Ministry
for Economic Affairs and Energy through an EXIST grant. From 2016-2018, it was developed
by Benjamin Buchfink as an independent researcher. From 2013-2015, the initial version
was developed by Benjamin Buchfink at the Huson lab, University of Tübingen, Germany.

\[[:email:Email](mailto:buchfink@gmail.com)\]
\[[X](https://x.com/bbuchfink)\]
\[[Bluesky](https://bsky.app/profile/bbuchfink.bsky.social)\]
\[[LinkedIn](https://www.linkedin.com/in/benjamin-j-buchfink-875692105/)\]
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
