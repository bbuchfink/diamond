# Command line options

## Commands

Commands are issued as the first parameter on the command line and set
the task to be run by the program.

  - `makedb`  
      
    Create a DIAMOND formatted reference database from a FASTA input
    file.

  - `blastp`  
      
    Align protein query sequences against a protein reference database.

  - `blastx`  
      
    Align translated DNA query sequences against a protein reference
    database.

  - `view`  
      
    Generate formatted output from DAA files.

  - `version`  
      
    Print version information.

  - `dbinfo`  
      
    Print information about a database file.

  - `help`  
      
    Print help message.
    
  - `test`
  
    Run a series of test cases and verify the output against reference
    hashes. This command will exit with code `0` if all
    tests have passed and `1` otherwise. Running this command requires write
    access to the current working directory.

## Makedb options

  - `--in <file>`  
      
    Path to the input protein reference database file in FASTA format
    (may be gzip compressed). If this parameter is omitted, the input
    will be read from `stdin`.

  - `--db/-d <file>`  
      
    Path to the output DIAMOND database file.

  - `--taxonmap <file>`  
      
    Path to mapping file that maps NCBI protein accession numbers to
    taxon ids (gzip compressed). This parameter is optional and needs to
    be supplied in order to provide taxonomy features. The file can be
    downloaded from NCBI: <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz>.
    
    A custom file following the same format may be supplied here. Note
    that the first line of this file is assumed to contain headings and
    will be ignored.

  - `--taxonnodes <file>`  
      
    Path to the `nodes.dmp` file from the NCBI taxonomy. This parameter
    is optional and needs to be supplied in order to provide taxonomy
    features. The file is contained within this archive downloadable at
    NCBI: <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip>.
    
  - `--taxonnames <file>`
    
    Path to the `names.dmp` file from the NCBI taxonomy. This parameter
    is optional and needs to be supplied in order to provide taxonomy
    features. The file is contained within this archive downloadable at
    NCBI: <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip>.

## General options

  - `--threads/-p #`  
      
    Number of CPU threads. By default, the program will auto-detect and
    use all available virtual cores on the machine.
    
  - `--quiet`
  
    Disable all terminal output.
    
  - `--verbose/-v`
  
    Enable more verbose terminal output.
    
  - `--log`
  
    Enable even more verbose terminal output, which is also written to
    a file named `diamond.log` is the current working directory.

## Input options

  - `--db/-d <file>`  
      
    Path to the DIAMOND database file.

  - `--query/-q <file>`  
      
    Path to the query input file in FASTA or FASTQ format (may be gzip
    compressed). If this parameter is omitted, the input will be read
    from `stdin`.

  - `--taxonlist <list>`  
      
    Comma-separated list of NCBI taxonomic IDs to filter the database
    by. Any taxonomic rank can be used, and only reference sequences
    matching one of the specified taxon ids will be searched against.
    Using this option requires setting the `--taxonmap` and `--taxonnodes`
    parameters for `makedb`.

  - `--query-gencode #`  
      
    Genetic code used for translation of query in BLASTX mode. A list of
    possible values can be found at the [NCBI
    website](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
    By default, the Standard Code is used. Note: changing the genetic
    code is currently not fully supported for the DAA format.

  - `--strand {both, plus, minus}`  
      
    Set strand of query to align for translated searches. By default
    both strands are searched.

  - `--min-orf/-l #`  
      
    Ignore translated sequences that do not contain an open reading
    frame of at least this length. By default this feature is disabled
    for sequences of length below 30, set to 20 for sequences of length
    below 100, and set to 40 otherwise. Setting this option to `1` will
    disable this feature.

## Alignment options

  - `--mid-sensitive`
  
    Enable the mid-sensitive mode which is between the default (fast) mode
    and the sensitive mode in sensitivity. *Option supported since v2.0.3*
  
    Without using any sensitivity option, the default (fast) mode will run which
    is designed for finding hits of >70% identity and short read alignment.

  - `--sensitive`
  
    Enable the sensitive mode designed for full sensitivity for hits >40% identity.
    
  - `--more-sensitive`  
      
    This mode is slightly more sensitive than the `--sensitive` mode.
    
  - `--very-sensitive`
  
    Enable the very-sensitive mode designed for best sensitivity including the
    twilight zone range of <40% identity.
    
    **Available since version 2.0.0.**
    
  - `--ultra-sensitive`
    
    Enable the ultra-sensitive mode which is yet more sensitive than the
    `--very-sensitive` mode.
    
    **Available since version 2.0.0.**

  - `--frameshift/-F #`  
      
    Penalty for frameshifts in DNA-vs-protein alignments. Values around
    15 are reasonable for this parameter. Enabling this feature will
    have the aligner tolerate missing bases in DNA sequences and is most
    recommended for long, error-prone sequences like MinION reads.  
    In the pairwise output format, frameshifts will be indicated by `\`
    and `/` for a shift by +1 and -1 nucleotide in the direction of
    translation respectively.
    Note that this feature is disabled by default.

  - `--gapopen #`  
      
    Gap open penalty.

  - `--gapextend #`  
      
    Gap extension penalty.

  - `--matrix <matrix name>`  
      
    Scoring matrix. The following matrices are supported, with the
    default being
    BLOSUM62.  
      
    
    | Matrix   | Supported values for (gap open)/(gap extend) | Default gap penalties |
    | :------- | :------------------------------------------- | :-------------------- |
    | BLOSUM45 | (10-13)/3; (12-16)/2; (16-19)/1              | 14/2                  |
    | BLOSUM50 | (9-13)/3; (12-16)/2; (15-19)/1               | 13/2                  |
    | BLOSUM62 | (6-11)/2; (9-13)/1                           | 11/1                  |
    | BLOSUM80 | (6-9)/2; 13/2; 25/2; (9-11)/1                | 10/1                  |
    | BLOSUM90 | (6-9)/2; (9-11)/1                            | 10/1                  |
    | PAM250   | (11-15)/3; (13-17)/2; (17-21)/1              | 14/2                  |
    | PAM70    | (6-8)/2; (9-11)/1                            | 10/1                  |
    | PAM30    | (5-7)/2; (8-10)/1                            | 9/1                   |
    
  - `--masking (0,1)`
  
    DIAMOND by default applies the tantan repeat masking algorithm to the
    query and target sequences as described in [1]. This masking procedure
    increases the specificity of alignments and serves to filter out
    spurious hits. If this is not desired, repeat masking can be disabled
    using `--masking 0`.
  
  - `--comp-based-stats (0,1)`  
      
    Compositional bias correction of alignment scores. 0 means no score
    correction, 1 means compositional bias correction as described in 
    [2]. Compositionally biased sequences often cause false
    positive matches, which are effectively filtered by this algorithm
    in a way similar to the composition based statistics used by BLAST .

  - `--algo (0,1)`  
      
    Algorithm for seed search. 0 means double-indexed and 1 means
    query-indexed. The double-indexed algorithm is the program’s main
    algorithm, but it is inefficient for very small query files, where
    the query-indexed algorithm should be used instead.  
    By default, the program will automatically choose one of the
    algorithms based on the size of the query and database files. The
    algorithm used will be displayed at program startup.  
    Note that while the two algorithms are configured to provide roughly
    the same sensitivity for the respective modes, results will not be
    exactly identical to each other.

## Output options

  - `--out/-o <file>`  
      
    Path to the output file. If this parameter is omitted, the results
    will be written to the standard output and all other program output
    will be suppressed.

  - `--outfmt/-f #`  
      
    Format of the output file. The following values are accepted:
    
      - `0`  
        BLAST pairwise format.
    
      - `5`  
        BLAST XML format.
    
      - `6`  
        BLAST tabular format (default). This format can be customized,
        the `6` may be followed by a space-separated list of the
        following keywords, each specifying a field of the output.
        
          - `qseqid`
            Query Seq - id
        
          - `qlen`
            Query sequence length
        
          - `sseqid`
            Subject Seq - id
        
          - `sallseqid`
            All subject Seq - id(s), separated by a ’;’
        
          - `slen`
            Subject sequence length
        
          - `qstart`
            Start of alignment in query*
        
          - `qend`
            End of alignment in query*
        
          - `sstart`
            Start of alignment in subject*
        
          - `send`
            End of alignment in subject*
        
          - `qseq`
            Aligned part of query sequence*
            
          - `full_qseq`
            Full query sequence
        
          - `sseq`
            Aligned part of subject sequence*
        
          - `full_sseq`
            Full subject sequence
        
          - `evalue`
            Expect value
        
          - `bitscore`
            Bit score
        
          - `score`
            Raw score
        
          - `length`
            Alignment length*
        
          - `pident`
            Percentage of identical matches*
        
          - `nident`
            Number of identical matches*
        
          - `mismatch`
            Number of mismatches*
        
          - `positive`
            Number of positive - scoring matches*
        
          - `gapopen`
            Number of gap openings*
        
          - `gaps`
            Total number of gaps*
        
          - `ppos`
            Percentage of positive - scoring matches*
        
          - `qframe`
            Query frame
        
          - `btop`
            Blast traceback operations(BTOP)*
            
          - `cigar`
            CIGAR string*
        
          - `staxids`
            Unique Subject Taxonomy ID(s), separated by a ’;’ (in
            numerical order). This field requires setting the
            `--taxonmap` parameter for `makedb`.
            
          - `sscinames`
            Unique Subject Scientific Name(s), separated by a ';'. This
            field requires setting the `--taxonmap` and `--taxonnames`
            parameters for `makedb`.
          
          - `sskingdoms`
            Unique Subject Super Kingdom(s), separated by a ';'. This
            field requires setting the `--taxonmap`, `--taxonnodes`
            and `--taxonnames` parameters for `makedb`.
            
          - `skingdoms`
            Unique Subject Kingdom(s), separated by a ';'. This
            field requires setting the `--taxonmap`, `--taxonnodes`
            and `--taxonnames` parameters for `makedb`.
            
          - `sphylums`
            Unique Subject Phylums(s), separated by a ';'. This
            field requires setting the `--taxonmap`, `--taxonnodes`
            and `--taxonnames` parameters for `makedb`.
        
          - `stitle`
            Subject Title
        
          - `salltitles`
            All Subject Title(s), separated by a ’\<\>’
        
          - `qcovhsp`
            Query Coverage Per HSP*
            
          - `scovhsp`
            Subject Coverage Per HSP*
        
          - `qtitle`
            Query title
            
          - `qqual`
            Query quality values for the aligned part of the query*
            
          - `full_qqual`
            Query quality values
            
          - `qstrand`
            Query strand
        
        By default, there are 12 preconfigured fields: `qseqid sseqid
        pident length mismatch gapopen qstart qend sstart send evalue
        bitscore`.
        
        *These fields require alignment traceback. If none of these fields
        are selected, traceback computation will be disabled, which improves
        performance and reduces use of temporary disk space.
    
      - `100`  
        DIAMOND alignment archive (DAA). The DAA format is a proprietary
        binary format that can subsequently be used to generate other
        output formats using the `view` command. It is also supported by
        MEGAN and allows a quick import of results. *Note that this format
        cannot be written to a stream output like stdout.*
    
      - `101`  
        SAM format.
    
      - `102`  
        Taxonomic classification. This format will not print alignments
        but only a taxonomic classification for each query using the LCA
        algorithm. The output lines consist of 3 tab-delimited fields:
        
        1.  Query ID
        
        2.  NCBI taxonomy ID (0 if unclassified)
        
        3.  E-value of the best alignment with a known taxonomic ID
            found for the query (0 if unclassified)
        
        The score range for the LCA algorithm is set by the `--top`
        parameter. The default value is 10 which means that all
        alignments whose score is at most 10% lower than the best score
        are considered for the LCA computation.
        
        Using this format requires setting the `--taxonmap` and
        `--taxonnodes` parameters for `makedb`.
    
      - `103`  
        [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).
        The custom fields in the format are AS (bit score), ZR (raw
        score) and ZE (e-value).

  - `--salltitles`  
      
    Include full length subject titles into the DAA format. By default,
    DAA files contain only the shortened sequence id (up to the first
    blank character).

  - `--sallseqid`  
      
    Include all subject ids into the DAA file. By default only the first
    id of each subject is included. As the subject ids are much shorter
    than the full titles this option will save space compared to the
    `--salltitles` option.

  - `--compress (0,1)`  
      
    Enable compression of the output file. 0 (default) means no
    compression, 1 means gzip compression.

  - `--max-target-seqs/-k #`  
      
    The maximum number of target sequences per query to report
    alignments for (default=25). Setting this to `-k0` will report all
    targets for which alignments were found.
    
    *Note that this parameter does not only affect the reporting, but
    also the algorithm as it is taken into account for heuristics that
    eliminate hits prior to full gapped extension.*

  - `--top #`  
      
    Report alignments within the given percentage range of the top
    alignment score for a query (overrides `--max-target-seqs` option).
    For example, setting `--top 10` will report all alignments whose
    score is at most 10% lower than the best alignment score for a
    query.
    
    *Note that this parameter does not only affect the reporting, but
    also the algorithm as it is taken into account for heuristics that
    eliminate hits prior to full gapped extension.*
        
  - `--max-hsps #`  
      
    The maximum number of HSPs (High-Scoring Segment Pairs) per target
    sequence to report for each query. The default policy is to report
    only the highest-scoring HSP for each target, while disregarding
    alternative, lower-scoring HSPs that are contained in the same target.
    **This is not to be confused with the `--max-target-seqs` option.**
    
    This parameter can be increased to report an alternative HSP if 
    its query and subject ranges are not enveloped by a higher scoring
    HSP and if it meets the e-value threshold. Setting this option to
    `--max-hsps 0` will report all alternative HSPs.

  - `--range-culling`  

    Enable hit culling with respect to the query range. This feature is
    designed for long query DNA sequences that may span several genes.
    In these cases, reporting the overall top N hits can cause hits to
    a lower-scoring gene to be superseded by a higher-scoring gene.
    Using this option, hit culling will be performed locally with respect
    to a hit's query range, thus reporting the locally top N hits while
    allowing more hits that span a different region of the query.

    Using this feature along with `-k 25` (default), a hit will only be
    deleted if at least 50% of its query range is spanned by at least 25
    higher or equal scoring hits.
    
    Using this feature along with `--top 10`, a hit will only be deleted
    if its score is more than 10% lower than that of a higher scoring
    hit over at least 50% of its query range.
    
    The overlap percentage is configurable using `--range-cover`. Note that
    this feature is currently only available in frameshift alignment mode.

  - `--evalue/-e #`  
      
    Maximum expected value to report an alignment (default=0.001).

  - `--min-score #`  
      
    Minimum bit score to report an alignment. Setting this option will
    override the `--evalue` parameter.

  - `--id #`  
      
    Report only alignments above the given percentage of sequence
    identity.
    
    *Note that using this option reduces performance.*

  - `--query-cover #`  
      
    Report only alignments above the given percentage of query cover.
    
    *Note that using this option reduces performance.*

  - `--subject-cover #`  
      
    Report only alignments above the given percentage of subject cover.
    
    *Note that using this option reduces performance.*

  - `--unal (0,1)`  
      
    Report unaligned queries (0=no, 1=yes). By default, unaligned
    queries are reported for the BLAST pairwise, BLAST XML and SAM
    format.

  - `--no-self-hits`  
      
    Suppress reporting of identical self-hits between sequences.
    The FASTA sequence identifier as well as the sequence of
    query and target need to be identical for a hit to be deleted.

## Memory & performance options

  - `--block-size/-b #`  
      
    Block size in billions of sequence letters to be processed at a
    time. This is the main parameter for controlling the program’s
    memory and disk space usage. Bigger numbers will increase the use
    of memory and temporary disk space, but also improve performance.
    The program can be expected to use roughly six times this number
    of memory (in GB).
    
    The default value is `-b2.0`. The parameter can be decreased for
    reducing memory use, as well as increased for better performance
    (values of >20 are not recommended).
    
    The very-sensitive and ultra-sensitive modes use `-b0.4` as a
    default. Note that these two modes benefit only slightly
    from increasing this parameter.
    
    Note that this parameter affects the algorithm and results will not
    be completely identical for different values of the block size.

  - `--tmpdir/-t <directory>`  
      
    Directory to be used for temporary storage. This is set to the
    output directory by default. The amount of disk space that will be
    used depends on the program’s settings and the input data. As a general
    rule it should be ensured that 100 GB of disk space are available here.
    
    If the program is being run in a cluster environment, and disk space is
    mounted over a network based file system, it is recommended to set this
    parameter to a fast local disk or to `/dev/shm` to avoid any I/O
    bottlenecks.

  - `--index-chunks/-c #`  
      
    The number of chunks for processing the seed index. This
    option can be additionally used to tune the performance. The default
    value is `-c4`, while setting this parameter to `-c1` instead will
    improve the performance at the cost of increased memory use. Note
    that the very-sensitive and ultra-sensitive modes use `-c1` by
    default.

## View options

  - `--daa/-a <file>`  
      
    Path to input file in DAA format.

  - `--out/-o <file>`  
      
    Path to output file. If this parameter is omitted, the results will
    be written to the standard output and all other program output will
    be suppressed.

These aligner parameters apply to the view command as well and work in
the same way: `--outfmt`, `--compress`, `--max-target-seqs`, `--top`. Note
that taxonomy features are currently not available for the DAA format.

## Advanced options

  - `--freq-sd #`  
      
    During the seed search, seeds with very high frequency in the
    queries or the database are ignored. This option sets the number of
    standard deviations above the average frequency for ignoring a seed.
    The default values are 50 in default mode, 20 in sensitive mode,
    200 in more-sensitive mode, 15 in very-sensitive mode and 20 in
    ultra-sensitive mode.
    
  - `--ext (banded-fast,banded-slow)`
  
    This option determines how band sizes are setup for banded Smith-Waterman
    extension. The `banded-slow` setting is slightly slower but more accurate.
    The defaults are `banded-fast` for the default and sensitive modes, and
    `banded-slow` for the more-sensitive, very-sensitive and ultra-sensitive
    modes.

  - `--band #`  
      
    Set a fixed band size for banded Smith Waterman extension. Setting this
    option overrides the preconfigured defaults set by the `--ext` option.
    Note that the band size is first determined by chaining, while this
    parameter provides an additional band margin in both directions.
    
  - `--no-ranking`
  
    Disable ranking heuristic. Ranking refers to a heuristic that eliminates
    hits prior to full gapped extension. It works by establishing a tentative order
    on the target sequences that were hit with respect to a single query. This
    order is determined by ungapped extension scores at seed hits. The aligner
    will compute gapped Smith-Waterman extensions for chunks of targets (see
    `--ext-chunk-size`) according to the ranking order until no alignment
    was found for the current chunk that meets the user-specified reporting
    criteria.
    
    Setting this option disables the ranking heuristic, causing gapped
    extensions to be computed for all seed hits.
    
    *Option supported since v2.0.3.*
    
  - `--ext-chunk-size #`
  
    The chunk size of target sequences used by the ranking heuristic (see
    above). Higher numbers improve the accuracy of the algorithm. The default
    is the number specified by the `--max-target-seqs/-k` option, rounded up
    to the next multiple of 32, with a minimum of 128 and a maximum of 400.
    When using `--top`, the default value is 128.
  
    *Option supported since v2.0.3.*

  - `--xml-blord-format`  
      
    Use gnlBL\_ORD\_ID style format for hit IDs in XML output.

  - `--range-cover #`  
      
    The minimum percentage of a hit's query range that needs to be spanned by
    higher scoring hits for a hit to be deleted in range culling mode
    (default=`50.0`).
    
  - `--culling-overlap #`
  
    The minimum percentage of a hit's query or target range that needs to be
    spanned by a higher scoring hit against the same target for a hit
    to be deleted. (default=`50.0`)
    
  - `--bin #`
  
    The number of bins for storing seed hits. Higher numbers will lead to the
    creation of more temporary files, reduce memory usage of the extension
    stage, but also slightly reduce efficiency. The default value is 16,
    except for the ultra-sensitive mode which uses `--bin 64` by default.
    
---

[1] Martin C. Frith; A new repeat-masking method enables specific detection
of homologous sequences, Nucleic Acids Research, Volume 39, Issue 4, 1 March
2011, Page e23, <https://doi.org/10.1093/nar/gkq1212>

[2] Maria Hauser, Martin Steinegger, Johannes Söding; MMseqs software suite
for fast and deep clustering and searching of large protein sequence sets,
Bioinformatics, Volume 32, Issue 9, 1 May 2016, Pages 1323–1330
