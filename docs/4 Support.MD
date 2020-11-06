
# Support & FAQ

Support is provided at the [forums](index.php?forums/), by
[email](mailto:buchfink@gmail.com) or the
[GitHub issues page](http://github.com/bbuchfink/diamond/issues/).
For support requests it is usually helpful if you provide details about
your data (size and type of sequences) and your execution environment
(CPU cores, RAM, operating system).

Community feedback and feature requests are also appreciated.

1.  **How to run the program on multiple input files?**
    
    It is not recommended to run more than one instance of DIAMOND
    simultaneously on the same machine, because it is more efficient if
    you allocate more ressources to a single task by increasing the
    block size parameter `-b` and using lower numbers for `-c`. For this
    reason, input files should be processed consecutively. Pipes can
    also be used process multiple input files. For example, this command
    will build a database from all `fasta.gz` files in the current
    directory:
    
        zcat *.fasta.gz | diamond makedb -d diamond_db

2.  **The program runs for a very long time and the output file remains
    empty.**
    
    This is normal. Larger databases are split into chunks and the
    output is kept in temporary files at first and only merged in the
    end.

3.  **There are no temporary files.**
    
    The temporary files are not visible using standard shell commands
    because `unlink` is invoked on them after creation, which ensures
    that the files will be removed in case of ungraceful termination of
    the program.

4.  **More than one alignment is reported for a single query and subject
    sequence.**
    
    Due to the definition of HSPs, which are pairs of segments that form
    a linear alignment, any number of these may be contained in a single
    query/subject pair. Prior to version 2.0.0, the program reported
    additional aligned domains by default as they may contain valuable
    information. *Note that the default behaviour has changed since version
    2.0.0*

5.  **The program does not appear to be sensitive enough / it does not
    find an alignment I expect it to find or that another aligner
    reported.**
    
    There are a number of reasons that can cause hits to not be found or
    reported. In particular, even identical self-alignments of protein
    sequences may not be reported for some of the reasons listed below.
    
    1.  The hit may not be deemed significant, i.e. its e-value is above
        the configured cutoff. This is possible even for short 100%
        identity matches. Also note that e-value computations are not
        compatible for different aligners most of the time.
    
    2.  Repeat masking and compositional bias correction may
        cause hits to be filtered. This is usually desirable as these
        methods are highly effective in eliminating false positives.
        Repeat masking can be disabled using `--masking 0` if required.
    
    3.  There are limitations in aligning very short sequences of length
        <20, and this task requires some tuning of parameters.
    
    4.  Seeds that occur too frequently in the input are ignored, which
        may cause problems if your data is extremely repetitive. This
        behaviour can be adjusted using the `--freq-sd` option.
    
    5.  The program may not be sensitive enough for your application.
        You should try the more sensitive modes (options `--sensitive`
        up to `--ultra-sensitive`). Still, perfect results cannot be
        guaranteed by heuristical algorithms, but only by a rigorous
        Smith-Waterman computation.

6.  **The program terminates with an unclear or no error message.**
    
    The task may have been killed by the system for running out of
    memory. Use `dmesg` to check for a respective kernel message. See
    [instructions on memory usage](index.php?pages/command_line_options/#memory-performance-options).
    If you run the program on a cluster, you may need to allocate more memory to
    the task.

7.  **Are results consistent for multiple calls of the program?**
    
    The program has a basic consistency guarantee that the same version
    run on the same input with the same parameters will produce
    identical results, regardless of the hardware, operating system,
    compiler, or number of CPU threads. If this was not the case, I
    would consider it a bug. This is also ensured using the built-in
    regression tests (`diamond test`).
    
    Beyond that, no guarantees of consistency can be made. DIAMOND like
    all fast aligners uses numerous heuristics that can affect results in
    unexpected ways.
