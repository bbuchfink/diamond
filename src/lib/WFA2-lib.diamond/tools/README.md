# WFA TOOLS

## WHAT IS THIS?

The WFA2-lib offers tools that exploits the different libraries features and serve as examples on how to use and exploit the library functions. Also, these tools serve for testing the library and benchmarking its performance against other state-of-the-art tools and libraries for pairwise alignment. If you are interested in benchmarking WFA2-lib and other algorithms (internal or external), checkout the branch `benchmark`.

* [Generate Dataset](#tool.generate)
* [Align Benchmark](#tool.align)

## <a name="tool.generate"></a> 1. GENERATE DATASET TOOL

The *generate-dataset* tool allows generating synthetic random datasets for testing and benchmarking purposes. This tool produces a simple output format (i.e., each pair of sequences in 2 lines) containing the pairs of sequences to be aligned using the *align-benchmark* tool. For example, the following command generates a dataset named 'sample.dataset.seq' of 5M pairs of 100 bases with an alignment error of 5% (i.e., 5 mismatches, insertions, or deletions per alignment).

```
$> ./bin/generate_dataset -n 5000000 -l 100 -e 0.05 -o sample.dataset.seq
```

### Command-line Options 

```
        --output|o <File>
          Filename/Path to the output dataset.
          
        --num-patterns|n <Integer>
          Total number of pairs pattern-text to generate.
          
        --length|l <Integer>
          Length of the generated pattern (ie., query or sequence) and text (i.e., target or reference)
          
        --pattern-length|P <Integer>
          Length of the generated pattern (ie., query or sequence)
        
        --text-length|T <Integer>
          Length of the generated text (i.e., target or reference)
          
        --error|e <Float>
          Total error-rate between the pattern and the text (allowing single-base mismatches, 
          insertions and deletions). This parameter may modify the final length of the text.
          
        --debug|g
          Output debug information.
          
        --help|h
          Outputs a succinct manual for the tool.
```

## <a name="tool.align"></a> 2. ALIGNMENT BENCHMARK TOOL

### Introduction to Alignment Benchmarking

The WFA2-lib includes the benchmarking tool *align-benchmark* to test and compare the performance of various pairwise alignment implementations. This tool takes as input a dataset containing pairs of sequences (i.e., pattern and text) to align. Patterns are preceded by the '>' symbol and texts by the '<' symbol. Example:

```
>ATTGGAAAATAGGATTGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTCGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTAGCTCGAAGCCCA
<GATTGGAAAATAGGATGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTGCTCGAAGCCCA
>CCGTAGAGTTAGACACTCGACCGTGGTGAATCCGCGACCACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCAGTGATTAAAC
<CCTAGAGTTAGACACTCGACCGTGGTGAATCCGCGATCTACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCGAGTGATTAAAC
[...]
```

Once you have the dataset ready, you can run the *align-benchmark* tool to benchmark the performance of a specific pairwise alignment method. For example, the WFA algorithm:

```
$> ./bin/align_benchmark -i sample.dataset.seq -a gap-affine-wfa
...processed 10000 reads (benchmark=125804.398 reads/s;alignment=188049.469 reads/s)
...processed 20000 reads (benchmark=117722.406 reads/s;alignment=180925.031 reads/s)
[...]
...processed 5000000 reads (benchmark=113844.039 reads/s;alignment=177325.281 reads/s)
[Benchmark]
=> Total.reads            5000000
=> Time.Benchmark        43.92 s  (    1   call,  43.92  s/call {min43.92s,Max43.92s})
  => Time.Alignment      28.20 s  ( 64.20 %) (    5 Mcalls,   5.64 us/call {min438ns,Max47.05ms})
```

The *align-benchmark* tool will finish and report overall benchmark time (including reading the input, setup, checking, etc.) and the time taken by the algorithm (i.e., *Time.Alignment*). If you want to measure the accuracy of the alignment method, you can add the option `--check` and all the alignments will be verified. 

```
$> ./bin/align_benchmark -i sample.dataset.seq -a gap-affine-wfa --check
...processed 10000 reads (benchmark=14596.232 reads/s;alignment=201373.984 reads/s)
...processed 20000 reads (benchmark=13807.268 reads/s;alignment=194224.922 reads/s)
[...]
...processed 5000000 reads (benchmark=10625.568 reads/s;alignment=131371.703 reads/s)
[Benchmark]
=> Total.reads            5000000
=> Time.Benchmark         7.84 m  (    1   call, 470.56  s/call {min470.56s,Max470.56s})
  => Time.Alignment      28.06 s  (  5.9 %) (    5 Mcalls,   5.61 us/call {min424ns,Max73.61ms})
[Accuracy]
 => Alignments.Correct        5.00 Malg        (100.00 %) (samples=5M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
 => Score.Correct             5.00 Malg        (100.00 %) (samples=5M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => Score.Total           147.01 Mscore uds.            (samples=5M{mean29.40,min0.00,Max40.00,Var37.00,StdDev6.00)}
     => Score.Diff            0.00 score uds.  (  0.00 %) (samples=0,--n/a--)}
 => CIGAR.Correct             0.00 alg         (  0.00 %) (samples=0,--n/a--)}
   => CIGAR.Matches         484.76 Mbases      ( 96.95 %) (samples=484M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Mismatches        7.77 Mbases      (  1.55 %) (samples=7M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Insertions        7.47 Mbases      (  1.49 %) (samples=7M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}
   => CIGAR.Deletions         7.47 Mbases      (  1.49 %) (samples=7M{mean1.00,min1.00,Max1.00,Var0.00,StdDev0.00)}

```

Using the `--check` option, the tool will report *Alignments.Correct* (i.e., total alignments that are correct, not necessarily optimal), and *Score.Correct* (i.e., total alignments that have the optimal score). Note that the overall benchmark time will increase due to the overhead introduced by the checking routine, however the *Time.Alignment* should remain the same.

### Algorithms & Implementations Summary

Summary of algorithms/methods implemented within the benchmarking tool. If you are interested 
in benchmarking WFA with other algorithms implemented or integrated into the WFA library, 
checkout branch `benchmark`.

|          Algorithm Name           |       Code-name       | Distance        |  Output   |    Library     |
|-----------------------------------|:---------------------:|:---------------:|:---------:|:--------------:|
| WFA Indel (LCS)                   | indel-wfa             | Indel           | Alignment | WFA2-lib       |
| Bit-Parallel-Myers (BPM)          | edit-bpm              | Edit            | Alignment | WFA2-lib       |
| DP Edit                           | edit-dp               | Edit            | Alignment | WFA2-lib       |
| DP Edit (Banded)                  | edit-dp-banded        | Edit            | Alignment | WFA2-lib       |
| WFA Edit                          | edit-wfa              | Edit            | Alignment | WFA2-lib       |
| DP Gap-linear (NW)                | gap-linear-nw         | Gap-linear      | Alignment | WFA2-lib       |
| WFA Gap-linear                    | gap-linear-wfa        | Gap-linear      | Alignment | WFA2-lib       |
| DP Gap-affine (SWG)               | gap-affine-swg        | Gap-affine      | Alignment | WFA2-lib       |
| DP Gap-affine Banded (SWG Banded) | gap-affine-swg-banded | Gap-affine      | Alignment | WFA2-lib       |
| WFA Gap-affine                    | gap-affine-wfa        | Gap-affine      | Alignment | WFA2-lib       |
| DP Gap-affine Dual-Cost           | gap-affine2p-dp       | Gap-affine (2P) | Alignment | WFA2-lib       |
| WFA Gap-affine Dual-Cost          | gap-affine2p-wfa      | Gap-affine (2P) | Alignment | WFA2-lib       |

* DP: Dynamic Programming
* SWG: Smith-Waterman-Gotoh
* NW: Needleman-Wunsch

### Command-line Options

#### - Input

```
          --algorithm|a <algorithm-code-name> 
            Selects pair-wise alignment algorithm/implementation.
                                                       
          --input|i <File>
            Filename/path to the input SEQ file. That is, file containing the sequence pairs to
            align. Sequences are stored one per line, grouped by pairs where the pattern is 
            preceded by '>' and text by '<'.
            
          --output|o <File>
            Filename/path of the output file containing a brief report of the alignment. Each line
            corresponds to the alignment of one input pair with the following tab-separated fields:
            <SCORE>  <CIGAR>
          
          --output-full <File> 
            Filename/path of the output file containing a report of the alignment. Each line
            corresponds to the alignment of one input pair with the following tab-separated fields:
            <PATTERN-LENGTH>  <TEXT-LENGTH>  <SCORE>  <PATTERN>  <TEXT>  <CIGAR>
```
                                     
#### - Penalties & Span

```                                                  
          --lineal-penalties|p M,X,I
            Selects gap-lineal penalties for those alignment algorithms that use this penalty model.
            Example: --lineal-penalties="0,1,2"
              M - Match penalty
              X - Mismatch penalty
              I - Indel penalty
                
          --affine-penalties|g M,X,O,E
            Selects gap-affine penalties for those alignment algorithms that use this penalty model.
            Example: --affine-penalties="0,4,2,6" 
              M - Match penalty
              X - Mismatch penalty
              O - Open penalty
              E - Extend penalty
          
          --affine2p-penalties M,X,O1,E1,O2,E2 
            Select gap-affine dual-cost penalties for those alignment algorithms that use this  
            penalty model. Example: --affine2p-penalties="0,4,2,6,20,2" 
              M  - Match penalty
              X  - Mismatch penalty
              O1 - Open penalty (gap 1)
              E1 - Extend penalty (gap 1)
              O2 - Open penalty (gap 2)
              E2 - Extend penalty (gap 2)
          
          --ends-free P0,Pf,T0,Tf
            Determines the maximum ends length to allow for free in the ends-free alignment mode.
            Example: --ends-free="100,100,0,0"
              P0 - Pattern begin (for free)
              Pf - Pattern end (for free)
              T0 - Text begin (for free)
              Tf - Text end (for free)
          
```
                         
#### - Wavefront parameters

```                                                                  
          --minimum-wavefront-length <INT>
            Selects the minimum wavefront length to trigger the WFA-Adapt reduction method.
            
          --maximum-difference-distance <INT>
            Selects the maximum difference distance for the WFA-Adapt reduction method.  
```
                   
#### - Others

```           
          --bandwidth <INT>
            Selects the bandwidth size for those algorithms that use bandwidth strategy. 
            
          --num-threads|t <INT>
            Sets the number of threads to use to align the sequences. If the multithreaded mode is
            enabled, the input is read in batches of multiple sequences and aligned using the number
            of threads configured. 
            
          --batch-size <INT>
            Selects the number of pairs of sequences to read per batch in the multithreaded mode.
            
          --check|c 'correct'|'score'|'alignment'                    
            Activates the verification of the alignment results. 
          
          --check-distance 'edit'|'gap-lineal'|'gap-affine'
            Select the alignment-model to use for verification of the results.
          
          --check-bandwidth <INT>
            Sets a bandwidth for the simple verification functions.

          --help|h
            Outputs a succinct manual for the tool.
```

## AUTHORS

  Santiago Marco-Sola \- santiagomsola@gmail.com     

## REPORTING BUGS

Feedback and bug reporting it's highly appreciated. Please report any issue or suggestion on github or by email to the main developer (santiagomsola@gmail.com).

## LICENSE

WFA Tools are distributed under MIT licence.
