# Experimental floating-point composition correction

`blastp --swipe --float-scores --comp-based-stats 1` uses floating-point
Smith–Waterman alignment with affine gap penalties. It consumes
the existing floating-point Hauser corrections without converting each one
to an integer. The dynamic-programming sums use `double`; the corrections
retain their existing `float` representation. Traceback follows stored
decisions, without floating-point equality comparisons.

The first pass scores targets together in SIMD lanes (four doubles with
AVX2, two with SSE2 or AArch64 NEON, with a scalar fallback elsewhere).
Targets of similar length are grouped to reduce padding. Substitution
profiles are cached per query residue within a batch. This pass needs only
rolling score arrays, with no per-cell traceback storage.

Score filtering and, when possible, top/k culling happen before traceback.
Coordinate/identity filters defer culling until the required tracebacks have
been checked. If the requested fields only need scores, traceback is skipped
entirely. `--float-score-kernel scalar` selects the independent scalar
reference for comparison, including calculation of its own optimal path.

This is an opt-in reference mode for investigating the effect of integer
correction rounding. The default SIMD search and its output are unchanged.
It is not a claim that the integer implementation is incorrect.

## Example

```
diamond blastp --query src/test/floating-scores.faa \
  --db src/test/floating-scores.faa --swipe --float-scores \
  --comp-based-stats 1 --masking 0 --max-target-seqs 0 \
  --outfmt 6 qseqid sseqid score bitscore --out fractional.tsv
```

Repeat without `--float-scores` for the integer reference. Floating raw scores
are used in bit-score conversion, significance filtering, ranking, and `--top`
culling. In this mode, the `score` and `bitscore` tabular fields use 17
significant digits rather than the usual compact output formatting.

## Limitations

- This initial reference implementation supports protein queries, CBS mode 1,
  one HSP per target, tabular output, and a database that fits in one block.
  Normalization by integer self-alignment scores is unsupported.
  Unsupported combinations fail with an explanation. Increase `--block-size`
  if necessary. DAA and intermediate binary records retain integer scores
  and are intentionally not supported here.
- This is exhaustive alignment, not the normal seeded search. Both integer
  and floating comparisons should use `--swipe`. It can be substantially
  slower than the normal integer SIMD kernels. The score pass uses linear
  memory in the longest target of each batch. When a traceback is required,
  it uses one byte per DP cell plus rolling double score arrays, per concurrent
  alignment. Use small inputs and few threads first; this mode is not yet
  intended for a full production database.
- Existing matrix statistical parameters are reused for bit scores and
  E-values. They have **not** been recalibrated for this fractional scoring
  scheme; E-values should be treated as experimental estimates.
- The model used to estimate `approx_pident` was trained on integer scores.
  Since this mode computes traceback, it supplies the exact identity
  percentage for that field and for `--approx-id` filtering instead.
- A one-residue query has no neighbouring residues, so its correction is zero.

## Validation

`ctest -R floating` checks fractional score preservation, affine gaps in both
directions, local restarts, empty input, and 200 small random matrices against
exhaustive enumeration of every positive local alignment path. End-to-end
tests cover the two protein examples, unrounded output, filtering, and
unsupported options. They also compare the scalar and SIMD results for mixed
lengths, partially filled batches, indels, unknown residues, score-only and
traceback output, and coordinate/identity filters, using two threads.

For an illustrative local timing comparison, run:

```
python3 src/test/benchmark_floating_scores.py build/diamond
```

It compares three runs of integer SIMD, fractional scalar, and fractional SIMD
on 20 queries against 200 targets, each 300 residues, on one thread with
masking disabled and exhaustive search. It also requires byte-identical
fractional scalar/SIMD output. The generated data are synthetic and the
timings include startup and input/output overhead.

An illustrative Linux x86_64/AVX2 Release run gave these medians:

| Mode | Time |
| --- | ---: |
| Integer SIMD | 0.130 s |
| Fractional scalar | 5.838 s |
| Fractional SIMD | 0.791 s |

The fractional SIMD mode was about 7.4 times faster than the scalar reference
and about 6.1 times slower than integer SIMD on this small comparison.
