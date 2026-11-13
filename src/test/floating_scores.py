# SPDX-License-Identifier: GPL-3.0-or-later
"""End-to-end checks for the opt-in fractional scoring path."""
import math
import random
from pathlib import Path
import subprocess
import sys
import tempfile


diamond = str(Path(sys.argv[1]).resolve())
fasta = str(Path(sys.argv[2]).resolve())


def run(*extra, floating=True, success=True):
    args = [diamond, "blastp", "-q", fasta, "-d", fasta, "--swipe",
            "--comp-based-stats", "1", "--masking", "0", "-p", "1",
            "--min-score", "0.01", "-f", "6", "qseqid", "sseqid", "score", "bitscore"]
    if floating:
        args.append("--float-scores")
    if "--top" not in extra:
        args.extend(["-k", "0"])
    args.extend(extra)
    result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if success:
        assert result.returncode == 0, result.stderr
        return {(q, t): (float(raw), float(bits))
                for q, t, raw, bits in (line.split("\t") for line in result.stdout.splitlines())}
    assert result.returncode != 0, "Unsupported combination unexpectedly succeeded"
    assert "--float-scores" in result.stderr, result.stderr


floating = run()
assert run("--float-score-kernel", "scalar") == floating
integer = run(floating=False)
expected = {"S0C0_1": (249.752639093498, 251, 101),
            "S0C10_1": (581.339136152734, 584, 229)}
for name, (raw_expected, integer_raw, integer_bits) in expected.items():
    raw, bits = floating[name, name]
    # HauserCorrection stores floats; allow the small float/double difference.
    assert math.isclose(raw, raw_expected, abs_tol=2e-5), (name, raw)
    assert math.isclose(bits, (0.267 * raw - math.log(0.041)) / math.log(2), abs_tol=1e-10)
    assert integer[name, name] == (integer_raw, integer_bits), integer[name, name]
    assert (name, name) in run("--min-score", str(bits - 1e-6))
    assert (name, name) not in run("--min-score", str(bits + 1e-6))

for percent in (0, 10, 100):
    selected = run("--top", str(percent))
    for query in expected:
        best = max(bits for (q, _), (_, bits) in floating.items() if q == query)
        wanted = {pair for pair, (_, bits) in floating.items()
                  if pair[0] == query and bits >= max(1, best * (1 - percent / 100))}
        assert {pair for pair in selected if pair[0] == query} == wanted

for options in (("--comp-based-stats", "0"), ("--comp-based-stats", "2"),
                ("--max-hsps", "2"), ("-f", "100"), ("--iterate",),
                ("--float-score-kernel", "invalid"), ("-f", "6", "normalized_bitscore"),
                ("--block-size", "0.00000005")):
    run(*options, success=False)

missing_swipe = subprocess.run(
    [diamond, "blastp", "-q", fasta, "-d", fasta, "--float-scores"],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
assert missing_swipe.returncode != 0 and "requires blastp --swipe" in missing_swipe.stderr

# The single-residue correction is zero (there are no neighbours).
with tempfile.TemporaryDirectory() as directory:
    tiny = Path(directory) / "one.faa"
    tiny.write_text(">one\nW\n")
    result = run("-q", str(tiny), "-d", str(tiny))
    assert result["one", "one"][0] == 11, result

# Compare SIMD with the independent scalar reference across mixed target
# lengths, incomplete SIMD batches, indels, unknown residues and filtering.
with tempfile.TemporaryDirectory() as directory:
    rng = random.Random(1003)
    alphabet = "ARNDCQEGHILKMFPSTWYV"
    sequences = ["W", "AXWY", "AAAAAAAAAA", "WWWWWWWWWW"]
    sequences += ["".join(rng.choices(alphabet, k=n)) for n in (13, 27, 42, 69, 91)]
    queries = Path(directory) / "queries.faa"
    targets = Path(directory) / "targets.faa"
    queries.write_text("".join(">q{}\n{}\n".format(i, s) for i, s in enumerate(sequences)))
    target_sequences = sequences + [s[:len(s)//2] + "WQEP" + s[len(s)//2:] for s in sequences]
    target_sequences += ["".join(rng.choices(alphabet, k=n)) for n in (3, 8, 18, 31, 48)]
    targets.write_text("".join(">t{}\n{}\n".format(i, s) for i, s in enumerate(target_sequences)))
    for fields in (("score", "bitscore"),
                   ("score", "bitscore", "qstart", "qend", "sstart", "send", "cigar", "pident", "gapopen")):
        for options in (("--top", "100"), ("-k", "2"),
                        ("--id", "70", "--query-cover", "50", "-k", "1")):
            outputs = []
            for kernel in ("scalar", "simd"):
                args = [diamond, "blastp", "-q", str(queries), "-d", str(targets),
                        "--swipe", "--float-scores", "--float-score-kernel", kernel,
                        "--comp-based-stats", "1", "--masking", "0", "--min-score", "0.01",
                        "-p", "2", "-f", "6", "qseqid", "sseqid", *fields, *options]
                result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                assert result.returncode == 0, result.stderr
                outputs.append(sorted(result.stdout.splitlines()))
            assert outputs[0] == outputs[1], (fields, options, outputs)

print("Fractional CLI checks passed, including SIMD/scalar equivalence and the integer baseline.")
