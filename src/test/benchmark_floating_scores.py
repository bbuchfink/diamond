# SPDX-License-Identifier: GPL-3.0-or-later
"""Small reproducible timing comparison; not a production search benchmark."""
from pathlib import Path
import random
import statistics
import subprocess
import sys
import tempfile
import time

diamond = str(Path(sys.argv[1]).resolve())
rng = random.Random(1003)
alphabet = "ARNDCQEGHILKMFPSTWYV"
targets = ["".join(rng.choices(alphabet, k=300)) for _ in range(200)]
with tempfile.TemporaryDirectory() as directory:
    for name, records in (("reference", targets), ("queries", targets[:20])):
        Path(directory, name + ".faa").write_text("".join(
            ">seq_{}\n{}\n".format(i, s) for i, s in enumerate(records)))
    outputs = {}
    for name, flags in (("integer SIMD", []),
                        ("fractional scalar", ["--float-scores", "--float-score-kernel", "scalar"]),
                        ("fractional SIMD", ["--float-scores"])):
        timings = []
        for _ in range(3):
            args = [diamond, "blastp", "-q", str(Path(directory, "queries.faa")),
                    "-d", str(Path(directory, "reference.faa")), "--swipe",
                    "--comp-based-stats", "1", "--masking", "0", "-p", "1",
                    "-f", "6", "qseqid", "sseqid", "score", "bitscore", "--quiet", *flags]
            start = time.perf_counter()
            result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            timings.append(time.perf_counter() - start)
            outputs[name] = result.stdout
        print("{}: median {:.3f}s; runs {}".format(name, statistics.median(timings), timings), flush=True)
    assert outputs["fractional scalar"] == outputs["fractional SIMD"], "Fractional results differ"
