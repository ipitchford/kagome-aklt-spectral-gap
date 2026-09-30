# Historical v0.3 development receipt

The record below describes the supplied Linux development run. It is not the current macOS publication replay; see REPLAY_RECEIPT.md.

# v0.3 release receipt — 30 September 2026

The frozen-source replay completed with **PASS_FULL**. All 18 proof-source files match the hashes recorded at the end of that execution.

| Item | Final execution |
|---|---:|
| Hub-angle tests | 50 / 50 pass |
| Weighted-block-angle tests | 4 / 4 pass |
| Largest hub matrix | 11,524 |
| Largest block matrix | 856 |
| Wall-clock time reported by the driver | 458.85 seconds |
| `/usr/bin/time` wall-clock time | 7 minutes 39.54 seconds |
| Reported peak resident set | 2,585,656 KiB (approximately 2.47 GiB) |
| Maximum hub error-budget / shift ratio | 0.062645977 |
| Original input-manifest files verified | 85 |

The exact star inertia, exact bowtie reconstruction and ordering, structural/geometry checks, and adversarial validator tests all pass. Selected checks are also replayed under optimized Python. The final driver rejects failed receipts or missing coverage; merely running Cholesky or loading stored successful scalars is not acceptance.

The theorem expression is **11(2−√3)(59−3√313)/3360 > 0.0051**, in the total-spin-four edge-projector normalization. The scope and exclusions are specified in the paper and `CLAIMS.json`. These are internal execution results; external specialist verification and an official research-assessment decision are not claimed.

Use `python tools/verify_manifest.py` to check release integrity before running the computational replay in a separate working copy. Manifest verification is not a mathematical verification of the theorem.
