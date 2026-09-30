# Publication replay — v0.3.1

Date: 30 September 2026. Result: **PASS_FULL**, all 20 driver stages; 15 pair, 35 triple and four block positive-definiteness checks. Scientific Python/C++ source hashes in evidence/new/replay_receipt.json match the final package. The v0.3.1 changes are exposition, references and packaging only.

Fresh reconstruction on macOS 27 ARM64, Python 3.14.7, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0, Apple Accelerate BLAS/LAPACK. Driver elapsed: 385.3066 seconds. Backend version is not exposed by NumPy. Compiler: Apple clang through g++; GMP headers/libraries supplied by Homebrew. Full environment output: logs/replay/environment.log.

The outer BSD time utility returned nonzero after the successful child run because the sandbox denied its sysctl resource query. This did not negate the driver's PASS_FULL or any of the 20 successful subprocess exit codes. Peak memory was not measured on this run; the supplied Linux audit reports about 2.47 GiB. That Linux record is reported evidence, not a new measurement here.

The same driver reconstructs inputs, rejects stale/missing stage receipts, exercises adversarial matrices and an all-rational 83-dimensional sector, and repeats selected checks under optimized Python. Quick mode is not full proof verification. Exact output bytes may vary by backend or timing; inequalities must pass each run.

A supplied separately written small-system checker, importing none of the candidate code, was rerun privately and passed. It independently assembles the exact two-spin projector and triangle spectrum, checks constants and geometry, and numerically corroborates star/bowtie formulas. These floating checks are not exact spectral-ordering certificates; this is implementation diversity, not an unaffiliated full reimplementation of all 54 tests. Private audit files are not redistributed without an explicit licence.

Trust boundaries: written tensor-to-model and sector/gluing arguments; exact arithmetic implementation; the documented conventional binary64/BLAS error model; code, compiler and hardware. No formal proof, external specialist review or worldwide priority determination is established. CI results are available separately on GitHub and are producer-coordinated replays, not peer review.
