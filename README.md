# Kagome AKLT: deterministic gap certificate — v0.3.1

**Research candidate for independent review.** The main paper contains a deterministic computer-assisted proof, under the stated arithmetic contract, of a uniform spectral gap for the undecorated spin-2 kagome AKLT Hamiltonian. It has not received external specialist verification or an official REF assessment.

## Main result

With each edge term the total-spin-four orthogonal projector, the proved lower expression is

**Δ = 11(2 − √3)(59 − 3√313) / 3360 = 0.0051971299023116… > 0.0051.**

The scope is every rectangular kagome torus with both sides at least three, and finite open patches whose edges are precisely unions of complete down-star edge sets. This is not a theorem about arbitrary induced open patches. It does not settle the square-lattice gap problem.

The former proposed constant **0.0407 is not certified here**. The former random-start filter is not a dependency of this proof. The new route uses 15 pair-sector, 35 triple-sector and four block-angle deterministic tests. The largest matrix has order 11,524.

## Start here

`paper/paper.pdf` is the complete revised manuscript; `paper/paper.tex` is its source. `RESEARCH_REPORT.md` explains the advance and its prospective quality. `REVIEW_RESPONSE_v0.3.1.md` answers the current supplied report; `REVIEW_RESPONSE.md` is historical. [AI_INDEX.md](AI_INDEX.md) provides a claim-level agent handoff. `ASSURANCE.md` specifies what the computations do and do not establish. `CLAIMS.json` is the machine-readable claim/status ledger.

Authoritative execution results are in `evidence/new/`, especially `replay_receipt.json`, `validated_all.json`, `block_gap_validated.json` and `validator_tests.json`. Logs of the frozen-source full replay are in `logs/replay/`; its resource and backend boundary is in `REPLAY_RECEIPT.md`. Historical supplied receipts are retained privately; preliminary numerical tests are segregated under exploratory/. They are not substitutes for the authoritative run.

## Full replay

The tested environment is Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0, GNU C++ and GMP, on x86-64 Linux with OpenBLAS 0.3.30. Install the three Python packages from `requirements.txt`; a C++ compiler and GMP development headers/libraries must also be available. The driver invokes `g++` with `-lgmpxx -lgmp`.

From the root of an unpacked working copy:

```sh
python -m pip install -r requirements.txt
python code/replay.py
```

The default is the **full** proof replay, not an audit of stored scalar inequalities. It rebuilds the integer/rational tensors and projector factors, runs all 54 positive-definiteness tests, reconstructs the exact star and bowtie results, checks symmetry/geometry/normalization, executes negative controls, and repeats selected checks under Python optimization. A rejected stage terminates the run. The final driver also checks the status and coverage of every authoritative receipt. Only `PASS_FULL` denotes a complete successful run.

On the development machine, a complete replay takes approximately eight minutes with two BLAS threads. Treat this as an observation, not a cross-platform runtime guarantee. Provide at least 4 GiB of free RAM; additional headroom is preferable. Actual final-run resource counters are supplied. No GPU, cloud account or python-flint dependency is needed. Publication replay also passed on macOS ARM64 with Accelerate in 385.31 seconds (see REPLAY_RECEIPT.md): the documented binary64 arithmetic contract must still hold, and the run must pass its own bounds. Low-precision or approximate matrix multiplication is not an admissible substitution.

```sh
python code/replay.py --quick
```

Quick mode intentionally omits the 50 main-sector tests and ends with `PASS_QUICK_NOT_FULL_PROOF`. It is **not** a complete verification of the kagome theorem.

The driver writes generated caches to `.work/`, overwrites execution receipts in `evidence/new/`, and writes new logs. Generated pickle caches are not distributed as proof inputs and are rebuilt before use. To preserve the release's byte-level integrity record, verify `MANIFEST.sha256` first and run the replay in a separate working copy. Numerical bounds may vary slightly by conforming backend; timings and regenerated receipts need not be byte-identical.

## Paper build

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
```

The manuscript uses standard LaTeX mathematics, TikZ and longtable. The checked PDF is included.

## Trust boundary

The exact part trusts Python/SymPy integer and rational arithmetic, the elementary GMP fraction-free elimination, and the correspondence between the explicitly constructed tensors and the stated model. The validated part additionally trusts conventional binary64 BLAS dot products and the outward scalar bounds. A Cholesky factor is only a candidate: all entries of its residual are checked. There is no exceptional random-start set and no probability-of-truth claim. The new pipeline is not a proof-assistant formalization and has not been externally audited.

## Provenance and licensing

Historical input archives and private review material are retained privately, not redistributed without clear licensing. SOURCE_INPUTS.json records their hashes. Original code is MIT; original prose/data are CC0-1.0. Dependencies and fonts retain their terms. No human mathematical contribution is inferred. See PROVENANCE.md and LICENSES.md.

## Immutable citation

Anonymous. Version 0.3.1-candidate. https://doi.org/10.5281/zenodo.23070704

https://github.com/ipitchford/kagome-aklt-spectral-gap/releases/tag/v0.3.1-candidate
