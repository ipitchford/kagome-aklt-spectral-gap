# AI index — kagome AKLT spectral-gap candidate

Status: unrefereed computer-assisted proof candidate, v0.3.1. Creator: Anonymous. Original text/data CC0-1.0, original code MIT; see [licences](LICENSES.md) and [provenance](PROVENANCE.md).

## Exact claim and scope

For spin two on the undecorated kagome lattice, with each edge term the orthogonal projector onto total spin four, the candidate proves

`gap(H) >= 11*(2-sqrt(3))*(59-3*sqrt(313))/3360 > 51/10000`.

Volumes: every rectangular torus with both periods at least three and every nonempty finite union of complete down-star edge sets. NOT arbitrary induced patches; NOT an unqualified infinite-volume theorem; NOT a square-lattice gap, LTQO or perturbation-stability theorem. The old 0.0407 bound is not certified. The current argument does not depend on a random-start filter or an implementation-level failure probability.

## Read in this order

1. [Manuscript](paper/paper.pdf), [standalone LaTeX](paper/paper.tex): definitions, exact theorem, grouped proof, worked 83-dimensional sector, all-sector table and arithmetic appendix.
2. [CLAIMS.json](CLAIMS.json): exact claim IDs and dependency graph. `uniform_kagome_gap` depends on `three_block_gap`, `weighted_block_gap` and exact ground-space construction.
3. [ASSURANCE.md](ASSURANCE.md): semantic and arithmetic trust boundary, sector-dependent grouping, rank-lift argument and nonclaims.
4. [SOURCES.md](SOURCES.md): dated bounded antecedent comparison; not global novelty clearance.
5. [REVIEW_RESPONSE_v0.3.1.md](REVIEW_RESPONSE_v0.3.1.md): disposition of the supplied review. The older response file is historical.

## Reproduce, do not just read receipts

Verify [MANIFEST.sha256](MANIFEST.sha256) with `python tools/verify_manifest.py` before running anything. Copy the extraction to a working directory because replay replaces generated receipts. Install [requirements.txt](requirements.txt), a C++17 compiler and GMP development libraries; reserve at least 4 GiB free RAM. Run `python code/replay.py`. Only `PASS_FULL` is full proof replay. `--quick` omits all 50 hub-angle checks and is explicitly not sufficient.

[code/replay.py](code/replay.py) reconstructs inputs rather than trusting stored inequalities. [rational_hub.py](code/certified/rational_hub.py) and [gram_exact.py](code/certified/gram_exact.py) construct exact tensors and rational frames. [certify_hub.py](code/certified/certify_hub.py) checks 15 pair and 35 triple sectors. [validated.py](code/certified/validated.py) propagates exact-input and floating errors and checks every factor-residual entry. [block_gap.py](code/certified/block_gap.py) supplies four local block tests. [validator_tests.py](code/certified/validator_tests.py) contains semantic negative controls and the small all-rational comparison.

## Evidence and safe downstream use

[REPLAY_RECEIPT.md](REPLAY_RECEIPT.md) records the publication replay and backend. [evidence/new/replay_receipt.json](evidence/new/replay_receipt.json) binds actual scientific source hashes to all driver stages. [validated_all.json](evidence/new/validated_all.json) contains all 50 hub tests; [block_gap_validated.json](evidence/new/block_gap_validated.json) contains four block tests. Logs are under `logs/replay/`.

Trust ordinary round-to-nearest binary64 and conventional bounded-length BLAS dot products under the stated error model, exact integer/rational constructors, and the written mapping between tensors and the physical model. Numerical lift columns and Cholesky factors are candidates, not trusted eigenvectors or exact factors. The exact common-range rank and positive-definiteness inequality are load-bearing.

Producer replay, implementation diversity on small examples, external reproduction, formal verification, human review and novelty are separate. Only the first two are recorded here; no unaffiliated full reproduction is established. Cite the exact release DOI/commit from [CITATION.cff](CITATION.cff). Do not silently broaden volume scope, improve the constant, turn a local matrix validation margin into a physical gap or treat an archived candidate as a peer-reviewed theorem.
