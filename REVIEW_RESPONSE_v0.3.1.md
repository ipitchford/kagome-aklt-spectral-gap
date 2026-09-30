# Review response — 30 September 2026

Decision: publish the corrected v0.3.1 unrefereed candidate. This is a response to supplied model-assisted review, not an external specialist endorsement or an official REF assessment. All required exposition, positioning and archival requests are addressed; optional new theorems and optimisation are not asserted.

| Review point | Action and evidence |
| --- | --- |
| Explain the sector decomposition through one concrete case | New worked `(mm,mm)|mm` subsection derives external character, private character, 83-dimensional zero-weight space and rank-seven lift. The table remains complete for all 15 pair and 35 triple types. |
| Explain partners, permutations and sector-dependent grouping | The worked subsection states the exact symmetry intertwiners and why the chosen pair can differ across invariant sectors. No fixed global grouping is claimed. |
| Separate exact input, numerical candidates, error propagation and acceptance | Expanded verification-dependency appendix separates all four interfaces. Candidate lift columns and Cholesky factors are not trusted as exact eigenvectors or exact factorizations. |
| Strengthen implementation diversity and state its limits | Full unchanged scientific-code replay passed on macOS ARM64/Accelerate in 385.31 seconds, versus the supplied Linux/OpenBLAS record. This is producer-coordinated cross-platform replay, not unaffiliated reproduction. The supplied separately written small-system checker is rerun privately; it does not verify all 54 large tests. |
| Position weighted finite-size and hierarchy antecedents | Added Lemm–Sandvik–Wang (2020) and Rai et al. (2026), distinguishing inherited finite-size/local-positivity ideas from the model-specific grouped certificate and explaining the failed pair-only threshold. No benchmark superiority is claimed. |
| Explain infinite-volume and stability scope | Retained exact finite-volume theorem: rectangular tori with both periods at least three and complete down-star edge unions. Added explicit nonclaims for arbitrary induced patches, infinite-volume representation-dependent statements, LTQO and perturbation stability. Related stability literature is context, not an imported conclusion. |
| Remove notation collision | Triangle projector is now indexed by `T_0`; nine-site down-star projector retains `D`. |
| Present the paper as mathematics, not a revision memo | Early v0.2 change log moved to version-history appendix; scientific statement, proof and references remain primary. |
| Square obstruction threshold | Retained as a separate restricted obstruction and moved the minor corrected-threshold history to the version appendix. No square-lattice gap claim. |
| Archive all load-bearing components and retain evidence | Public package supplies manuscript, exact constructors, checker, full receipts/logs, pinned requirements, manifest, claim ledger and substantive AI index. Private supplied reviews and historical input archives are retained outside the public payload with hashes, rather than redistributed without a clear licence. |

## Checks and unresolved boundaries

The main executable mathematics was not changed for this revision. Full replay reconstructs exact tensors, checks 15+35+4 positive-definiteness targets, exact local ordering, geometry, adversarial controls and selected optimized-Python runs. Package gates bind those source hashes to the final archive. The mathematical semantic bridge and stated binary64/BLAS error model remain trusted; no proof-assistant verification, worldwide novelty proof, external human review, optimal gap constant or quantitative stability result is supplied.

The old `REVIEW_RESPONSE.md` concerns v0.2-to-v0.3 development and is retained as historical context, not the response for this publication revision.
