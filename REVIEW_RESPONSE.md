# Response to the supplied referee report

This table describes the v0.3 changes. “Addressed” means addressed by the revised proof and executed internal checks; it does not mean that an independent referee has accepted the response.

| Review issue | Action and evidence | Status |
|---|---|---|
| Random-start probabilities do not cover the floating implementation | Removed the large random filter from the proof. Replaced it with exact rational projector definitions, rank-lift inertia and all-entry validated shifted-Cholesky residuals. No implementation-level failure probability is asserted. | Addressed by a new proof route |
| Wrong Schwinger-polynomial normalization | Corrected the map to `|k>/sqrt(binomial(4,k))`. The exact two-site regression accepts the corrected singlet state and rejects the printed v0.2 convention. | Addressed |
| Unverified exact-to-numerical projector interface | Reconstruct integer polynomial maps; exact Grams and rational frames; check exact local inverse-factor equations and global magnetization selection; bound each rational conversion and subsequent contraction. No numerical rank or approximate inverse defines the target projector. | Addressed, subject to the stated software contract |
| Fragile sampled deflation checks | Eliminated the need to validate the numerical deflation subspace. An arbitrary dyadic rank-k lift is sufficient once the lifted exact form is proved positive. The common-space count remains exact. | Addressed analytically |
| Star root ordering loses rigor at final float comparisons | Replaced the minimality argument with independent reconstruction of an integer 381-by-381 threshold matrix and exact inertia `(336,45,0)`. Together with 44 exact zero modes and the invariant quadratic witness, this proves the minimum without supplied factors. | Addressed; no claim about the next positive star eigenvalue retained |
| Bowtie root ordering loses certified radii | Reconstruct all highest-weight/leaf-parity blocks over Q, verify full dimension coverage and metric self-adjointness, factor exactly and use rational Sturm counts. Integer roots are excluded algebraically. | Addressed |
| Separate existence from the optional numerical constant | Establish μ₃ ≥ (73−√313)/80 > 2/3 first. A separate four-test argument establishes γ₀ ≥ 11(2−√3)/42. The old numerical constant is withdrawn from the main theorem. | Addressed; explicit new bound >0.0051 |
| Independent structural checks of the three-block operator | Exact dihedral/S₃ symmetry checks, explicit rational contraction formulas, exact common counts, complete sector list, an alternative normalized-tensor exploration and an all-rational 83-dimensional triple certificate. | Addressed internally; full external replication remains outstanding |
| “Common kernel” incorrectly used for bowtie common ground space | Replaced by common range: 19 at magnetization zero; common kernel is 190. Exact embedded ranks are 105 each. | Addressed |
| Square threshold comparison near N=40 | Record t₃₄(2)/2 = 33/2246 < 0.015. No finite-patch gap is inferred from this threshold calculation. | Addressed |
| Geometry and presentation | Added a labelled 21-site shared-hub graph, identified each nine-site block, tabulated every certification sector, separated exact/validated/exploratory claims, and added the theorem-to-code map. | Addressed |
| Square obstruction must retain restricted scope | The theorem concerns cross-supported translation densities after disjoint-edge products are removed. Larger supports, other decompositions and square-lattice gap existence are not ruled out. | Scope preserved |
| Priority and quality calibration | Bounded current discovery found no competing undecorated kagome proof, but no exhaustive worldwide priority claim is made. The quality assessment remains prospective and non-official. | Explicit limitation |

## Tractable extensions completed

The rectangular-torus statement now covers both side lengths independently, provided each is at least three. The same local certificates establish a uniform gap for any finite union of complete down-star edge sets, including disconnected unions. The proof explicitly handles hubs containing one, two or three selected blocks. This is not extended to arbitrary induced site patches.

The weighted-block calculation is new: it compares the product of three disjoint up-triangle ground spaces to the central down-triangle ground space. Four parity sectors suffice, the largest of order 856. It turns the existence theorem into an explicit quantitative result without restoring the old block filter.

## Developments during the internal audit

Before the frozen replay, the implementation was strengthened to include the cross term in the metric-radius product, global exact magnetization selection for local projector factors, finite/symmetry guards, underflow allowances, and fail-closed block/replay acceptance. The final source hashes in `replay_receipt.json` identify the version actually used for the frozen full replay. Preliminary development receipts are not authoritative.
