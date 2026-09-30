# Research development report — v0.3

## Conclusion

The tractable route was not to make the original 106,935,090-dimensional randomized test more elaborate. It was to replace the spectral question with a weaker sufficient inequality that could be certified completely.

The resulting candidate contains a deterministic computer-assisted proof of a uniform spin-2 kagome AKLT gap, with explicit lower expression

**11(2 − √3)(59 − 3√313) / 3360 = 0.00519712990231160577… > 0.0051.**

Every required finite inequality has been executed in the supplied implementation, starting from freshly reconstructed exact rational tensors. The output remains a research candidate requiring independent specialist assessment. The new proof/checker has not been externally reviewed, and this report is not an official REF determination.

## 1. Why this is a substantive improvement

The supplied review identified a concentrated assurance failure: a small exceptional-start probability was being attached to a floating implementation whose operator construction, deflation and rounding errors had not been bounded. It also identified a false printed polynomial convention and incomplete exact spectral-ordering decisions. Those were not cosmetic matters.

The revision changes the proof dependency graph. Exact rational tensors now define the local projectors. A rank-lift inertia lemma replaces certified numerical deflation. Complete matrix residual checks, with explicit rounding bounds, replace random-start filtering. The relevant common-ground dimensions are supplied by exact polynomial algebra rather than numerical rank estimates.

The reduction is not an assertion that an 11,524-dimensional matrix has the full spectrum of the former 106-million-dimensional operator. It is a new, deliberately weaker **sufficient bound** on its gap. The loss in the constant is the price of a far smaller verification problem.

The original proposed 0.0407 bound is not retained. A smaller rigorously supported bound is more valuable here than a sharper unvalidated one: the central result is a positive lower bound independent of lattice size.

## 2. The new mathematical route

The original down-star cover and three-block hub counting survive. If γ₀ is the weighted nine-site block gap and μ₃ the three-block gap, the global lower bound is γ₀(3μ₃−2). Existence only requires μ₃ > 2/3.

Instead of testing the three-block operator directly, group two block ground projectors first. If their nontrivial overlap is at most η, their two-block excited sum dominates (1−η) times the complement of their common ground space. Compare that common space with the third block, with overlap θ. The weighted two-projection formula then yields

`μ₃ ≥ [2 − η − sqrt(η² + 4(1−η)θ²)] / 2`.

The certified inputs are η < 7/40 and θ² < 1/20. The latter is checked with a suitable grouping on each exact symmetry sector, not asserted for an untested fixed global grouping. They imply

**μ₃ ≥ (73−√313)/80 = 0.691352424838… > 2/3.**

Each private six-site branch has independent dihedral symmetry. Its exact decomposition gives five representative types, with dimensions 78, 66, 21, 15 and 72; the last has two equivalent partners. Together with branch permutations and magnetization restriction, this leaves 15 pair types and 35 triple types. The largest proof matrix has order 11,524, with common-ground dimension 1,085.

## 3. An additional quantitative block result

The nine-site block itself admits a different small reduction. Its three up triangles are disjoint. Compare their product ground space with the ground space of the central down triangle. The triangle gap is exactly 11/14. Four parity-sector tests prove that the squared overlap in this comparison is below 2/3; their orders are 856, 462, 252 and 141.

This yields

**γ₀ ≥ 11(2−√3)/42 = 0.070177169446…**

and hence the explicit uniform constant above. The old nine-site numerical filter is not needed, either for existence or for this new quantitative statement.

## 4. What makes the new acceptance deterministic

In exact compressed coordinates, B is the positive metric and C the contraction of an orthogonal projector. The generalized eigenvalue one has a known common-ground multiplicity k. For a threshold q < 1, the form qB−C already has at least k negative directions.

The checker supplies any dyadic matrix Y with k columns and proves qB−C+YYᵀ positive definite. A rank-k perturbation cannot remove more than k nonpositive directions. Therefore the original form has exactly the known k negative eigenvalues and no further nonpositive eigenvalue. All nontrivial squared overlaps are below q.

This matters because Y does not have to be certified as the exact common-ground projector. Its approximate construction is only a way of proposing a useful lift. The proof comes from the complete positivity check.

Every rational-to-floating conversion used as operator input is bounded against the exact rational value. Contraction errors, metric-product cross terms, symmetrization, dyadic rank-lift products and final residual computations are bounded. The final test computes every matrix entry of a shifted-Cholesky residual. Cholesky success by itself is not treated as a certificate.

For the largest test, the shift is 0.00006103515625 and the total error budget is below 0.000002976. The resulting positive margin exceeds 0.000058. These are congruence-matrix validation quantities, not a second physical gap estimate.

## 5. Executed verification

The full driver reconstructs its exact inputs before any main validation. It executes all 50 hub-angle and four block-angle tests, all of which pass. Its final receipt checks stage status and coverage rather than relying only on successful process termination. Resource exhaustion or a rejected inequality cannot produce `PASS_FULL`.

The smaller exact calculations were also strengthened. The star minimum is proved by an independently built 381-by-381 integer threshold matrix with inertia 336 positive, 45 negative and no zero eigenvalues. Its 44 exact zero-energy states and the explicit quadratic eigenvalue identify the smallest positive star eigenvalue without trusting supplied factors.

The bowtie calculation reconstructs exact highest-weight and leaf-parity blocks of order at most 22. Their dimensions cover all 3,125 physical states. Rational Sturm counts establish the largest noninteger sum-projector eigenvalue, and therefore the exact overlap (47+√4417)/138. The common-range/common-kernel terminology and counts are corrected.

An all-rational 83-dimensional triple test provides a further internal cross-check. At q=1/20 its lifted form has 83 positive eigenvalues. At q=3/100 it has 19 negative eigenvalues, so the deliberately stronger threshold is rejected. The floating checker independently rejects that same overstrong target. Additional negative controls reject fake Cholesky factors, indefinite matrices, insufficient-rank lifts and asymmetric inputs. Matrix-product error enclosures are compared entry by entry with rational arithmetic. Selected checks pass in optimized Python as well.

The original 85-file manifest also verifies. Full provenance is retained. Preliminary exploration and development logs are explicitly separated from the final proof receipts.

## 6. Extensions that were completed

The proof now covers rectangular tori with both side lengths at least three, rather than only equal-sided tori. Exact geometry checks cover all 36 pairs of side lengths from three to eight and support the general coordinate argument.

It also proves the same uniform bound for any finite union of down-star **edge sets**. Partial hubs have one, two or three selected blocks and retain the necessary local gap bounds. The proof accommodates disconnected patches and degenerate ground spaces. Checks cover 511 nonempty selections from a nine-block test family.

This extension must not be advertised as a theorem for arbitrary induced open patches. Adding missing edges can shrink the kernel, and positivity alone does not preserve the relevant positive spectral gap under that change. Perturbative stability and local topological quantum order are not claimed.

The square-lattice result remains exactly delimited: a positive density on the minimal cross cannot represent the specified translation-summed expression once the disjoint-edge products have been removed. It is not a square-lattice gap proof or an obstruction to all larger-support methods.

## 7. Prospective four-star assessment

My internal assessment is that this is now a **credible four-star research candidate**, not merely the corrected conditional/numerical paper envisaged by the review. The reason is the change in what is established: the revised route supplies an explicit, deterministic uniform-gap argument for the undecorated spin-2 kagome target identified by GPW, rather than favorable numerical evidence alone.

The strongest case for originality and significance is the particular two-plus-one/symmetry/certification combination and the resulting theorem. The largest gain in rigor is the completed exact-to-computed contract and exhaustive spectral acceptance, together with removal of the incorrect probability language and repaired local algebra. Better documentation alone would not have produced this change.

This assessment remains conditional on independent confirmation of the new mathematics and checker, and on fuller priority review. All adversarial tests here are internal. They do not warrant claiming that an independent expert has verified the full result. The bounded literature search found no competing proof for the same undecorated kagome model, but it was not an exhaustive MathSciNet/zbMATH/citation-network review. A REF grade is not guaranteed by solving a named problem, by the size of the computation, or by the number of passing tests.

The highest-value remaining step is independent adversarial review of this **new** proof package: specifically the sector-dependent gluing, exact polynomial/projector correspondence and rounding-error validator. Repeating the original large randomized experiment is no longer the missing proof obligation.

## Sources and assessment context

### Mathematical comparator and verification methodology

Guo, W., Pomata, N., & Wei, T.-C. (2021). Nonzero spectral gap in several uniformly spin-2 and hybrid spin-1 and spin-2 AKLT models. *Physical Review Research, 3*, 013255. https://doi.org/10.1103/PhysRevResearch.3.013255

Rump, S. M. (2011). Verified bounds for singular values, in particular for the spectral norm of a matrix and its inverse. *BIT Numerical Mathematics, 51*(2), 367–384. Author-hosted full text: https://www.tuhh.de/ti3/paper/rump/Ru10a.pdf

The full manuscript also cites the original AKLT construction, Knabe counting and the earlier two-dimensional AKLT overlap literature. Verified Cholesky and local-to-global gap arguments are credited as established methods, not presented as inventions of this revision.

### REF context

Research Excellence Framework. (2026). *Section 4—Contributions to Knowledge and Understanding guidance*. REF 2029; last updated 22 January 2026. https://2029.ref.ac.uk/guidance/section-4-contributions-to-knowledge-and-understanding-cku-guidance/

The official page consulted on 30 September 2026 still states that full panel working methods and criteria will be developed. The rating above is therefore a prospective scholarly calibration, not an official assessment or a determination of institutional submission eligibility.
