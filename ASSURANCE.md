# Assurance and proof dependencies

## Status of the main claim

The main output is a **deterministic computer-assisted proof candidate with a fully executed internal replay**, not a randomized numerical conjecture and not an externally certified theorem. “Deterministic” describes the acceptance argument: it does not claim infallibility of the implementation or hardware.

The claim is a positive uniform spectral gap in the specified spin-2 kagome model, normalized by rank-nine two-site edge projectors. The proof uses no numerical lower bound imported from the v0.2 filter. The explicit final expression is positive by elementary rational square comparisons.

## Exact input contract

The polynomial monomial corresponds to the orthonormal spin state divided by the square root of the binomial coefficient. The wrong reciprocal convention is an explicit rejected regression case. All ground spaces are images of multiplication by products of coprime singlet polynomials. This supplies injectivity, dimensions and spin content.

Integer Gram matrices are reconstructed from those maps. A mathematical upper bound on every signed accumulation is checked before int64 arithmetic. Branch symmetry maps use integer orbit sums. Private and external rational coordinate frames are invertible exact congruences. Positive local pivots and reconstruction identities are checked. The private metric is diagonal with entries in [1,4); this is not a claimed condition-number bound on every larger matrix.

The local projector factors are defined over the rationals. Every magnetization block solve is checked exactly, and all off-block entries of the Gram and factor matrices are checked to vanish. Every rational-to-binary64 conversion used as proof input is compared to the exact rational value; its maximum error is retained. Floating zero is never used to discard a nonzero exact tensor coefficient.

## Spectral contract

For a compressed exact projector, B is the positive coordinate Gram form and C is its positive contraction, with 0 ≤ C ≤ B. The generalized eigenvalue one occurs on the exact common ground space. The common-space dimension is obtained from injective polynomial boundary coordinates, including all mixed symmetry partners.

For a threshold q < 1, the form qB − C has at least k negative eigenvalues. The checker constructs an arbitrary dyadic matrix Y with k columns and proves qB − C + YYᵀ positive definite. The rank bound then proves that qB − C has exactly k negative eigenvalues, no zero eigenvalues and no further nonpositive direction. Hence every nontrivial squared principal cosine is strictly below q.

The computed columns of Y do **not** need certified containment, orthogonality or accurate normalization. Their only trusted structural property is their number. They suggest a useful rank lift; positive definiteness of the lifted exact form is what is verified.

The triple grouping is sector-dependent. Fifteen unordered pair types and 35 unordered triple types cover the exact branch decomposition, with an explicit chosen pair on each triple sector. A single fixed grouping is not asserted to satisfy the same bound globally without its additional ordered tests.

## Floating arithmetic contract

The implementation uses conventional binary64, round to nearest, bounded-length standard BLAS dot products, and no overflow. It does not permit Strassen-like cancellation patterns outside the stated error analysis, approximate GPU multipliers, or reduced-precision accumulation. Fused multiply-add is acceptable. Backend identity and environment are recorded.

A factor 4ku, where u = 2^-53, bounds length-k dot-product roundoff in the allowed regime. Positive bound computations use outward nextafter operations. Exact input radii, both first-order contributions and cross terms, are propagated. A deliberately large absolute underflow allowance is included. Symmetrization is part of the error budget, not an assumed identity for an unsymmetrized matrix.

For the final target T, the checker bounds ||T − A||₂. It forms a dyadically shifted center and obtains a Cholesky factor L as a candidate. It computes every entry of the factor residual, in row blocks, with dot-product, subtraction and norm-accumulation errors bounded. Acceptance requires the shift to exceed the sum of the input, shift and full-residual bounds. Successful Cholesky without the residual test is never sufficient.

For the largest matrix, order 11,524, the shift is 2^-14. The total certified error is below 2.976 × 10^-6, leaving a positive congruence-matrix margin above 5.80 × 10^-5. That is a validation margin, not the physical lattice-gap constant.

## Executed checks

The authoritative full replay includes 50 hub-angle positive-definiteness checks and four weighted-block-angle checks. It also reconstructs the star threshold matrix and proves its inertia with integer arithmetic; reconstructs the bowtie highest-weight blocks and orders their roots with rational Sturm counts; checks all six branch representatives/partners and full dimension coverage; checks 36 rectangular torus instances and 511 finite block selections; and verifies the normalization regression.

Adversarial tests reject an indefinite target, an insufficient-rank lift, an asymmetric input, a fake Cholesky factor and an overstrong triple threshold. Every entry of a collection of matrix-product enclosures is checked against exact rational multiplication, including cancellation and underflow examples. A separate 83-dimensional triple target is independently assembled entry by entry over the rationals. Integer inertia is (83 positive, 0 negative, 0 zero) at q = 1/20, but (64 positive, 19 negative, 0 zero) at q = 3/100. Selected checkers are run with Python optimization enabled.

A full replay is accepted only after all stages return successfully and their authoritative JSON statuses and expected sector counts are checked. Existing authoritative receipts are removed at the start, preventing a stale success record from standing in for a missing calculation. A resource failure is a failed execution, not mathematical evidence against the claimed gap and never an acceptance.

## Remaining assurance work

Independent specialist review of the tensor/geometry correspondence, the sector-dependent gluing argument and the roundoff checker has not occurred. Neither a second independently written full-size validated implementation nor proof-assistant formalization is supplied. The all-rational cross-check is one small sector, not an exact-arithmetic replacement of all 54 main tests. Agreement among implementations within this workflow is internal corroboration, not external replication.

The original 106,935,090-dimensional filter was not rerun. The old 0.0407 constant, its claimed implementation-level failure probabilities, a square-lattice gap, arbitrary induced-patch gaps, perturbative stability and local topological quantum order are not certified by this release.
