# SVD, Krylov, and randomized-SVD evaluation for skrCCA

## Conclusion

The reviewer is correct for exactly two cell types: after the PC scores have
been whitened, skrCCA is an ordinary singular-vector problem in the small
PC-space cross-operator. One dense SVD returns all CC axes exactly. We have
therefore added this as a production special case.

The same reduction does **not** hold for three or more cell types. CoPro uses a
separate unit-variance constraint for every cell type. An eigendecomposition of
one stacked block operator replaces these constraints by a single global norm,
so it is a relaxation rather than the same estimator. On colon D3, the block
eigenvector had highly unequal cell-type block norms (0.236--0.688; coefficient
of variation 0.486) and, after block normalization, did not satisfy the CoPro
stationarity equations (relative residual 0.0388 versus (3.1\times10^{-9})
for the CoPro solution).

Krylov/Lanczos SVD is useful only if the retained PC dimension becomes much
larger than it is now. With the current 40 PCs, the exact dense SVD takes about
0.20 ms and operator construction takes about 3.35 ms. `irlba` reduces the SVD
itself to about 0.10 ms, but changes the form-and-solve time only from 3.55 to
3.45 ms. Randomized SVD is slower than `irlba` at every tested dimension and is
approximate, so it provides no benefit here.

## Algebra

For cell type (i), sample (s), let (X_{is}\in\mathbb{R}^{n_{is}\times p})
be its PC score matrix and let
(K_{ij,s}\in\mathbb{R}^{n_{is}\times n_{js}}) be the spatial kernel between
types (i) and (j). CoPro precomputes the small (p\times p) operators

\[
C_{ij}=\sum_s X_{is}^{\mathsf T}K_{ij,s}X_{js}.
\]

Its first-axis objective is

\[
\max_{w_1,\ldots,w_m}\ \sum_{i<j} w_i^{\mathsf T}C_{ij}w_j,
\qquad w_i^{\mathsf T}D_iw_i=1\quad\text{for every }i,
\]

where (D_i=I) for the recommended whitened-PC analysis and
(D_i=\operatorname{diag}(\mathrm{sdev}_i^2)) otherwise.

### Two cell types

For (m=2), define

\[
A=D_1^{-1/2}C_{12}D_2^{-1/2}=U\Sigma V^{\mathsf T}.
\]

Then all axes are obtained exactly as

\[
W_1=D_1^{-1/2}U_{[:,1:r]},\qquad
W_2=D_2^{-1/2}V_{[:,1:r]}.
\]

Thus the reviewer's (Y=KX) intuition is right, with the small clarification
that cross-cell-type CoPro generally has different cell sets on the two sides:
the operator is (X_1^{\mathsf T}K_{12}X_2), not (X^{\mathsf T}KX).

### Multiple samples

The reviewer's suggested absorption of samples into the kernel is also exact.
Let

\[
\widetilde X_i=\begin{bmatrix}X_{i1}\\ \vdots\\ X_{iS}\end{bmatrix},
\qquad
\widetilde K_{ij}=\operatorname{blockdiag}
  (K_{ij,1},\ldots,K_{ij,S}).
\]

Then

\[
\widetilde X_i^{\mathsf T}\widetilde K_{ij}\widetilde X_j
=\sum_s X_{is}^{\mathsf T}K_{ij,s}X_{js}=C_{ij}.
\]

Constructing the large block-diagonal kernel is unnecessary: summing the
(p\times p) operators is the same calculation with less memory. Numerical
checks agreed to (1.14\times10^{-13}). Because projection/deflation is linear
when weights are shared across samples, deflating each sample and then summing
also agreed with summing first and deflating once to
(4.97\times10^{-14}). The implementation now follows the latter route.

#### Memory behavior

The production implementation never constructs either the stacked cell-level
matrix or the block-diagonal kernel above. `compute_Y_multi_slide()` computes
each sample's (p\times p) contribution and reduces those contributions by
summation, one cell-type pair at a time. Its persistent operator storage is
(O(m(m-1)p^2)), independent of the number of cells and samples; parallel
construction temporarily holds at most (O(Sp^2)) values for the current
cell-type pair. For the evaluated case (three cell types, 13 samples, 40 PCs),
these numeric payloads are approximately 75 KiB persistent plus 163 KiB
temporary, excluding small R object overhead.

The pre-change first-axis implementation already used this memory-efficient
summation. The change made here is for later axes: the old path retained the
full pairwise operator list for every sample, requiring
(O(Sm(m-1)p^2)) storage, whereas the new path reuses the single summed
operator. Thus the sample aggregation change reduces memory; it does not trade
memory for speed. The large block matrices in the benchmark script are used
only to test the reviewer's proposed algebra and are never called by CoPro.

### Three or more cell types

One can form a symmetric block matrix (B) whose off-diagonal blocks are
(C_{ij}). A leading eigenvector of (B) solves

\[
\max_z z^{\mathsf T}Bz\quad\text{subject to}\quad \lVert z\rVert=1.
\]

CoPro instead requires one norm constraint per block. Its stationarity
conditions are

\[
\sum_{j\ne i}C_{ij}w_j=\lambda_iD_iw_i,
\]

with a generally different multiplier (\lambda_i) for each cell type. A
single eigenproblem imposes one multiplier and lets signal-rich cell types take
more of the global norm. This is why multiple cell types cannot be absorbed
into an ordinary two-view CCA without changing the estimator. This is the
SUMCOR form of multiset CCA; multiple extensions of two-set CCA are known to be
non-equivalent ([Kettenring, 1971](https://doi.org/10.1093/biomet/58.3.433)).

## Numerical evaluation

All real-data results use the colon D3 vignette data, 40 PCs, spatial bandwidth
0.01, and four axes. The reproducible script is
[`evaluate_svd_krylov_randomized.R`](evaluate_svd_krylov_randomized.R).

### Exact two-type result and runtime

The pre-change coordinate solver and dense SVD produced absolute cosines of
1.000000 for both cell types on every one of the first four axes, and their
objectives equaled the corresponding singular values.

| Method | Median time | Speedup over pre-change |
|---|---:|---:|
| Pre-change coordinate + sequential deflation | 35.00 ms | 1.0x |
| Exported two-step API with exact SVD | 7.00 ms | 5.0x |
| Production path: form (C_{12}) once + dense SVD | 3.55 ms | 9.86x |
| Form once + `irlba` | 3.45 ms | 10.14x |
| Dense SVD only | 0.20 ms | -- |
| `irlba` only | 0.10 ms | -- |
| Randomized SVD only | 0.30 ms | -- |

The production improvement is real within skrCCA, but it should be stated in
proportion: on the earlier full colon D3 pipeline, distance construction,
kernel construction, and downstream scoring took seconds, whereas all skrCCA
optimization over five bandwidths took 0.325 s. SVD improves the optimizer; it
does not remove the cell-level kernel construction cost.

### Krylov and randomized SVD crossover

We compared exact base/LAPACK SVD, augmented implicitly restarted Lanczos
bidiagonalization (`irlba`), and randomized range finding with oversampling 10
and two power iterations. Synthetic square operators had slowly decaying
singular values (\sigma_j=j^{-0.6}), which is deliberately less favorable to
randomized approximation.

| (p) | Exact SVD | `irlba` | Randomized | `irlba` max subspace error | Randomized max subspace error |
|---:|---:|---:|---:|---:|---:|
| 40 | 0.20 ms | 0.10 ms | 0.25 ms | (2.0\times10^{-7}) | 0.0075 |
| 100 | 0.95 ms | 0.175 ms | 0.40 ms | (1.3\times10^{-8}) | 0.0148 |
| 250 | 7.18 ms | 0.65 ms | 0.75 ms | (1.2\times10^{-9}) | 0.0150 |
| 500 | 36.0 ms | 1.55 ms | 2.35 ms | (9.2\times10^{-10}) | 0.0076 |

Lanczos bidiagonalization is designed for a few singular triplets of large
sparse matrices ([Baglama and Reichel, 2005](https://doi.org/10.1137/04060593X)),
and randomized decompositions are most compelling for matrices that are large,
sparse, distributed, or expensive to pass over
([Halko, Martinsson, and Tropp, 2011](https://doi.org/10.1137/090771806)).
Those are not the current (40\times40) CoPro operators. If CoPro later retains
hundreds of PCs, `irlba` is the appropriate adaptive alternative; randomized
SVD is not indicated by these results.

### Multi-type block eigensolver

| Method | Objective | KKT residual | Median solver time |
|---|---:|---:|---:|
| CoPro block-coordinate solver | 2171.928 | (3.08\times10^{-9}) | 41 ms |
| Block eigenvector, then normalize each block | 2170.116 | 0.0388 | 2 ms |
| Block eigenvector, then refine with CoPro | 2171.928 | (1.59\times10^{-9}) | 26 ms |

The block eigenvector is not a valid final solution. It is a useful symmetric
initializer and reduced this small solve from 41 to 26 ms after refinement, but
the absolute gain is only 15 ms and has so far been checked on one real data
set. We therefore do not make it the default.

## Higher CC axes

For two cell types, higher axes require no special algorithm: the same SVD
returns them all simultaneously and exactly.

For more than two cell types, the previous pairwise rank-one subtraction is not
equivalent to orthogonal projection because a multi-set stationary vector is
not a singular vector of every pairwise (C_{ij}). On colon D3, the maximum
within-cell-type orthogonality error across four axes was 0.423. Replacing it
with

\[
C_{ij}\leftarrow (I-W_iW_i^{\mathsf T})C_{ij}
                    (I-W_jW_j^{\mathsf T})
\]

reduced that error to (6.7\times10^{-16}). The new projected axes remained
very similar to the old axes for CC2 (absolute cosines 0.996--1.000 across cell
types); differences increased for the weaker CC4 (0.877--0.978), where the old
directions contained the most lower-axis leakage. Projection is now the default
for unweighted analyses with more than two cell types. The recommended
`scalePCs = TRUE` path is unweighted; the less-used weighted multi-type path
retains its previous deflation until a weighted projection is formally added.

We also tested a simultaneous product-Stiefel update for (r) axes,

\[
\max_{W_i^{\mathsf T}W_i=I_r}
\sum_{i<j}\operatorname{tr}(W_i^{\mathsf T}C_{ij}W_j),
\]

using polar-factor block updates. It was fast (18 ms), orthogonal to machine
precision, and achieved a slightly higher total four-axis objective than the
old sequential method (4391.823 versus 4391.797). However, it trades a small
amount of the leading-axis objective for total subspace fit (CC1 2170.793
versus 2171.928) and changes the ordered, sequential null used for higher-axis
testing. It is promising for future work but is not a drop-in replacement.

## Implemented changes

- `runSkrCCA()` now detects the ordinary two-cell-type, non-transferred case,
  forms the aggregated PC-space operator once, and returns all requested axes
  from one exact dense SVD.
- The public first/subsequent-component optimizers use the same exact two-type
  solution. A transferred first axis retains conditional sequential deflation.
- Multi-slide optimization retains the existing memory-efficient operator sum
  and now reuses that sum for higher-axis deflation; no stacked matrix or
  block-diagonal sample kernel is materialized.
- Unweighted analyses with more than two cell types use full projected
  deflation, restoring within-cell-type orthogonality.
- The cached conditional-permutation CC1 path uses the same exact SVD as the
  uncached path.
- New tests cover exact SVD equivalence, weighted constraints, multi-sample
  aggregation, and multi-type higher-axis orthogonality.

Machine-readable results are in the CSV files in this directory, with the most
compact combined output in [`summary.txt`](summary.txt).

## Suggested Methods replacement

For cell type (i) in sample (s), let (X_{is}) denote the whitened PC-score
matrix and (K_{ij,s}) the spatial kernel linking cells of types (i) and
(j). We first computed the PC-space cross-operators
(C_{ij}=\sum_sX_{is}^{\mathsf T}K_{ij,s}X_{js}); summing these operators is
equivalent to stacking samples and using a block-diagonal spatial kernel. We
then maximized
(\sum_{i<j}w_i^{\mathsf T}C_{ij}w_j) subject to
(\lVert w_i\rVert_2=1) for every cell type. With two cell types this problem
reduces exactly to the SVD of (C_{12}), and all requested canonical axes were
obtained in one decomposition. With three or more cell types the separate norm
constraints yield a multiset (SUMCOR) CCA problem with a distinct Lagrange
multiplier for each cell type, so a single stacked eigendecomposition is not
equivalent; we used block-coordinate updates
(w_i\leftarrow\operatorname{normalize}(\sum_{j\ne i}C_{ij}w_j)). For later
axes, the operators were projected onto the orthogonal complements of the
previous weight vectors for each cell type before applying the same updates.

## Suggested reviewer response

We thank the reviewer for this helpful suggestion. We agree that the
two-cell-type case has a direct SVD solution. Writing
(C_{12}=\sum_sX_{1s}^{\mathsf T}K_{12,s}X_{2s}), skrCCA with whitened PCs
maximizes (w_1^{\mathsf T}C_{12}w_2) subject to
(\lVert w_1\rVert=\lVert w_2\rVert=1); hence (w_1) and (w_2) are the left
and right singular vectors of (C_{12}). We have revised the implementation to
use one exact dense SVD for all requested axes in this case. On the colon D3
data (40 PCs, four axes), the SVD weights were identical to the iterative
weights (absolute cosine 1.000000 for every axis) and reduced the skrCCA solve
from 35.0 to 3.55 ms (9.9-fold).

We also agree algebraically that samples can be absorbed into the spatial
operator. Stacking
the sample-specific (X_{is}) and defining
(\widetilde K_{ij}=\operatorname{blockdiag}(K_{ij,1},\ldots,K_{ij,S})) gives
(\widetilde X_i^{\mathsf T}\widetilde K_{ij}\widetilde X_j
=\sum_sX_{is}^{\mathsf T}K_{ij,s}X_{js}). We now perform the equivalent, more
memory-efficient operation of summing the small PC-space operators. This sum
was already used for the first axis; we now also reuse it for higher-axis
deflation rather than storing sample-specific residual operators. At no point
do we materialize the stacked cell-level matrix or block-diagonal kernel.

The reason we retain a different strategy for three or more cell types is that
the problem is then multiset rather than ordinary two-set CCA. CoPro constrains
each cell type separately,
(w_i^{\mathsf T}w_i=1), leading to stationary equations
(\sum_{j\ne i}C_{ij}w_j=\lambda_iw_i) with cell-type-specific multipliers.
A single eigendecomposition of a stacked block matrix instead imposes one
global norm and is not equivalent. Empirically, its cell-type block norms
ranged from 0.236 to 0.688 on colon D3 and the block-normalized vector failed
the CoPro stationarity equations (relative residual 0.0388, versus
(3.1\times10^{-9}) for the block-coordinate solution). A block eigenvector
was useful only as an initializer.

Finally, we evaluated truncated Krylov and randomized SVD. At the current
(40\times40) operator size, exact SVD required 0.20 ms; Krylov reduced this by
only 0.10 ms and did not materially change the operator-plus-solve time, while
randomized SVD was slower and introduced approximation error. Krylov methods
became advantageous for synthetic operators with hundreds of PCs, so they are
a sensible future adaptive path if that regime is supported, but they do not
improve the present analysis. We have added these clarifications and the direct
two-type SVD formulation to Methods.
