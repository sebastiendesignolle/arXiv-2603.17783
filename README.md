# Local model for arbitrary qubit POVMs

This directory contains the measurement constructions and verification code for a local model of the two-qubit Werner state under arbitrary qubit positive operator-valued measurements (POVMs).
The finite Bell scenario has 2,306 measurements for Alice, 121 for Bob, and four output labels per measurement.
Its symmetry reduction has 32,712 coordinates.
The verification combines certified measurement approximations with an exact reconstruction of a saved Frank–Wolfe iterate at visibility $v_0=69/100$.

The resulting all-POVM visibility is

$$
 v_{\mathrm{cert}}=\eta_A\eta_B\nu v_0
 =\frac{153579203590380337120294913823}{259479667522670300000000000000}
 \sqrt{\frac{620+185\sqrt5+\sqrt{30(12905+5701\sqrt5)}}{2245}}
 \approx 0.546064019665139.
$$

Alice's shrinking factor is the exact radical in this expression.
Bob's certified factor and the finite-model correction are

$$
 \eta_B=\frac{860491571453}{10^{12}},\qquad
 \nu=\frac{2586644231183439}{2594796675226703}.
$$

The final bound is algebraic; the finite local-model reconstruction and its error bounds use rational arithmetic.

## Reproducing the result

Run the following commands from the BellPolytopes repository root:

```sh
julia --startup-file=no --project=qubit_povm/s2306_121 -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=qubit_povm/s2306_121 qubit_povm/s2306_121/reproduce.jl
```

The script verifies both shrinking factors, reconstructs the saved local model, and writes `result.toml`.
The report contains the exact visibility expression, rational coefficients of its square in a specified number-field basis, certified rational enclosures, reconstruction parameters, and SHA-256 hashes of the proof inputs.
The enclosure bounds describe the displayed certified visibility, not an interval containing the optimal nonlocality threshold.
The report is an output and is never used as a proof input.
No optimization is required to reproduce the certificate.

The option `--skip-shrinking` repeats only the finite-model reconstruction after the coverage proof has been checked separately.
The option `--integer` selects the independent integer reconstruction backend instead of the default exact integer-valued BLAS backend.
The report records whether shrinking-factor verification was performed in that run.
All proof-data paths are relative to the source files, independently of the working directory.

Oscar provides rigorous real-ball and exact algebraic arithmetic for Bob's coverage proof and is loaded only when that proof is requested.
Measurement construction and active-set reconstruction do not require loading Oscar.
The saved `.dat` file uses Julia `Serialization` with MultiFloats weights; the verified environment uses Julia 1.13.2 and MultiFloats 3.3.2.
Retaining the environment versions is recommended when archiving this serialized model.
The project resolves BellPolytopes from the repository root through a relative path.

## Numerical and exact APIs

The constructors evaluate the algebraic measurement data directly in the requested floating-point type:

```julia
include("qubit_povm/s2306_121/Measurements.jl")
using .SingletPOVMs, BellPolytopes, MultiFloats

T = Float64x2
A = measurement_set(; T)
B = template_set(; T)
s = build_scenario(A, B)  # p, o, q, deflate, inflate
```

`Float64`, `Float64x2`, `Float64x3`, and `BigFloat` are covered by the regression tests.
Set the working precision before constructing `BigFloat` measurements.
Floating-point coordinates are evaluated from exact formulas rather than obtained by promoting rounded Float64 data.
The two parties may use different numeric types; the scenario construction promotes them to a common type.

For a new Frank–Wolfe calculation:

```julia
data = certified_povm_data(A, B)
ball = povm_radius(data)
radius = SingletPOVMs.downward_convert(T, ball.radius)

result = bell_frank_wolfe(s.p; o=s.o, prob=true, sym=true,
    deflate=s.deflate, inflate=s.inflate, radius,
    v0=T(69)/100, epsilon=parse(T, "1e-14"), mode_last=-1)

cert = active_set_certificate(result, data; v0=69//100, q=s.q)
bound = exact_alice_visibility(cert)
bound.expression  # Exact radical formula for the certified visibility
bound.squared     # Exact number-field element; the visibility is its nonnegative square root
(bound.lower, bound.upper)  # Rigorous rational enclosure
```

`active_set_certificate` accepts a solver result, an active set, or an `ActiveSetStorageMulti` with weights admitting exact conversion to `Rational{BigInt}`.
MultiFloats limbs are converted exactly, and stored BigFloat precision is preserved.
Use `69//100` to specify exactly 0.69; a floating-point argument certifies its represented value.
The input storage is not modified.

`cert.visibility_lower` retains a rational-only fallback using `922602194543/10^12` for Alice.
The final analytical result uses `exact_alice_visibility(cert)` instead.
Its `prefactor` is the rational number $\eta_B\nu v_0$, and its `alice_squared` field is the exact squared inradius.
`reproduce()` returns the same exact representation as `cert.visibility_exact` and reports it as the main result.
The exact factor needs no additional rounding correction because the residual calculation already encloses the ideal algebraic probabilities.
Call `certify_shrinking()` to verify the measurement-coverage assumptions independently of any finite-model calculation.

## Files and proof inputs

| File | Purpose |
|---|---|
| [Measurements.jl](Measurements.jl) | Module, exact measurement constructions, and numerical constructors. |
| [exact_arithmetic.jl](exact_arithmetic.jl) | Number-field arithmetic, rational intervals, exact linear solves, and support enumeration. |
| [symmetry.jl](symmetry.jl) | Group actions, identical-effect symmetries, dimension count, orbit tensor, and probabilities. |
| [shrinking.jl](shrinking.jl) | Alice's exact inradius, complete coverage entry point, and algebraic effect enclosures. |
| [outer_cover.jl](outer_cover.jl) | Independent verification of Bob's outer approximation and simulation witnesses. |
| [certification.jl](certification.jl) | Local radii, forbidden-output repair, and exact active-set post-processing. |
| [reproduce.jl](reproduce.jl) | Complete verification of the supplied finite model and all-POVM bound. |
| [runtests.jl](runtests.jl) | Arithmetic, numeric-type, ordering, symmetry, and certificate regression tests. |
| [outer_cover.toml](outer_cover.toml) | Algebraic construction recipes, proposed geometry and group actions, and seed metadata. |
| [centers_2_targets.tsv](centers_2_targets.tsv) | The 42,547 canonical outer-target support representatives. |
| [centers_2_witnesses.tsv](centers_2_witnesses.tsv) | Candidate simulation bases and mixtures for those targets. |
| [s2306_121_v6900.dat](s2306_121_v6900.dat) | The 7,421 saved deterministic strategies and their numerical mixture weights. |

The TOML geometry, both TSV files, and the saved active set are required proof inputs.
Their geometric and algebraic claims are checked by the verifier; small numerical LP residuals are not accepted as certificates.
The verification uses no Python code, LP solver, or data outside this directory.

## Measurement approximation

Write each qubit effect as $M_a=w_a I+\vec z_a\cdot\vec\sigma$, with $w_a\geq\|\vec z_a\|_2$, $\sum_a w_a=1$, and $\sum_a\vec z_a=0$.
Depolarization replaces $\vec z_a$ by $\eta\vec z_a$ while preserving the scalar weight.
A finite measurement set has shrinking factor at least $\eta$ if every depolarized POVM can be simulated by a random choice from that set followed by classical output processing.

For a finite vector list $V=(\vec v_i)_i$, define

$$
 W(V)=\left\{\lambda\geq0:\sum_i\lambda_i=1,\quad\sum_i\lambda_i\vec v_i=0\right\}.
$$

Its vertices are supported on at most four indices and define tuples $Q_i=\lambda_i(I+\vec v_i\cdot\vec\sigma)$ summing to the identity.
If the vectors are unit and their convex hull contains the ball of radius $\eta$, these tuples are physical POVMs that simulate all qubit POVMs depolarized by $\eta$.
If instead their convex hull contains the unit ball, the tuples may have nonpositive effects, but simulating every depolarized tuple suffices to simulate every depolarized physical POVM.
These statements follow by decomposing each effect's normalized Bloch vector into the relevant hull, decomposing the resulting weight vector into vertices of $W(V)$, and applying classical output processing.
They apply to arbitrary finite numbers of requested outcomes.

### Alice's measurements and exact shrinking factor

Alice uses the unit directions of a regular icosahedron and its aligned dual dodecahedron.
With $\phi=(1+\sqrt5)/2$, these are the normalized cyclic permutations of $(0,\pm1,\pm\phi)$ and $(0,\pm\phi,\pm\phi^{-1})$, together with the normalized vectors $(\pm1,\pm1,\pm1)$.
The 32 directions give 2,306 vertices of the weight polytope: 16 binary, 120 ternary, and 2,170 quaternary POVMs.

The implementation uses the exact field $K=\mathbb Q(s,t)$, where $s=\sqrt5>0$ and $t=\sqrt{(15-3s)/10}>0$, with basis $(1,s,t,st)$.
The first twelve raw vectors are $t(0,\pm1,\pm\phi)$ and cyclic permutations; the remaining twenty are the dodecahedral directions before normalization.
Every raw vector has squared norm three, so division by $\sqrt3$ gives the unit directions.

For each support of size two, three, or four, the verifier solves the normalization and zero-vector-sum equations exactly and retains full-column-rank systems with strictly positive weights.
Such solutions are exactly the vertices of the weight polytope.
The enumeration considers supports modulo the verified symmetry group and then expands their complete orbits.
Measurements are ordered first by support size and then lexicographically.

The hull verifier exhausts all supporting planes modulo symmetry and identifies 60 triangular facets.
For a raw supporting plane $\vec n\cdot\vec x=d$, its squared distance in unit-vector coordinates is $d^2/(3\vec n\cdot\vec n)$.
Every facet has exactly the squared distance

$$
 \eta_A^2=\frac{620+185s+(505+205s)t}{2245}
 =\frac{620+185\sqrt5+\sqrt{30(12905+5701\sqrt5)}}{2245}.
$$

This is the squared inradius $\eta_{16}^2$ of Designolle et al., *Physical Review Research* **5**, 043059 (2023), arXiv:2302.04721.
The code also verifies the equation $1-40X+530X^2-2480X^3+2245X^4=0$ exactly.
The positive square root gives $\eta_A\approx0.9226021945439895$ and is retained in the final visibility.
`alice_inradius_squared()` returns its exact square as a `Field` element.

### Bob's measurements

Bob's set consists of the distinct images of seven seeds under the full icosahedral group, identifying measurements that differ only by an output permutation.
The entries below give raw directions $\vec d_a$ and scalar weights $w_a$; the effects are $w_a(I+\vec d_a\cdot\vec\sigma/\|\vec d_a\|_2)$.

| Seed | Directions and weights | Orbit size |
|---|---|---:|
| Icosahedral-axis PVM | $\pm(0,1,\phi)$, weights $1/2$. | 6 |
| Coordinate-axis PVM | $\pm(1,0,0)$, weights $1/2$. | 15 |
| Threefold-axis PVM | $\pm(1,1,1)$, weights $1/2$. | 10 |
| Trine | Cyclic permutations of $(7+3\sqrt5,-2,-5-3\sqrt5)$, weights $1/3$. | 20 |
| Regular tetrahedron | $(1,1,1)$, $(1,-1,-1)$, $(-1,1,-1)$, $(-1,-1,1)$, weights $1/4$. | 10 |
| Four-outcome seed | Parameters $(u,v)=(1/4,1/2)$ in the formula below. | 30 |
| Four-outcome seed | Parameters $(u,v)=(1/8,3/4)$ in the formula below. | 30 |

For the final two seeds, set

$$
 x=\frac{2u}{1+u^2},\quad z=\frac{1-u^2}{1+u^2},\quad
 y=\frac{2v}{1+v^2},\quad t=\frac{1-v^2}{1+v^2}.
$$

Their unit directions are $(x,0,z)$, $(-x,0,z)$, $(0,y,-t)$, and $(0,-y,-t)$, with weights $t/[2(z+t)]$ on the first pair and $z/[2(z+t)]$ on the second.
The weights sum to one and their weighted vector sum is zero.
The seven orbit sizes sum to 121.

The order is seed order, group order, and then first occurrence modulo output permutations.
Shared ray directions are identified by exact positive collinearity; rounded keys only locate candidates for that exact check.
Every seed's rank, positivity, normalization, and balance are verified algebraically.

### Output representation and symmetry reduction

The first three icosahedron rays are independent.
Enumerating ordered triples with the same Gram matrix determines all 120 orthogonal symmetries, including reflections.
The code checks orthogonality, the image of every ray, bijectivity, faithfulness, and closure of the permutation group.
Reflections preserve single-qubit effect positivity and classical simulations, while simultaneous reflections preserve the singlet dot products.
The measurement-set symmetry does not require every orthogonal transformation to be a unitary channel.

All measurement arrays use four output labels.
For a binary measurement, both effects are split into equal halves, with supports `[i,j,i,j]` and weights `1/4`.
For a ternary measurement with a unique largest weight, that effect is split into two equal effects; otherwise, a zero fourth effect is appended.
Coarsening these refinements recovers the original measurements and preserves their shrinking statements.
Alice has no zero effects, whereas Bob's twenty trines each have one.

In addition to the geometric symmetries, the reduction includes independent permutations of identical effects within each measurement.
Burnside's lemma gives the invariant dimension as `sum_g fixed_A(g)*fixed_B(g)/120 = 32712`.
The tensor `q` assigns orbit labels in Julia linear index order for `(a,b,x,y)`.
Equal numerical probabilities alone are not used to identify symmetry orbits.
The deflation and inflation maps preserve the full probability-tensor Euclidean norm on invariant tensors.

## Bob's coverage certificate

Start with the twelve unit icosahedron vertices and add the normalized centers of all triangular facets of their convex hull, obtaining 32 unit directions.
Repeat this operation on the new hull to obtain a set $N$ of 92 unit normals.
The outer polytope is

$$
 P=\{\vec r:\vec n\cdot\vec r\leq1\text{ for every }\vec n\in N\}.
$$

Every defining plane is tangent to the unit Bloch sphere, so $P$ contains the entire Bloch ball.
The 180 vertices of $P$ remain at their original positions outside the sphere and are never normalized.
The geometric inradius of $\operatorname{conv}(N)$ is not an additional factor in Bob's certified shrinking.

### Geometry and completeness

`outer_cover.toml` specifies integer normalized-sum recipes rather than rounded coordinates.
The verifier proves that the two construction rounds contain all 20 and all 60 facet centers, respectively.
It matches the proposed action on the initial twelve vertices to the exact orthogonal group and verifies equivariance of every subsequent recipe.

For each proposed final hull face with vertices $a,b,c$, set $n=(b-a)\times(c-a)$ and $d=n\cdot a$, with positive orientation.
Every other normal must lie strictly inside the supporting halfspace.
There are $2\cdot92-4=180$ distinct supporting triangles, proving completeness of the hull description.
Their polar vertices are $r=n/d$; this also determines the induced group action on the outer vertices.
Rigorous real-ball bounds prove strict inequalities, and exact algebraic calculations resolve the required incidences.

The weight polytope on the 180 outer vertices has 90 binary, 600 ternary, and 5,059,530 quaternary extreme tuples, totaling 5,060,220.
To establish these counts, use the $m=90$ antipodal axes and the complete list of central planes containing $k\geq3$ distinct axes:

```
N2 = m
N3 = 2 * sum_planes binomial(k, 3)
N4 = 2 * (binomial(m, 4)
          - sum_planes (binomial(k, 3) * (m-k) + binomial(k, 4)))
```

Each triple of distinct coplanar axes admits two positive sign choices.
Four axes with every triple independent have one dependence with nonzero coefficients and likewise admit two positive sign choices.
A four-axis set containing a dependent triple determines a unique offending plane, which accounts for the subtraction in `N4`.

Plane candidates are grouped using certified bins of the independent entries of the normal projector $nn^T/\|n\|^2$.
Equal planes necessarily occupy the same bin; equality within each candidate bin is then proved by a reflection fixing its rays or by exact algebraic determinant identities.
The target file contains 42,547 distinct canonical supports whose group-orbit sizes sum to the independently established total.
Every support is also checked for the required rank and strictly positive normalized weights.
These checks establish completeness of the target list.

### Simulation witnesses

`centers_2_targets.tsv` has one sorted, one-based support per line.
The corresponding line of `centers_2_witnesses.tsv` has the format

```
k  i1 ... ik  t  nt  [seed group out1 out2 out3 out4 lambda] ...
```

There are `nt` seven-entry terms, without literal brackets.
The seed labels refer to the seven unrefined seeds, and the group labels refer to the `V_actions` order in `outer_cover.toml`.
The verifier matches that order to the exact measurement group.
The four output labels describe a deterministic map to the target's `k` outcomes; unused source outcomes have zero effects.
Decimal coefficients are parsed as Float64 and then interpreted exactly as dyadic rationals.
They select candidate solutions and do not define feasibility tolerances.

The verifier reconstructs the source effects from their exact definitions with rigorous 192-bit real-ball enclosures.
For 41,574 targets, interval Gaussian elimination certifies an invertible simulation basis and strictly positive solution coefficients, including a certified interval for the achievable shrinking factor.
For the remaining 973 targets, it normalizes the proposed nonnegative mixture coefficients exactly and repairs the residual constructively.

For this repair, let $S$ be the simulated POVM and $Q_a=w_a(I+\vec r_a\cdot\vec\sigma)$ the outer target.
Relative to $\mathcal D_t(Q)$, write the scalar and vector errors as $dw_a$ and $d\vec z_a$ and choose

$$
 e\geq\max\left(0,\max_a\frac{dw_a+\|d\vec z_a\|_1}{w_a}\right).
$$

The correction effects $C_a=(e w_a-dw_a)I-d\vec z_a\cdot\vec\sigma$ sum to $eI$ and have scalar part at least the vector's $\ell_1$ norm.
Consequently, $C/e$ is simulable using the three coordinate projective measurements already present in Bob's set.
Explicitly, for coefficients $(s_a,\vec z_a)$ summing to $(1,0)$, choose coordinate measurement $j$ with probability $\sum_a|z_{aj}|$ and map its signs to outcome $a$ with probabilities $(|z_{aj}|\pm z_{aj})/\sum_a|z_{aj}|$.
The remaining scalar mass is produced by a constant output, and zero-probability branches are omitted.
The identity $(S+C)/(1+e)=\mathcal D_{t/(1+e)}(Q)$ proves the repaired bound; no correction is needed when $e=0$.

The minimum certified factor over all target orbits, rounded downward to twelve decimal places, is $860491571453/10^{12}$.
Smaller factors are obtained by mixing with the trivial tuple $w_a I$.
The outer-cover decomposition therefore proves this shrinking factor for every qubit POVM.

## Local radius and finite-model certification

For the singlet, the ideal probability and white-noise tensors are

$$
 p(a,b|x,y)=w^A_{a|x}w^B_{b|y}-\vec z^A_{a|x}\cdot\vec z^B_{b|y},\qquad
 o(a,b|x,y)=w^A_{a|x}w^B_{b|y}.
$$

State white noise therefore need not produce uniform outcome probabilities.
All norms below are Euclidean norms of the full probability tensor.
Zero effects restrict the relevant ball to the normalized, no-signaling affine space where their probabilities remain zero.

### Radius construction

For party $n$ and setting $x$, let $k_{nx}$ be the number of allowed outputs and let $P_{nx}$ be the projector onto zero-sum vectors on that support.
Let $u_n$ denote the single-party behavior uniform on each allowed support and define

```
K_n = max_x k_nx
beta_n = sum_x 1/k_nx
alpha_n = min_{x,a allowed} k_nx * w_n[a,x]
H = 16*(K_A-1)^2*(K_B-1)^2
    + K_A*(K_A-1)/beta_B + K_B*(K_B-1)/beta_A
R_face = alpha_A * alpha_B / sqrt(H)
```

The binary correlation local polytope contains the unit Frobenius ball: for a coefficient matrix $M$ and independent uniform signs $b$, $\mathbb E\|Mb\|_2^2=\|M\|_F^2$, while $\max_a a^TMb=\|Mb\|_1\geq\|Mb\|_2$.
Its support function is therefore at least $\|M\|_F$.
A simultaneous sign flip makes the local marginals zero without changing the correlations.

Introduce binary settings indexed by allowed outcome-setting pairs and use the local response maps $u_n+P_ns/[2(K_n-1)]$.
These are probability distributions because each projected sign coordinate has magnitude at most $2(k_{nx}-1)/k_{nx}$.
They map the binary correlation ball to a local ball of zero-marginal perturbations with radius $1/[4(K_A-1)(K_B-1)]$ about $u_A\otimes u_B$.

Every normalized no-signaling perturbation decomposes orthogonally into a zero-marginal component and the two marginal components $\delta_A\otimes u_B$ and $u_A\otimes\delta_B$.
The single-party probability simplices admit radius $1/\sqrt{K_n(K_n-1)}$.
The marginal components therefore admit radii $\sqrt{\beta_B/[K_A(K_A-1)]}$ and $\sqrt{\beta_A/[K_B(K_B-1)]}$.
Convexly combining the three constructions and applying Cauchy–Schwarz gives a ball of radius $1/\sqrt H$ at the uniform product behavior.

For the physical weights, write $w_n=\alpha_n u_n+(1-\alpha_n)v_n$ with a valid single-party behavior $v_n$.
The product $w_A\otimes w_B$ contains $u_A\otimes u_B$ with mixture weight $\alpha_A\alpha_B$.
Scaling and translating the uniform ball gives the stated face radius.
The code uses rational lower bounds on positive weights and rounds the final square root downward.

An independent covariance construction is also available through `povm_radius(...; method=:covariance)`.
With $L_n=\sum_x k_{nx}(k_{nx}-1)$, its squared denominator is $L_A/\beta_B+L_B/\beta_A+L_A L_B$.
For independent uniform deterministic responses, the centered single-party covariance is $C_n=\bigoplus_x P_{nx}/k_{nx}$, and the product covariance is $G=C_A\otimes C_B+C_A\otimes u_Bu_B^T+u_Au_A^T\otimes C_B$.
The three ranges are orthogonal, and direct inversion gives $\|G^+(d-u)\|^2=L_A/\beta_B+L_B/\beta_A+L_A L_B$ for each deterministic vertex $d$.
Reweighting the uniform deterministic distribution by $1+\langle G^+(d-u),\delta\rangle$ reconstructs $u+\delta$ and remains nonnegative within the reciprocal-square-root radius.
The default `method=:best` selects the stronger of the two certified constructions; here it selects the binary-correlation construction.

### Forbidden outputs and solver radius

Sending each forbidden output uniformly to the allowed outputs is a local channel that fixes both the target and noise behaviors.
Its operator norm is $\sqrt{4/\min_x k_{nx}}$ for each party.
The post-processor applies this channel exactly and uses the face radius with the repaired residual.
A solver correction based on an unrepaired residual must instead use the face radius divided by the product of these norms.
Here Alice always has four allowed outputs and Bob has three on twenty settings, giving

| Quantity | Certified value or decimal display |
|---|---:|
| Face radius $R$ | $2586644231183439/576460752303423488$ |
| Face radius, decimal | $0.004487112471833891$ |
| Safe solver constant `ball.radius` | $0.003885953390246$ |
| Squared repair norm | $4/3$ |

Splitting binary effects reduces the worst repair norm and supplies identical-effect permutations for the symmetry reduction.
The radius and residual both depend on the output representation, so the face radius must not be substituted into a correction using an unrepaired residual.

### Dyadic reconstruction and analytical correction

A dyadic rational has an integer numerator and a power-of-two denominator.
The verifier replaces the stored mixture weights by nonnegative dyadic counts summing exactly to $2^{52}$ and reconstructs the local behavior by integer sums.
The BLAS backend multiplies nonnegative integer-valued arrays with every intermediate sum bounded by $2^{52}$, within the exact integer range of Float64 arithmetic.
Subsequent orbit accumulation uses integers, so this backend gives an exact reconstruction without rounding the input weights to Float64.
The integer backend provides an independent implementation.

Orbit averaging is a convex mixture of locality-preserving relabelings.
Uniform forbidden-output repair commutes with these symmetries and preserves locality.
The target probabilities are enclosed directly from the ideal algebraic effects using rational intervals.
Every residual contribution is rounded outward before its squared norm is accumulated, and an integer square root supplies a certified upper bound.

For the supplied active set, the result is

$$
 \epsilon=\frac{15549553}{1099511627776}\approx1.414223606844\times10^{-5},\qquad
 \nu=\frac{R}{R+\epsilon}=\frac{2586644231183439}{2594796675226703}.
$$

Here $\epsilon$ is an upper bound on the distance from the ideal target to the reconstructed local point, whereas $R$ is a lower bound on the admissible local radius.
The certificate uses 7,421 atoms, a 52-bit mixture grid, 40-bit effect and residual grids, and a 60-bit radius grid.
The solver's cached iterate and numerical objective are not proof inputs.

If $x$ is the reconstructed local point and $p_{v_0}$ is the target, then $y=o+\nu(p_{v_0}-x)/(1-\nu)$ lies in the certified local ball whenever $\epsilon>0$.
Thus $p_{\nu v_0}=\nu x+(1-\nu)y$ is local; if $\epsilon=0$, the target is already local.
Applying the measurement simulations on both sides gives the all-POVM visibility $\eta_A\eta_B\nu v_0$ with Alice's exact algebraic factor.

## Tests

```sh
julia --startup-file=no --check-bounds=yes --project=qubit_povm/s2306_121 qubit_povm/s2306_121/runtests.jl
```

The suite checks exact arithmetic, arbitrary-type construction, measurement and orbit ordering, probability normalization, complete coverage, rejection of malformed witnesses, and the saved local-model certificate.
It also compares the BLAS and integer reconstruction backends on a smaller instance and verifies the exact algebraic visibility representation.
