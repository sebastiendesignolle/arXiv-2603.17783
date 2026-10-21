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
 \eta_B=\frac{860491571453}{10^{12}}\qquad\text{and}\qquad\nu=\frac{2586644231183439}{2594796675226703}.
$$

The final bound is algebraic; the finite local-model reconstruction and its error bounds use rational arithmetic.

## Reproducing the result

Run the following commands from the repository root:

```sh
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project reproduce.jl
```

The script verifies both shrinking factors, reconstructs the saved local model, and writes `result.toml`.
The report contains the exact visibility expression, rational coefficients of its square in a specified number-field basis, certified rational enclosures, reconstruction parameters, and SHA-256 hashes of the proof inputs.
The enclosure bounds describe the displayed certified visibility, not an interval containing the optimal nonlocality threshold.
The report is an output and is never used as a proof input.
No optimization is required to reproduce the certificate.

The option `--skip-shrinking` repeats only the finite-model reconstruction after the coverage proof has been checked separately.
The option `--integer` selects the independent integer reconstruction backend instead of the default exact integer-valued BLAS backend.
The report records whether shrinking-factor verification was performed in that run.

Oscar provides rigorous real-ball and exact algebraic arithmetic for Bob's coverage proof and is loaded only when that proof is requested.
Measurement construction and active-set reconstruction do not require loading Oscar.
The saved `.dat` file uses Julia `Serialization` with MultiFloats weights; the verified environment uses Julia 1.13.2 and MultiFloats 3.3.2.

## Numerical and exact APIs

The constructors evaluate the algebraic measurement data directly in the requested floating-point type:

```julia
include("Measurements.jl")
using .SingletPOVMs, BellPolytopes, MultiFloats

T = Float64x2
A = measurement_set(; T)
B = template_set(; T)
s = build_scenario(A, B)  # p, o, q, deflate, inflate
```

For a new Frank–Wolfe calculation:

```julia
data = certified_povm_data(A, B)
ball = povm_radius(data)
radius = SingletPOVMs.downward_convert(T, ball.radius)

result = bell_frank_wolfe(s.p; o=s.o, prob=true, sym=true,
    deflate=s.deflate, inflate=s.inflate, radius,
    v0=T(69)/100, epsilon=T(1e-10), mode_last=-1)

cert = active_set_certificate(result, data; v0=69//100, q=s.q)
bound = exact_alice_visibility(cert)
bound.expression  # Exact radical formula for the certified visibility
bound.squared     # Exact number-field element; the visibility is its nonnegative square root
(bound.lower, bound.upper)  # Rigorous rational enclosure
```

`active_set_certificate` accepts a solver result, an active set, or an `ActiveSetStorageMulti` with weights admitting exact conversion to `Rational{BigInt}`.
Use `69//100` to specify exactly 0.69; a floating-point argument certifies its represented value.

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
| [outer_cover.toml](outer_cover.toml) | Algebraic construction recipes, proposed geometry and group actions, and seed metadata. |
| [centers_2_targets.tsv](centers_2_targets.tsv) | The 42,547 canonical outer-target support representatives. |
| [centers_2_witnesses.tsv](centers_2_witnesses.tsv) | Candidate simulation bases and mixtures for those targets. |
| [s2306_121_v6900.dat](s2306_121_v6900.dat) | The 7,421 saved deterministic strategies and their numerical mixture weights. |
