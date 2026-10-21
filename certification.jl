# Exact active-set reconstruction. Generic dyadic integer BLAS and rounding
# already live in BellPolytopes; only the POVM-specific face repair and
# algebraic residual enclosure belong here. No Python or parent-folder files.
"""
    povm_radius(data; bits=60, method=:best)
    povm_radius(A, B=A; coordinate_bits=40, bits=60)

A proven rational Euclidean radius about the POVM state-white-noise tensor
`o[a,b,x,y]=wA[a,x]*wB[b,y]`. The norm is over the FULL probability tensor.
Return rational `radius`, `face_radius`, their squared lower bounds,
`uniform_fraction`, `projection_norm_squared`, `full_support`, and `data`.

When `full_support=false`, the actual ball is confined to the normalized
no-signaling FACE where padded zero outcomes stay zero. `face_radius` is its
radius. The returned `radius` includes a further contraction for a local
projection onto that face. It is a valid residual-correction constant for
BellPolytopes even when its iterate uses padded outcomes, provided target
and noise lie in the face. It is NOT a ball in the larger affine space.
The post-processor projects the orbit-averaged model onto the face using
exact rational arithmetic and uses `face_radius`.
By default `method=:best` selects the larger of the covariance and binary-lift
bounds. Use `:covariance` to reproduce the original radius.
The constructor radius is conservative; use the exact post-processor for a
final bound, even when displaying floating solver corrections.
"""
function povm_radius(data; bits::Int = 60, method::Symbol = :best)
    method in (:best, :covariance, :binary_lift) || throw(ArgumentError("unknown radius method"))
    bits > 0 || throw(ArgumentError("bits must be positive"))
    D = BigInt(data.denominator)
    alpha = Rational{BigInt}[]
    beta = Rational{BigInt}[]
    L = BigInt[]
    projection2 = big(1) // 1
    kmax = Int[]
    full_support = true
    for (lo, hi) in zip(data.lo, data.hi)
        a = big(1) // 1
        b = big(0) // 1
        ell = big(0)
        minimum_k = size(lo, 2)
        maximum_k = 0
        for x in axes(lo, 3)
            positive = findall(>(0), lo[1, :, x])
            k = length(positive)
            k > 0 || error("increase coordinate_bits to resolve positive weights")
            all(lo[1, j, x] > 0 || (lo[1, j, x] == hi[1, j, x] == 0) for j in axes(lo, 2)) || error("increase coordinate_bits to separate zero and positive effects")
            sum(lo[1, :, x]) <= D <= sum(hi[1, :, x]) || error("weights cannot sum to one")
            a = min(a, k * minimum(lo[1, positive, x]) // D)
            b += 1 // big(k)
            ell += k * (k - 1)
            full_support &= k == size(lo, 2)
            minimum_k = min(minimum_k, k)
            maximum_k = max(maximum_k, k)
        end
        push!(alpha, a)
        push!(beta, b)
        push!(L, ell)
        projection2 *= size(lo, 2) // big(minimum_k)
        push!(kmax, maximum_k)
    end
    H = L[1] / beta[2] + L[2] / beta[1] + L[1] * L[2]
    H > 0 || error("the affine space has dimension zero")
    t = prod(alpha)
    covariance2 = t * t / H
    # Lift the unit Frobenius ball of binary correlations through local
    # response maps u + P*s/(2*(kmax-1)), and add the orthogonal marginals.
    Hlift = 16big(prod(kmax .- 1))^2 + kmax[1] * (kmax[1] - 1) / beta[2] + kmax[2] * (kmax[2] - 1) / beta[1]
    binary_lift2 = t * t / Hlift
    chosen = method == :best ? (binary_lift2 >= covariance2 ? :binary_lift : :covariance) : method
    face2 = chosen == :binary_lift ? binary_lift2 : covariance2
    r2 = face2 / projection2
    scale = big(1) << bits
    r = isqrt(fld(numerator(r2) * scale^2, denominator(r2))) // scale
    r > 0 || error("increase bits to obtain a positive radius")
    face_radius = isqrt(fld(numerator(face2) * scale^2, denominator(face2))) // scale
    return (;
        radius = r,
        radius_squared = r2,
        face_radius,
        face_radius_squared = face2,
        uniform_fraction = t,
        projection_norm_squared = projection2,
        full_support,
        method = chosen,
        covariance_radius_squared = covariance2 / projection2,
        binary_lift_radius_squared = binary_lift2 / projection2,
        data
    )
end
povm_radius(A::MeasurementSet, B::MeasurementSet = A; coordinate_bits = 40, bits = 60, method = :best) = povm_radius(certified_povm_data(A, B; coordinate_bits); bits, method)

"""
    primal_target(ball; relative_loss=1//100)

Exact solver residual/primal targets for a fractional certification loss δ.
Real inputs (including MultiFloats and BigFloat) are converted to their exact
represented rational values; use e.g. `1//100` for an exact decimal fraction.
Uses the safe solver constant: epsilon <= radius*δ/(1-δ), primal=epsilon²/2.
The final certificate also includes geometry/grid rounding; leave some margin.
"""
function primal_target(ball; relative_loss = 1 // 100)
    delta = exact_real(relative_loss, "relative_loss")
    0 < delta < 1 || throw(ArgumentError("relative_loss must lie in (0,1)"))
    residual = ball.radius * delta / (1 - delta)
    return (; relative_loss = delta, residual_target = residual, primal_target = residual^2 / 2)
end

_storage(s::BellPolytopes.ActiveSetStorageMulti) = s
_storage(s::Tuple) = length(s) == 8 ? _storage(s[5]) : throw(ArgumentError("expected a bell_frank_wolfe result"))
_storage(s) = _storage(BellPolytopes.ActiveSetStorage(s))

# Apply uniform local repair AFTER orbit averaging. Since the physical
# relabeling group preserves the support, this channel commutes with its
# Reynolds projection and the result is still constant on each q orbit.
function uniform_face_projection(orbit, data, q)
    all(all(>(0), h[1, :, :]) for h in data.hi) && return orbit
    mult = orbit.multiplicities
    values = [BigInt(s) // big(m) for (s, m) in zip(orbit.numerators, mult)]
    indices = CartesianIndices(orbit.dims)
    linear = LinearIndices(orbit.dims)
    reps = zeros(Int, length(mult))
    for i in 1:length(indices)
        j = q === nothing ? i : q[i]
        reps[j] == 0 && (reps[j] = i)
    end
    value(a, b, x, y) =
        let i = linear[a, b, x, y]
            values[q === nothing ? i : q[i]]
        end
    projected = map(reps) do rep
        a, b, x, y = Tuple(indices[rep])
        ha, hb = data.hi
        (ha[1, a, x] == 0 || hb[1, b, y] == 0) && return big(0) // 1
        badA = findall(iszero, ha[1, :, x])
        badB = findall(iszero, hb[1, :, y])
        ka, kb = 4 - length(badA), 4 - length(badB)
        z = value(a, b, x, y)
        for aa in badA
            z += value(aa, b, x, y) / ka
        end
        for bb in badB
            z += value(a, bb, x, y) / kb
        end
        for aa in badA, bb in badB
            z += value(aa, bb, x, y) / (ka * kb)
        end
        z
    end
    scale = foldl(lcm, denominator.(projected); init = big(1))
    numerators = [numerator(v) * div(scale, denominator(v)) * m for (v, m) in zip(projected, mult)]
    return merge(orbit, (; numerators, denominator = BigInt(orbit.denominator) * scale))
end

# Integer interval products; T is selected once, outside the tensor loop.
function product_bounds(::Type{T}, al, ah, bl, bh) where {T}
    p, q, r, s = T(al) * bl, T(al) * bh, T(ah) * bl, T(ah) * bh
    return min(p, q, r, s), max(p, q, r, s)
end
function residual_squares(::Type{T}, data, orbit, q, visibility, residual_bits) where {T}
    al, bl = data.lo
    ah, bh = data.hi
    F = T(1) << residual_bits
    vn, vd = T(numerator(visibility)), T(denominator(visibility))
    divisor = vd * (T(1) << (2data.bits - residual_bits))
    local_floor = [T(fld(BigInt(s) * BigInt(F), BigInt(orbit.denominator) * m)) for (s, m) in zip(orbit.numerators, orbit.multiplicities)]
    H = zero(T)
    for (linear, i) in enumerate(CartesianIndices(orbit.dims))
        a, b, x, y = Tuple(i)
        j = q === nothing ? linear : q[linear]
        if ah[1, a, x] == 0 || bh[1, b, y] == 0
            orbit.numerators[j] == 0 || error("q does not preserve the zero-outcome face")
            continue
        end
        wl, wh = product_bounds(T, al[1, a, x], ah[1, a, x], bl[1, b, y], bh[1, b, y])
        zl = zero(T)
        zh = zero(T)
        for c in 2:4
            l, h = product_bounds(T, al[c, a, x], ah[c, a, x], bl[c, b, y], bh[c, b, y])
            zl += l
            zh += h
        end
        lower = fld(vd * wl - vn * zh, divisor)
        upper = fld(vd * wh - vn * zl, divisor)
        e = max(abs(lower - local_floor[j]), abs(upper - local_floor[j])) + 1
        H += e * e
    end
    return BigInt(H)
end

"""
    active_set_certificate(storage_or_result, data; v0, q=nothing,
        weight_bits=52, residual_bits=40, radius_bits=60, backend=:blas,
        row_block=512, column_block=512, atom_block=4096)

Produce a RATIONAL singlet white-noise visibility bound from a probability
active set, without rebuilding its quadratic cache. Stored weights and `v0`
may use any real type supporting exact `Rational{BigInt}` conversion, including
Float64x2, Float64x3, BigFloat, and rationals. They are never narrowed to Float64.
`v0` must describe the finite visibility used for that saved model. A floating
`v0` certifies its exact represented value; use `69//100` for exactly 0.69.
The output bounds remain Rational{BigInt}, independently of the input type.
Use exact_alice_visibility(cert) to retain Alice's exact algebraic inradius
in the final all-POVM visibility instead of its rational fallback.
The 52-bit weight grid and integer-valued Float64 BLAS are deliberate exact
arithmetic stages, not conversions of the high-precision input weights.
`data=certified_povm_data(A,B)` fixes the order of settings and outcomes.
Pass the SAME genuine orbit tensor `q` used by the run; omitting it is valid
but usually gives a worse residual. The caller must ensure that q describes
locality-preserving relabelings, e.g. build_scenario(A,B).q.

1. Round and normalize stored weights to dyadic counts.
2. Reconstruct and orbit-average using exact blocked integer-valued GEMM.
3. Apply uniform local repair of zero outcomes, using exact rational orbit sums.
   `repair=:deterministic` retains the previous pre-averaging response remap.
4. Bound every full-tensor residual entry with algebraic/dyadic intervals.
5. Apply nu=r/(r+epsilon) and both certified shrinking factors.

Returns `visibility_lower=etaA*etaB*v0*nu`, `finite_visibility`, `radius`,
`residual_norm_upper`, exact squared residual upper bound, and replay data.
No floating target, cached iterate, or numerical norm is trusted. The input
storage is not modified. Algebraic rounding error is included in the bound.
"""
function active_set_certificate(result, data; v0, q = nothing, weight_bits::Int = 52, residual_bits::Int = 40, radius_bits::Int = 60, backend::Symbol = :blas, row_block::Int = 512, column_block::Int = 512, atom_block::Int = 4096, repair::Symbol = :uniform)
    repair in (:uniform, :deterministic) || throw(ArgumentError("unknown face repair"))
    visibility = exact_real(v0, "v0")
    0 <= visibility <= 1 || throw(ArgumentError("v0 must lie in [0,1]"))
    1 <= residual_bits <= 52 && residual_bits <= 2data.bits || throw(ArgumentError("require 1 <= residual_bits <= min(52,2coordinate_bits)"))
    storage = _storage(result)
    length(storage.ax) == 2 && storage.o == [4, 4] || throw(ArgumentError("expected a bipartite four-output probability storage"))
    ax = copy.(storage.ax)
    remapped = 0
    for n in 1:2
        size(ax[n], 2) == size(data.lo[n], 3) || throw(DimensionMismatch("stored settings and POVMs disagree"))
        for x in axes(ax[n], 2)
            replacement = findfirst(>(0), data.lo[n][1, :, x])
            replacement === nothing && error("no positive effect")
            for k in axes(ax[n], 1)
                a = ax[n][k, x]
                1 <= a <= 4 || error("invalid stored response")
                if data.hi[n][1, a, x] == 0
                    repair == :deterministic && (ax[n][k, x] = replacement)
                    remapped += 1
                end
            end
        end
    end
    weights = dyadic_weights([exact_real(w, "active-set weight") for w in storage.weights]; bits = weight_bits)
    orbit = dyadic_probability_orbit_sums(ax, weights.numerators, storage.o, q; backend, row_block, column_block, atom_block)
    raw_orbit = orbit
    repair == :uniform && (orbit = uniform_face_projection(orbit, data, q))
    ball = povm_radius(data; bits = radius_bits)
    # Bound all products and the total squared sum before choosing Int128.
    F = big(1) << residual_bits
    D = big(1) << data.bits
    safe = max(4denominator(visibility) * D^2, prod(BigInt.(orbit.dims)) * (5F + 1)^2) <= typemax(Int128)
    H = residual_squares(safe ? Int128 : BigInt, data, orbit, q, visibility, residual_bits)
    root = isqrt(H)
    root^2 < H && (root += 1)
    epsilon = root // F
    nu = ball.face_radius / (ball.face_radius + epsilon)
    finite_visibility = visibility * nu
    return (;
        visibility_lower = prod(data.eta) * finite_visibility,
        finite_visibility,
        nu,
        v0 = visibility,
        radius = ball.face_radius,
        residual_norm_upper = epsilon,
        residual_squared_upper = H // F^2,
        eta = data.eta,
        remapped_responses = remapped,
        weights,
        orbit,
        ax,
        q,
        data,
        backend,
        weight_bits,
        residual_bits,
        radius_bits,
        repair,
        raw_orbit
    )
end
