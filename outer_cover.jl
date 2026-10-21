"""
Independent verifier for the 92-normal cover of Bob's 121-measurement family.

Loaded lazily into SingletPOVMs by certify_shrinking(). Every geometry entry and
simulation column is checked; the stored floating data are proposals only.
Rigorous 192-bit Arb enclosures prove strict inequalities and feasible bases.
Oscar's exact algebraic real numbers resolve coplanar incidences when necessary.
"""
module OuterCover92
import Oscar as O
using LinearAlgebra, Combinatorics, TOML
const P = parentmodule(@__MODULE__)
const Q = Rational{BigInt}
const B = O.ArbField(192)
const A = O.algebraic_closure(O.QQ)
const ROOT = @__DIR__
dot3(a, b) = sum(a[i] * b[i] for i in 1:3)
cross3(a, b) = [a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1]]
normalize3(v) = v ./ sqrt(dot3(v, v))
lower(x) = Q(BigFloat(x, RoundDown))
upper(x) = Q(BigFloat(x, RoundUp))
fraction(x) = string(numerator(x), "/", denominator(x))
ball(x::P.Field) = sum(B(c) * z for (c, z) in zip(x.c, (B(1), sqrt(B(5)), sqrt((15 - 3sqrt(B(5))) / 10), sqrt(B(5)) * sqrt((15 - 3sqrt(B(5))) / 10))))
logmsg(args...) = (println(args...); flush(stdout))

mutable struct Geometry
    mesh::Dict
    V::Vector
    W::Vector
    va::Vector
    wa::Vector
    anti::Int
    exactV::Dict{Int,Vector}
    exactW::Dict{Int,Vector}
end

function exact_vertex(g::Geometry, i; outer = false)
    cache = outer ? g.exactW : g.exactV
    haskey(cache, i) && return cache[i]
    v = if outer
        a, b, c = [exact_vertex(g, j) for j in g.mesh["faces"][i]]
        n = cross3(b - a, c - a)
        n ./ dot3(n, a)
    elseif i <= 12
        a, b, k = g.mesh["raw"][i]
        normalize3(circshift([A(0), A(a), A(b) * (1 + sqrt(A(5))) / 2], k))
    else
        terms = g.mesh["operations"][i-12]
        normalize3([sum(A(c) * exact_vertex(g, j)[k] for (j, c) in terms) for k in 1:3])
    end
    cache[i] = v
end

"Verify exact group action from the construction DAG, without algebraic expansion."
function geometry(mesh, G)
    V = Vector{typeof(B(0))}[]
    for (a, b, k) in mesh["raw"]
        push!(V, normalize3(circshift([B(0), B(a), B(b) * (1 + sqrt(B(5))) / 2], k)))
    end
    opkey(terms) = Tuple(sort([(Int(i), Int(c)) for (i, c) in terms]))
    lookup = Dict(opkey(t) => j + 12 for (j, t) in enumerate(mesh["operations"]))
    @assert length(lookup) == length(mesh["operations"])
    for terms in mesh["operations"]
        @assert all(t -> 1 <= t[1] <= length(V) && t[2] > 0, terms)
        push!(V, normalize3([sum(B(c) * V[j][k] for (j, c) in terms) for k in 1:3]))
    end
    va = [Int.(g) for g in mesh["V_actions"]]
    order = Dict(Tuple(a[1:12]) => i for (i, a) in enumerate(G.actions))
    perm = [order[Tuple(a[1:12])] for a in va]
    @assert sort(perm) == collect(1:120)
    for action in va
        @assert sort(action) == collect(eachindex(V))
        for (i, terms) in enumerate(mesh["operations"])
            @assert lookup[opkey([(action[j], c) for (j, c) in terms])] == action[i+12]
        end
    end
    anti = only(findall(i -> G.matrices[perm[i]] == -Matrix{P.Field}(I, 3, 3), 1:120))
    faces = [sort(Int.(f)) for f in mesh["faces"]]
    @assert length(Set(Tuple.(faces))) == length(faces) == 2length(V) - 4
    flookup = Dict(Tuple(f) => i for (i, f) in enumerate(faces))
    wa = [[flookup[Tuple(sort(g[f]))] for f in faces] for g in va]
    @assert wa == mesh["W_actions"]
    seen = Set{Int}()
    reps = Int[]
    radii = typeof(B(0))[]
    for (i, face) in enumerate(faces)
        i in seen && continue
        push!(reps, i)
        foreach(g -> push!(seen, g[i]), wa)
        a, b, c = V[face]
        n = cross3(b - a, c - a)
        d = dot3(n, a)
        if d < 0
            n = -n
            d = -d
        end
        @assert d > 0
        # Every nonincident vertex is STRICTLY inside the supporting plane.
        @assert all(j -> j in face || dot3(n, V[j]) < d, eachindex(V))
        push!(radii, d * d / dot3(n, n))
    end
    @assert length(seen) == length(faces)
    # All points have unit norm by construction; strict supporting triangles
    # and F=2V-4 prove that these are ALL facets, not merely a selection.
    W = [
        begin
            a, b, c = V[f]
            n = cross3(b - a, c - a)
            n ./ dot3(n, a)
        end for f in faces
    ]
    g = Geometry(mesh, V, W, va, wa, anti, Dict{Int,Vector}(), Dict{Int,Vector}())
    best = argmin(Float64.(radii))
    chosen = reps[best]
    exactradius(i) = begin
        a, b, c = [exact_vertex(g, j) for j in faces[i]]
        n = cross3(b - a, c - a)
        dot3(n, a)^2 / dot3(n, n)
    end
    for j in eachindex(reps)
        j == best && continue
        radii[j] > radii[best] && continue
        @assert exactradius(reps[j]) >= exactradius(chosen)
    end
    eta = sqrt(radii[best])
    info = Dict(
        "vertices" => length(V),
        "W_vertices" => length(W),
        "facet_orbits" => length(reps),
        "eta_decimal" => Float64(eta),
        "eta_lower" => fraction(lower(eta)),
        "eta_upper" => fraction(upper(eta)),
        "limiting_face" => faces[chosen],
        "algebraic_definition" => "sqrt(dot(cross(v_b-v_a,v_c-v_a),v_a)^2 / dot(cross(v_b-v_a,v_c-v_a),cross(v_b-v_a,v_c-v_a)))"
    )
    g, info, perm, () -> exactradius(chosen)
end

"""
Exact counts without enumerating C(n,4) supports. For m antipodal axes and
each distinct central plane containing k>=3 axes:
N2=m, N3=2sum C(k,3), N4=2[C(m,4)-sum(C(k,3)(m-k)+C(k,4))].
For four axes with every triple independent, the one linear dependence has
four nonzero coefficients and exactly two sign choices give positive weights.
A bad four-set belongs to a unique plane class containing a dependent triple.

Plane classes are grouped by certified bins of the six entries of nn'/||n||².
This removes sign/pivot ambiguities. Equal planes MUST receive the same bin;
coplanarity within every multi-pair bin is proved exactly, modulo symmetry.
"""
function count_povms(g; outer = false, deadline = Inf, verbose = true)
    X = outer ? g.W : g.V
    actions = outer ? g.wa : g.va
    axes = Int.(g.mesh[outer ? "W_axes" : "V_axes"])
    m = length(axes)
    @assert length(X) == 2m && sort(vcat(axes, actions[g.anti][axes])) == collect(eachindex(X))
    pairbins = Dict{NTuple{6,Int},Vector{Tuple{Int,Int}}}()
    # A transcendental offset cannot place an algebraic projector entry
    # exactly on a bin boundary (a rational offset can do so).
    D = B(10^12 + 39)
    shift = O.const_pi(B) / 10
    for i in 1:(m-1), j in (i+1):m
        n = cross3(X[axes[i]], X[axes[j]])
        s = dot3(n, n)
        @assert s > 0
        key = Tuple(Int(floor(BigInt, D * n[a] * n[b] / s + shift)) for a in 1:3 for b in a:3)
        push!(get!(pairbins, key, Tuple{Int,Int}[]), (i, j))
    end
    planes = Tuple[]
    for pairs in values(pairbins)
        length(pairs) == 1 && continue
        ids = sort!(unique(vcat([[i, j] for (i, j) in pairs]...)))
        @assert length(pairs) == binomial(length(ids), 2)
        push!(planes, Tuple(ids))
    end
    allplanes = Set(planes)
    seen = Set{Tuple}()
    planeorbits = 0
    axisof = zeros(Int, length(X))
    for (i, j) in enumerate(axes)
        axisof[j] = axisof[actions[g.anti][j]] = i
    end
    axactions = [axisof[a[axes]] for a in actions]
    for ids in sort(planes)
        ids in seen && continue
        time() < deadline || error("plane certificate time limit")
        planeorbits += 1
        # A nontrivial orthogonal map fixing >=3 distinct rays pointwise
        # fixes a plane. This often proves the largest incidences immediately.
        geometric = [axes[i] for i in ids]
        reflection = any(a -> a != collect(eachindex(X)) && all(i -> a[i] == i, geometric), actions)
        if !reflection
            a, b = [exact_vertex(g, axes[i]; outer) for i in ids[1:2]]
            n = cross3(a, b)
            @assert all(i -> iszero(dot3(n, exact_vertex(g, axes[i]; outer))), ids[3:end])
        end
        for a in axactions
            image = Tuple(sort([a[i] for i in ids]))
            @assert image in allplanes
            push!(seen, image)
        end
        verbose && planeorbits % 25 == 0 && logmsg("    ", outer ? "W" : "V", " plane orbits certified: ", planeorbits)
    end
    @assert seen == allplanes
    triangles = sum(binomial(BigInt(length(p)), 3) for p in planes; init = big(0))
    bad = sum(binomial(BigInt(length(p)), 3) * (m - length(p)) + binomial(BigInt(length(p)), 4) for p in planes; init = big(0))
    counts = [BigInt(m), 2triangles, 2(binomial(BigInt(m), 4) - bad)]
    @assert counts == g.mesh[outer ? "W_counts" : "V_counts"]
    (; counts, total = sum(counts), planes, axes, axisof, planeorbits)
end

function target_weights(g, ids, pc)
    k = length(ids)
    X = g.W[ids]
    c = if k == 2
        @assert g.wa[g.anti][ids[1]] == ids[2]
        [B(1), B(1)]
    elseif k == 3
        ax = sort(pc.axisof[ids])
        @assert any(p -> all(i -> i in p, ax), pc.planes)
        n = cross3(X[1], X[2])
        rows = setdiff(1:3, [argmax(abs.(Float64.(n)))])
        cross3([x[rows[1]] for x in X], [x[rows[2]] for x in X])
    else
        [(-1)^(i - 1) * dot3(X[setdiff(1:4, [i])[1]], cross3(X[setdiff(1:4, [i])[2]], X[setdiff(1:4, [i])[3]])) for i in 1:4]
    end
    w = c ./ sum(c)
    @assert all(x -> x > 0, w)
    w
end

"Interval Gaussian elimination; nonzero pivot intervals certify nonsingularity."
function solve_ball(M, b)
    n = length(b)
    C = hcat(M, b)
    for j in 1:n
        p = j - 1 + argmax(abs.(Float64.(C[j:n, j])))
        C[j, :], C[p, :] = copy(C[p, :]), copy(C[j, :])
        pivot = C[j, j]
        (pivot > 0 || pivot < 0) || return nothing
        C[j, j:(n+1)] .= C[j, j:(n+1)] ./ pivot
        for i in (j+1):n
            v = C[i, j]
            C[i, j:(n+1)] .= C[i, j:(n+1)] .- v .* C[j, j:(n+1)]
        end
    end
    x = fill(B(0), n)
    for j in n:-1:1
        x[j] = C[j, end] - sum(C[j, l] * x[l] for l in (j+1):n; init = B(0))
    end
    x
end

function simulators(G, perm)
    seeds = P.seed_measurements(G.matrices)
    @assert seeds[2].V == [P.Field[1, 0, 0], P.Field[-1, 0, 0]] && seeds[2].w == fill(P.Field(1 // 2), 2)
    @assert all(j -> any(R -> P.matvec(R, P.Field[1, 0, 0]) == P.Field.(collect(1:3) .== j), G.matrices), 1:3)
    stored = TOML.parsefile(joinpath(ROOT, "outer_cover.toml"))["seeds"]
    @assert length(stored) == length(seeds) == 7
    decode(c) = P.Field(Tuple(P.rat.(c)))
    for (s, e) in zip(seeds, stored)
        @assert s.V == [decode.(v) for v in e["vectors"]]
        @assert s.w == decode.(e["weights"]) && s.N == decode(e["norm"])
    end
    E = Array{typeof(B(0)),4}(undef, 4, 4, 7, 120)
    fill!(E, B(0))
    for j in 1:7, h in 1:120
        seed = seeds[j]
        scale = sqrt(ball(seed.N))
        for (a, (v, w)) in enumerate(zip(seed.V, seed.w))
            E[1, a, j, h] = ball(w)
            for (c, z) in enumerate(P.matvec(G.matrices[perm[h]], v))
                E[c+1, a, j, h] = ball(w * z) / scale
            end
        end
    end
    E
end

function certify_simulations(g, pc, E; path = joinpath(ROOT, "centers_2_witnesses.tsv"), verbose = true)
    name = g.mesh["name"]
    targets = [parse.(Int, split(line)) for line in eachline(joinpath(ROOT, name * "_targets.tsv"))]
    @assert length(Set(Tuple.(targets))) == length(targets)
    total = big(0)
    for ids in targets
        orbit = Set(Tuple(sort(a[ids])) for a in g.wa)
        @assert Tuple(ids) == minimum(orbit)
        total += length(orbit)
    end
    @assert total == pc.total
    lines = readlines(path)
    @assert length(lines) == length(targets) "one witness is required for every target orbit"
    best = one(Q)
    numericbest = 1.0
    basis_count = 0
    repair_count = 0
    worst = 0
    for (index, (line, ids)) in enumerate(zip(lines, targets))
        s = split(line)
        k = parse(Int, s[1])
        @assert k == length(ids) && parse.(Int, s[2:(k+1)]) == ids
        t = Q(parse(Float64, s[k+2]))
        nt = parse(Int, s[k+3])
        @assert 0 < t <= 1 && length(s) == k + 3 + 7nt
        weights = target_weights(g, ids, pc)
        columns = Vector{typeof(B(0))}[]
        coeff = Q[]
        for j in 0:(nt-1)
            seed, group, mapping, lambda = parse_term(s[(k+4+7j):(k+10+7j)], k)
            F = fill(B(0), 4, k)
            for a in 1:4
                F[:, mapping[a]] .+= E[:, a, seed, group]
            end
            push!(columns, vcat(vec(F[:, 1:(k-1)]), B(1)))
            push!(coeff, lambda)
        end
        n = 4(k - 1) + 1
        rhs = fill(B(0), n)
        rhs[end] = B(1)
        direction = fill(B(0), n)
        for a in 1:(k-1)
            rhs[4a-3] = weights[a]
            direction[(4a-2):4a] = weights[a] .* g.W[ids[a]]
        end
        bound = nothing
        if nt + 1 == n
            z = solve_ball(hcat(columns..., -direction), rhs)
            if z !== nothing && all(x -> x > 0, z) && z[end] < 1
                bound = lower(z[end])
                basis_count += 1
            end
        end
        if bound === nothing
            # Same constructive coordinate-PVM repair as the old certificate.
            coeff ./= sum(coeff)
            sim = sum(B(c) .* v for (c, v) in zip(coeff, columns))
            S = fill(B(0), 4, k)
            S[:, 1:(k-1)] = reshape(sim[1:(n-1)], 4, k - 1)
            S[:, k] = [B(1), B(0), B(0), B(0)] - vec(sum(S[:, 1:(k-1)]; dims = 2))
            eps = zero(Q)
            for a in 1:k
                dw = S[1, a] - weights[a]
                dz = [S[c+1, a] - B(t) * weights[a] * g.W[ids[a]][c] for c in 1:3]
                eps = max(eps, (upper(dw) + sum(max(abs(lower(x)), abs(upper(x))) for x in dz)) / lower(weights[a]))
            end
            bound = t / (1 + eps)
            repair_count += 1
        end
        if bound < best
            best = bound
            worst = index
        end
        numericbest = min(numericbest, Float64(t))
        verbose && index % 1000 == 0 && logmsg("    certified ", index, "/", length(targets), "; bound ", Float64(best))
    end
    decimal = fld(numerator(best) * big(10)^12, denominator(best)) // big(10)^12
    Dict(
        "status" => "certified",
        "orbits" => length(targets),
        "eta_rational" => fraction(decimal),
        "eta_decimal" => Float64(decimal),
        "numerical_minimum" => numericbest,
        "worst_target" => worst,
        "basis_certificates" => basis_count,
        "repair_certificates" => repair_count
    )
end

"Validate a single stored simulator term; floating weights are interpreted exactly."
function parse_term(term, k)
    @assert length(term) == 7 && 2 <= k <= 4
    seed, group = parse.(Int, term[1:2])
    mapping = parse.(Int, term[3:6])
    lambda = Q(parse(Float64, term[7]))
    @assert 1 <= seed <= 7 && 1 <= group <= 120 && all(a -> 1 <= a <= k, mapping) && lambda >= 0
    seed, group, mapping, lambda
end

"Recompute the entire 92-normal coverage proof; stored reports are not inputs."
function certify(; witnesses = joinpath(ROOT, "centers_2_witnesses.tsv"), verbose = true)
    mesh = TOML.parsefile(joinpath(ROOT, "outer_cover.toml"))["mesh"]
    @assert mesh["name"] == "centers_2" && length(mesh["operations"]) == 80
    G = P.exact_group()
    verbose && logmsg("Certifying the 92-normal outer cover...")
    g, info, perm, _ = geometry(mesh, G)
    @assert length(g.V) == 92 && length(g.W) == 180
    # Verify that the construction consists of two complete rounds of facet
    # centers (20 then 60), not an arbitrary list of normalized positive sums.
    for (n, range) in ((12, 13:32), (32, 33:92))
        @assert length(range) == 2n - 4
        for i in range
            terms = mesh["operations"][i-12]
            @assert length(terms) == 3 && all(t -> 1 <= t[1] <= n && t[2] == 1, terms)
            ids = first.(terms)
            @assert length(unique(ids)) == 3
            a, b, c = g.V[ids]
            normal = cross3(b - a, c - a)
            d = dot3(normal, a)
            if d < 0
                normal = -normal
                d = -d
            end
            @assert d > 0 && all(j -> j in ids || dot3(normal, g.V[j]) < d, 1:n)
        end
    end
    pc = count_povms(g; outer = true, verbose)
    @assert pc.counts == [90, 600, 5059530] && pc.total == 5060220
    result = certify_simulations(g, pc, simulators(G, perm); path = witnesses, verbose)
    @assert result["status"] == "certified" && result["orbits"] == 42547
    result
end
end # module OuterCover92
