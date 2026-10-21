# Genuine orbit reduction, preserving the original q labels exactly.
# Outcome flags are pairs (outcome, input). Each geometric transformation acts
# on input and outcome simultaneously. Local permutations of identical effects
# are then quotiented independently for each input and party. For A != B there
# is no party exchange. Burnside counts avoid allocating q to predict its size.
# Orbit labels are assigned at the first occurrence in Julia tensor order.
# Equal numerical probabilities alone never justify identifying coordinates.
support_key(ids) = (u = sort(unique(filter(!iszero, ids))); Tuple(vcat(u, zeros(Int, 4 - length(u)))))
function flag_actions(M::MeasurementSet)
    for x in 1:length(M), a in 1:4, b in 1:4
        if M.supports[a, x] == M.supports[b, x] && M.exact.weights[a, x] != M.exact.weights[b, x]
            throw(ArgumentError("repeated rays must have equal weights for the local swap reduction"))
        end
    end
    lookup = Dict(support_key(M.supports[:, x]) => x for x in 1:length(M))
    actions = Vector{Int}[]
    for g in M.group
        action = zeros(Int, 4length(M))
        for x in 1:length(M)
            ids = sort(unique(filter(!iszero, M.supports[:, x])))
            key = support_key(vcat(g[ids], zeros(Int, 4 - length(ids))))
            haskey(lookup, key) || error("selected measurements are not group-closed")
            y = lookup[key]
            for a in 1:4
                id = M.supports[a, x]
                target = id == 0 ? 0 : g[id]
                source = findall(==(id), M.supports[:, x])
                dest = findall(==(target), M.supports[:, y])
                length(source) == length(dest) || error("refinement is not covariant")
                b = dest[findfirst(==(a), source)]
                action[a+4(x-1)] = b + 4(y - 1)
                M.exact.weights[a, x] == M.exact.weights[b, y] || error("invalid exact outcome action")
            end
        end
        push!(actions, action)
    end
    return actions
end

# Quotient identical effects (including zeros) locally. These are genuine
# outcome permutations, independent at each input and at each party.
function collapsed_flags(M)
    return [a + 4(x - 1) for x in 1:length(M) for a in 1:4 if findfirst(==(M.supports[a, x]), M.supports[:, x]) == a]
end

function same_measurements(A, B)
    return A.supports == B.supports && A.vertices == B.vertices && A.weights == B.weights
end

"""Exact integer Burnside count, before allocating any dense tensors."""
function orbit_dimension(A::MeasurementSet, B::MeasurementSet = A; exchange::Bool = same_measurements(A, B))
    exchange && !same_measurements(A, B) && throw(ArgumentError("party exchange requires identical sets"))
    actsA, actsB = flag_actions(A), flag_actions(B)
    flagsA, flagsB = collapsed_flags(A), collapsed_flags(B)
    total = 0
    for (a, b) in zip(actsA, actsB)
        total += count(i -> a[i] == i, flagsA) * count(i -> b[i] == i, flagsB)
        exchange && (total += count(i -> a[a[i]] == i, flagsA))
    end
    denominator = length(actsA) * (exchange ? 2 : 1)
    total % denominator == 0 || error("noninteger Burnside count")
    return total ÷ denominator
end

function outcome_range(M, flag)
    a = mod1(flag, 4)
    x = cld(flag, 4)
    return findall(==(M.supports[a, x]), M.supports[:, x])
end

function orbit_tensor(A, B, exchange, dim)
    actsA, actsB = flag_actions(A), flag_actions(B)
    rangesA = [outcome_range(A, i) for i in 1:4length(A)]
    rangesB = [outcome_range(B, i) for i in 1:4length(B)]
    q = zeros(UInt32, 4, 4, length(A), length(B))
    label = UInt32(0)
    for y in 1:length(B), x in 1:length(A), b in 1:4, a in 1:4
        q[a, b, x, y] == 0 || continue
        label += 1
        for (ga, gb) in zip(actsA, actsB)
            fa, fb = ga[a+4(x-1)], gb[b+4(y-1)]
            xx, yy = cld(fa, 4), cld(fb, 4)
            for aa in rangesA[fa], bb in rangesB[fb]
                q[aa, bb, xx, yy] = label
                exchange && (q[bb, aa, yy, xx] = label)
            end
        end
    end
    Int(label) == dim || error("orbit tensor disagrees with Burnside count")
    return q
end

"""Return singlet probabilities `p[a,b,x,y]` and physical white-noise `o`.
The latter is w[a,x]*w[b,y], generally NOT the uniform distribution 1/16.
"""
function probability_tensors(A::MeasurementSet, B::MeasurementSet = A)
    T = promote_type(eltype(A.vertices), eltype(B.vertices))
    p = zeros(T, 4, 4, length(A), length(B))
    o = similar(p)
    dots = Matrix{T}(A.vertices) * Matrix{T}(B.vertices)'
    for y in 1:length(B), x in 1:length(A), b in 1:4, a in 1:4
        u, v = A.supports[a, x], B.supports[b, y]
        o[a, b, x, y] = A.weights[a, x] * B.weights[b, y]
        p[a, b, x, y] = u == 0 || v == 0 ? zero(T) : o[a, b, x, y] * (1 - dots[u, v])
    end
    return p, o
end

"""
    build_scenario(A, B=A; orbit_limit=100_000, exchange=...)

Check the reduced dimension before allocating q, p or o. The default cap is
100,000; pass a larger `orbit_limit` or `nothing` to explore complete sets. q contains group orbits, not classes of equal probability values.
Use `bell_frank_wolfe(s.p; o=s.o, prob=true, sym=true,
                     deflate=s.deflate, inflate=s.inflate, ...)`.
"""
function build_scenario(A::MeasurementSet, B::MeasurementSet = A; orbit_limit::Union{Nothing,Int} = 100_000, exchange::Bool = same_measurements(A, B))
    orbit_limit === nothing || orbit_limit > 0 || throw(ArgumentError("orbit_limit must be positive or nothing"))
    dimension = orbit_dimension(A, B; exchange)
    (orbit_limit === nothing || dimension <= orbit_limit) || throw(ArgumentError("$dimension orbits exceed the $orbit_limit cap; choose smaller sets, or explicitly increase orbit_limit for a complete-set experiment"))
    q = orbit_tensor(A, B, exchange, dimension)
    p, o = probability_tensors(A, B)
    deflate, inflate = BellPolytopes.build_deflate_inflate_q(eltype(p), q)
    return (; A, B, p, o, q, deflate, inflate, dimension, eta_A = A.eta, eta_B = B.eta, shrinking_product = downward_convert(eltype(p), exact_real(A.eta, "eta") * exact_real(B.eta, "eta")), exchange)
end
