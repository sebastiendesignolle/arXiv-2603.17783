module SingletPOVMs
using BellPolytopes, Combinatorics, LinearAlgebra, TOML, SHA
export measurement_set, template_set, build_scenario, scenario, certified_povm_data, certify_shrinking, MeasurementSet, povm_radius, active_set_certificate, primal_target, alice_inradius_squared, exact_alice_visibility
include("exact_arithmetic.jl")

# Rational fallback for APIs returning rational bounds; the final analytical
# visibility uses sqrt(alice_inradius_squared()) through exact_alice_visibility.
const ETA_A = big(922602194543) // 10^12
const ETA_B = big(860491571453) // 10^12

"""
The 32 directions: 12 icosahedron vertices and 20 aligned dual vertices.
Raw vectors are scaled by sqrt(3), so they lie in Field and all have norm²=3.
The order is the legacy Julia order: a,b,k loops for the first orbit, then
lexicographically ordered dual vertices. Changing it invalidates saved atoms.
"""
function alice_rays()
    ico = [RT .* circshift(Field[0, a, b*PHI], k) for a in (-1, 1) for b in (-1, 1) for k in 0:2]
    dual = vcat([Field[a, b, c] for a in (-1, 1) for b in (-1, 1) for c in (-1, 1)], [circshift(Field[0, a*PHI, b/PHI], k) for a in (-1, 1) for b in (-1, 1) for k in 0:2])
    sort!(dual; by = Tuple)
    V = vcat(ico, dual)
    @assert all(v -> edot(v, v) == Field(3), V)
    V
end

"""
Enumerate all 120 orthogonal symmetries in their original order. The ordered
image of the first three independent icosahedron vertices determines a unique
linear transformation. Retain exactly the triples with the same Gram matrix.
Verify orthogonality, all vertex images, faithfulness and permutation closure.
Reflections act on Bloch vectors and preserve effects and singlet probabilities;
they need not be implemented as unitary conjugations for this argument.
"""
function exact_group(V = alice_rays())
    B = hcat(V[1:3]...)
    gram = B' * B
    matrices = Matrix{Field}[]
    actions = Vector{Int}[]
    lookup = Dict(Tuple(v) => i for (i, v) in enumerate(V))
    for ids in permutations(1:12, 3)
        C = hcat(V[ids]...)
        C' * C == gram || continue
        R = permutedims(hcat([exact_solve(B', C[i, :]) for i in 1:3]...))
        @assert R * R' == Matrix{Field}(I, 3, 3)
        g = [lookup[Tuple(matvec(R, v))] for v in V]
        push!(matrices, R)
        push!(actions, g)
    end
    verify_group(actions)
    (; matrices, actions)
end
function verify_group(actions)
    permutations = Set(Tuple(g) for g in actions)
    @assert length(actions) == length(permutations) == 120
    @assert all(g -> sort(g) == collect(1:length(g)), actions)
    @assert all(Tuple(g[h]) in permutations for g in actions for h in actions)
    nothing
end

"""
The seven seeds of the 121-setting family, before covariant outcome refinement.
Return (raw directions, scalar weights, common squared norm).

1--3: PVM axes (0,1,phi), (1,0,0), (1,1,1).
4: a planar trine orthogonal to (1,1,1); a mirror fixes its orientation.
5: a regular tetrahedron, weights 1/4.
6--7: directions (+/-x,0,z), (0,+/-y,-t), with x=2u/(1+u²),
   z=(1-u²)/(1+u²), y=2v/(1+v²), t=(1-v²)/(1+v²).
   Weights are t/[2(z+t)] twice and z/[2(z+t)] twice. The parameters are
   (u,v)=(1/4,1/2) and (1/8,3/4). Thus the two opposite horizontal pairs
   balance the vertical component despite their unequal weights.
"""
function seed_measurements(G)
    seeds = NamedTuple[]
    function add(V, w)
        N = edot(V[1], V[1])
        @assert positive(N) && all(v -> edot(v, v) == N, V)
        @assert sum(w) == Field(1) && all(positive, w)
        @assert all(iszero, sum((a .* v for (a, v) in zip(w, V))))
        @assert exact_weights(V, collect(1:length(V))) == w
        push!(seeds, (; V, w, N))
    end
    for v in (Field[0, 1, PHI], Field[1, 0, 0], Field[1, 1, 1])
        add([v, -v], fill(Field(1 // 2), 2))
    end
    U = Field[1, 1, 1]
    R = first(R for R in G if matvec(R, U) == U && edot(R[1, :], ecross(R[2, :], R[3, :])) == Field(-1))
    v = Field[1, 0, 0] + R[:, 1]
    v = v .- edot(v, U) / 3
    add([v, circshift(v, 1), circshift(v, 2)], fill(Field(1 // 3), 3))
    add([Field[a, b, c] for a in (-1, 1) for b in (-1, 1) for c in (-1, 1) if a * b * c == 1], fill(Field(1 // 4), 4))
    for (u, v) in ((Field(1 // 4), Field(1 // 2)), (Field(1 // 8), Field(3 // 4)))
        x = 2u / (1 + u * u)
        z = (1 - u * u) / (1 + u * u)
        y = 2v / (1 + v * v)
        t = (1 - v * v) / (1 + v * v)
        add([Field[x, 0, z], Field[-x, 0, z], Field[0, y, -t], Field[0, -y, -t]], [t / (2(z + t)), t / (2(z + t)), z / (2(z + t)), z / (2(z + t))])
    end
    seeds
end

struct ExactMeasurements
    rays::Vector{Vector{Field}}
    norms::Vector{Field}
    supports::Matrix{Int}
    weights::Matrix{Field}
    group::Vector{Vector{Int}}
end
"Floating effects in T, with their exact algebraic construction attached."
struct MeasurementSet{T<:AbstractFloat}
    name::Symbol
    vertices::Matrix{T}
    supports::Matrix{Int}
    weights::Matrix{T}
    group::Vector{Vector{Int}}
    eta::T
    exact::ExactMeasurements
end
Base.length(M::MeasurementSet) = size(M.supports, 2)

# Hash only for candidate lookup; the positive-collinearity check proves actual
# equality of normalized rays. This handles shared rays with different norms.
function ray_key(v, N)
    # 256-bit evaluation is independent of the caller's arithmetic and precision.
    setprecision(BigFloat, 256) do
        Tuple(round(Int, 10^12 * floating(BigFloat, x) / sqrt(floating(BigFloat, N))) for x in v)
    end
end

function refine_exact(rays, norms, S, W, G)
    S = copy(S)
    W = copy(W)
    for x in axes(S, 2)
        k = count(!iszero, S[:, x])
        if k == 2
            @assert W[1, x] == W[2, x] == Field(1 // 2)
            S[3:4, x] = S[1:2, x]
            W[:, x] .= Field(1 // 4)
        elseif k == 3
            order = sortperm(W[1:3, x]; rev = true)
            if positive(W[order[1], x] - W[order[2], x])
                a = order[1]
                S[4, x] = S[a, x]
                W[a, x] /= 2
                W[4, x] = W[a, x]
            end
        end
    end
    ExactMeasurements(rays, norms, S, W, G)
end

const CONSTRUCTIONS = Ref{Any}(nothing)
function exact_constructions()
    CONSTRUCTIONS[] !== nothing && return CONSTRUCTIONS[]
    V = alice_rays()
    group = exact_group(V)
    G = group.matrices
    cat = positive_catalogue(V, group.actions)
    ids = sort!(collect(keys(cat.valid)); by = s -> (length(s), s))
    @assert length(ids) == 2306
    S = hcat([vcat(collect(s), zeros(Int, 4 - length(s))) for s in ids]...)
    W = hcat([vcat(cat.valid[s], fill(Field(0), 4 - length(s))) for s in ids]...)
    A = refine_exact(V, fill(Field(3), length(V)), S, W, group.actions)

    rays = Vector{Field}[]
    norms = Field[]
    lookup = Dict{NTuple{3,Int},Int}()
    supports = Vector{Int}[]
    weights = Vector{Field}[]
    models = Dict{NTuple{4,Int},Int}()
    function rayid(v, N)
        key = ray_key(v, N)
        if haskey(lookup, key)
            i = lookup[key]
            @assert all(iszero, ecross(v, rays[i])) && positive(edot(v, rays[i]))
            return i
        end
        push!(rays, v)
        push!(norms, N)
        lookup[key] = length(rays)
    end
    seeds = seed_measurements(G)
    seed_counts = Int[]
    for seed in seeds
        before = length(supports)
        for R in G
            ids = [rayid(matvec(R, v), seed.N) for v in seed.V]
            @assert length(unique(ids)) == length(ids)
            order = sortperm(ids)
            s = vcat(ids[order], zeros(Int, 4 - length(ids)))
            w = vcat(seed.w[order], fill(Field(0), 4 - length(ids)))
            key = Tuple(s)
            if haskey(models, key)
                @assert weights[models[key]] == w
            else
                push!(supports, s)
                push!(weights, w)
                models[key] = length(supports)
            end
        end
        push!(seed_counts, length(supports) - before)
    end
    @assert length(supports) == 121
    actions = Vector{Int}[]
    for R in G
        g = Int[]
        for (v, N) in zip(rays, norms)
            z = matvec(R, v)
            i = lookup[ray_key(z, N)]
            @assert all(iszero, ecross(z, rays[i])) && positive(edot(z, rays[i]))
            push!(g, i)
        end
        push!(actions, g)
    end
    verify_group(actions)
    B = refine_exact(rays, norms, hcat(supports...), hcat(weights...), actions)
    CONSTRUCTIONS[] = (; A, B, G, seeds, seed_counts)
end

function exact_real(x::Real, name)
    isfinite(x) || throw(ArgumentError("$name must be finite"))
    QQ(x)
end
function exact_real(x::BigFloat, name)
    isfinite(x) || throw(ArgumentError("$name must be finite"))
    setprecision(BigFloat, max(precision(BigFloat), precision(x))) do
        QQ(x)
    end
end
function downward_convert(::Type{T}, bound::Rational) where {T<:AbstractFloat}
    x = T(bound)
    while exact_real(x, "bound") > bound
        x = prevfloat(x)
    end
    x
end
function numerical(E, name, eta, ::Type{T}) where {T<:AbstractFloat}
    V = reduce(vcat, permutedims.([floating.(T, v) ./ sqrt(floating(T, N)) for (v, N) in zip(E.rays, E.norms)]))
    MeasurementSet(name, V, copy(E.supports), floating.(T, E.weights), deepcopy(E.group), downward_convert(T, eta), E)
end
"All 2306 positive circuits of the 32-ray mesh, in the saved active-set order."
measurement_set(; T::Type{<:AbstractFloat} = Float64) = numerical(exact_constructions().A, :alice, ETA_A, T)
"The seven seed orbits (121 POVMs), in the saved active-set order."
template_set(; T::Type{<:AbstractFloat} = Float64) = numerical(exact_constructions().B, :bob, ETA_B, T)

include("symmetry.jl")
include("shrinking.jl")
include("certification.jl")
scenario(; T::Type{<:AbstractFloat} = Float64, kwargs...) = build_scenario(measurement_set(; T), template_set(; T); kwargs...)
end
