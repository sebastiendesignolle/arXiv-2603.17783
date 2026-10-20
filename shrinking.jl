# Independent coverage certification. The witness file supplies trial convex
# decompositions only. Neither its claimed geometry nor floating LP tolerances
# are trusted. Alice uses Field and rational intervals; the lazily loaded Bob
# verifier uses rigorous real balls and exact algebraic numbers.

"Exact squared inradius eta_16² from DIB+23, in the real embedding of Field."
alice_inradius_squared() = (620+185S5+(505+205S5)*RT)/2245

const ALICE_SHRINKING_EXPRESSION =
    "sqrt((620+185*sqrt(5)+sqrt(30*(12905+5701*sqrt(5))))/2245)"

"""
    exact_alice_visibility(cert)

Retain Alice's exact algebraic shrinking factor in a finite-model certificate
returned by active_set_certificate for this construction.
The certified visibility is the NONNEGATIVE square root of `result.squared`,
an exact Field element in Q(sqrt(5),sqrt((15-3sqrt(5))/10)).
Equivalently it is `result.prefactor * sqrt(result.alice_squared)`.
`expression` gives this exact radical formula; `lower` and `upper` are rigorous
rational enclosures of its value, not bounds bracketing the optimal threshold.
No floating-point approximation or Oscar initialization enters this conversion.

The coverage proof must be checked with certify_shrinking().
The dyadic residual already encloses the ideal algebraic measurement behavior,
so substituting the exact inradius requires no new residual or radius estimate.
The rational-only `cert.visibility_lower` remains available as a weaker bound.
"""
function exact_alice_visibility(cert)
    cert.eta==(ETA_A,ETA_B) ||
        throw(ArgumentError("expected shrinking metadata for the Alice/Bob construction"))
    prefactor=cert.eta[2]*cert.finite_visibility
    prefactor>=0 || throw(ArgumentError("finite visibility must be nonnegative"))
    alice_squared=alice_inradius_squared()
    squared=prefactor^2*alice_squared
    enclosure=iszero(prefactor) ? Interval(0) : sqrt_interval(Interval(squared))
    expression="($(numerator(prefactor))/$(denominator(prefactor)))*"*ALICE_SHRINKING_EXPRESSION
    (;prefactor,alice_squared,squared,lower=enclosure.lo,upper=enclosure.hi,expression)
end

function mesh_inradius(V,actions)
    @assert all(any(u->u==-v,V) for v in V)
    @assert exact_solve(hcat(V[1:3]...),[1,0,0])!==nothing
    seen=Set{UInt64}();lower=rat(1);faces=Set{Tuple}()
    for ids in combinations(1:length(V),3)
        support_mask(ids) in seen && continue
        for g in actions;push!(seen,support_mask(g[ids]));end
        a,b,c=V[ids];n=ecross(b-a,c-a);all(iszero,n) && continue
        d=edot(n,a);delta=[edot(n,v)-d for v in V]
        any(positive,delta)&&any(x->positive(-x),delta) && continue
        @assert !iszero(d)
        face=findall(iszero,delta)
        for g in actions;push!(faces,Tuple(sort(g[face])));end
        r2=d*d/(3edot(n,n));@assert positive(r2)
        # Every facet has this distance. This connects the exhaustive geometry
        # check to the analytical expression, rather than recognizing decimals.
        @assert r2==alice_inradius_squared()
        @assert 1-40r2+530r2^2-2480r2^3+2245r2^4==Field(0)
        lower=min(lower,Interval(r2).lo)
    end
    @assert length(faces)==60
    sqrt_interval(Interval(lower)).lo
end

"""
    certify_shrinking(; witnesses=joinpath(@__DIR__,"centers_2_witnesses.tsv"), verbose=true)

Prove Alice's exact inradius and eta_B >= 0.860491571453.
Alice: exhaust the supporting facets of the 32-ray hull and identify every
squared facet distance with alice_inradius_squared() in the exact field.
The `eta` tuple retains rational lower bounds for rational-only consumers;
`alice_squared` specifies the exact shrinking factor by its positive square root.
Bob: verify the 92 unit normals, 180 polar vertices, complete catalogue of
5,060,220 outer generators, and all 42,547 simulation-orbit witnesses.
The Bob verifier uses Oscar's rigorous real balls and exact algebraic numbers;
Oscar is loaded only when this function is called. No numerical feasibility
tolerance, Python, or LP solver is used. See README.md for the proof and formats.
"""
function certify_shrinking(;witnesses=joinpath(@__DIR__,"centers_2_witnesses.tsv"),verbose=true)
    c=exact_constructions()
    verbose && println("Certifying the 32-ray hull inradius...")
    etaA=mesh_inradius(c.A.rays,c.A.group)
    @assert etaA>=ETA_A
    if !isdefined(@__MODULE__,:OuterCover92)
        Base.include(@__MODULE__,joinpath(@__DIR__,"outer_cover.jl"))
    end
    # Both binding lookup and invocation must occur in the latest world:
    # Julia 1.12+ also tracks the world age of newly introduced module bindings.
    bob=Base.invokelatest() do
        getfield(@__MODULE__,:OuterCover92).certify(;witnesses,verbose)
    end
    etaB=rat(bob["eta_rational"])
    @assert etaB>=ETA_B
    verbose && println("Verified 2306 Alice measurements, 121 Bob measurements, and 42547 outer-generator certificates.")
    (;eta=(ETA_A,ETA_B),alice_squared=alice_inradius_squared(),
      computed_bounds=(etaA,etaB),outer_generators=5060220,
      outer_orbits=42547,basis_certificates=bob["basis_certificates"],
      repair_certificates=bob["repair_certificates"],
      witness_sha256=bytes2hex(sha256(read(witnesses))))
end

"""
    certified_povm_data(A, B; coordinate_bits=40)

Enclose ideal effects w*I+z⋅sigma directly in Julia. Numeric T is not a proof
input: the exact construction attached to each set supplies every coefficient.
The returned rational coverage bounds must be checked with certify_shrinking
before claiming the final all-POVM theorem. No full probability tensor needed.
"""
function certified_povm_data(A::MeasurementSet,B::MeasurementSet;coordinate_bits::Int=40)
    8<=coordinate_bits<=50 || throw(ArgumentError("coordinate_bits must lie in 8:50"))
    D=Int64(1)<<coordinate_bits
    sides=map((A,B)) do M
        E=M.exact
        M.supports==E.supports && M.group==E.group || error("measurement order/group was modified")
        lo=zeros(Int64,4,4,length(M));hi=similar(lo)
        scales=sqrt_interval.(Interval.(E.norms))
        for x in 1:length(M),a in 1:4
            i=E.supports[a,x];w=E.weights[a,x]
            coefficients=i==0 ? fill(Interval(0),4) :
                [Interval(w);[Interval(w*z)/scales[i] for z in E.rays[i]]]
            for c in 1:4
                lo[c,a,x]=floor(BigInt,coefficients[c].lo*D)
                hi[c,a,x]=ceil(BigInt,coefficients[c].hi*D)
            end
        end
        (;lo,hi,eta=M.name==:alice ? ETA_A : M.name==:bob ? ETA_B : error("unknown measurement set"))
    end
    (;lo=Tuple(s.lo for s in sides),hi=Tuple(s.hi for s in sides),
      eta=Tuple(s.eta for s in sides),bits=coordinate_bits,denominator=D)
end
