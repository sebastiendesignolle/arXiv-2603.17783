# Internal exact arithmetic. The real embedding is s=sqrt(5)>0 and
# t=sqrt((15-3s)/10)>0; Alice's raw rays all have squared norm 3.
# No floating feasibility test enters any certificate.
const QQ = Rational{BigInt}
rat(x::Integer) = BigInt(x)//big(1)
rat(x::Rational) = QQ(x)
function rat(x::AbstractString)
    p=split(x,'/'); length(p)<=2 || error("invalid rational")
    parse(BigInt,p[1]) // (length(p)==1 ? big(1) : parse(BigInt,p[2]))
end
struct Field <: Number
    c::NTuple{4,QQ}
end
Field(x::Union{Integer,Rational}) = Field((rat(x),rat(0),rat(0),rat(0)))
Field(x::Field)=x
Base.zero(::Type{Field})=Field(0)
Base.one(::Type{Field})=Field(1)
Base.zero(::Field)=Field(0)
Base.one(::Field)=Field(1)
Base.conj(x::Field)=x
Base.transpose(x::Field)=x
Base.adjoint(x::Field)=x
Base.real(x::Field)=x
Base.abs2(x::Field)=x*x
Base.iszero(x::Field)=all(iszero,x.c)
Base.:(==)(x::Field,y::Field)=x.c==y.c
Base.hash(x::Field,h::UInt)=hash(x.c,h)
Base.convert(::Type{Field},x::Union{Integer,Rational})=Field(x)
Base.promote_rule(::Type{Field},::Type{T}) where {T<:Union{Integer,Rational}}=Field
Base.:+(x::Field,y::Field)=Field(x.c.+y.c)
Base.:-(x::Field)=Field(.-x.c)
Base.:-(x::Field,y::Field)=x+-y
base_mul(a,b)=(a[1]*b[1]+5a[2]*b[2],a[1]*b[2]+a[2]*b[1])
function Base.:*(x::Field,y::Field)
    a,b,c,d=x.c;e,f,g,h=y.c
    uu=base_mul((a,b),(e,f));vv=base_mul(base_mul((c,d),(g,h)),(3//2,-3//10))
    uv=base_mul((a,b),(g,h));vu=base_mul((c,d),(e,f))
    Field((uu[1]+vv[1],uu[2]+vv[2],uv[1]+vu[1],uv[2]+vu[2]))
end
function Base.inv(x::Field)
    iszero(x) && throw(DivideError())
    a,b,c,d=x.c
    if iszero(b)&&iszero(c)&&iszero(d);return Field(inv(a));end
    uu=base_mul((a,b),(a,b));vv=base_mul(base_mul((c,d),(c,d)),(3//2,-3//10))
    e,f=uu.-vv;den=e*e-5f*f
    u=base_mul((a,b),(e/den,-f/den));v=base_mul((-c,-d),(e/den,-f/den))
    Field((u...,v...))
end
Base.:/(x::Field,y::Field)=x*inv(y)

"Closed rational interval; every operation is outward and independent of MPFR."
struct Interval
    lo::QQ
    hi::QQ
    function Interval(a::QQ,b::QQ)
        a<=b || error("reversed interval");new(a,b)
    end
end
Interval(x::Union{Integer,Rational})=Interval(rat(x),rat(x))
Base.:+(a::Interval,b::Interval)=Interval(a.lo+b.lo,a.hi+b.hi)
Base.:-(a::Interval)=Interval(-a.hi,-a.lo)
Base.:-(a::Interval,b::Interval)=a+-b
function Base.:*(a::Interval,b::Interval)
    z=(a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi);Interval(minimum(z),maximum(z))
end
function Base.:/(a::Interval,b::Interval)
    b.lo*b.hi>0 || error("interval divisor contains zero")
    a*Interval(inv(b.hi),inv(b.lo))
end
abs_upper(a::Interval)=max(abs(a.lo),abs(a.hi))
function sqrt_interval(a::Interval;bits=120)
    a.lo>=0 || error("negative square root")
    D=big(1)<<bits
    lower=isqrt(fld(numerator(a.lo)*D^2,denominator(a.lo)))
    upper=isqrt(fld(numerator(a.hi)*D^2,denominator(a.hi)))+1
    Interval(lower//D,upper//D)
end
const SI=sqrt_interval(Interval(5))
const TI=sqrt_interval((Interval(15)-Interval(3)*SI)/Interval(10))
const STI=SI*TI
function Interval(x::Field)
    z=Interval(x.c[1])
    for (c,b) in zip(x.c[2:4],(SI,TI,STI));z=z+Interval(c)*b;end
    z
end
function positive(x::Field)
    iszero(x) && return false
    z=Interval(x)
    z.lo<=0<=z.hi && error("increase algebraic root enclosure precision")
    z.lo>0
end
Base.isless(x::Field,y::Field)=positive(y-x)
const S5=Field((rat(0),rat(1),rat(0),rat(0)))
const RT=Field((rat(0),rat(0),rat(1),rat(0)))
const PHI=(1+S5)/2

function floating(::Type{T},x::Field) where {T<:AbstractFloat}
    a,b,c,d=T.(x.c);s=sqrt(T(5));t=sqrt((15-3s)/10)
    a+b*s+(c+d*s)*t
end
edot(a,b)=sum(x*y for (x,y) in zip(a,b);init=Field(0))
ecross(a,b)=[a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1]]
matvec(A,v)=[edot(row,v) for row in eachrow(A)]

"Exact full-column-rank solve, including consistency of overdetermined rows."
function exact_solve(A,b)
    m,n=size(A);m>=n || return nothing
    B=hcat(Field.(A),Field.(b))
    for j in 1:n
        pivot=findfirst(i->!iszero(B[i,j]),j:m)
        pivot===nothing && return nothing
        i=j+pivot-1;B[j,:],B[i,:]=copy(B[i,:]),copy(B[j,:])
        a=B[j,j];B[j,:].=B[j,:]./a
        for i in 1:m
            i==j && continue
            a=B[i,j];iszero(a) || (B[i,:].=B[i,:].-a.*B[j,:])
        end
    end
    any(!iszero,B[n+1:m,end]) && return nothing
    B[1:n,end]
end
function exact_weights(V,ids)
    A=vcat(fill(Field(1),1,length(ids)),hcat(V[ids]...))
    w=exact_solve(A,[1,0,0,0])
    w!==nothing && all(positive,w) ? w : nothing
end
support_mask(ids)=foldl((a,i)->a|(UInt64(1)<<(i-1)),ids;init=UInt64(0))

"Exhaust all supports of size 2--4 modulo a verified orthogonal group."
function positive_catalogue(V,group)
    length(V)<=64 || error("support mask supports at most 64 rays")
    seen=Set{UInt64}();valid=Dict{Tuple,Vector{Field}}();representatives=Tuple[]
    for k in 2:4, ids in combinations(1:length(V),k)
        support_mask(ids) in seen && continue
        for g in group;push!(seen,support_mask(g[ids]));end
        w=exact_weights(V,ids);w===nothing && continue
        push!(representatives,Tuple(ids))
        for g in group
            p=sortperm(g[ids]);key=Tuple(g[ids][p]);value=w[p]
            haskey(valid,key) && valid[key]!=value && error("inconsistent group orbit")
            valid[key]=value
        end
    end
    (;valid,representatives)
end
