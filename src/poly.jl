abstract type AbstractUnivariateMonomial <: AbstractMonomial end
struct Uninomial{V} <: AbstractUnivariateMonomial
    v::V
    d::Int
end
Uninomial{T}() where T = Uninomial()
const DO_NOT_CHECK_VAR = typemax(Int)
Uninomial() = Uninomial(DO_NOT_CHECK_VAR, 0)
isunivariate(x::Uninomial) = true

#Base.promote_rule(u::Type{Uninomial{T}}, ::Type{Uninomial{Nothing}}) where {T} = u

degree(m::Uninomial) = m.d
Base.isless(x::Uninomial, y::Uninomial) = isless(degree(x), degree(y))
function Base.:(==)(x::Uninomial, y::Uninomial)
    mx = metadata(x)
    my = metadata(y)
    if mx === nothing || mx === DO_NOT_CHECK_VAR || my === nothing || my === DO_NOT_CHECK_VAR
        return degree(x) == degree(y)
    else
        return mx == my && degree(x) == degree(y)
    end
end
var(m::Uninomial) = m.v

#coeff(m::Uninomial) = 1
Base.one(::Type{<:Uninomial}) = Uninomial(nothing, 0)
Base.show(io::IO, m::Uninomial) = print_single_monomial(io, var(m), degree(m))
dropmetadata(t::Type{<:Uninomial}) = Uninomial{Nothing}
dropmetadata(t::Uninomial) = Uninomial(nothing, degree(t))
addmetadata(t::Uninomial, meta) = Uninomial(meta, degree(t))
Base.convert(::Type{<:Uninomial{Nothing}}, m::Uninomial) = dropmetadata(m)
Base.convert(::Type{<:Term{T,<:Uninomial{Nothing}}}, m::Term{T,<:Uninomial}) where T = dropmetadata(m)
Base.promote_rule(t::Type{<:Uninomial{Nothing}}, ::Type{<:Uninomial}) = t

function Base.:(*)(x::Uninomial, y::Uninomial)
    v = checkmetadata(x, y)
    Uninomial(v, degree(x) + degree(y))
end

#struct Uniterm{T,M<:AbstractMonomial} <: AbstractTerm
#    coeff::T
#    uninomial::M
#end
#Base.convert(::Type{<:Uniterm}, t::Uninomial) = Uniterm(t)
#Base.convert(::Type{<:Uniterm}, t::CoeffType) = Uniterm(t, Uninomial())
#Uniterm(t::Uninomial) = Uniterm(coeff(t), t)

#coeff(t::Uniterm) = t.coeff
#monomial(t::Uniterm) = t.uninomial

#Base.zero(t::Uniterm) = Uniterm(zero(coeff(t)), one(monomial(t)))
#Base.one(t::Uniterm) = Uniterm(one(coeff(t)), one(monomial(t)))

#=
struct SPoly{T} <: AbstractPolynomial
    coeffs::Vector{T}
    exps::Vector{UInt}
    v::IDType

    SPoly(coeffs::Vector{T}, exps::Vector{UInt}, v::IDType) where T = SPoly{T}(coeffs, exps, v)
    function SPoly{T}(coeffs::Vector{T}, exps::Vector{UInt}, v::IDType) where T
        length(coeffs) == length(exps) || error("coeffs and exps' length must match!")
        new{T}(coeffs, exps, v)
    end
end
Base.convert(::Type{<:SPoly}, m::Uninomial) = SPoly(Uniterm(m))
Base.convert(::Type{<:SPoly}, t::Uniterm) = SPoly(t)
#Base.convert(::Type{<:SPoly}, t::CoeffType) = SPoly(convert(Uniterm, t))
SPoly(t::Uniterm) = SPoly([coeff(t)], [degree(t)], var(t))
SPoly(c::Number, v) = SPoly([c], [zero(UInt)], v)
const EMPTY_EXPS = UInt[]
SPoly(t::Uniterm, id::IDType) = SPoly(typeof(coeff(t))[], EMPTY_EXPS, id)
Base.similar(p::SPoly) = SPoly(similar(coeffs(p)), similar(p.exps), var(p))

var(p::SPoly) = p.v
coeffs(p::SPoly) = p.coeffs
coeff(p::SPoly, i) = coeffs(p)[i]
term(p::SPoly, i) = Uniterm(p.coeffs[i], Uninomial(var(p), p.exps[i]))
terms(p::SPoly) = (term(p, i) for i in eachindex(coeffs(p)))
lt(p::SPoly) = term(p, 1)
lc(p::SPoly) = coeff(p, 1)
degree(p::SPoly) = iszero(p) ? -1 : p.exps[1]
Base.iszero(p::SPoly) = isempty(p.exps) || iszero(lt(p))

Base.zero(t::SPoly) = zero(lt(t))
Base.one(t::SPoly) = one(lt(t))

Base.deleteat!(p::SPoly, i::Int) = (deleteat!(p.coeffs, i); deleteat!(p.exps, i); nothing)
Base.copy(p::SPoly) = SPoly(copy(coeffs(p)), copy(p.exps), var(p))
Base.copy(p::SPoly{<:AbstractPolynomial}) = SPoly(map(copy, coeffs(p)), copy(p.exps), var(p))
=#

struct SPoly{T,M} <: AbstractUnivariatePolynomial{T,M}
    terms::T
    metadata::M
end
Base.similar(p::SPoly) = SPoly(similar(terms(p)), metadata(p))
function Base.:(==)(p::SPoly, q::SPoly)
    p === q && return true
    terms(p) == terms(q)
    #TODO: check symbols?
end

function checkmetadata(a::Union{SPoly,Uninomial}, b::Union{SPoly,Uninomial})
    x, y = metadata(a), metadata(b)
    if x === nothing || x === DO_NOT_CHECK_VAR
        return y
    elseif y === nothing || y === DO_NOT_CHECK_VAR
        return x
    else
        x == y || error("$a and $b contain different variables!")
        return x
    end
end

function print_poly_term(io, t)
    need_par = !isone(monomial(t))
    need_par && print(io, '(')
    show(io, coeff(t))
    need_par && print(io, ')')
    need_par && show(io, monomial(t))
    nothing
end

function Base.show(io::IO, p::SPoly{<:AbstractPolynomial})
    ts = terms(p)
    if isempty(ts)
        print(io, 0)
        return
    end
    n = length(ts)
    t1, tr = Iterators.peel(ts)
    print_poly_term(io, t1)
    for t in tr
        print(io, " + ")
        print_poly_term(io, t)
    end
end

#########################
# Arithmetic/Algorithms #
#########################

Base.:(^)(x::Uninomial, n::Integer) = Uninomial(var(x), degree(x)+n-1)

#=
function Base.:(+)(x::Uninomial, y::Uninomial)
    v = checkmetadata(x, y)
    if x < y
        x, y = y, x
    end
    SPoly([], v)
end
function Base.:(+)(x::Uniterm, y::Uniterm)
    v = checkmetadata(x, y)
    if ismatch(x, y)
        c = x.coeff + y.coeff
        return iszero(c) ? SPoly(x, var(x)) : SPoly(Uniterm(c, monomial(x)))
    else
        if x < y
            x, y = y, x
        end
        return SPoly([coeff(x), coeff(y)], [degree(x), degree(y)], var(x))
    end
end
Base.:(*)(x::CoeffType, y::Uninomial) = Uniterm(x, y)
Base.:(*)(x::Uninomial, y::CoeffType) = y * x

Base.:(*)(x::Uniterm, y::Uniterm) = Uniterm(coeff(x) * coeff(y), monomial(x) * monomial(y))

function Base.:*(p::SPoly, x::Uniterm)
    check_poly(p, x)
    if iszero(x)
        return SPoly(x)
    elseif isone(x)
        return p
    else
        cfs = similar(coeffs(p))
        exps = similar(p.exps)
        for (i, t) in enumerate(terms(p))
            nt = t * x
            cfs[i] = coeff(nt)
            exps[i] = degree(nt)
        end
        return SPoly(cfs, exps, var(x))
    end
end
Base.:*(x::CoeffType, y::Uniterm) = Uniterm(x * coeff(y), monomial(y))
Base.:*(x::Uniterm, y::CoeffType) = y * x

Base.:*(x::SPoly, y::SPoly) = sum(t->x * t, terms(y))
smul(x, y::SPoly) = SPoly(x * coeffs(y), y.exps, var(y))
Base.:*(x::CoeffType, y::SPoly) = smul(x, y)
Base.:*(x::SPoly, y::CoeffType) = y * x
Base.:*(x::T, y::SPoly{T}) where {T<:AbstractPolynomial} = smul(x, y)
Base.:*(x::SPoly{T}, y::T) where {T<:AbstractPolynomial} = y * x

addcoef(x::Uniterm, c) = (c += coeff(x); return iszero(c), Uniterm(c, monomial(x)))
addcoef(x::Uniterm, c::Uniterm) = addcoef(x, coeff(c))
Base.:(-)(x::Uniterm) = Uniterm(-coeff(x), monomial(x))
sub!(p::SPoly, x::Uniterm) = add!(p, -x)
function add!(p::SPoly, x::Uniterm)
    iszero(x) && return p
    for (i, t) in enumerate(terms(p))
        if ismatch(t, x)
            iz, t = addcoef(t, x)
            if iz
                deleteat!(p, i)
            else
                p.coeffs[i] = coeff(t)
            end
            return p
        elseif t < x
            insert!(p.coeffs, i, coeff(x))
            insert!(p.exps, i, degree(x))
            return p
        end
    end
    push!(p.coeffs, coeff(x))
    push!(p.exps, degree(x))
    return p
end

Base.div(x::Uninomial, y::Uninomial) = (x/y)[1]
function Base.:(/)(x::Uninomial, y::Uninomial)
    check_poly(x, y)
    xd, yd = degree(x), degree(y)
    xd >= yd ? (Uninomial(xd - yd, var(x)), false) : (x, true)
end
=#

# dest = (p * a).^n
function mulpow!(dest, p::SPoly, a, n::Integer)
    @assert n >= 0
    ts = terms(dest)
    pts = terms(p)
    resize!(ts, length(pts))
    for (i, t) in enumerate(pts)
        m = monomial(t)
        M = typeof(m)
        c, d = coeff(t), degree(m)
        ts[i] = Term(c * a, M(metadata(m), d + n))
        #dest.coeffs[i] = p.coeffs[i] * a
        #dest.exps[i] = p.exps[i] + n
    end
    return dest
end

function fnmadd!(p::SPoly, l, s)
    for li in terms(l)
        sub!(p, li * s)
    end
end
function mul!(p::SPoly, l)
    #cs = coeffs(p)
    ts = terms(p)
    for i in eachindex(ts)
        t = ts[i]
        ts[i] = Term(coeff(t) * l, monomial(t)) #cs[i] *= l
    end
end
#=function _copyto!(x::Vector, y::Vector)
    resize!(x, length(y))
    copyto!(x, y)
    nothing
end
function Base.copyto!(x::Vector{MPoly}, y::Vector{MPoly})
    tsx = terms(x)
    tsy = terms(y)
    resize!(tsx, length(tsy))
    for i in eachindex(tsy)
        tsx[i] = copy(tsy[i])
    end
    x
end
function Base.copyto!(x::MPoly, y::MPoly)
    _copyto!(terms(x), terms(y))
end
function Base.copyto!(x::SPoly, y::SPoly)
    _copyto!(coeffs(x), coeffs(y))
    _copyto!(x.exps, y.exps)
end =#
function pseudorem(p::SPoly, d::SPoly)
    checkmetadata(p, d)
    degree(p) < degree(d) && return p
    k = degree(p) - degree(d) + 1
    l = lc(d)
    dd = similar(d)
    # pp = similar(p)
    while !iszero(p) && degree(p) >= degree(d)
        s = mulpow!(dd, d, lc(p), degree(p) - degree(d))
        p = copy(p)
        mul!(p, l)
        sub!(p, s)
        k -= 1
    end
    return l^k * p
end

function content(a::SPoly)
    ts = terms(a)
    length(ts) == 1 && return lc(a)
    g = gcd(coeff(ts[1]), coeff(ts[2]))
    # TODO: short circuit and split to small terms
    for i in 3:length(ts)
        g = gcd(g, coeff(ts[i]))
    end
    g
end

function univariate_gcd(x::AbstractPolynomial, y::AbstractPolynomial)
    while !iszero(y)
        x = pseudorem(x, y)
        x, y = y, x
    end
    return x#, a / x, b / x
end

#function divexact(a::SPoly{T}, b::T) where {T<:AbstractPolynomial}
function divexact(a::SPoly, b)
    SPoly([Term(divexact(coeff(t), b), monomial(t)) for t in terms(a)], metadata(a))
end

divexact(c, b) = ((d, r) = divrem(c, b); @assert iszero(r); d;)
#function divexact(a::SPoly, b)
#    SPoly([Term(divexact(coeff(t), b), monomial(t)) for t in terms(a)], metadata(a))
#    #cfs = [divexact(c, b) for c in coeffs(a)]
#    #SPoly(cfs, a.exps, var(a))
#end

function Base.gcd(x::SPoly, y::SPoly)
    if degree(y) > degree(x)
        x, y = y, x
    end
    if iszero(y)
        return x
    elseif isone(y)
        return y
    end
    c1, x = contprim(x)
    c2, y = contprim(y)
    c = gcd(c1, c2)
    CT = monomialtype(eltype(terms(x)))
    g = one(CT)
    h = one(CT)
    while true
        r = pseudorem(x, y)
        iszero(r) && break
        degree(r) == 0 && return SPoly(c, var(x))

        d = degree(x) - degree(y)
        x, y = y, divexact(r, g*h^d)
        g = lc(x)
        h = d > 1 ? divexact(g^d, h^(d - 1)) : h^(1 - d)*g^d
    end
    return c*primpart(y)
end

primpart(p) = divexact(p, content(p))
contprim(p) = (c = content(p); (c, divexact(p, c)))

##############
# Conversion #
##############

# mpoly2poly(5x⁴ + 2x³ + x² + xy + y^2, 0x00001 [y])
# ->
#              ---------------- constants
#              |
# y^2 + xy + (5x⁴ + 2x³ + x²)
#
# M2P(a, v)
function emplace_back!(poly, ts, perm, chunk_start_idx, idx, olddegree)
    v = metadata(poly)
    if olddegree > 0
        t = ts[perm[chunk_start_idx]]
        coeff = MPoly(term2polycoeff(t, v))
        for i in chunk_start_idx+1:idx-1
            t = ts[perm[i]]
            add!(coeff, term2polycoeff(t, v))
        end
    else
        coeff = MPoly(ts[perm[chunk_start_idx]])
        for i in chunk_start_idx+1:idx-1
            add!(coeff, ts[perm[i]])
        end
    end
    push!(poly.terms, Term(coeff, Uninomial(nothing, olddegree)))
    #push!(poly.exps, olddegree)
    #push!(poly.coeffs, coeff)
    nothing
end
term2polycoeff(t::Term{<:Any,<:Monomial}, v::IDType) = Term(coeff(t), Monomial(filter(!isequal(v), monomial(t).ids)))
function SPoly(p::MPoly, v::IDType)
    ts = terms(p)
    pows = map(ts) do t
        count(isequal(v), monomial(t).ids)
    end
    perm = sortperm(pows, rev=true)
    sterms = Term{typeof(p),Uninomial{Nothing}}[]
    #exps = UInt[]
    poly = SPoly(sterms, v)

    olddegree = pows[perm[1]]
    chunk_start_idx = idx = 1
    while idx <= length(perm)
        degree = pows[perm[idx]]
        if olddegree != degree # new chunk
            emplace_back!(poly, ts, perm, chunk_start_idx, idx, olddegree)
            chunk_start_idx = idx
            olddegree = degree
        end
        idx += 1
    end
    if !isempty(chunk_start_idx:idx-1) # remember to handle the last one
        emplace_back!(poly, ts, perm, chunk_start_idx, idx, olddegree)
    end
    return poly
end

function univariate_to_multivariate(g::SPoly)
    #cfs = coeffs(g)
    #eps = g.exps
    ts = terms(g)
    v = var(g)
    @assert !isempty(ts)
    #sum(zip(cfs, eps)) do (c, e)
    #    c * Monomial(fill(v, e))
    #end
    s = coeff(ts[1]) * Monomial(fill(v, degree(ts[1])))
    for i in 2:length(ts)
        t = ts[i]
        c = coeff(t)
        e = degree(t)
        add!(s, c * Monomial(fill(v, e)))
    end
    return s
end
