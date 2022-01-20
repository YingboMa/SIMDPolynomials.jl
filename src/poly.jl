struct Uninomial <: AbstractMonomial
    v::IDType
    d::UInt
end

Base.isless(x::Uninomial, y::Uninomial) = isless(degree(x), degree(y))
degree(m::Uninomial) = m.d
var(m::Uninomial) = m.v

coeff(m::Uninomial) = 1
Base.one(m::Uninomial) = Uninomial(m.v, 0)
Base.show(io::IO, m::Uninomial) = print_single_monomial(io, var(m), degree(m))

struct Uniterm{T,M<:AbstractMonomial} <: AbstractTerm
    coeff::T
    uninomial::M
end
Base.convert(::Type{<:Uniterm}, t::Uninomial) = Uniterm(t)
#Base.convert(::Type{<:Uniterm}, t::CoeffType) = Uniterm(t, Uninomial())
Uniterm(t::Uninomial) = Uniterm(coeff(t), t)

coeff(t::Uniterm) = t.coeff
monomial(t::Uniterm) = t.uninomial
var(t::Uniterm) = var(monomial(t))

struct SparsePoly{T} <: AbstractPolynomial
    coeffs::Vector{T}
    exps::Vector{UInt}
    v::IDType

    SparsePoly(coeffs::Vector{T}, exps::Vector{UInt}, v::IDType) where T = SparsePoly{T}(coeffs, exps, v)
    function SparsePoly{T}(coeffs::Vector{T}, exps::Vector{UInt}, v::IDType) where T
        length(coeffs) == length(exps) || error("coeffs and exps' length must match!")
        new{T}(coeffs, exps, v)
    end
end
Base.convert(::Type{<:SparsePoly}, m::Uninomial) = SparsePoly(Uniterm(m))
Base.convert(::Type{<:SparsePoly}, t::Uniterm) = SparsePoly(t)
#Base.convert(::Type{<:SparsePoly}, t::CoeffType) = SparsePoly(convert(Uniterm, t))
SparsePoly(t::Uniterm) = SparsePoly([coeff(t)], UInt[degree(t)], var(t))
SparsePoly(c::Number, v) = SparsePoly([c], [zero(UInt)], v)
const EMPTY_EXPS = UInt[]
#SparsePoly(t::Uniterm, id::IDType) = SparsePoly(typeof(coeff(t))[], EMPTY_EXPS, id)
Base.similar(p::SparsePoly) = SparsePoly(similar(coeffs(p)), similar(p.exps), var(p))

var(p::SparsePoly) = p.v
coeffs(p::SparsePoly) = p.coeffs
coeff(p::SparsePoly, i) = coeffs(p)[i]
term(p::SparsePoly, i) = Uniterm(p.coeffs[i], Uninomial(var(p), p.exps[i]))
terms(p::SparsePoly) = (term(p, i) for i in eachindex(coeffs(p)))
lt(p::SparsePoly) = term(p, 1)
lc(p::SparsePoly) = coeff(p, 1)
degree(p::SparsePoly) = iszero(p) ? -1 : p.exps[1]
Base.iszero(p::SparsePoly) = isempty(p.exps) || iszero(lt(p))

Base.zero(t::SparsePoly) = zero(lt(t))
Base.one(t::SparsePoly) = one(lt(t))

Base.deleteat!(p::SparsePoly, i::Int) = (deleteat!(p.coeffs, i); deleteat!(p.exps, i); nothing)
Base.copy(p::SparsePoly) = SparsePoly(copy(coeffs(p)), copy(p.exps), var(p))
Base.copy(p::SparsePoly{<:AbstractPolynomial}) = SparsePoly(map(copy, coeffs(p)), copy(p.exps), var(p))
Base.:(==)(p::SparsePoly, q::SparsePoly) = all(x->x[1] == x[2], zip(terms(p), terms(q)))

function check_poly(x, y)
    v = var(x)
    v == var(y) || error("$x and $y contain different variables!")
    nothing
end

function print_poly_term(io, t)
    need_par = !isone(monomial(t))
    need_par && print(io, '(')
    show(io, coeff(t))
    need_par && print(io, ')')
    need_par && show(io, monomial(t))
    nothing
end

function Base.show(io::IO, p::SparsePoly{<:AbstractPolynomial})
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

function Base.:(+)(x::Uninomial, y::Uninomial)
    check_poly(x, y)
    if x < y
        x, y = y, x
    end
    SparsePoly([1, 1], [degree(x), degree(y)], var(x))
end
function Base.:(+)(x::Uniterm, y::Uniterm)
    check_poly(x, y)
    if ismatch(x, y)
        c = x.coeff + y.coeff
        return iszero(c) ? SparsePoly(typeof(x)[], EMPTY_EXPS, var(x)) : SparsePoly(Uniterm(c, monomial(x)))
    else
        if iszero(x) && iszero(y)
            c = coeff(y)
            return SparsePoly(typeof(c)[], EMPTY_EXPS, var(x))
        elseif iszero(x)
            exp = degree(y)
            c = coeff(y)
            return SparsePoly([c], UInt[exp], var(x))
        elseif iszero(y)
            exp = degree(x)
            c = coeff(x)
            return SparsePoly([c], UInt[exp], var(x))
        else
            if x < y
                x, y = y, x
            end
            return SparsePoly([coeff(x), coeff(y)], UInt[degree(x), degree(y)], var(x))
        end
    end
end
Base.:(*)(x::CoeffType, y::Uninomial) = Uniterm(x, y)
Base.:(*)(x::Uninomial, y::CoeffType) = y * x

function Base.:(*)(x::Uninomial, y::Uninomial)
    check_poly(x, y)
    Uninomial(var(x), degree(x) + degree(y))
end
Base.:(*)(x::Uniterm, y::Uniterm) = Uniterm(coeff(x) * coeff(y), monomial(x) * monomial(y))

function Base.:*(p::SparsePoly, x::Uniterm)
    check_poly(p, x)
    if iszero(x)
        return SparsePoly(x)
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
        return SparsePoly(cfs, exps, var(x))
    end
end
Base.:*(x::CoeffType, y::Uniterm) = Uniterm(x * coeff(y), monomial(y))
Base.:*(x::Uniterm, y::CoeffType) = y * x

Base.:*(x::SparsePoly, y::SparsePoly) = sum(t->x * t, terms(y))
smul(x, y::SparsePoly) = SparsePoly(x * coeffs(y), y.exps, var(y))
Base.:*(x::CoeffType, y::SparsePoly) = smul(x, y)
Base.:*(x::SparsePoly, y::CoeffType) = y * x
Base.:*(x::T, y::SparsePoly{T}) where {T<:AbstractPolynomial} = smul(x, y)
Base.:*(x::SparsePoly{T}, y::T) where {T<:AbstractPolynomial} = y * x

addcoef(x::Uniterm, c) = (c += coeff(x); return iszero(c), Uniterm(c, monomial(x)))
addcoef(x::Uniterm, c::Uniterm) = addcoef(x, coeff(c))
Base.:(-)(x::Uniterm) = Uniterm(-coeff(x), monomial(x))
sub!(p::SparsePoly, x::Uniterm) = add!(p, -x)
function add!(p::SparsePoly, x::Uniterm)
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

# dest = (p * a).^n
function mulpow!(dest, p::SparsePoly, a, n::Integer)
    @assert n >= 0
    for i in eachindex(p.exps)
        dest.coeffs[i] = p.coeffs[i] * a
        dest.exps[i] = p.exps[i] + n
    end
    return dest
end

function fnmadd!(p::SparsePoly, l, s)
    for li in terms(l)
        sub!(p, li * s)
    end
end
function mul!(p::SparsePoly, l)
    cs = coeffs(p)
    for i in eachindex(cs)
        cs[i] *= l
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
function Base.copyto!(x::SparsePoly, y::SparsePoly)
    _copyto!(coeffs(x), coeffs(y))
    _copyto!(x.exps, y.exps)
end =#
function pseudorem(p::SparsePoly, d::SparsePoly)
    check_poly(p, d)
    degree(p) < degree(d) && return p
    k = degree(p) - degree(d) + 1
    l = lc(d)
    dd = similar(d)
    # pp = similar(p)
    while !iszero(p) && degree(p) >= degree(d)
        s = mulpow!(dd, d, lc(p), degree(p) - degree(d))
        # copyto!(pp, p);
        # p, pp = pp, p
        p = copy(p)
        mul!(p, l)
        sub!(p, s)
        k -= 1
    end
    return l^k * p
end

function termwise_content(a::MPoly)
    ts = terms(a)
    length(ts) == 1 && return ts[1]
    g = gcd(ts[1], ts[2])
    isone(g) || for i in 3:length(ts)
        g = gcd(g, ts[i])
        isone(g) && break
    end
    g
end

function content(a::SparsePoly)
    cfs = coeffs(a)
    length(cfs) == 1 && return first(cfs)
    # If any coeff here is a multivariate polynomial with only one term, then we
    # just need to compute the termwize content.
    MP = eltype(cfs)
    if MP <: AbstractPolynomial
        for i in eachindex(cfs)
            ts = terms(cfs[i])
            if length(ts) == 1
                g = gcd(termwise_content(cfs[1]), termwise_content(cfs[2]))
                isone(g) || for i in 3:length(cfs)
                    g = gcd(g, termwise_content(cfs[i]))
                    isone(g) && break
                end
                return MP(g)
            end
        end
    end

    g = gcd(cfs[1], cfs[2])
    isone(g) || for i in 3:length(cfs)
        g = gcd(g, cfs[i])
        isone(g) && break
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

function divexact(a::SparsePoly{T}, b::T) where {T<:AbstractPolynomial}
    cfs = [divexact(c, b) for c in coeffs(a)]
    SparsePoly(cfs, a.exps, var(a))
end

divexact(c, b) = ((d, r) = divrem(c, b); @assert iszero(r); d;)
function divexact(a::SparsePoly, b)
    cfs = [divexact(c, b) for c in coeffs(a)]
    SparsePoly(cfs, a.exps, var(a))
end

function Base.gcd(x::SparsePoly, y::SparsePoly)
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
    g = one(eltype(coeffs(x)))
    h = one(eltype(coeffs(x)))
    while true
        r = pseudorem(x, y)
        iszero(r) && break
        degree(r) == 0 && return SparsePoly(c, var(x))

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
    v = var(poly)
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
    push!(poly.exps, olddegree)
    push!(poly.coeffs, coeff)
    nothing
end
term2polycoeff(t::Term, v::IDType) = Term(coeff(t), rmid(monomial(t), v))
function SparsePoly(p::MPoly, v::IDType)
    ts = terms(p)
    pows = map(ts) do t
        degree(monomial(t), v)
    end
    perm = sortperm(pows, rev=true)
    coeffs = typeof(p)[]
    exps = UInt[]
    poly = SparsePoly(coeffs, exps, v)

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

univariate_to_multivariate(g::SparsePoly{<:AbstractPolynomial}) = univariate_to_multivariate(g, monomialtype(lc(g)))
function univariate_to_multivariate(g, ::Type{<:Monomial})
    cfs = coeffs(g)
    eps = g.exps
    v = var(g)
    @assert !isempty(eps)
    #sum(zip(cfs, eps)) do (c, e)
    #    c * Monomial(fill(v, e))
    #end
    s = cfs[1] * Monomial(fill(v, eps[1]))
    for i in 2:length(cfs)
        c = cfs[i]
        e = eps[i]
        add!(s, c * Monomial(fill(v, e)))
    end
    return s
end

function univariate_to_multivariate(g, P::Type{<:PackedMonomial})
    cfs = coeffs(g)
    eps = g.exps
    v = var(g)
    @assert !isempty(eps)
    s = cfs[1] * P(v)^eps[1]
    for i in 2:length(cfs)
        c = cfs[i]
        e = eps[i]
        add!(s, c * P(v)^e)
    end
    return s
end
