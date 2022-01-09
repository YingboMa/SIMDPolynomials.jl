# coeff, length, shift_left

struct Uninomial <: AbstractMonomial
    v::IDType
    d::UInt
end
#Uninomial() = Uninomial(nothing, 0)
Base.isless(x::Uninomial, y::Uninomial) = isless(degree(x), degree(y))
degree(m::Uninomial) = m.d
var(m::Uninomial) = m.v

coeff(m::Uninomial) = 1
Base.isone(m::Uninomial) = iszero(degree(m))
Base.one(m::Uninomial) = Uninomial(m.v, 0)
Base.show(io::IO, m::Uninomial) = print_single_monomial(io, var(m), degree(m))

struct Uniterm{T} <: AbstractTerm
    coeff::T
    uninomial::Uninomial
end
Base.convert(::Type{<:Uniterm}, t::Uninomial) = Uniterm(t)
#Base.convert(::Type{<:Uniterm}, t::Rat) = Uniterm(t, Uninomial())
Uniterm(t::Uninomial) = Uniterm(coeff(t), t)

coeff(t::Uniterm) = t.coeff
monomial(t::Uniterm) = t.uninomial
var(t::Uniterm) = var(monomial(t))

degree(t::Uniterm) = degree(monomial(t))
Base.isless(x::Uniterm, y::Uniterm) = isless(monomial(x), monomial(y))

Base.iszero(t::Uniterm) = iszero(coeff(t))
Base.isone(x::Uniterm) = isone(coeff(x)) && isone(monomial(x))
Base.zero(t::Uniterm) = Uniterm(zero(coeff(t)), one(monomial(t)))
Base.one(t::Uniterm) = Uniterm(one(coeff(t)), one(monomial(t)))

struct SparsePoly{T} <: AbstractPoly
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
#Base.convert(::Type{<:SparsePoly}, t::Rat) = SparsePoly(convert(Uniterm, t))
SparsePoly(t::Uniterm) = SparsePoly([coeff(t)], [degree(t)], var(t))

var(p::SparsePoly) = p.v
coeffs(p::SparsePoly) = p.coeffs
coeff(p::SparsePoly, i) = coeffs(p)[i]
term(p::SparsePoly, i) = Uniterm(p.coeffs[i], Uninomial(var(p), p.exps[i]))
terms(p::SparsePoly) = (term(p, i) for i in eachindex(coeffs(p)))
lt(p::SparsePoly) = term(p, lastindex(p.coeffs))
degree(p::SparsePoly) = iszero(p) ? -1 : p.exps[end]
Base.iszero(p::SparsePoly) = isempty(p.exps)

Base.zero(t::SparsePoly) = zero(lt(t))
Base.one(t::SparsePoly) = one(lt(t))
Base.deleteat!(p::SparsePoly, i::Int) = (deleteat!(p.coeff, i); deleteat!(p.exps, i);nothing)
Base.copy(p::SparsePoly) = SparsePoly(deepcopy(coeffs(p)), copy(p.exps), var(p))

#########################
# Arithmetic/Algorithms #
#########################

Base.:(^)(x::Uninomial, n::Integer) = Uninomial(var(x), degree(x)+n-1)

function Base.:(+)(x::Uninomial, y::Uninomial)
    v = var(x)
    v == var(y) || error("$x + $y contains different variables!")
    if x < y
        x, y = y, x
    end
    SparsePoly([1, 1], [degree(x), degree(y)], v)
end
function Base.:(+)(x::Uniterm, y::Uniterm)
    v = var(x)
    v == var(y) || error("$x + $y contains different variables!")
    if ismatch(x, y)
        c = x.coeff + y.coeff
        return iszero(c) ? SparsePoly(zero(x)) : SparsePoly(Uniterm(c, monomial(x)))
    else
        if x < y
            x, y = y, x
        end
        return SparsePoly([coeff(x), coeff(y)], [degree(x), degree(y)], v)
    end
end
Base.:(*)(x::Number, y::Uninomial) = Uniterm(x, y)
Base.:(*)(x::Uninomial, y::Number) = y * x

function Base.:(*)(x::Uninomial, y::Uninomial)
    v = var(x)
    v == var(y) || error("$x * $y contains different variables!")
    Uninomial(v, degree(x) + degree(y))
end
Base.:(*)(x::Uniterm, y::Uniterm) = Uniterm(coeff(x) * coeff(y), monomial(x) * monomial(y))

function Base.:*(p::SparsePoly, x::Uniterm)
    v = var(x)
    v == var(p) || error("$p * $x contains different variables!")
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
        return SparsePoly(cfs, exps, v)
    end
end
Base.:*(x::SparsePoly, y::SparsePoly) = sum(t->x * t, terms(y))

addcoef(x::Uniterm, c) = (c += coeff(x); return iszero(c), Uniterm(c, monomial(x)))
addcoef(x::Uniterm, c::Uniterm) = addcoef(x, coeff(c))
function addterm!(p::SparsePoly, x)
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
    v = var(x)
    @assert v == var(y)
    xd, yd = degree(x), degree(y)
    xd >= yd ? (Uninomial(xd - yd, v), false) : (x, true)
end

#=
function pseudorem(p::SparsePoly, d::SparsePoly)
    p = copy(p)
    while !isempty(terms(p))
        nx, fail = lt(p) / lt(d)
        if fail
            return p
        else
            p -= d * nx
        end
    end
    return p
end
=#

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
term2polycoeff(t::Term, v::IDType) = Term(t.coeff, Monomial(filter(!isequal(v), monomial(t).ids)))
function SparsePoly(p::MPoly, v::IDType)
    ts = terms(p)
    pows = map(ts) do t
        count(isequal(v), monomial(t).ids)
    end
    perm = sortperm(pows, rev=true)
    coeffs = Term[]
    exps = UInt[]

    sorted_pows = pows[perm]
    olddegree = sorted_pows[1]
    chunk_start_idx = idx = 1
    while idx <= length(sorted_pows)
        degree = sorted_pows[idx]
        if olddegree != degree # new chunk
            if olddegree > 0
                t = ts[perm[chunk_start_idx]]
                coeff = term2polycoeff(t, v)
                for i in chunk_start_idx+1:idx-1
                    t = ts[perm[i]] # TODO
                    coeff += term2polycoeff(t, v)
                end
            else
                coeff = ts[perm[chunk_start_idx]]
                for i in chunk_start_idx+1:idx-1
                    coeff += ts[perm[i]] # TODO
                end
            end
            push!(exps, olddegree)
            push!(coeffs, coeff)

            chunk_start_idx = idx
            olddegree = degree
        end
        idx += 1
    end
    if !isempty(chunk_start_idx:idx-1) # remember to handle the last one
        if olddegree > 0
            t = ts[perm[chunk_start_idx]]
            coeff = term2polycoeff(t, v)
            for i in chunk_start_idx+1:idx-1
                t = ts[perm[i]] # TODO
                coeff += term2polycoeff(t, v)
            end
        else
            coeff = ts[perm[chunk_start_idx]]
            for i in chunk_start_idx+1:idx-1
                coeff += ts[perm[i]] # TODO
            end
        end
        push!(coeffs, coeff)
        push!(exps, olddegree)
    end
    return SparsePoly(coeffs, exps, v)
end
