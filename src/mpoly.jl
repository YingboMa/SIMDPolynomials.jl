const Rat = Union{Rational{<:Integer},Integer}

const IDType = UInt32
const EMPTY_IDS = IDType[]

struct Monomial <: Number
    ids::Vector{IDType}
end
Monomial() = Monomial(EMPTY_IDS)

const VARNAME_DICT = Dict{IDType, String}(0 => "x", 1 => "y", 2 => "z")
const PRETTY_PRINT = Ref(true)
const SUPERSCRIPTS = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']

function int2superscript(x)
    mapreduce(d->SUPERSCRIPTS[d + 1], *, Iterators.reverse(digits(x)))
end
function print_single_monomial(io, v, o, star=false)
    star && (PRETTY_PRINT[] || print(io, '*'))
    print(io, VARNAME_DICT[v])
    if o > 1
        if PRETTY_PRINT[]
            print(io, int2superscript(o))
        else
            print(io, '^', o)
        end
    elseif o == 1 || o == 0
    else
        error("unreachable")
    end
end
function Base.show(io::IO, x::Monomial)
    isempty(x.ids) && (print(io, '1'); return)

    ids = x.ids
    v = first(ids)
    count = 0
    star = false
    for id in ids
        if id == v
            count += 1
        else
            print_single_monomial(io, v, count, star)
            star = true
            v = id
            count = 1
        end
    end
    if count > 0
        print_single_monomial(io, v, count, star)
    end
    return nothing
end

Base.:-(y::Monomial) = Term(-1, y)
Base.:+(y::Monomial, x::Monomial) = Term(y) + Term(x)
function Base.:*(y::Monomial, x::Monomial)
    ids = y.ids
    i = j = 1
    n0 = length(ids)
    n1 = length(x.ids)
    r = Monomial(similar(ids, n0 + n1))
    T = eltype(ids)
    M = typemax(T)
    for k in 1:(n0 + n1)
        a = i <= n0 ? ids[i] : M
        b = j <= n1 ? x.ids[j] : M
        if a < b
            i += 1
            r.ids[k] = a
        else
            j += 1
            r.ids[k] = b
        end
    end
    return r
end

function Base.:(/)(y::Monomial, x::Monomial)
    n = Monomial(similar(y.ids, 0))
    i = j = 1
    ids = y.ids
    n0 = length(ids)
    n1 = length(x.ids)
    T = eltype(ids)
    M = typemax(T)
    while (i + j - 2) < (n0 + n1)
        a = i <= n0 ? ids[i] : M
        b = j <= n1 ? x.ids[j] : M
        if a < b
            push!(n.ids, a)
            i += 1
        elseif a == b
            i += 1
            j += 1
        else
            return n, true
        end
    end
    return n, false
end

Base.:(==)(x::Monomial, y::Monomial) = (x === y) || (x.ids == y.ids)

Base.isone(x::Monomial) = isempty(x.ids)
degree(x::Monomial) = length(x.ids)
# graded lex
function Base.isless(x::Monomial, y::Monomial)
    dx = degree(x)
    dy = degree(y)
    dx < dy && return true
    dx > dy && return false
    for i in 1:dx
        x.ids[i] > y.ids[i] && return true
    end
    return false
end

Base.:*(c::Rat, m::Monomial) = Term(c, m)
Base.:*(m::Monomial, c::Rat) = c * m

abstract type AbstractTerm <: Number end
struct Term <: AbstractTerm
    coeff::Rational{Int}
    monomial::Monomial
end
monomial(x::Term) = x.monomial
Term(x) = Term(x, Monomial())
Term(m::Monomial) = Term(1, m)
Base.Base.promote_rule(::Type{Term}, ::Type{Monomial}) = Term
print_coeff(io::IO, coeff) = isinteger(coeff) ? print(io, Integer(coeff)) : print(io, coeff)
function Base.show(io::IO, x::Term)
    printed = false
    if x.coeff != 1
        if x.coeff == -1
            print(io, '-')
        else
            print_coeff(io, x.coeff)
            printed = true
            PRETTY_PRINT[] || print(io, '*')
        end
    end
    m = monomial(x)
    if !(printed && isone(m))
        show(io, m)
    end
end

ismatch(x::Term, y::Term) = monomial(x) == monomial(y)
Base.isless(x::Term, y::Term) = isless(monomial(x), monomial(y))
Base.iszero(x::Term) = iszero(x.coeff)
Base.isone(x::Term) = isone(x.coeff) && isone(monomial(x))
addcoef(x::Term, c) = (c += x.coeff; return iszero(c), Term(c, monomial(x)))
addcoef(x::Term, c::Term) = addcoef(x, c.coeff)
subcoef(x::Term, c) = (c = x.coeff - c; return iszero(c), Term(c, monomial(x)))
subcoef(x::Term, c::Term) = subcoef(x, c.coeff)

Base.isinteger(x::Term) = isinteger(x.coeff)
Base.:*(x::Term, y::Term) = Term(x.coeff * y.coeff, monomial(x) * monomial(y))
function Base.:+(x::Term, y::Term)
    if ismatch(x, y)
        c = x.coeff + y.coeff
        return iszero(c) ? Poly(Term(0)) : Poly(Term(c, monomial(x)))
    else
        return x < y ? Poly(Term[y, x]) : Poly(Term[x, y])
    end
end

Base.:-(x::Term) = Term(-x.coeff, monomial(x))
function Base.:-(x::Term, y::Term)
    if ismatch(x, y)
        c = x.coeff - y.coeff
        return iszero(c) ? Poly(Term(0)) : Poly(Term(c, x.coeff))
    else
        y = -y
        return x < y ? Poly(Term[y, x]) : Poly(Term[x, y])
    end
end

function Base.:(/)(y::AbstractTerm, x::AbstractTerm)
    m, fail = monomial(y) / monomial(x)
    if fail
        Term(y.coeff, m), fail
    else
        Term(y.coeff/x.coeff, m), fail
    end
end
function Base.:(/)(y::Uniterm, x::Uniterm)
    m, fail = monomial(y) / monomial(x)
    if fail
        Uniterm(y.coeff, m), fail
    else
        d, r = divrem(y.coeff, x.coeff)
        if iszero(r)
            Uniterm(d, m), fail
        else
            Uniterm(d, m), true
        end
    end
end

const EMPTY_TERMS = Term[]
abstract type AbstractPoly <: Number end
struct Poly <: AbstractPoly
    terms::Vector{Term}
end
Poly() = Poly(EMPTY_TERMS)
Poly(x::Term) = Poly([x])
Poly(x::Union{Rat,Monomial}) = Poly(Term(x))
Base.Base.promote_rule(::Type{Poly}, ::Type{Monomial}) = Poly
Base.Base.promote_rule(::Type{Poly}, ::Type{Term}) = Poly
Base.Base.promote_rule(::Type{Poly}, ::Type{<:Rat}) = Poly
terms(x::Poly) = x.terms
Base.empty!(x::Poly) = (empty!(terms(x)); x)
function Base.show(io::IO, p::Poly)
    ts = terms(p)
    if isempty(ts)
        print(io, 0)
        return
    end
    n = length(ts)
    show(io, ts[1])
    for i in 2:n
        t = ts[i]
        if t.coeff < 0
            print(io, " - ")
            t = -t
        else
            print(io, " + ")
        end
        show(io, t)
    end
end

Base.copy(x::Poly) = Poly(copy(terms(x)))
#function largercopy(x::Poly, i::Int)
#    n = length(terms(x))
#    terms = similar(terms(x), i + n)
#    copyto!(terms, 1, terms(x), 1, n)
#    Poly(terms)
#end
Base.:(==)(x::Poly, y::Poly) = x === y || (terms(x) == terms(y))
Base.iszero(x::Poly) = isempty(terms(x))

Base.:*(x::Term, p::Poly) = p * x
function Base.:*(p::Poly, x::Term)
    iszero(x) ? Poly() :
        isone(x) ? p : Poly(Term[t * x for t in terms(p)])
end
function Base.:*(p::Poly, x::Poly)
    sum(t->p * t, terms(x))
end
Base.:+(p::Poly, x::Term) = addterm!(copy(p), x)
Base.:-(p::Poly, x::Term) = subterm!(copy(p), x)

subterm!(p::Poly, x) = addterm!(p, -x)
function addterm!(p::Poly, x)
    iszero(x) && return p
    ts = terms(p)
    for (i, t) in enumerate(ts)
        if ismatch(t, x)
            iz, t = addcoef(t, x)
            if iz
                deleteat!(ts, i)
            else
                ts[i] = t
            end
            return p
        elseif t < x
            insert!(ts, i, x)
            return p
        end
    end
    push!(ts, x)
    return p
end

function Base.:+(p::Poly, x::Poly)
    s = copy(p)
    for t in terms(x)
        addterm!(s, t)
    end
    return s
end
function Base.:-(p::Poly, x::Poly)
    s = copy(p)
    for t in terms(x)
        subterm!(s, t)
    end
    return s
end

lt(p::AbstractPoly) = first(terms(p))
function rmlt!(p::Poly)
    ts = terms(p)
    popfirst!(ts)
    return p
end
function takelt!(p::Poly, x::Poly)
    addterm!(p, lt(x))
    rmlt!(x)
    return p
end

function Base.divrem(p::Poly, d::Poly)
    p = copy(p)
    q = Poly(similar(terms(p), 0))
    r = Poly(similar(terms(p), 0))
    while !isempty(terms(p))
        nx, fail = lt(p) / lt(d)
        if fail
            takelt!(r, p)
        else
            p -= d * nx
            q += nx
        end
    end
    return q, r
end

function Base.rem(p::AbstractPoly, d::AbstractPoly)
    p = copy(p)
    r = Poly(similar(terms(p), 0))
    while !isempty(terms(p))
        nx, fail = lt(p) / lt(d)
        if fail
            takelt!(r, p)
        else
            p -= d * nx
        end
    end
    return r
end

function Base.gcd(x::Monomial, y::Monomial)
    g = Monomial(similar(x.ids, 0))
    a = Monomial(similar(x.ids, 0))
    b = Monomial(similar(x.ids, 0))

    i = j = 1
    n0 = length(x.ids)
    n1 = length(y.ids)
    n = n0 + n1 + 2
    T = eltype(x.ids)
    M = typemax(T)
    while i + j < n
        xk = i <= n0 ? x.ids[i] : M
        yk = j <= n1 ? y.ids[j] : M
        if xk < yk
            push!(a.ids, xk)
            i += 1
        elseif xk > yk
            push!(b.ids, yk)
            j += 1
        else
            push!(g.ids, xk)
            i += 1
            j += 1
        end
    end
    return g, a, b
end

function Base.gcd(x::Term, y::Term)
    g, a, b = gcd(x, y)
    gr = gcd(x.coeff, y.coeff)
    return Term(gr, g), Term(x.coeff / gr, a), Term(y.coeff / gr, b)
end

function Base.:(/)(x::Poly, y::Poly)
    d, r = divrem(x, y)
    @assert iszero(r)
    d
end

function univariate_gcd(a::AbstractPoly, b::AbstractPoly)
    x = copy(a)
    y = copy(b)
    while !iszero(y)
        x = rem(x, y)
        x, y = y, x
    end
    return x#, a / x, b / x
end

# mpoly2poly(5x⁴ + 2x³ + x² + xy + y^2, 0x00001 [y])
# ->
#              ---------------- constants
#              |
# y^2 + xy + (5x⁴ + 2x³ + x²)
#
# M2P(a, v)
struct Uninomial
    d::Int
    v::IDType
end
struct Uniterm <: AbstractTerm
    coeff::Poly
    uninomial::Uninomial
end
struct UniPoly <: AbstractPoly
    #dc::Vector{Tuple{UInt32, Vector{Term}}}
    terms::Vector{Uniterm}
end
terms(p::UniPoly) = p.terms
Base.iszero(p::UniPoly) = isempty(terms(p))
monomial(p::Uniterm) = p.uninomial

function UniPoly(p::Poly, v::IDType)
    ts = terms(p)
    pows = map(ts) do t
        count(isequal(v), monomial(t).ids)
    end
    perm = sortperm(pows, rev=true)
    sorted_pows = pows[perm]
    olddegree = sorted_pows[1]
    chunk_start_idx = idx = 1
    sth = Uniterm[]
    while idx <= length(sorted_pows)
        degree = sorted_pows[idx]
        if olddegree != degree # new chunk
            coeff = ts[perm[chunk_start_idx:idx-1]]
            if olddegree > 0
                coeff = map(coeff) do t
                    Term(t.coeff, Monomial(filter(!isequal(v), monomial(t).ids)))
                end
            end
            push!(sth, Uniterm(sum(coeff), Uninomial(olddegree, v)))
            chunk_start_idx = idx
            olddegree = degree
        end
        idx += 1
    end
    coeff = ts[perm[chunk_start_idx:idx-1]]
    if olddegree > 0
        coeff = map(coeff) do t
            Term(t.coeff, Monomial(filter(!isequal(v), monomial(t).ids)))
        end
    end
    push!(sth, Uniterm(sum(coeff), Uninomial(olddegree, v)))
    return UniPoly(sth)
end

function Base.divrem(p::UniPoly, d::UniPoly)
    p = copy(p)
    q = UniPoly(similar(terms(p), 0))
    while !isempty(terms(p))
        nx, fail = lt(p) / lt(d)
        if fail
            return q, p
        else
            p -= d * nx
            q += nx
        end
    end
    return q, p
end

function Base.rem(p::UniPoly, d::UniPoly)
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

function Base.:(/)(x::Uninomial, y::Uninomial)
    # d::Int
    # v::IDType
    @assert x.v == y.v
    x.d >= y.d ? (Uninomial(x.v, x.d - y.d), false) : (x, true)
end
