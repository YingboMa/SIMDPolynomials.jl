const Rat = Union{Rational{<:Integer},Integer}

const IDType = UInt32
const NOT_A_VAR = typemax(IDType)
const EMPTY_IDS = IDType[]

abstract type AbstractMonomial <: Number end

struct Monomial <: AbstractMonomial
    ids::Vector{IDType}
end
Monomial() = Monomial(EMPTY_IDS)
Base.promote_rule(::Type{Monomial}, ::Type{<:Rat}) = Term

const VARNAME_DICT = Dict{IDType, String}(0 => "x", 1 => "y", 2 => "z")
const PRETTY_PRINT = Ref(true)
const SUPERSCRIPTS = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']

function int2superscript(x)
    mapreduce(d->SUPERSCRIPTS[d + 1], *, Iterators.reverse(digits(x)))
end
function print_single_monomial(io, v, o, star=false)
    iszero(o) && return print(io, 1)
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
Base.:^(y::Monomial, n::Integer) = iszero(n) ? Monomial() : Base.power_by_squaring(y, n)
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

Base.isone(x::Monomial) = iszero(degree(x))
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
coeff(x::Term) = x.coeff
degree(x::Term) = degree(monomial(x))

Base.promote_rule(t::Type{<:AbstractTerm}, ::Type{<:AbstractMonomial}) = t
Base.promote_rule(t::Type{<:AbstractTerm}, ::Type{<:Rat}) = t
Base.:(==)(x::AbstractTerm, y::AbstractTerm) = x === y || (coeff(x) == coeff(y) && monomial(x) == monomial(y))

print_coeff(io::IO, coeff) = isinteger(coeff) ? print(io, Integer(coeff)) : print(io, coeff)
function Base.show(io::IO, x::AbstractTerm)
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

ismatch(x::AbstractTerm, y::AbstractTerm) = monomial(x) == monomial(y)
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
        return iszero(c) ? MPoly(Term(0)) : MPoly(Term(c, monomial(x)))
    else
        return x < y ? MPoly(Term[y, x]) : MPoly(Term[x, y])
    end
end

Base.:-(x::Term) = Term(-x.coeff, monomial(x))
function Base.:-(x::Term, y::Term)
    if ismatch(x, y)
        c = x.coeff - y.coeff
        return iszero(c) ? MPoly(Term(0)) : MPoly(Term(c, x.coeff))
    else
        y = -y
        return x < y ? MPoly(Term[y, x]) : MPoly(Term[x, y])
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

const EMPTY_TERMS = Term[]
abstract type AbstractPoly <: Number end
struct MPoly <: AbstractPoly
    terms::Vector{Term}
end
MPoly() = MPoly(EMPTY_TERMS)
MPoly(x::Term) = MPoly([x])
MPoly(x::Union{Rat,Monomial}) = MPoly(Term(x))
Base.promote_rule(p::Type{<:AbstractPoly}, ::Type{<:AbstractMonomial}) = p
Base.promote_rule(p::Type{<:AbstractPoly}, ::Type{<:AbstractTerm}) = p
Base.promote_rule(p::Type{<:AbstractPoly}, ::Type{<:Rat}) = p
terms(x::MPoly) = x.terms
Base.empty!(x::MPoly) = (empty!(terms(x)); x)
function Base.show(io::IO, p::AbstractPoly)
    ts = terms(p)
    if isempty(ts)
        print(io, 0)
        return
    end
    n = length(ts)
    t1, tr = Iterators.peel(ts)
    show(io, t1)
    for t in tr
        if t.coeff < 0
            print(io, " - ")
            t = -t
        else
            print(io, " + ")
        end
        show(io, t)
    end
end

Base.copy(x::MPoly) = MPoly(copy(terms(x)))
#function largercopy(x::MPoly, i::Int)
#    n = length(terms(x))
#    terms = similar(terms(x), i + n)
#    copyto!(terms, 1, terms(x), 1, n)
#    MPoly(terms)
#end
Base.:(==)(x::MPoly, y::MPoly) = x === y || (terms(x) == terms(y))
Base.iszero(x::MPoly) = isempty(terms(x))
Base.isone(x::MPoly) = (ts = terms(x); length(ts) == 1 && isone(only(ts)))

Base.:*(x::Term, p::MPoly) = p * x
function Base.:*(p::MPoly, x::Term)
    if iszero(x)
        return MPoly()
    elseif isone(x)
        return p
    else
        ts = Term[t * x for t in terms(p) if !iszero(t)]
        MPoly(ts)
    end
end
function Base.:*(p::MPoly, x::MPoly)
    sum(t->p * t, terms(x))
end
Base.:+(p::MPoly, x::Term) = addterm!(copy(p), x)
Base.:-(p::MPoly, x::Term) = subterm!(copy(p), x)

Base.:-(p::MPoly) = -1 * p
subterm!(p::MPoly, x) = addterm!(p, -x)
function addterm!(p::MPoly, x)
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

function Base.:+(p::AbstractPoly, x::AbstractPoly)
    s = copy(p)
    for t in terms(x)
        addterm!(s, t)
    end
    return s
end
function Base.:-(p::AbstractPoly, x::AbstractPoly)
    s = copy(p)
    for t in terms(x)
        subterm!(s, t)
    end
    return s
end

lt(p::AbstractPoly) = first(terms(p))
function rmlt!(p::MPoly)
    ts = terms(p)
    popfirst!(ts)
    return p
end
function takelt!(p::MPoly, x::MPoly)
    addterm!(p, lt(x))
    rmlt!(x)
    return p
end

function Base.divrem(p::MPoly, d::MPoly)
    p = copy(p)
    q = MPoly(similar(terms(p), 0))
    r = MPoly(similar(terms(p), 0))
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
    r = MPoly(similar(terms(p), 0))
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
    return g#, a, b
end

function Base.gcd(x::Term, y::Term)
    #g, a, b = gcd(monomial(x), monomial(y))
    g = gcd(monomial(x), monomial(y))
    gr = gcd(x.coeff, y.coeff)
    return Term(gr, g)#, Term(x.coeff / gr, a), Term(y.coeff / gr, b)
end

function divexact(x::MPoly, y::MPoly)
    d, r = divrem(x, y)
    @assert iszero(r)
    d
end

Base.:(/)(x::MPoly, y::MPoly) = divexact(x, y)

Base.gcd(x::Term, y::MPoly) = gcd(MPoly(x), y)
Base.gcd(x::MPoly, y::Term) = gcd(y, x)
function Base.gcd(x::MPoly, y::MPoly)
    # trival case
    if iszero(x) || isone(y)
        return y
    elseif iszero(y) || isone(x)
        return x
    end
    x == y && return x

    v1, p1 = to_univariate(x)
    v2, p2 = to_univariate(y)
    if v1 < v2
        x, y = y, x
        v1, v2 = v2, v1
        p1, p2 = p2, p1
    end
    # v2 < v1
    # both are constants
    (v1 == NOT_A_VAR && v2 == NOT_A_VAR) && return MPoly(gcd(lt(x), lt(y)))
    if v2 < v1
        # `v2` in p2 doesn't exist in `x`, so the gcd at this level is 1 and we
        # just move on to the next level
        return gcd(x, content(p2))
    end
    v1 == v2 || error("unreachable")

    g = gcd(p1, p2)
    return univariate_to_multivariate(g)
end

function pick_var(x::MPoly)
    ts = terms(x)
    v = NOT_A_VAR
    for (i, t) in enumerate(ts)
        if degree(t) > 0
            m = monomial(t)
            vv = m.ids[1] # get the minvar
            if vv < v
                v = vv
            end
        end
    end
    return v
end

function to_univariate(x::MPoly)
    v = pick_var(x)
    v, (v == NOT_A_VAR ? nothing : SparsePoly(x, v))
end

function univariate_to_multivariate(g::SparsePoly{<:AbstractPoly})
    cfs = coeffs(g)
    eps = g.exps
    v = var(g)
    # TODO
    sum(zip(cfs, eps)) do (c, e)
        c * Monomial(fill(v, e))
    end
end
