# default polynomial type
#=
struct MPoly{T} <: AbstractPolynomial{T}
    terms::T
end
=#
#MPoly() = MPoly(EMPTY_TERMS)
#MPoly(x::AbstractTerm) = MPoly([x])
#MPoly(x::M) where {M<:AbstractMonomial} = MPoly(Term(x))
MPoly{T}(x::AbstractTerm) where T = MPoly([x])
MPoly{T}(x::M) where {T,M<:AbstractMonomial} = MPoly(Term(x))
MPoly{T}(x::CoeffType) where {T} = MPoly(eltype(T)(x))
terms(x::MPoly) = x.terms
function nvariables(p::MPoly)
    if monomialtype(p) <: PackedMonomial
        nvariables(monomial(lt(p)))
    else
        error()
    end
end
coeffs(x::MPoly) = (coeff(t) for t in terms(x))

Base.copy(x::MPoly) = MPoly(map(copy, terms(x)))
Base.one(::Type{<:MPoly{T}}) where T = MPoly(eltype(T)(1))

Base.:*(x::AbstractTerm, p::MPoly) = p * x
function Base.:*(p::MPoly, x::T) where {T<:AbstractTerm}
    if iszero(p)
        return p
    elseif iszero(x)
        return MPoly(emptyterm(T))
    elseif isone(x)
        return p
    else
        MPoly(T[t * x for t in terms(p) if !iszero(t)])
    end
end
function Base.:*(p::MPoly, x::MPoly)
    ts = terms(x)
    out = similar(terms(p), length(ts) * length(terms(p)))
    current_len = 0
    s = MPoly(out)
    for tp in terms(p)
        start = 1
        for tx in terms(x)
            _, (current_len, start) = _add!(s, tp * tx, current_len, start)
        end
    end
    resize!(out, current_len)
    (debugmode() && !issorted(terms(s), rev=true)) && throw("Polynomial not sorted!")

    return s
end
Base.:+(p::MPoly, x::AbstractTerm) = add!(copy(p), x)
Base.:-(p::MPoly, x::AbstractTerm) = sub!(copy(p), x)

Base.:-(p::MPoly) = -1 * p
sub!(p::MPoly, x::AbstractTerm) = add!(p, -x)

addcoef(x::T, c) where {T<:AbstractTerm} = (c += coeff(x); return iszero(c), T(c, monomial(x)))
addcoef(x::T, c::T) where {T<:AbstractTerm} = addcoef(x, c.coeff)
subcoef(x::T, c) where {T<:AbstractTerm} = (c = coeff(x) - c; return iszero(c), T(c, monomial(x)))
subcoef(x::T, c::T) where {T<:AbstractTerm} = subcoef(x, c.coeff)

function add!(p::MPoly, x::AbstractTerm)
    q, (current_len, _) = _add!(p, x, length(terms(p)), 1);
    resize!(terms(q), current_len)
    (debugmode() && !issorted(terms(q), rev=true)) && throw("Polynomial not sorted!")
    q
end

@inline function _add!(p::MPoly, x::AbstractTerm, current_len, start)
    iszero(x) && return p, (current_len, 1)
    ts = terms(p)
    i = start
    for i in start:current_len
        t = ts[i]
        if t < x
            current_len += 1
            if current_len > length(ts)
                insert!(ts, i, x)
            else
                @inbounds @simd for j in current_len:-1:i+1
                    ts[j] = ts[j-1]
                end
                ts[i] = x
            end
            return p, (current_len, i)
        elseif ismatch(t, x)
            return p, _add_term_matched(ts, i, t, x, current_len)
        end
    end
    current_len += 1
    if current_len > length(ts)
        push!(ts, x)
    else
        ts[current_len] = x
    end
    return p, (current_len, i)
end

@inline function _add_term_matched(ts, i, t, x, current_len)
    iz, t = addcoef(t, x)
    if iz
        # deleteat is somehow faster
        deleteat!(ts, i)
        return current_len-1, max(1, i-1)
        #current_len -= 1
        #@inbounds @simd for j in i:current_len
        #    ts[j] = ts[j+1]
        #end
        #return current_len, max(1, i-1)
    else
        ts[i] = t
        return current_len, i
    end
end

@inline function _add_rev!(p::MPoly, x::AbstractTerm, _end=1)
    iszero(x) && return p, 1
    ts = terms(p)
    i = _end
    for i in length(ts):-1:_end
        t = ts[i]
        if t > x
            insert!(ts, i+1, x)
            return p, i
        elseif ismatch(t, x)
            return p, _add_term_matched_rev(ts, i, t, x)
        end
    end
    push!(ts, x)
    return p, i
end

function _add_term_matched_rev(ts, i, t, x)
    iz, t = addcoef(t, x)
    if iz
        deleteat!(ts, i)
        return max(1, i-1)
    else
        ts[i] = t
        return i
    end
end

function add!(p::AbstractPolynomial, x::AbstractPolynomial)
    tp = terms(p)
    tx = terms(x)
    if tp isa Vector
        sizehint!(tp, length(tp) + length(tx))
    end
    for t in tx
        add!(p, t)
    end
    return p
end
function sub!(p::AbstractPolynomial, x::AbstractPolynomial)
    tp = terms(p)
    tx = terms(x)
    if tp isa Vector
        sizehint!(tp, length(tp) + length(tx))
    end
    for t in tx
        sub!(p, t)
    end
    return p
end
Base.:+(p::AbstractPolynomial, x::AbstractPolynomial) = add!(copy(p), x)
Base.:-(p::AbstractPolynomial, x::AbstractPolynomial) = sub!(copy(p), x)

lt(p::AbstractPolynomial) = first(terms(p))
lc(p::AbstractPolynomial) = coeff(lt(p))
function rmlt!(p::MPoly)
    ts = terms(p)
    popfirst!(ts)
    return p
end
function takelt!(p::MPoly, x::MPoly)
    add!(p, lt(x))
    rmlt!(x)
    return p
end

function fnmadd!(p::MPoly, l, s)
    for li in terms(l)
        sub!(p, li * s)
    end
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
            #p -= d * nx
            #q += nx
            fnmadd!(p, d, nx)
            add!(q, nx)
        end
    end
    return q, r
end

function divexact(p::MPoly, d::MPoly)
    p = copy(p)
    ts = similar(terms(p), 0)
    sizehint!(ts, length(terms(d))-1)
    q = MPoly(ts)
    _end = 1
    while !isempty(terms(p))
        nx, fail = lt(p) / lt(d)
        @assert !fail
        fnmadd!(p, d, nx)
        _, _end = _add_rev!(q, nx, _end)
    end
    (debugmode() && !issorted(terms(q), rev=true)) && throw("Polynomial not sorted!")
    return q
end

function Base.rem(p::AbstractPolynomial, d::AbstractPolynomial)
    p = copy(p)
    r = MPoly(similar(terms(p), 0))
    while !isempty(terms(p))
        nx, fail = lt(p) / lt(d)
        if fail
            takelt!(r, p)
        else
            sub!(p, d * nx)
        end
    end
    return r
end

# TODO
Base.:(/)(x::MPoly, y::MPoly) = divexact(x, y)

Base.gcd(x::Union{AbstractTerm,AbstractMonomial}, y::MPoly) = gcd(MPoly(x), y)
Base.gcd(x::MPoly, y::Union{AbstractTerm,AbstractMonomial}) = gcd(y, x)
function Base.gcd(x::MPoly, y::MPoly)
    # trival case
    if iszero(x) || isone(y)
        return y
    elseif iszero(y) || isone(x) || x == y
        return x
    end

    v1, p1 = to_univariate(x)
    v2, p2 = to_univariate(y)
    if v1 < v2
        x, y = y, x
        v1, v2 = v2, v1
        p1, p2 = p2, p1
    end
    # v2 < v1
    # both are constants
    v2 == NOT_A_VAR && return MPoly(gcd(lt(x), lt(y)))
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
            vv = firstid(m) # get the minvar
            if vv < v
                v = vv
            end
        end
    end
    return IDType(v)
end

function to_univariate(x::MPoly)
    v = pick_var(x)
    v, (v == NOT_A_VAR ? nothing : SparsePoly(x, v))
end
