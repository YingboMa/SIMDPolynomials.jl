const IDType = UInt32
const NOT_A_VAR = typemax(IDType)
const EMPTY_IDS = IDType[]

struct Monomial <: AbstractMonomial
    ids::Vector{IDType}
end
Monomial() = Monomial(EMPTY_IDS)
degree(x::Monomial) = length(x.ids)
Base.copy(x::Monomial) = Monomial(copy(x.ids))

const VARNAME_DICT = Dict{IDType, String}(0 => "x", 1 => "y", 2 => "z")
const SUPERSCRIPTS = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']
const LOWERSCRIPTS = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']

function int2superscript(x)
    mapreduce(d->SUPERSCRIPTS[d + 1], *, Iterators.reverse(digits(x)))
end
function int2lowerscript(x)
    mapreduce(d->LOWERSCRIPTS[d + 1], *, Iterators.reverse(digits(x)))
end
function print_single_monomial(io, v, o, star=false)
    iszero(o) && return print(io, 1)
    star && (PRETTY_PRINT[] || print(io, '*'))
    if v == DO_NOT_CHECK_VAR || v === nothing
        print(io, "◌")
    else
        print(io, VARNAME_DICT[v])
    end
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

# graded lex
function Base.isless(x::Monomial, y::Monomial)
    dx = degree(x)
    dy = degree(y)
    dx < dy && return true
    dx > dy && return false
    for i in 1:dx
        vx = x.ids[i]
        vy = y.ids[i]
        vx != vy && return vx > vy
    end
    return false
end

function Base.gcd(x::Monomial, y::Monomial)
    g = Monomial(similar(x.ids, 0))
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
            i += 1
        elseif xk > yk
            j += 1
        else
            push!(g.ids, xk)
            i += 1
            j += 1
        end
    end
    return g
end
