###
### Contract: mutation without explicit copy is undefined behavior.
###

const PRETTY_PRINT = Ref(true)
const CoeffType = Union{Rational{<:Integer},Integer}

###
### AbstractMonomial: isless, degree
###
abstract type AbstractMonomial <: Number end

Base.isone(x::AbstractMonomial) = iszero(degree(x))
Base.one(::Type{<:T}) where {T<:AbstractMonomial} = T()

Base.:-(y::AbstractMonomial) = Term(-1, y)
Base.:+(y::T, x::T) where {T <: AbstractMonomial} = Term(y) + Term(x)
Base.:*(c::CoeffType, m::AbstractMonomial) = Term(c, m)
Base.:*(m::AbstractMonomial, c::CoeffType) = c * m
Base.:^(y::T, n::Integer) where {T <: AbstractMonomial} = iszero(n) ? T() : Base.power_by_squaring(y, n)
monomialtype(x::AbstractMonomial) = typeof(x)

###
### AbstractTerm: `coeff` and `monomial`
###
abstract type AbstractTerm <: Number end
Base.promote_rule(::Type{C}, ::Type{M}) where {C<:CoeffType,M<:AbstractMonomial} = Term{C,M}
Base.promote_rule(t::Type{<:AbstractTerm}, ::Type{<:AbstractMonomial}) = t
Base.promote_rule(t::Type{<:AbstractTerm}, ::Type{<:CoeffType}) = t
@generated emptyterm(T::Type{<:AbstractTerm}) = T[]

degree(x::AbstractTerm) = degree(monomial(x))
ismatch(x::T, y::T) where {T<:AbstractTerm} = monomial(x) == monomial(y)
Base.:(==)(x::AbstractTerm, y::AbstractTerm) = x === y || (coeff(x) == coeff(y) && monomial(x) == monomial(y))
Base.copy(x::T) where {T<:AbstractTerm} = T(copy(coeff(x)), copy(monomial(x)))
Base.isless(x::T, y::T) where {T<:AbstractTerm} = isless(monomial(x), monomial(y))
Base.iszero(x::AbstractTerm) = iszero(coeff(x))
Base.isone(x::AbstractTerm) = isone(coeff(x)) && isone(monomial(x))
Base.zero(t::T) where {T<:AbstractTerm} = parameterless_type(T)(zero(coeff(t)), one(monomial(t)))
Base.one(t::T) where {T<:AbstractTerm} = parameterless_type(T)(one(coeff(t)), one(monomial(t)))
Base.isinteger(x::AbstractTerm) = isinteger(coeff(x))

monomialtype(x::AbstractTerm) = monomialtype(monomial(x))

Base.:*(x::T, y::T) where {T<:AbstractTerm} = T(coeff(x) * coeff(y), monomial(x) * monomial(y))
# TODO
function Base.:(/)(y::T, x::T) where {T<:AbstractTerm}
    m, fail = monomial(y) / monomial(x)
    if fail
        T(coeff(y), m), fail
    else
        T(coeff(y)/coeff(x), m), fail
    end
end

function Base.:+(x::T, y::T) where {T<:AbstractTerm}
    if ismatch(x, y)
        c = coeff(x) + coeff(y)
        return iszero(c) ? MPoly(emptyterm(T)) : MPoly(T(c, monomial(x)))
    else
        if x < y
            x, y = y, x
        end
        return MPoly(T[x, y])
    end
end

Base.:-(x::T) where {T<:AbstractTerm} = T(-coeff(x), monomial(x))
function Base.:-(x::T, y::T) where {T<:AbstractTerm}
    if ismatch(x, y)
        c = coeff(x) - coeff(y)
        return iszero(c) ? MPoly(T[]) : MPoly(T(c, monomial(x)))
    else
        y = -y
        if x < y
            x, y = y, x
        end
        return MPoly(T[x, y])
    end
end

function Base.gcd(x::T, y::T) where {T<:AbstractTerm}
    #g, a, b = gcd(monomial(x), monomial(y))
    g = gcd(monomial(x), monomial(y))
    gr = gcd(x.coeff, y.coeff)
    return T(gr, g)#, Term(x.coeff / gr, a), Term(y.coeff / gr, b)
end

print_coeff(io::IO, coeff) = isinteger(coeff) ? print(io, Integer(coeff)) : print(io, coeff)
function Base.show(io::IO, x::AbstractTerm)
    printed = false
    if coeff(x) != 1
        if coeff(x) == -1
            print(io, '-')
        else
            print_coeff(io, coeff(x))
            printed = true
            PRETTY_PRINT[] || print(io, '*')
        end
    end
    m = monomial(x)
    if !(printed && isone(m))
        show(io, m)
    end
end

# default term type
struct Term{C,M<:AbstractMonomial} <: AbstractTerm
    coeff::C
    monomial::M
end
coeff(x::Term) = x.coeff
monomial(x::Term) = x.monomial
Term{M}(x) where {M<:AbstractMonomial} = Term(x, M())
Term(x::M) where {M<:AbstractMonomial} = Term(1, x)
Term{A,B}(x::M) where {A,B,M<:AbstractMonomial} = Term(1, x)
Term{A,B}(x) where {A,B} = Term(x, B())
#const EMPTY_TERM = Term[]

###
### AbstractPolynomial: terms, copy
###
abstract type AbstractPolynomial <: Number end

Base.promote_rule(p::Type{<:AbstractPolynomial}, ::Type{<:AbstractMonomial}) = p
Base.promote_rule(p::Type{<:AbstractPolynomial}, ::Type{<:AbstractTerm}) = p
Base.promote_rule(p::Type{<:AbstractPolynomial}, ::Type{<:CoeffType}) = p

Base.iszero(x::AbstractPolynomial) = isempty(terms(x))
Base.isone(x::AbstractPolynomial) = (ts = terms(x); length(ts) == 1 && isone(only(ts)))
Base.:(==)(x::T, y::T) where {T<:AbstractPolynomial} = x === y || (terms(x) == terms(y))
monomialtype(p::AbstractPolynomial) = monomialtype(lt(p))

function Base.show(io::IO, p::AbstractPolynomial)
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
