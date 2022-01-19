###
### Contract: mutation without explicit copy is undefined behavior.
###

const PRETTY_PRINT = Ref(true)
const CoeffType = Union{Rational{<:Integer},Integer}

###
### AbstractMonomial: isless, degree, isunivariate
###
abstract type AbstractPolynomialLike <: Number end
abstract type AbstractMonomial <: AbstractPolynomialLike end

Base.isone(x::AbstractMonomial) = iszero(degree(x))

Base.:-(y::AbstractMonomial) = Term(-1, y)
Base.:+(y::T, x::T) where {T <: AbstractMonomial} = Term(y) + Term(x)
Base.:*(c::CoeffType, m::AbstractMonomial) = Term(c, m)
Base.:*(m::AbstractMonomial, c::CoeffType) = c * m
Base.:^(y::T, n::Integer) where {T <: AbstractMonomial} = iszero(n) ? T() : Base.power_by_squaring(y, n)

# only used for the univariate case
var(t) = metadata(t)
metadata(x) = nothing
isunivariate(x::AbstractMonomial) = false
dropmetadata(t) = t

###
### AbstractTerm: `coeff` and `monomial`
###
abstract type AbstractTerm <: AbstractPolynomialLike end
Base.promote_rule(::Type{C}, ::Type{M}) where {C<:CoeffType,M<:AbstractMonomial} = Term{C,M}
Base.promote_rule(t::Type{<:AbstractTerm}, ::Type{<:AbstractMonomial}) = t
Base.promote_rule(t::Type{<:AbstractTerm}, ::Type{<:CoeffType}) = t
@generated emptyterm(T::Type{<:AbstractTerm}) = T[]

metadata(x::AbstractTerm) = metadata(monomial(x))
degree(x::AbstractTerm) = degree(monomial(x))
ismatch(x::AbstractTerm, y::AbstractTerm) = monomial(x) == monomial(y)
Base.:(==)(x::AbstractTerm, y::AbstractTerm) = x === y || (coeff(x) == coeff(y) && monomial(x) == monomial(y))
Base.copy(x::T) where {T<:AbstractTerm} = T(copy(coeff(x)), copy(monomial(x)))
Base.isless(x::AbstractTerm, y::AbstractTerm) = isless(monomial(x), monomial(y))
Base.iszero(x::AbstractTerm) = iszero(coeff(x))
Base.isone(x::AbstractTerm) = isone(coeff(x)) && isone(monomial(x))
Base.isinteger(x::AbstractTerm) = isinteger(coeff(x))

dropmetadata(t::T) where {T<:AbstractTerm} = parameterless_type(T)(coeff(t), dropmetadata(monomial(t)))
addmetadata(t::T, meta) where {T<:AbstractTerm} = parameterless_type(T)(coeff(t), addmetadata(monomial(t), meta))

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

function add_or_sub(f::Union{typeof(+), typeof(-)}, x::T1, y::T1) where {T1<:AbstractTerm}
    T = dropmetadata(T1)
    v = checkmetadata(x, y)
    x = dropmetadata(x)
    y = dropmetadata(y)
    P = typeof(monomial(x)) <: Uninomial ? SPoly : MPoly
    if ismatch(x, y)
        c = f(coeff(x), coeff(y))
        return iszero(c) ? P(emptyterm(T), v) : P(T(c, monomial(x)), v)
    else
        y = f(y)
        if x < y
            x, y = y, x
        end
        return P(T[x, y], v)
    end
end
Base.:+(x::T, y::T) where {T<:AbstractTerm} = add_or_sub(+, x, y)
Base.:-(x::T, y::T) where {T<:AbstractTerm} = add_or_sub(-, x, y)
Base.:-(x::T) where {T<:AbstractTerm} = T(-coeff(x), monomial(x))

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
coefftype(::Type{<:Term{C}}) where {C} = C
monomialtype(::Type{<:Term{C,M}}) where {C,M} = M
# TODO
Term{C,M}(x) where {C,M<:AbstractMonomial} = Term(x, M())
#Term(x::M) where {M<:AbstractMonomial} = Term(1, M())
Term{C,M1}(x::M) where {C,M<:AbstractMonomial,M1<:AbstractMonomial} = Term(one(C), x)
Term(x::M) where {M<:AbstractMonomial} = Term(1, x)
Base.promote_rule(::Type{<:Term{C1,M1}}, ::Type{<:Term{C2,M2}}) where {C1,C2,M1,M2} = Term{promote_type(C1,C2), promote_type(M1,M2)}
#Term{C,M}(x::Term) where {C,M} = Term(covert(C, coeff(x)), covert(M, monomial(x)))
dropmetadata(T::Type{<:Term{C,M}}) where {C,M} = parameterless_type(T){C,dropmetadata(M)}

###
### AbstractPolynomial: terms, copy
###
abstract type AbstractPolynomial <: AbstractPolynomialLike end

Base.promote_rule(p::Type{<:AbstractPolynomial}, ::Type{<:AbstractMonomial}) = p
Base.promote_rule(p::Type{<:AbstractPolynomial}, ::Type{<:AbstractTerm}) = p
Base.promote_rule(p::Type{<:AbstractPolynomial}, ::Type{<:CoeffType}) = p

Base.iszero(x::AbstractPolynomial) = isempty(terms(x))
Base.isone(x::AbstractPolynomial) = (ts = terms(x); length(ts) == 1 && isone(only(ts)))
Base.:(==)(x::T, y::T) where {T<:AbstractPolynomial} = x === y || (terms(x) == terms(y))

lt(p::AbstractPolynomial) = first(terms(p))
lc(p::AbstractPolynomial) = coeff(lt(p))
degree(p::AbstractPolynomial) = degree(lt(p))

checkmetadata(x::AbstractPolynomialLike, y::AbstractPolynomialLike) = metadata(x)

function Base.show(io::IO, p::AbstractPolynomial)
    ts = terms(p)
    if isempty(ts)
        print(io, 0)
        return
    end
    md = metadata(p)
    n = length(ts)
    t1, tr = Iterators.peel(ts)
    show(io, addmetadata(t1, md))
    for t′ in tr
        t = addmetadata(t′, md)
        if t.coeff < 0
            print(io, " - ")
            t = -t
        else
            print(io, " + ")
        end
        show(io, t)
    end
end
