# E = exponent bits = 7 = 1 byte
# L = local var bits = 8
# global: x0, x1, ..., x(2^16 - 1)
# compression: global x -> local l_0, ..., l_(L-1)
# the max exp is 2^E-1
# each variable takes (E) = 8 bits
# all variabes take (L+1)*(E+1) = 256 bytes
# Monomial with capacity L+1: [tot degree bits; l_{L-1} exp bits, ..., l_1 exp bits, l_0 exp bits]
# Polynomial: [metadata, coeff n, monomial n, coeff n-1, monomial n-1, ..., coeff 0, monomial 0]
# Metadata: x_i <-> l_i bijection
# local to global: bitdict; Vector( global id ) ->
# global to local: hash map? searchsortedfirst the local to global?
# searchsorted is log2(2^L) = L; given that L is bounded, likely faster
# than hashmap (and takes up less space, as we reuse same mem)
# x_i < x_j <=> l_i < l_j
# (L + 1) * (E + 1)
# cld(ans, 64) # implement via + 63 and then >> 6
# align to 8 UInt64s
var_per_UInt64(::Val{E}) where {E} = 64 ÷ (new_E(Val(E)) + 1)
function calc_K(::Val{L},::Val{E}) where {L,E}
    @assert E <= 63
    # step 1: figure out how many E+1 we can fit inside a UInt64
    needed = cld((L + 1), var_per_UInt64(Val(E)))
    if needed <= 2
        needed
    elseif needed <= 4
        4
    else
        (needed + 7) & -8
    end
end
new_E(::Val{E}) where {E} = max(Base._nextpow2(1+E)-1, 7)

"""
    PackedMonomial{L,E}

Bit packed monomial with maximum of L variables and E bits of exponents.
"""
struct PackedMonomial{L,E,K} <: AbstractMonomial
    bits::NTuple{K,UInt64}
    function PackedMonomial{L,E}() where {L,E}
        EN = new_E(Val(E))
        K = calc_K(Val(L),Val(EN))
        PackedMonomial{L,EN,K}()
    end
    function PackedMonomial{L,E,K}() where {L,E,K}
        PackedMonomial{L,E,K}(ntuple(Ret(zero(UInt64)), Val(K)))
    end
    PackedMonomial{L,E}(bits::NTuple{K,UInt64}) where {L,E,K} = new{L,new_E(Val(E)),K}(bits)
    PackedMonomial{L,E,K}(bits::NTuple{K,UInt64}) where {L,E,K} = new{L,new_E(Val(E)),K}(bits)
end
Base.copy(m::PackedMonomial) = m
nvariables(x::PackedMonomial{L}) where L = L
function MultivariatePolynomials.exponents(m::PackedMonomial{L,E,K}) where {L,E,K}
    tup = m.bits
    k = 0
    vpu = var_per_UInt64(Val(E))
    es = ntuple(Ret(0), L)
    isone(m) && (return es)
    for i in eachindex(tup)
        int = tup[i] << ((i == 1) * (E+1))
        for j in 1:vpu - (i == 1)
            exponent = int >> ((E+1)*(vpu-1))
            if exponent > 0
                es = Base.setindex(es, Int(exponent), k+1)
            end
            k += 1
            int <<= (E+1)
            L == k && break
        end
    end
    return es
end

PackedMonomial{L,E,K}(i::Integer) where {L,E,K} = PackedMonomial{L,E}(i)
function PackedMonomial{L,OE}(i::Integer) where {L,OE} # x_i
    @assert i < L
    E = new_E(Val(OE))
    t = PackedMonomial{L,E}().bits
    vpu = var_per_UInt64(Val(E))
    d, r = divrem(i+1, vpu)
    # o = [ 1, 0, 0, 0]
    o = one(UInt64) << ((E+1) * (vpu-1)) # `1` in `0`th typewriter (most significant) chunk
    b = o >> (r*(E+1)) # sets the bit for `i`th variable
    if d == 0 # degree and exponent on `i` are in same `UInt64`
        b |= o # set the total degree bit in same `UInt64` as `i`th monomial
    else
        t = Base.setindex(t, o, 1) # set the tot deg bit in first `UInt64`
    end
    t = Base.setindex(t, b, d+1)
    PackedMonomial{L,E}(t)
end

function degree(p::PackedMonomial{L,E}) where {L,E}
    (p.bits[1] >> ((E+1) * (var_per_UInt64(Val(E))-1))) % Int
end

function uint_type(::Val{E}) where {E}
    if E < 8
        return UInt8
    elseif E < 16
        return UInt16
    elseif E < 32
        return UInt32
    else
        return UInt64
    end
end
using SIMD

function calc_degree(p::PackedMonomial{L,E,K}) where {L,E,K}
    U = uint_type(Val(E))
    vpu = var_per_UInt64(Val(E))
    KS = vpu * K
    vp = reinterpret(Vec{KS,U},Vec{K,UInt64}(p.bits))
    sum(p) - p[vpu]
end
# we have 4 slots
# [ , , , ]
# we only need two of them: 1 for total order, 1 for a chunk
# var_per_UInt64 == 4
# L = 2
#
# [total degree, chunk, chunk, __]
# [__, total_degree, chunk, chunk] >> E+1
# Iterate right to left:
# [ a, b, c, d] mask top 3 => d; iter by shift right 1
# [ 0, a, b, c] mask top 3 => c; iter by shift right 1
# [ 0, 0, a, b] mask top 3 => b; iter by shift right 1
# [ 0, 0, 0, a] mask top 3 => a; iter by shift right 1
# number == 0 -> terminate
# Iterate left to right (typewriter)
# [ a, b, c, d] shift by 3 => a; iter by shift left 1
# [ b, c, d, 0] shift by 3 => b; iter by shift left 1
# [ c, d, 0, 0] shift by 3 => c; iter by shift left 1
# [ d, 0, 0, 0] shift by 3 => d; iter by shift left 1
# number == 0 -> terminate
function Base.gcd(p::T, q::T) where {L,E,K,T<:PackedMonomial{L,E,K}}
    U = uint_type(Val(E))
    vpu = var_per_UInt64(Val(E))
    KS = vpu * K
    vp = reinterpret(Vec{KS,U},Vec{K,UInt64}(p.bits))
    vq = reinterpret(Vec{KS,U},Vec{K,UInt64}(q.bits))
    g = min(vp, vq)
    total_degree = sum(g) - g[vpu]
    g = Base.setindex(g, total_degree, vpu)
    T(Tuple(reinterpret(Vec{K,UInt64}, g)))
end

# the first chunk is the total degree
Base.isone(x::PackedMonomial) = iszero(first(x.bits))
function MultivariatePolynomials.grlex(x::T, y::T) where {T<:PackedMonomial}
    x.bits < y.bits ? -1 : (x == y ? 0 : 1)
end
Base.isless(x::T, y::T) where {T<:PackedMonomial} = x.bits < y.bits
@noinline overflowed_error() = throw(Base.OverflowError("Try increasing E."))
@inline function check_zero_mask(::Val{E}) where {E}
    shifted = E
    x = one(UInt64) << E
    while shifted < 64
        x <<= 1
        x |= one(UInt64)
        x <<= E
        shifted += E
    end
    return x
end
@inline function zero_bits(x::UInt64, ::Val{E}) where {E}
    # every E+1st bit should be zero, so we create a mask, and then check for
    # zero
    (check_zero_mask(Val(E)) & x)
end
@inline function check_zero_bits(x::UInt64, ::Val{E}) where {E}
    # every E+1st bit should be zero, so we create a mask, and then check for
    # zero
    zero_bits(x, Val(E)) != zero(UInt64)
end

function Base.:*(x::T, y::T) where {L,E,T<:PackedMonomial{L,E}}
    xys = _fmap(+, x.bits, y.bits)
    o = reduce_tup(|, _fmap(Base.Fix2(zero_bits, Val(E)), xys))

    o != zero(UInt64) && overflowed_error()
    return T(xys)
end

function MultivariatePolynomials.divides(y::T, x::T) where {L,E,T<:PackedMonomial{L,E}}
    xys = _fmap(-, x.bits, y.bits)
    o = reduce_tup(|, _fmap(Base.Fix2(zero_bits, Val(E)), xys))

    #o != zero(UInt64) && overflowed_error()
    #return T(xys), o != zero(UInt64)
    return o == zero(UInt64)
end

MultivariatePolynomials.constantmonomial(::M) where {M<:PackedMonomial} = constantmonomial(M)
function MultivariatePolynomials.constantmonomial(M::Type{<:PackedMonomial})
    M()
end

# TODO: make it fast
function MultivariatePolynomials.mapexponents!(f::F, m1::M, m2::M) where {F<:Function, L, E, M<:PackedMonomial{L,E}}
    MultivariatePolynomials.mapexponents(f, m1, m2)
end
function MultivariatePolynomials.mapexponents(f::F, m1::M, m2::M) where {F<:Function, L, E, M<:PackedMonomial{L,E}}
    e1 = exponents(m1)
    e2 = exponents(m2)
    ne = map(f, e1, e2)
    prod(((i, e),)->PackedMonomial{L,E}(i-1)^e, enumerate(ne))
end

function MultivariatePolynomials._div(x::T, y::T) where {L,E,T<:PackedMonomial{L,E}}
    xys = _fmap(-, x.bits, y.bits)
    o = reduce_tup(|, _fmap(Base.Fix2(zero_bits, Val(E)), xys))

    #o != zero(UInt64) && overflowed_error()
    #return T(xys), o != zero(UInt64)
    return T(xys)#, o != zero(UInt64)
end
function Base.:^(x::T, y::Integer) where {L,E,T<:PackedMonomial{L,E}}
    xys = _fmap(*, x.bits, y)
    o = reduce_tup(|, _fmap(Base.Fix2(zero_bits, Val(E)), xys))

    o != zero(UInt64) && overflowed_error()
    return T(xys)
end
zero_nondegreemask(::Val{7 }) = 0xff00000000000000
zero_nondegreemask(::Val{15}) = 0xffff000000000000
zero_nondegreemask(::Val{31}) = 0xffffffff00000000
zero_nondegreemask(::Val{63}) = 0xffffffffffffffff
function firstid(m::PackedMonomial{L,E,K}) where {L,E,K}
    b = m.bits[1] & (~zero_nondegreemask(Val(E)))
    if (b != zero(b))
        return (leading_zeros(b) ÷ (E+1)) - 1
    end
    vpu = var_per_UInt64(Val(E))
    b = vpu - 1
    for k in 2:K
        bk = m.bits[k]
        if (bk != zero(bk))
            return (leading_zeros(bk) ÷ (E+1)) + b
        end
        b += vpu
    end
    error("firstTermID should only be called if degree > 0.");
    return 0
end

function degree(m::PackedMonomial{L,E,K}, id) where {L,E,K}
    vpu = var_per_UInt64(Val(E))
    d = K == 1 ? 0 : (id + 1) ÷ vpu
    r = (id + 1) - (d * vpu)
    b = m.bits[d+1] << (r * (E + 1))
    return b >> ((E + 1) * (vpu - 1))
end

function rmid(x::T, id) where {L,E,K,T<:PackedMonomial{L,E,K}}
    bits = x.bits
    vpu = var_per_UInt64(Val(E))
    d = K == 1 ? 0 : (id + 1) ÷ vpu
    r = (id + 1) - (d * vpu)
    m = zero_nondegreemask(Val(E)) >> (r*(E+1))
    oldBits = bits[d+1]
    b = oldBits & (~m)
    remDegree = (oldBits & m) >> ((E+1)*(vpu-1-r))

    o = remDegree << ((E + 1) * (vpu - 1))
    if d != zero(d)
        bits = Base.setindex(bits, b, d+1)
        bits = Base.setindex(bits, bits[1] - o, 1)
    else
        bits = Base.setindex(bits, b - o, 1)
    end
    return T(bits)
end

function Base.show(io::IO, m::PackedMonomial{L,E}) where {L,E}
    tup = m.bits
    k = 0
    vpu = var_per_UInt64(Val(E))
    isone(m) && (print(io, '1'); return)
    for i in eachindex(tup)
        int = tup[i] << ((i == 1) * (E+1))
        for j in 1:vpu - (i == 1)
            exponent = int >> ((E+1)*(vpu-1))
            if exponent > 0
                print(io, 'x')
                print(io, int2lowerscript(k))
                exponent > 1 && print(io, int2superscript(exponent))
            end
            k += 1
            int <<= (E+1)
            L == k && break
        end
    end
end

struct Variable{L,E} <: AbstractVariable
    id::UInt
end

function Base.:(==)(v1::V, v2::V) where {V<:Variable}
    v1 === v2
end
function Base.isless(v1::V, v2::V) where {V<:Variable}
    isless(v1.id, v2.id)
end
function MultivariatePolynomials.name_base_indices(v::Variable)
    "x", (v.id,)
end
MultivariatePolynomials.name(v::Variable) = string("x", int2lowerscript(v.id))
MultivariatePolynomials.variable_union_type(::Type{<:PackedMonomial{L,E}}) where {L,E} = Variable{L,E}

MultivariatePolynomials.variables(::Type{<:PackedMonomial{L,E}}) where {L, E} = ntuple(i->Variable{L,E}(i-1), Val(L))
MultivariatePolynomials.variables(::M) where {M<:PackedMonomial} = variables(M)

MultivariatePolynomials.termtype(M::Type{<:PackedMonomial}, ::Type{T}) where {T} = Term{T, M}
MultivariatePolynomials.polynomialtype(::Type{Term{T, M}}) where {T, M<:PackedMonomial} = Polynomial{T, Term{T, M}, Vector{Term{T, M}}}

MultivariatePolynomials.nvariables(p::Polynomial{T, Term{T, M}}) where {T, L, M<:PackedMonomial{L}} = L
MultivariatePolynomials.variables(p::Polynomial{T, Term{T, M}}) where {T, M<:PackedMonomial} = variables(M)
Base.:(^)(x::Variable{L,E}, p::Integer) where {L,E} = PackedMonomial{L,E}(x.id) ^ p
function Base.:(==)(m1::T, m2::T) where {T<:PackedMonomial}
    m1 === m2
end

function MultivariatePolynomials.substitute(::MultivariatePolynomials.Subs, v::V, p::Pair{V, Int64}) where {L,E,V<:Variable{L, E}}
    v == p[1] ? p[2] : v
end
