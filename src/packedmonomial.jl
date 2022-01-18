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
function calc_K(::Val{L},::Val{E}) where {L,E}
    needed = (((L+1)*(E+1)+63) >> 6)
    if needed <= 2
        needed
    elseif needed <= 4
        4
    else
        (needed + 7) & -8
    end
end
struct PackedMonomial{L,E,K}
    bits::NTuple{K,UInt64}
    function PackedMonomial{L,E}() where {L,E}
        PackedMonomial{L,E}(ntuple(Returns(zero(UInt64)), calc_K(Val(L),Val(E))))
    end
    PackedMonomial{L,E}(bits::NTuple{K,UInt64}) where {L,E,K} = new{L,E,K}(bits)
end

Base.isless(x::T, y::T) where {T<:PackedMonomial} = x.bits < y.bits
@noinline overflowed_error() = throw("Overflowed. Try increasing E.")
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

@generated function reduce_tup(f::F, inds::Tuple{Vararg{Any,N}}) where {F,N}
    q = Expr(:block, Expr(:meta, :inline, :propagate_inbounds))
    if N == 1
        push!(q.args, :(inds[1]))
        return q
    end
    syms = Vector{Symbol}(undef, N)
    i = 0
    for n ∈ 1:N
        syms[n] = iₙ = Symbol(:i_, (i += 1))
        push!(q.args, Expr(:(=), iₙ, Expr(:ref, :inds, n)))
    end
    W =  1 << (8sizeof(N) - 2 - leading_zeros(N))
    while W > 0
        _N = length(syms)
        for _ ∈ 2W:W:_N
            for w ∈ 1:W
                new_sym = Symbol(:i_, (i += 1))
                push!(q.args, Expr(:(=), new_sym, Expr(:call, :f, syms[w], syms[w+W])))
                syms[w] = new_sym
            end
            deleteat!(syms, 1+W:2W)
        end
        W >>>= 1
    end
    q
end

struct GetWithShift{F,T}
    f::F
    tup::T
    shift::Int
end
(g::GetWithShift)(i) = g.f(g.tup[i + g.shift])

struct GetWithShift2{F,T}
    f::F
    tup1::T
    tup2::T
end
(g::GetWithShift2)(i) = g.f(g.tup1[i], g.tup2[i])

function _fmap(f::F, x::Tuple{Vararg{Any,K}}) where {K,F}
    ntuple(GetWithShift(f, x, 0), Val(K))
end
function _fmap(f::F, x::Tuple{Vararg{Any,K}}, y::Tuple{Vararg{Any,K}}) where {K,F}
    ntuple(GetWithShift2(f, x, y), Val(K))
end
function Base.:*(x::T, y::T) where {L,E,T<:PackedMonomial{L,E}}
    xys = _fmap(+, x.bits, y.bits)
    o = reduce_tup(|, _fmap(Base.Fix2(zero_bits, Val(E)), xys))

    o != zero(UInt64) && overflowed_error()
    return xys
end
