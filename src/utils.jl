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
struct GetWithShift3{F,T,S}
    f::F
    tup::T
    s::S
end
(g::GetWithShift3)(i) = g.f(g.tup[i], g.s)

function _fmap(f::F, x::Tuple{Vararg{Any,K}}) where {K,F}
    ntuple(GetWithShift(f, x, 0), Val(K))
end
function _fmap(f::F, x::Tuple{Vararg{Any,K}}, y::Tuple{Vararg{Any,K}}) where {K,F}
    ntuple(GetWithShift2(f, x, y), Val(K))
end
function _fmap(f::F, x::Tuple{Vararg{Any,K}}, y) where {K,F}
    ntuple(GetWithShift3(f, x, y), Val(K))
end
