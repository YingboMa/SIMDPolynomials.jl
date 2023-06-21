const IDType = UInt32
const NOT_A_VAR = typemax(IDType)
const EMPTY_IDS = IDType[]

struct Variable <: AbstractVariable
    id::IDType
end

function Base.:(^)(x::Variable, p::Integer)
    return Monomial([x.id for _ in 1:p])
end

function MP.monomial_type(::Type{Variable})
    return Monomial
end

struct Monomial <: MP.AbstractMonomial
    ids::Vector{IDType}
end
Monomial() = Monomial(EMPTY_IDS)
MP.constant_monomial(::Union{Monomial,Type{Monomial}}) = Monomial()
degree(x::Monomial) = length(x.ids)
Base.copy(x::Monomial) = Monomial(copy(x.ids))
#TODO: optimize
function nvariables(x::Monomial)
    nvar = 0
    lastvar::Int = -1
    for id in x.ids
        if id != lastvar
            lastvar = id
            nvar += 1
        end
    end
    return nvar
end

firstid(m::Monomial) = m.ids[1]
degree(m::Monomial, id) = count(isequal(id), m.ids)
rmid(m::Monomial, id) = Monomial(filter(!isequal(id), m.ids))

function _p(y::Monomial, x::Monomial)
    ids = y.ids
    i = j = 1
    n0 = length(ids)
    n1 = length(x.ids)
    #r = Monomial(similar(ids, n0 + n1))
    r = Monomial(zeros(IDType, n0 + n1))
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
function MP.compare(x::Monomial, y::Monomial)
    if x < y
        return -1
    elseif x == y
        return 0
    else
        return 1
    end
end
MP.isconstant(x::Monomial) = isempty(x.ids)
MP.variables(x::Monomial) = Variable.(unique(x.ids))
function MP.variables(p::Polynomial{T, Term{T, Monomial}}) where {T}
    v = Variable[]
    for t in MP.terms(p)
        append!(v, MP.variables(MP.monomial(t)))
    end
    sort!(v, rev = true)
    unique!(v)
    return v
end
function MP.exponents(x::Monomial)
    exps = Int[]
    cur_id = nothing
    for i in x.ids
        if !isnothing(cur_id) && cur_id == i
            exps[end] += 1
        else
            cur_id = i
            push!(exps, 1)
        end
    end
    return exps
end

function _degree(m::Monomial, offset)
    k = 1
    while (offset + k) in eachindex(m.ids) && m.ids[offset] == m.ids[offset + k]
        k += 1
    end
    return k
end

function MP.map_exponents(op, m1::Monomial, m2::Monomial)
    g = Monomial(similar(m1.ids, 0))
    i = firstindex(m1.ids)
    j = firstindex(m2.ids)
    while i in eachindex(m1.ids) || j in eachindex(m2.ids)
        if i in eachindex(m1.ids) && j in eachindex(m2.ids) && m1.ids[i] == m2.ids[j]
            v = m1.ids[i]
            di = _degree(m1, i)
            i += di
            dj = _degree(m2, j)
            j += dj
        elseif i in eachindex(m1.ids) && (!(j in eachindex(m2.ids)) || (m1.ids[i] < m2.ids[j]))
            v = m1.ids[i]
            di = _degree(m1, i)
            dj = 0
            i += di
        else
            v = m2.ids[j]
            di = 0
            dj = _degree(m2, j)
            j += dj
        end
        for _ in 1:op(di, dj)
            push!(g.ids, v)
        end
    end
    return g
end

function MP.map_exponents!(op, m1::Monomial, m2::Monomial)
    return MP.map_exponents(op, m1, m2)
end

function MP.divides(m1::Monomial, m2::Monomial)
    i = firstindex(m1.ids)
    j = firstindex(m2.ids)
    while i in eachindex(m1.ids) || j in eachindex(m2.ids)
        if i in eachindex(m1.ids) && j in eachindex(m2.ids) && m1.ids[i] == m2.ids[j]
            v = m1.ids[i]
            di = _degree(m1, i)
            i += di
            dj = _degree(m2, j)
            j += dj
        elseif i in eachindex(m1.ids) && (!(j in eachindex(m2.ids)) || (m1.ids[i] < m2.ids[j]))
            v = m1.ids[i]
            di = _degree(m1, i)
            dj = 0
            i += di
        else
            v = m2.ids[j]
            di = 0
            dj = _degree(m2, j)
            j += dj
        end
        if di > dj
            return false
        end
    end
    return true
end


function MA.promote_operation(
    ::typeof(MP.substitute),
    ::Type{MP.Subs},
    ::Type{Monomial},
    ::Type{Pair{Variable,T}},
) where {T}
    U = MA.promote_operation(^, T, Int)
    return Term{U,Monomial}
end
