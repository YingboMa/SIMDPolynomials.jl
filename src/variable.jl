abstract type AbstractVariable <: MP.AbstractVariable end

function Base.:(==)(v1::V, v2::V) where {V<:AbstractVariable}
    v1 === v2
end
function Base.isless(v1::V, v2::V) where {V<:AbstractVariable}
    isless(v2.id, v1.id)
end
function MP.name_base_indices(v::AbstractVariable)
    "x", (v.id,)
end
MP.name(v::AbstractVariable) = string("x", MP.unicode_subscript(v.id))
function MP.substitute(::MP.Subs, v::V, p::Pair{V,Int64}) where {V<:AbstractVariable}
    v == p[1] ? p[2] : v
end
