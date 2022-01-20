module LoopPoly

export Uninomial, Uniterm, Poly
export Monomial, Term, MPoly
export terms, term, coeffs, coeff, lt, var, content, univariate_gcd
export mpoly2poly
export univariate_gcd, pseudorem
export PackedMonomial

struct Ret{V} <: Function
    value::V
    Ret{V}(value) where {V} = new{V}(value)
    Ret(value) = new{Core.Typeof(value)}(value)
end

(obj::Ret)(args...; kw...) = obj.value

debugmode() = false

include("interface.jl")
include("utils.jl")
include("monomial.jl")
include("packedmonomial.jl")
include("mpoly.jl")
include("poly.jl")

end
