module SIMDPolynomials

import MutableArithmetics as MA
import MultivariatePolynomials as MP
const Term = MP.Term
const Polynomial = MP.Polynomial
export Uninomial, Uniterm, Poly
export PackedMonomial, Monomial, Term, MPoly
export terms, coeffs, content, contprim
export pseudorem

struct Ret{V} <: Function
    value::V
    Ret{V}(value) where {V} = new{V}(value)
    Ret(value) = new{Core.Typeof(value)}(value)
end

(obj::Ret)(args...; kw...) = obj.value

debugmode() = false

include("variable.jl")
include("utils.jl")
include("monomial.jl")
include("packedmonomial.jl")

end
