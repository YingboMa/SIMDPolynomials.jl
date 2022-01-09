module LoopPoly

export Uninomial, Uniterm, Poly
export Monomial, Term, MPoly
export terms, term, coeffs, coeff, lt, var, content
export mpoly2poly
export univariate_gcd, pseudorem

include("mpoly.jl")
include("poly.jl")

end
