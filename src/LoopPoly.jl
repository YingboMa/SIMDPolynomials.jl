module LoopPoly

export Uninomial, Uniterm, Poly
export Monomial, Term, MPoly
export terms, term, coeffs, coeff, lt, var, content, univariate_gcd
export mpoly2poly
export univariate_gcd, pseudorem
export PackedMonomial

include("mpoly.jl")
include("packedmonomial.jl")
include("poly.jl")

end
