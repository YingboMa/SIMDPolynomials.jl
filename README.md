# SIMDPolynomials

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://YingboMa.github.io/SIMDPolynomials.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://YingboMa.github.io/SIMDPolynomials.jl/dev)
[![Build Status](https://github.com/YingboMa/SIMDPolynomials.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/YingboMa/SIMDPolynomials.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/YingboMa/SIMDPolynomials.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/YingboMa/SIMDPolynomials.jl)

---

SIMDPolynomials.jl uses bit-packed monomials so that most of the operations on
multivariate monomials can be done in a few CPU instructions. Currently, it is
only optimized for relatively small polynomials. Contributions, especially
optimizations for large polynomials, are welcome!

Examples:
```julia
julia> using SIMDPolynomials, BenchmarkTools

julia> x, y, z, t = [PackedMonomial{4,7}(i) for i in 0:3]; # PackedMonomial with maximum of 4 variables and 7 bits of exponents.

julia> p = x * y + 3 * (z * t)
x₀x₁ + 3x₂x₃

julia> q = (p + 1) * p
x₀²x₁² + 6x₀x₁x₂x₃ + 9x₂²x₃² + x₀x₁ + 3x₂x₃

julia> @btime gcd($p, $q)
  4.019 μs (94 allocations: 8.06 KiB)
x₀x₁ + 3x₂x₃

julia> begin
           c1 = 10*(x * z + x)
           c2 = 2*(x^2 + z)
           c3 = 2*(2 - z  )
           c4 = 20*(x * z^2)
           e1 = 0
           e2 = 5
           e3 = 7
           e4 = 10
           p = c1 * y^e1 + c2 * y^e2 + c3 * y^e3 + c4 * y^e4
           q = prod(i->p + i, 0:3)
       end;

julia> @btime for i in 0:3
           gcd($p + i, $q)
       end
  350.086 μs (1159 allocations: 588.06 KiB)
```

The same micro-benchmark using AbstractAlgebra:
```julia
julia> using AbstractAlgebra, BenchmarkTools

julia> R, (x, y, z, t) = PolynomialRing(AbstractAlgebra.Integers{Int}(), [:x, :y, :z, :t], ordering=:deglex);

julia> p = x * y + 3 * (z * t)
x*y + 3*z*t

julia> q = (p + 1) * p
x^2*y^2 + 6*x*y*z*t + 9*z^2*t^2 + x*y + 3*z*t

julia> @btime gcd($p, $q) # SIMDPolynomials.jl is 30x faster
  119.795 μs (1320 allocations: 89.17 KiB)
x*y + 3*z*t

julia> begin
           c1 = 10*(x * z + x)
           c2 = 2*(x^2 + z)
           c3 = 2*(2 - z  )
           c4 = 20*(x * z^2)
           e1 = 0
           e2 = 5
           e3 = 7
           e4 = 10
           p = c1 * y^e1 + c2 * y^e2 + c3 * y^e3 + c4 * y^e4
           q = prod(i->p + i, 0:3)
       end;

julia> @btime for i in 0:3 # SIMDPolynomials.jl is 14x faster
           gcd($p + i, $q)
       end
  4.934 ms (32235 allocations: 3.43 MiB)
```

The same micro-benchmark using DynamicPolynomials:
```julia
julia> using DynamicPolynomials, BenchmarkTools

julia> @polyvar x y z t;

julia> p = x * y + 3 * (z * t)
xy + 3zt

julia> q = (p + 1) * p
x²y² + 6xyzt + 9z²t² + xy + 3zt

julia> @btime gcd($p, $q)  # SIMDPolynomials.jl is 65x faster
  264.561 μs (4962 allocations: 298.19 KiB)
xy + 3zt

julia> begin
           c1 = 10*(x * z + x)
           c2 = 2*(x^2 + z)
           c3 = 2*(2 - z  )
           c4 = 20*(x * z^2)
           e1 = 0
           e2 = 5
           e3 = 7
           e4 = 10
           p = c1 * y^e1 + c2 * y^e2 + c3 * y^e3 + c4 * y^e4
           q = prod(i->p + i, 0:3)
       end;

julia> @btime for i in 0:3 # SIMDPolynomials.jl is 82x faster
           gcd($p + i, $q)
       end
  28.943 ms (529642 allocations: 31.20 MiB)
```
