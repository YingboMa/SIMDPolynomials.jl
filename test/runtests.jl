using LoopPoly
using Test

@testset "pseudorem" begin
    x = Uninomial(0, 1)
    p = 2x^10 + x^7 + 7*x^2 + x + 3x
    q = (p+one(p)) * (p+2one(p)) * (p+3one(p))
    @test pseudorem(q, p) == 12582912*one(p)
    q = x^7 + 20*one(x)
    @test pseudorem(q, p) == q
    @test pseudorem(p, q) == -40*x^3 + 7*x^2 + 4*x - 20*one(x)
    q = x^6 + 23*one(x)
    @test pseudorem(p, q) == -46*x^4 + 7*x^2 - 19*x
end

@testset "MPoly2Poly" begin
    x, y, z = [Monomial([i]) for i in 0:2]
    c1 = x * z + x
    c2 = x^2 + z
    c3 = 2 - z
    c4 = x * z^2 + x
    e1 = 0
    e2 = 5
    e3 = 7
    e4 = 10
    p = c1 * y^e1 + c2 * y^e2 + c3 * y^e3 + c4 * y^e4
    pp = LoopPoly.SparsePoly(p, y.ids[1]);
    @test var(pp) == y.ids[1]
    @test coeffs(pp) == [c4, c3, c2, c1]
    @test pp.exps == [e4, e3, e2, e1]
end
