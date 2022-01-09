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
