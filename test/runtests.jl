using LoopPoly
using Test

LoopPoly.debugmode() = true

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

@testset "GCD" begin
    x, y, z = [Monomial([i]) for i in 0:2]
    c1 = 10*(x * z + x)
    c2 = 2*(x^2 + z)
    c3 = 2*(2 - z  )
    c4 = 20*(x * z^2)
    e1 = 0
    e2 = 5
    e3 = 7
    e4 = 10
    p = c1 * y^e1 + c2 * y^e2 + c3 * y^e3 + c4 * y^e4
    pp = LoopPoly.SparsePoly(p, y.ids[1])
    @test var(pp) == y.ids[1]
    @test coeffs(pp) == [c4, c3, c2, c1]
    @test pp.exps == [e4, e3, e2, e1]
    q = prod(i->p + i, 0:3)
    for i in 0:3
        @test gcd(p + i, q) == p + i
    end

    k = y^2 + 1
    @test gcd(x*k, z*k) == k
    @test gcd(z*k, x*k) == k
    @test gcd(x*k, (z+1)*k) == k
    @test gcd((z+1)*k, x*k) == k
    @test gcd((z+1)*k, x*k) == k
    @test gcd(x*k, p*k) == k
    @test gcd(p*k, x*k) == k

    p = 2*x*y
    q = 2*x*y + x
    @test gcd(p, q) == x
    @test gcd(q, p) == x

    p = x*y + y
    q = -x*y - y
    @test gcd(p, q) == p
    @test gcd(q, p) == p
end
