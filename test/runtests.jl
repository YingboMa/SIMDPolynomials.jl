using SIMDPolynomials
#using SIMDPolynomials: var
import MultivariatePolynomials as MP
using Test

SIMDPolynomials.debugmode() = true

#@testset "Fuzz SparsePoly" begin
#    x = Uninomial(0, 1)
#    for i in 1:10000
#        p = sum(rand(-2:2)*x^i for i in 0:10)
#        q = sum(rand(-2:2)*x^i for i in 0:5)
#        try
#            g = gcd(p, q)
#            pdg = div(p, g)
#            qdg = div(q, g)
#            @test gcd(pdg, qdg) == one(x)
#            pq = p * q
#            @test div(pq, p) == q
#            @test div(pq, q) == p
#        catch
#            display(p)
#            display(q)
#            rethrow()
#        end
#    end
#end

@testset "Fuzz MPoly" begin
    for monos in ([PackedMonomial{4,8}(i) for i in 0:2], [Monomial([i]) for i in 0:2])
        x, y, z = monos
        for _ in 1:10000
            pterms = [rand(-2:2)*prod(rand([x, y, z])^i for i in 0:3) for j in 0:7]
            p = sum(pterms)
            q = sum(rand(-2:2)*prod(rand([x, y, z])^i for i in 0:2)  for j in 0:5)
            try
                g = gcd(p, q)
                pdg = div(p, g)
                qdg = div(q, g)
                @test gcd(pdg, qdg) == one(x)
            catch
                display(p)
                display(q)
                rethrow()
            end
        end
    end
end

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

function test_gcd(x, y)
    g1 = gcd(x, y)
    g2 = gcd(y, x)
    @test sign(MP.leading_coefficient(g1)) * g1 == sign(MP.leading_coefficient(g2)) * g2
    return g1
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
    pp = SIMDPolynomials.SparsePoly(p, y.ids[1])
    #@test var(pp) == y.ids[1]
    @test coeffs(pp) == [c4, c3, c2, c1]
    @test pp.exps == [e4, e3, e2, e1]
    q = prod(i->p + i, 0:3)
    @test length(terms(q)) == 262
    for i in 0:3
        @test test_gcd(p + i, q) == p + i
    end

    k = y^2 + 1
    @test test_gcd(x*k, z*k) == k
    @test test_gcd(x*k, (z+1)*k) == k
    @test test_gcd(x*k, p*k) == k

    p = 2*x*y
    q = 2*x*y + x
    @test test_gcd(q, p) == x

    p = x*y + y
    q = -x*y - y
    @test test_gcd(q, p) == p

    @test test_gcd((-x + 1) * (y^2+1), -x+1) == -x + 1
end

@testset "PackedMonomial" begin
    x, y, z, t = [PackedMonomial{4,8}(i) for i in 0:3]
    m = (x^10 * y^3 * z^2 * t^4)^6
    n = (x^9 * y^8 * z^3 * t^3)^6
    g = gcd(m, n)
    @test g == x^54 * y^18 *z^12 * t^18
    @test SIMDPolynomials.degree(g) == 54 + 18 + 12 + 18
    x, y, z, t = [PackedMonomial{4,7}(i) for i in 0:3]
    m = (x^10 * y^3 * z^2 * t^4)^6
    @test_throws Base.OverflowError (x^9 * y^8 * z^3 * t^3)^6
    g = gcd(m, m)
    @test g == m
    @test SIMDPolynomials.degree(g) == 60 + 18 + 12 + 24

    c1 = 10*(x * z + x)
    c2 = 2*(x^2 + z)
    c3 = 2*(2 - z  )
    c4 = 20*(x * z^2)
    e1 = 0
    e2 = 5
    e3 = 7
    e4 = 10
    p = c1 * y^e1 + c2 * y^e2 + c3 * y^e3 + c4 * y^e4
    q = prod(i->p + i, 0:3);
    @test length(terms(q)) == 262
    for i in 0:3
        @test test_gcd(p + i, q) == p + i
    end

    k = y^2 + 1
    @test test_gcd(x*k, z*k) == k
    @test test_gcd(x*k, (z+1)*k) == k
    @test test_gcd(x*k, p*k) == k

    p = 2*x*y
    q = 2*x*y + x
    @test test_gcd(q, p) == x

    p = x*y + y
    q = -x*y - y
    @test test_gcd(q, p) == p

    @test test_gcd((-x + 1) * (y^2+1), -x+1) == -x + 1

    p1 = 1 + x + y + z + t
    n = 15
    p = p1^n;
    @test test_gcd(p, p1) == p1
    k = p;
    for i in 1:n
        k = div(k, p1)
    end
    @test k == 1
    q = (p + 1) * p;
    @test gcd(q, p) == p
end
