using svd_IceSheetDEM:rsvd
using Test, TSVD, LinearAlgebra

@testset "test rsvd with m = $m, n = $n" for
    (m, n) in ((100, 100),
               (1000, 150),
               (50, 900))
    A = randn(m,n)
    U, S, V = svd(A)
    # for k in [2]
        Ur, Sr, Vr = rsvd(A, min(m,n))
        @test norm(Sr - S) < sqrt(eps(eltype(A)))
        @test norm(abs.(transpose(Ur)*U) - I) < sqrt(eps(eltype(A)))
        @test norm(abs.(transpose(Vr)*V) - I) < sqrt(eps(eltype(A)))
    # end
end
@testset "truncated rsvd" begin
    A = [1. 2 3; 4 5 6; 7 8 9]
    U, S, V = svd(A)
    for k = 2:3
        Ur, Sr, Vr = rsvd(A, k)
        @test norm(Sr - S[1:k]) < sqrt(eps(eltype(A)))
        @test norm(abs.(transpose(Ur)*U[:,1:k]) - I) < sqrt(eps(eltype(A)))
        @test norm(abs.(transpose(Vr)*V[:,1:k]) - I) < sqrt(eps(eltype(A)))
    end
end
