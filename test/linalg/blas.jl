# This file is a part of Julia. License is MIT: http://julialang.org/license

using Base.Test

@testset "BLAS" begin

    @testset "syr for eltype $elty" for elty in (Float32, Float64, Complex{Float32}, Complex{Float64})
        A = rand(elty, 5, 5)
        @test triu(A[1,:] * A[1,:].') ≈ BLAS.syr!('U', one(elty), A[1,:], zeros(elty, 5, 5))
        @test tril(A[1,:] * A[1,:].') ≈ BLAS.syr!('L', one(elty), A[1,:], zeros(elty, 5, 5))
        @test triu(A[1,:] * A[1,:].') ≈ BLAS.syr!('U', one(elty), view(A, 1, :), zeros(elty, 5, 5))
        @test tril(A[1,:] * A[1,:].') ≈ BLAS.syr!('L', one(elty), view(A, 1, :), zeros(elty, 5, 5))
    end

    @testset "her for eltype $elty" for elty in (Complex{Float32}, Complex{Float64})
        @test triu(A[1,:] * A[1,:]') ≈ BLAS.her!('U', one(real(elty)), A[1,:], zeros(elty, 5, 5))
        @test tril(A[1,:] * A[1,:]') ≈ BLAS.her!('L', one(real(elty)), A[1,:], zeros(elty, 5, 5))
        @test triu(A[1,:] * A[1,:]') ≈ BLAS.her!('U', one(real(elty)), view(A, 1, :), zeros(elty, 5, 5))
        @test tril(A[1,:] * A[1,:]') ≈ BLAS.her!('L', one(real(elty)), view(A, 1, :), zeros(elty, 5, 5))
    end
end
