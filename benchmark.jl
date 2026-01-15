using LinearAlgebra
#using BenchmarkTools

function time_block_method(K::Int, n::Int, M::AbstractMatrix)
    eig_KxK_diagblocks(K, n, M)  # warmup (JIT compiled language)
    t = @elapsed begin 
        eigs, vecs = eig_KxK_diagblocks(K, n, M)
    end

    return  (
        eigenvalues = eigs,
        eigenvectors = vecs,
        time = t
    )
end

function time_dense_eig(M::AbstractMatrix)
    eigen(M)
    t = @elapsed begin 
        E = eigen(M)
    end

    return (
        eigenvalues = E.values,
        eigenvectors = E.vectors,
        time = t
    )
end