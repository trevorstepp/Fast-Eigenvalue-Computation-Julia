function build_block_matrix(K::Int, n::Int; seed::Union{Int, Nothing}=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    # allocate matrix
    M = zeros(ComplexF64, K * n, K * n)

    # loop over diagonal positions inside each n x n block
    for l in 1:n
        B_l = randn(K, K)
        
        # add to each block using B_l
        for i in 1:K 
            for j in 1:K 
                M[(i - 1) * n + l, (j - 1) * n + l] = B_l[i, j]
            end
        end
    end

    return M
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Random
    using LinearAlgebra
    println(build_block_matrix(2, 2; seed=0))
end