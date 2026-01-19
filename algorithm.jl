function eig_KxK_diagblocks(K::Int, n::Int, matrix::AbstractMatrix)
    # hold eigenvalues
    eigs = zeros(ComplexF64, K * n)
    # hold eigenvectors
    vecs = zeros(ComplexF64, K * n, K * n)
    # column counter
    col = 1 

    # for each l = 1, ..., n
    for l in 1:n
        M_l = zeros(ComplexF64, K, K)
        for i in 1:K
            for j in 1:K
                M_l[i, j] = matrix[(i - 1) * n + l, (j - 1) * n + l]
            end
        end

        # get eigenvalues and eigenvectors
        E = eigen(M_l)
        eigvals = E.values
        eigvecs = E.vectors

        # store in eigs and vecs
        eigs[col:col + K - 1] = eigvals
        V = zeros(ComplexF64, K * n, K)
        for i in 1:K
            for j in 1:K
                V[(i - 1) * n + l, j] = eigvecs[i, j]
            end
        end
        vecs[:, col:col + K - 1] = V
        col += K
    end

    return eigs, vecs
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Example usage
    A = ComplexF64[
        1 0 3 0;
        0 2 0 4;
        5 0 7 0;
        0 6 0 8
    ]

    eigs, vecs = eig_KxK_diagblocks(2, 2, A)

    println("****KxK function results****\n")
    println("eigenvalues = $eigs")
    println("\neigenvectors = $vecs\n\n")

    println("****Regular Julia****\n")
    E = eigen(A)
    println("eigenvalues = $(E.values)")
    println("\neigenvectors = $(E.vectors)")
end