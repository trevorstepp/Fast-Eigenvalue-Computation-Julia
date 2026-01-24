const K_LARGEST = 10  # number of largest eigenvalues to compare

function verify_results(
    matrix::AbstractMatrix,
    alg_eigs::AbstractVector,
    alg_vecs::AbstractMatrix,
    reg_eigs::AbstractVector,
    reg_vecs::AbstractMatrix,
    atol::Real = 1e-8
)

    # sort eigenvalues
    alg_dom = dominant_eigs(alg_eigs, K_LARGEST)
    reg_dom = dominant_eigs(reg_eigs, K_LARGEST)

    # check equivalence
    eig_check = isapprox(alg_dom, reg_dom; atol=atol, rtol=0)

    # keep track of maximum (worst) residual and average residual
    max_res = 0.0
    res_sum = 0.0
    len = length(alg_eigs)

    # compute residual for each eigenvalue, eigenvector pair
    for i in eachindex(alg_eigs)
        eigval = alg_eigs[i]
        eigvec = alg_vecs[:, i]
        r = norm(matrix * eigvec - eigval * eigvec)
        res_sum += r
        max_res = max(max_res, r)
    end

    return (
        num_eigenvalues_compared = K_LARGEST,
        eigenvalues_match = eig_check,
        max_residual = max_res,
        mean_residual = res_sum / len
    )
end

# helper function only used in this file
function dominant_eigs(eigs::AbstractVector, k_largest::Int)
    index_arr = sortperm(abs.(eigs), rev=true)
    return eigs[index_arr][1:k_largest]
end