function verify_results(
    matrix::AbstractMatrix,
    alg_eigs::AbstractVector,
    alg_vecs::AbstractMatrix,
    reg_eigs::AbstractVector,
    reg_vecs::AbstractMatrix,
    atol::Real = 1e-8
)

    # sort eigenvalues
    alg_eigs_soted = sort(alg_eigs, by = x -> (real(x), imag(x)))
    reg_eigs_sorted = sort(reg_eigs, by = x -> (real(x), imag(x)))

    # check equivalence
    eig_check = isapprox(alg_eigs_soted, reg_eigs_sorted; atol=atol, rtol=0)

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
        eigenvalues_match = eig_check,
        max_residual = max_res,
        mean_residual = res_sum / len
    )
end