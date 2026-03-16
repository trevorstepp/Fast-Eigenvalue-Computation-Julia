include("BlockEig.jl")
using .BlockEig
using Printf
using Statistics
using LinearAlgebra
using BenchmarkTools

#BenchmarkTools.DEFAULT_PARAMETERS.samples = 5
#BenchmarkTools.DEFAULT_PARAMETERS.evals = 1

const N_VALUES = [100, 250, 500, 750, 1000, 1500, 2000]
const K_VAL = 3

function run_experiments()
    block_time = Float64[]
    dense_time = Float64[]
    block_residuals = Float64[]

    for n in N_VALUES
        println("\nRunning n = $n")
        # build matrix
        M = build_block_matrix(K_VAL, n; seed=0)

        # timing
        t_block = @benchmark eig_KxK_diagblocks($K_VAL, $n, $M)
        t_dense = @benchmark eigen($M)

        # get results
        block_eigs, block_vecs = eig_KxK_diagblocks(K_VAL, n, M)
        E = eigen(M)

        # store times 
        push!(block_time, median(t_block).time / 1e9)
        push!(dense_time, median(t_dense).time / 1e9)

        # verification
        check = verify_results(
            M,
            block_eigs,
            block_vecs,
            E.values,
            E.vectors
        )

        @printf(
            "Eigenvalues match (top %d): %s\n", 
            check.num_eigenvalues_compared,
            check.eigenvalues_match
        )
        @printf("Maximum residual: %.2e\n", check.max_residual)
        @printf("Mean residual: %.2e\n", check.mean_residual)

        # store residual
        push!(block_residuals, check.max_residual)
    end

    return block_time, dense_time, block_residuals
end

if abspath(PROGRAM_FILE) == @__FILE__
    block_time, dense_time, block_residuals = run_experiments()

    # save results for plotting
    using CSV, DataFrames
    df = DataFrame(
        K = K_VAL,
        n = N_VALUES,
        block_time = block_time,
        dense_time = dense_time,
        max_residual = block_residuals
    )
    CSV.write("julia_timings.csv", df)
end