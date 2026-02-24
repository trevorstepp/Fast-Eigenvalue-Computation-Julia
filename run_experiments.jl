include("BlockEig.jl")
using .BlockEig
using Printf
using Statistics
using LinearAlgebra

const N_VALUES = [100, 250, 500, 750, 1000, 1500, 2000]
const K_VAL = 3
const NUM_RUNS = 5

function run_experiments()
    block_time = Float64[]
    dense_time = Float64[]
    block_max_residuals = Float64[]

    for n in N_VALUES
        # build matrix
        M = build_block_matrix(K_VAL, n; seed=0)

        # hold run times of inner loop
        run_times_block = Float64[]
        run_times_dense = Float64[]

        # track largest residual 
        max_res_block = -Inf

        # warmup once per matrix size (JIT compiled language)
        eig_KxK_diagblocks(K_VAL, n, M)
        eigen(M)

        for i in 1:NUM_RUNS

            println("\nRunning n = $n, iteration $i")

            # timing
            block_result = time_block_method(K_VAL, n, M)
            dense_result = time_dense_eig(M)

            # store times 
            push!(run_times_block, block_result.time)
            push!(run_times_dense, dense_result.time)

            # verification
            check = verify_results(
                M,
                block_result.eigenvalues,
                block_result.eigenvectors,
                dense_result.eigenvalues,
                dense_result.eigenvectors
            )

            # update max residual
            max_res_block = max(max_res_block, check.max_residual)

            @printf(
                "Eigenvalues match (top %d): %s\n", 
                check.num_eigenvalues_compared,
                check.eigenvalues_match
            )
            @printf("Maximum residual: %.2e\n", check.max_residual)
            @printf("Mean residual: %.2e\n", check.mean_residual)
        end

        # store average times
        avg_time_block = mean(run_times_block)
        avg_time_dense = mean(run_times_dense)

        push!(block_time, avg_time_block)
        push!(dense_time, avg_time_dense)

        # store max residual
        push!(block_max_residuals, max_res_block)
    end

    return block_time, dense_time, block_max_residuals
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