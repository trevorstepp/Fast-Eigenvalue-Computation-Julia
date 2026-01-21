include("BlockEig.jl")
using .BlockEig
using Printf

const N_VALUES = [100, 250, 500, 750, 1000, 1500, 2000]
const K = 3

function run_experiments()
    block_time = Float64[]
    dense_time = Float64[]

    for n in N_VALUES
        println("\nRunning n = $n")

        # build matrix
        M = build_block_matrix(K, n; seed=0)

        # timing
        block_result = time_block_method(K, n, M)
        dense_result = time_dense_eig(M)

        # store times 
        push!(block_time, block_result.time)
        push!(dense_time, dense_result.time)

        # verification
        check = verify_results(
            M,
            block_result.eigenvalues,
            block_result.eigenvectors,
            dense_result.eigenvalues,
            dense_result.eigenvectors
        )

        println("Eigenvalues match: ", check.eigenvalues_match)
        @printf("Maximum residual: %.2e\n", check.max_residual)
        @printf("Mean residual: %.2e\n", check.mean_residual)
    end

    return block_time, dense_time
end

if abspath(PROGRAM_FILE) == @__FILE__
    block_time, dense_time = run_experiments()

    # save results for plotting
    using CSV, DataFrames
    df = DataFrame(
        n = N_VALUES,
        block_time = block_time,
        dense_time = dense_time
    )
    CSV.write("timings.csv", df)
end