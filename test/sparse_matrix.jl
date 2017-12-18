using PATHSolver
using Base.Test


# Verify that pulling data directly out of the SparseMatrixCSC form gives
# exactly what the previous hand-coded algorithm gave.
@testset "sparse matrix" begin
  function sparse_matrix_reference(A::AbstractSparseMatrix)
    m, n = size(A)
    @assert m==n

    col_start = Array{Int}(n)
    col_len = Array{Int}(n)
    row = Array{Int}(0)
    data = Array{Float64}(0)
    for j in 1:n
      if j==1
        col_start[j] = 1
      else
        col_start[j] = col_start[j-1] + col_len[j-1]
      end

      col_len[j] = 0
      for i in 1:n
        if A[i,j] != 0.0
          col_len[j] += 1
          push!(row, i)
          push!(data, A[i,j])
        end
      end
    end
    return col_start, col_len, row, data
  end

  sparse_matrix_reference(A::Matrix) = sparse_matrix_reference(sparse(A))

  srand(42)
  for i in 1:100
    n = rand(5:20)
    M = sprandn(n, n, 0.1)
    @test PATHSolver.sparse_matrix(M) == sparse_matrix_reference(M)

    Mf = full(M)
    @test PATHSolver.sparse_matrix(Mf) == sparse_matrix_reference(Mf)
    @test PATHSolver.sparse_matrix(Mf) == PATHSolver.sparse_matrix(M)
  end
end

