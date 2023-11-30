# naive implementation of qr without q
function qrR!(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    n <= m || throw(
        DimensionMismatch("A must have more rows than columns, received size $((m, n))"),
    ) # easier to reason about
    qr!(A)
    Av = view(A, 1:n, 1:n)
    triu!(Av)
    return Av
end

# i.e n-step prediction in Gauss-Markov process using naive qr without q
function lyapunov_recursion_qr(n, U0, Φ, QU, pre_array)
    m = LinearAlgebra.checksquare(U0)
    U = U0
    for _ = 1:n
        pre_array[1:m, 1:m] .= QU
        mul!(view(pre_array, m+1:2m, 1:m), U, Φ')
        U .= qrR!(pre_array)
    end
    return U
end

# i.e n-step prediction in Gauss-Markov process using non-allocating qr without q
function lyapunov_recursion_qrwoq(n, U0, Φ, QU, pre_array)
    m = LinearAlgebra.checksquare(U0)
    U = U0
    for _ = 1:n
        pre_array[1:m, 1:m] .= QU
        mul!(view(pre_array, m+1:2m, 1:m), U, Φ')
        U .= qrwoq!(pre_array)
    end
    return U
end

# make inout matrices to lyapunov_recursion
function make_test_problem(T, m)
    λ = T(1)
    A = λ * (I - T(2) * tril(ones(T, m, m)))
    Φ = exp(A)
    QU = diagm(0 => ones(T, m))
    U0 = diagm(0 => ones(T, m))
    return U0, Φ, QU
end

# benchmark qr algorithms
function lyapunov_benchmark(T, n, m)
    U0, Φ, QU = make_test_problem(T, m)
    pre_array = similar(U0, 2m, m)
    U01 = copy(U0)
    U02 = copy(U0)

    b1 = @benchmark lyapunov_recursion_qr($n, $U01, $Φ, $QU, $pre_array)
    b2 = @benchmark lyapunov_recursion_qrwoq($n, $U02, $Φ, $QU, $pre_array)
    return b1, b2
end