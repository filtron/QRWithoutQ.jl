module QRWithoutQ

using LinearAlgebra

include("reflectors.jl")
export clarfg!

include("unblocked.jl")
export qr_unblocked!, qrwoq_unblocked!


"""
    qrwoq!(A::AbstractMatrix{T}, A::AbstractMatrix{T}) where {T}

Computes the upper triangular qr factor of A in-place and returns a view into A giving the upper triangular factor.
The scalar factor the in
"""
qrwoq!(A::AbstractMatrix{T}) where {T} = qrwoq_unblocked!(A)
export qrwoq!


#=
function qrwoq_unblocked!(A::AbstractMatrix{T}) where {T}
    LinearAlgebra.require_one_based_indexing(A)
    m, n = size(A)
    n <= m || throw(
        DimensionMismatch("A must have more rows than columns, received shape $((m, n))"),
    ) # easier to reason about

    for k = 1:min(m - 1 + !(T <: Real), n)
        x = view(A, k:m, k)
        τk = LinearAlgebra.reflector!(x)
        LinearAlgebra.reflectorApply!(x, τk, view(A, k:m, k+1:n))
    end
    Av = view(A, 1:n, 1:n)
    triu!(Av)
    return Av
end
=#

end
