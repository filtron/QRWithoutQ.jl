"""
    qr_unblocked!(A::AbstractMatrix{T}, τ) where {T}

Computes the unblocked QR decomposition of A, in-place.
The vector τ is used for storing the scalar coefficients for the Householder reflectors.
"""
function qr_unblocked!(A::AbstractMatrix{T}, τ) where {T}
    LinearAlgebra.require_one_based_indexing(A)
    m, n = size(A)
    length(τ) != min(m, n) && throw(DimensionMismatch("τ must be of length equal to the smallest dimension of A."))
    for k = 1:min(m - 1 + !(T<:Real), n)
        x = view(A, k:m, k)
        τk = clarfg!(x)
        τ[k] = τk
        apply_reflector!(x, τk, view(A, k:m, k + 1:n))
    end
end

function qr_unblocked_wy!(A::AbstractMatrix{T}, W::AbstractMatrix{T}, u::AbstractVector{T}) where {T}
    LinearAlgebra.require_one_based_indexing(A)
    mA, nA = size(A)
    mW, nW = size(W)
    N = div(n, bz)+1
    for k = 1:min(m, n)
        x = view(A, k:m, k)
        τk = clarfg!(x)

        # form  - σ * w
        u[1] = one(T)
        @views u[2:mW] .= x[2:mA]
        u .= -τk .* u
        # TODO: finish this

    end
end



"""
    qrwoq_unblocked!(A::AbstractMatrix{T}) where {T}

Computes the R factor in the unblocked QR decomposition, in-place.
R is returned as a view of A and all Householder vectors are destroyed.
"""
function qrwoq_unblocked!(A::AbstractMatrix{T}) where {T}
    LinearAlgebra.require_one_based_indexing(A)
    m, n = size(A)
    for k = 1:min(m - 1 + !(T<:Real), n)
        x = view(A, k:m, k)
        τk = clarfg!(x)
        apply_reflector!(x, τk, view(A, k:m, k + 1:n))
    end
    Av = view(A, 1:min(m, n), 1:n)
    triu!(Av)
    return Av
end
