using LinearAlgebra

"""
    safemin(::Type{T})

Smallest element of type T for which it is safe to divide by.
Intended to be an implementation of LAPACKs slamch( 'S' ) / slamch( 'E' )
but I could be wrong...
"""
safemin(::Type{T}) where {T} = floatmin(T) / eps(T) # this only works for real numbers

"""
    _upscale_input!(xtail::AbstractMatrix{T}, α, β, sfmin, rsfmin) where {T}

Rescale the input to clarfg! if it is too small.
"""
function _upscale_input!(xtail::AbstractVector{T}, α, β, sfmin, rsfmin) where {T}
    knt = 0
    while abs(β) < sfmin && knt < 20
        knt = knt + 1
        @. xtail = rsfmin * xtail
        β = rsfmin * β
        α = rsfmin * α
    end
    return α, β, knt
end


"""
    _downcale_beta(β, knt, sfmin)

Undo the rescaling of β.
"""
function _downcale_beta(β, knt, sfmin)
    for _ in 1:knt
        β = sfmin * β
    end
    return β
end

"""
    clarfg!(x::AbstractVector{T}) where {T}

Intended to be a reasonably faithful implementation of clarg from LAPACK.
"""
@inline function clarfg!(x::AbstractVector{T}) where {T}
    LinearAlgebra.require_one_based_indexing(x)
    n = length(x)
    n == 0 && return zero(eltype(x))
    @inbounds begin
        α = x[1] # LAPACK input is (α, x[2:end])
        xnorm = norm(x)
        if iszero(xnorm)
            τ = zero(α/xnorm)
            return τ
        end
        β = -T(copysign(xnorm, real(α)))
        # rescale x to safe range if necesssary
        sfmin = safemin(eltype(x))
        rsfmin = one(eltype(x)) / sfmin
        α, β, knt = _upscale_input!(view(x, 2:length(x)), α, β, sfmin, rsfmin)
        τ = (β - α) / β
        α = one(eltype(x)) / (α - β)
        for i = 2:n
            x[i] = α * x[i]
        end
        β =  _downcale_beta(β, knt, rsfmin)
        x[1] = β
        return τ
    end
end

"""
    apply_reflector!(x::AbstractVector, τ, A::AbstractVecOrMat)

"""
@inline function apply_reflector!(x::AbstractVector, τ, A::AbstractVecOrMat)
    LinearAlgebra.require_one_based_indexing(x)
    LinearAlgebra.require_one_based_indexing(A)
    m, n = size(A)
    if length(x) != m
        throw(DimensionMismatch("reflector has length $(length(x)), which must match the first dimension of matrix A, $m"))
    end
    m == 0 && return A
    @inbounds for j = 1:n
        Aj, xj = view(A, 2:m, j), view(x, 2:m)
        vAj = conj(τ)*(A[1, j] + dot(xj, Aj))
        A[1, j] -= vAj
        axpy!(-vAj, xj, Aj)
    end
    return A
end

