abstract type AbstractReflector end

# SFMIN from LAPACK

"""
    safemin(::Type{T})

Computes the smallest element s::T such that one(T) / s does not overflow.
"""
safemin(::Type{T}) where {T} = floatmin(T) / eps(T) # this only works for real numbers

# follow LAPACK
function reflector!(x::AbstractVector{T}) where {T}
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
        knt = 0
        while abs(β) < sfmin && knt < 20
            knt = knt + 1
            @. x = rsfmin * x
            β = rsfmin * β
            α = rsfmin * α
        end

        τ = (β - α) / β
        α = one(eltype(x)) / (α - β)
        for i = 2:n
            x[i] = α * x[i]
        end

        for _ in 1:knt
            β = sfmin * β
        end
        x[1] = β

        return τ

    end
end

