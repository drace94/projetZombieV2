function f(p, alpha, uⁿ)
    return [-p*uⁿ[1]*uⁿ[2], p*uⁿ[1]*uⁿ[2] - alpha*uⁿ[2], alpha*uⁿ[2]]
end

function RK4(p, alpha, uⁿ, dt)
    k1 = f(p, alpha, uⁿ)
    k2 = f(p, alpha, uⁿ .+ (dt / 2.0) .* k1)
    k3 = f(p, alpha, uⁿ .+ (dt / 2.0) .* k2)
    k4 = f(p, alpha, uⁿ .+ dt .* k3)
    uⁿ⁺¹ = uⁿ .+ (dt/6.0) .* (k1 .+ k2./2.0 .+ k3./2.0 .+ k4)
    return uⁿ⁺¹
end

function iter_heat(A,B,uⁿ)
    return A \ (B * uⁿ)
end