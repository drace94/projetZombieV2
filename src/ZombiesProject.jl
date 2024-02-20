module ZombiesProject

using Printf, ToeplitzMatrices, CairoMakie

export run, animation

# ODE system: f(U) = dU/dt, U = [S,I,R]
# with: dS/dt = -β*S*I
#       dI/dt = β*S*I - γ*I
#       dR/dt = γ*I
function f(β::Float64, γ::Float64, Uⁿ::Vector{Float64})
    ∂ₜS = -β * Uⁿ[1] * Uⁿ[2]
    ∂ₜI =  β * Uⁿ[1] * Uⁿ[2] - γ * Uⁿ[2]
    ∂ₜR =  γ * Uⁿ[2]
    return [∂ₜS, ∂ₜI, ∂ₜR]
end

# RK4 method for time integration
function RK4(β::Float64, γ::Float64, Δt::Float64, Uⁿ::Vector{Float64})
    k₁ = f(β, γ, Uⁿ)
    k₂ = f(β, γ, Uⁿ .+ (Δt / 2.) .* k₁)
    k₃ = f(β, γ, Uⁿ .+ (Δt / 2.) .* k₂)
    k₄ = f(β, γ, Uⁿ .+ Δt .* k₃)
    Uⁿ⁺¹ = Uⁿ .+ (Δt/6.) .* (k₁ .+ k₂./2. .+ k₃./2. .+ k₄)
    return Uⁿ⁺¹
end

# solve SIR system at time tⁿ
function compute_sir(Nx²::Int64, β::Float64, γ::Float64, Δt::Float64, Uⁿ::Vector{Vector{Float64}})
    for i=1:Nx²
        Uⁿi = RK4(β, γ, Δt, [Uⁿ[1][i], Uⁿ[2][i], Uⁿ[3][i]])
        Uⁿ[1][i] = Uⁿi[1]
        Uⁿ[2][i] = Uⁿi[2]
        Uⁿ[3][i] = Uⁿi[3]
    end
    return Uⁿ
end

# heat equation dUⁿ/dt = α * d²Uⁿ/dx² associated FD linear system Auⁿ⁺¹ = Buⁿ
function centered_fd_system(Nx²::Int64, Δx::Float64, Δt::Float64, α::Float64)
    # create matrices A and B using Circulant() function
    Nx = Int(round(sqrt(Nx²)))
    # !!! TO BE FIXED !!!
    A₁ = zeros(Nx²)
    A₁[1     ] = 1. + 2. * α * Δt / (Δx^2)
    A₁[2     ] = - α * Δt / (2. * Δx^2)
    A₁[Nx+1  ] = - α * Δt / (2. * Δx^2)
    A₁[end-Nx] = - α * Δt / (2. * Δx^2)
    A₁[end   ] = - α * Δt / (2. * Δx^2)
    A = Circulant(A₁)
    B₁ = zeros(Nx²)
    B₁[1     ] = 1. - 2. * α * Δt / Δx^2
    B₁[2     ] = α * Δt / (2. * Δx^2)
    B₁[Nx+1  ] = α * Δt / (2. * Δx^2)
    B₁[end-Nx] = α * Δt / (2. * Δx^2)
    B₁[end   ] = α * Δt / (2. * Δx^2)
    B = Circulant(B₁)
    return A, B
end

# solve heat equation on each variables of Uⁿ at time tⁿ
function compute_heat(A::ToeplitzMatrices.Circulant{Float64}, B::ToeplitzMatrices.Circulant{Float64}, Uⁿ::Vector{Vector{Float64}})
    for i=1:3
        Uⁿ[i] = A \ (B * Uⁿ[i])
    end
    return Uⁿ
end

# solve global model over time: RK4 integration + heat equation at each step
function run(; scheme_param::Tuple{Int64,Float64,Float64}, model_param::Vector{Float64}, U₀::Vector{Vector{Float64}}, Δtₛ::Float64, tₑ::Float64)
    println("\nRun simulation...")
    start = time()
    tᵢ = 0.
    # unpack simulation parameters
    Nx, Δx, Δt = scheme_param
    α, β, γ = model_param
    Uⁿ = copy(U₀)
    Nx² = Int(Nx^2)
    # manual save of initial state into time series
    tₛ = [tᵢ]
    Sₛ = [U₀[1]]
    Iₛ = [U₀[2]]
    Rₛ = [U₀[3]]
    # iterate over time
    A, B = centered_fd_system(Nx², Δx, Δt, α)
    while tᵢ < tₑ
        # solve SIR model + heat equation
        Uⁿ = compute_sir(Nx², β, γ, Δt, Uⁿ)
        Uⁿ = compute_heat(A, B, Uⁿ)
        # save data to time series based on Δtₛ step
        if (abs(tᵢ - floor(tᵢ / Δtₛ) * Δtₛ) < Δt || tᵢ == tₑ)
            time_progress(tᵢ, tₑ)
            push!(tₛ, tᵢ)
            push!(Sₛ, Uⁿ[1])
            push!(Iₛ, Uⁿ[2])
            push!(Rₛ, Uⁿ[3])
        end
        # increment current simulation time tᵢ
        if tᵢ + Δt < tₑ
            tᵢ += Δt
        else
            Δt = tₑ - tᵢ
            tᵢ += Δt
        end
    end
    elapsed = time() - start
    @printf("\nCompleted. \n- computation time = %.4e s\n\n", elapsed)
    return [tₛ,Sₛ,Iₛ,Rₛ]
end

# custom progress bar function
function time_progress(tᵢ::Float64, tₑ::Float64; bar_size::Int64=30)
    progress = Int(floor(bar_size * tᵢ / tₑ))
    print("\r[")
    for i in 1:progress
        print("#")
    end
    for i in progress+1:bar_size
        print(" ")
    end
    @printf("]  tᵢ = %06.2f / %06.2f", tᵢ, tₑ)
    flush(stdout)
end

# animation from time series (heatmap)
function animation(L::Float64, tₛ::Vector{Float64}, Mₛ::Vector{Vector{Float64}}, folder_path::String, filename::String)
    @printf("\nCreate animation...")
    start = time()
    # time series size
    Ns = length(tₛ)
    Nx² = length(Mₛ[1])
    Nx = Int(round(sqrt(Nx²)))
    # mesh range
    x = range(0, stop=L, length=Nx)
    # frame handling
    function update_frame(ax, i)
        empty!(ax)
        # heatmap of Mᵢ=Mₛ(tᵢ)=Mₛ(tₛ(i)), reshaped from flattened format
        Mᵢ = Mₛ[i]
        M = reshape(Mᵢ, Nx, Nx)
        hm = heatmap!(ax, x, x, M, ratio=1, interpolation=true)
        # Colorbar(fig[1,2], hm)
        # plot current simulation time
        # t_string = @sprintf("t = %06.2f", tₛ[i])
        # ax.title("test")
        # text!(0.95*xhigh, 0.95*yhigh, text=t_string, align=(:right,:top))
    end
    # figure
    fig = Figure()
    ax = Axis(fig[1,1], xlabel=L"x", ylabel=L"y")
    # create and save the animation
    file_path = joinpath(folder_path, filename * ".mp4")
    record(fig, file_path, 1:Ns, framerate = 30) do i
        update_frame(ax, i)
    end
    elapsed = time() - start
    @printf("\nCompleted. \n- computation time = %.4e s\n", elapsed)
end

end # module ZombiesProject
