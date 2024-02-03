include("model_lib.jl")

# define spatial domain Ω = [0,L]x[0,L]
L = 1.0

# numerical scheme parameters
# 1. space discretization
Nx = 51
Δx = L / Nx
# 2. time step
Δt = 1e-2
scheme_param = (Nx,Δx,Δt)

# model parameters
# 1. heat equation thermal diffusivity α
α = 2e-5
# 2. SIR model infection rate β, recovery rate γ
β = 5e-1
γ = 1e-2
model_param = [α,β,γ]

# initialize SIR variables
Nx² = Nx^2    # flattened matrices size
S₀, I₀, R₀ = zeros(Nx²), zeros(Nx²), zeros(Nx²)
S₀[Nx² ÷ 2 + 1] = 80.
I₀[Nx² ÷ 2 + 1] = 20.
U₀ = [S₀,I₀,R₀]

# run simulation up to final time tₑ with saving rate Δtₛ to build time series
tₑ = 200.
Δtₛ = tₑ / 100.
tₛ,Sₛ,Iₛ,Rₛ = run(scheme_param=scheme_param, model_param=model_param, U₀=U₀, Δtₛ=Δtₛ, tₑ=tₑ)

# heatmap from time series
animation(L,tₛ,Sₛ,"anim_S")
animation(L,tₛ,Iₛ,"anim_I")