include("Projet_Zombie_lib.jl")
using LoopVectorization, CairoMakie, Printf, LinearAlgebra, ToeplitzMatrices, JLD

# Spatial domain
N_x = 20
N_y = 20

# Parameters of heat equation
dt = 1e-3
alpha = 1.0e-1
L = 1.0
domain_x = [-L/2.0,L/2.0]
domain_y = [-L/2.0,L/2.0]
dx = (domain_x[2] - domain_x[1]) / N_x
dy = dx
N_save = 1000

first_column_A = zeros(N_x * N_y)
first_column_A[1] = 1.0 + 2.0 * alpha * dt / (dx ^ 2)
first_column_A[2] = - alpha * dt / (2.0 * dx ^ 2)
first_column_A[N_x+1] = - alpha * dt / (2.0 * dx ^ 2)
first_column_A[end-N_x] = - alpha * dt / (2.0 * dx ^ 2)
first_column_A[end] = - alpha * dt  / (2.0 * dx ^ 2)
A = Circulant(first_column_A)

first_column_B = zeros(N_x * N_y)
first_column_B[1] = 1.0 - 2.0 * alpha * dt / (dx ^ 2)
first_column_B[2] = alpha * dt / (2.0 * dx ^ 2)
first_column_B[N_x+1] = alpha * dt / (2.0 * dx ^ 2)
first_column_B[end-N_x] = alpha * dt / (2.0 * dx ^ 2)
first_column_B[end] = alpha * dt  / (2.0 * dx ^ 2)
B = Circulant(first_column_B)

# Array to keep values
S = zeros(N_x * N_y, N_save + 1)
I = zeros(N_x * N_y, N_save + 1)
R = zeros(N_x * N_y, N_save + 1)

# Init
S₀ = zeros(N_x * N_y)
I₀ = zeros(N_x * N_y)
R₀ = zeros(N_x * N_y)
S₀[N_x * N_y ÷ 2] = 499.0
I₀[N_x * N_y ÷ 2] = 1.0

Sⁿ = copy(S₀)
Iⁿ = copy(I₀)
Rⁿ = copy(R₀)

# Intermediate values
S_ = zeros(N_x * N_y, 3)

# Time and save parameters
t = 0.0
t_save = 0.0
index_save = 1

t_final = 100.0
alpha = 5.0e-2
p = 1.0e-3

t_affichage = range(0.0,t_final,N_save+1)

@time while t < t_final
    
    # Splitting : S-I-R ODEs then heat
    for i=1:N_x * N_y
        global S_[i,:] = RK4(p, alpha, [Sⁿ[i],Iⁿ[i],Rⁿ[i]], dt)
    end

    for i=1:3
        global S_[:,i] = iter_heat(A,B,S_[:,i])
    end

    if t + dt >= t_final
       global dt = t_final - t
    end

    global t += dt
    global t_save += dt
    global Sⁿ = copy(S_[:,1])
    global Iⁿ = copy(S_[:,2])
    global Rⁿ = copy(S_[:,3])

    if t_save > t_final / 10.0
        println(((t / t_final) * 100),"% de la run effectuée (sauvegarde)")
        global index_save += 1
        global t_save = 0.0
        global S[:,index_save] = copy(Sⁿ)
        global I[:,index_save] = copy(Iⁿ)
        global R[:,index_save] = copy(Rⁿ)
    end
end

save("sauvegardes/sains.jld","S",S)
save("sauvegardes/infectes.jld","I",I)
save("sauvegardes/proteges.jld","R",R)

println("Matrice sauvegardées, vous pouvez passer sur le code postpro")