# Zombie Project (V2)

## About the project

Upgrading of an old research project at INSA Toulouse (2020), originally written in Python and here implemented in Julia.<br>
The objective of this project is to simulate the spread of an epidemic (Zombie for example), based on the SIR compartmental model, over a 2D grid:
```math
\begin{cases}
    &\displaystyle\frac{\partial S}{\partial t} &= &-\beta SI \\[1em]
    &\displaystyle\frac{\partial I}{\partial t} &= &\beta SI - \gamma I \\[1em]
    &\displaystyle\frac{\partial R}{\partial t} &= &\gamma I
\end{cases},
```
where $S$, $I$ and $R$ refer to *susceptible*, *infectious* and *removed* (immune) individuals respectively. The coefficient $\beta$ denotes the disease transmission rate (often introduced as $\beta/N$ in the literature, where $N$ is the total population), while $\gamma$ denotes the cure rate.<br>
We also introduce in our global model a notion of population displacement, here modeled by the heat equation:
```math
\begin{equation}
    \displaystyle\frac{\partial u}{\partial t} = \nu\Updelta u = \nu\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial x^2}\right)
\end{equation},
```
where $u(t,x,y)=\{S,I,R\}$ and $\nu$ is a positive coefficient that sets the diffusion speed.<br>
The initial condition $u(t=0)=u_0(x,y)$ will therefore be set by the initial population $S_0$, $I_0$ and $R_0$.<br>
We impose circular boundary conditions over the domain $\Omega=[0,L_x]\times[0,L_y]$ to solve the heat equation, which implies:
```math
\begin{equation}
    u(t,x=0,y) = u(t,x=L_x,y) \;\;\text{ and }\;\; u(t,x,y=0) = u(t,x,y=L_y)
\end{equation}
```

## Features

The `ZombiesProject`module instantiates an environment containing the required library dependencies for the proper functioning of scripts. It involves splitting solve the SIR model and the heat equation separately between instants $t^n$ and $t^{n+1}$, using the following numerical schemes:
- A Fouth-Order Runge-Kutta (RK4) method for the SIR model
- A Crank-Nicolson scheme in time and a centered finite difference scheme in space for the heat equation
This modules also offers some post-processing functionalities for analyzing simulation results.

## Usage

To run a test case defined in the `main.jl` script, start by activating the environment:

```bash
julia --project=/path/to/InOutProject
```

Or alternatively if you are already in the Julia REPL:

```julia
using Pkg
Pkg.activate("/path/to/InOutProject")
```

Then you juste need to execute the script as follows:

```julia
include("path/to/main.jl")
```
