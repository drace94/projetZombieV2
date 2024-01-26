# Projet Zombie V2

Ré-écriture pour l'exercice en Julia du code développé par Hugo Valayer & Gabriel Depaillat en 2020 (disponible sur un autre repo github) du projet supervisé par David Sanchez (INSA Toulouse). Résolution par méthode de splitting d'un système SIR et d'une équation de la chaleur.

L'équation de la chaleur résolue est la suivante :
```math
\begin{cases}
    \frac{\partial u}{\partial t}(t,x) - \alpha\Updelta u(t,x) = 0 \\
    u(t,x) = u_0(x) \\
    u(t,-L/2) = u(t,L/2)
\end{cases}
```
où $u\in\{S,I,R\}$. Cette équation est resolue à chaque itération en temps avec un schéma aux différences finies en espace et Crank-Nicholson en temps.

Le système SIR s'écrit de la façon suivante :
```math
\begin{cases}
    \frac{\partial S}{\partial t}(t,x) = -pS(t,x)I(t,x) \\
    \frac{\partial I}{\partial t}(t,x) = pS(t,x)I(t,x) - \alpha I(t,x) \\
    \frac{\partial R}{\partial t}(t,x) = \alpha I(t,x)
\end{cases}
```

Le schéma utilisé pour résoudre en temps est RK4. Ainsi, à chaque itération en temps, on résoud $N_x.N_y$ systèmes SIR et 3 équations de la chaleur. \
Pour lancer le code il faut écrire :
```
julia projet_zombie.jl
```