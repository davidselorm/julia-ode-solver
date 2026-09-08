# julia-ode-solver

High-performance adaptive Runge-Kutta numerical differential equation integration in Julia.

## Features
- **Classical RK4**: 4th-order fixed-step Runge-Kutta integrator.
- **Dormand-Prince (RK45)**: 5th-order embedded adaptive step controller with error budget estimation.
- **Symplectic Verlet**: Energy-conserving symplectic integrator for celestial mechanics and Hamiltonian dynamics.
- **Lorenz Attractor & Orbital Models**: Built-in chaotic and oscillatory dynamical systems.

## Usage
```julia
using ODESolver

f(t, y) = -2.0 * y
y0 = 1.0
dt = 0.05
y_next = rk4_step(f, 0.0, y0, dt)
```
