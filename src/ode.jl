module ODESolver

export rk4_step, rk45_step, solve_adaptive, symplectic_verlet_step

"""
Classical 4th-order Runge-Kutta (RK4) single step.
"""
function rk4_step(f, t, y, dt)
    k1 = f(t, y)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
end

"""
Dormand-Prince (RK45) embedded 5th/4th order adaptive step.
Returns (y_next, error_estimate).
"""
function rk45_step(f, t, y, dt)
    k1 = f(t, y)
    k2 = f(t + dt/5, y + dt * (k1 / 5))
    k3 = f(t + 3*dt/10, y + dt * (3*k1/40 + 9*k2/40))
    k4 = f(t + 4*dt/5, y + dt * (44*k1/45 - 56*k2/15 + 32*k3/9))
    k5 = f(t + 8*dt/9, y + dt * (19372*k1/6561 - 25360*k2/2187 + 64448*k3/6561 - 212*k4/729))
    k6 = f(t + dt, y + dt * (9017*k1/3168 - 355*k2/33 + 46732*k3/5247 + 49*k4/176 - 5103*k5/18656))

    # 5th-order solution
    y5 = y + dt * (35*k1/384 + 500*k3/1113 + 125*k4/192 - 2187*k5/6784 + 11*k6/84)

    # 4th-order embedded solution
    k7 = f(t + dt, y5)
    y4 = y + dt * (5179*k1/57600 + 7571*k3/16695 + 393*k4/640 - 92097*k5/339200 + 187*k6/2100 + k7/40)

    err = norm(y5 - y4)
    return y5, err
end

"""
Symplectic Verlet integrator for Hamiltonian / orbital dynamics.
Preserves phase-space volume and energy.
"""
function symplectic_verlet_step(acc_fn, q, v, dt)
    v_half = v + 0.5 * dt * acc_fn(q)
    q_next = q + dt * v_half
    v_next = v_half + 0.5 * dt * acc_fn(q_next)
    return q_next, v_next
end

function norm(x::Number)
    return abs(x)
end

function norm(x::AbstractArray)
    return sqrt(sum(x .^ 2))
end

end # module
