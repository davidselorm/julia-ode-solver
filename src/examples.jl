module ODEExamples

# 1. Harmonic Oscillator
function harmonic_oscillator(t, y)
    # y[1] = position, y[2] = velocity
    return [y[2], -y[1]]
end

# 2. Lorenz Chaotic Attractor
function lorenz_system(t, y; sigma=10.0, rho=28.0, beta=8.0/3.0)
    # y = [x, y, z]
    dx = sigma * (y[2] - y[1])
    dy = y[1] * (rho - y[3]) - y[2]
    dz = y[1] * y[2] - beta * y[3]
    return [dx, dy, dz]
end

end # module
