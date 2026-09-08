module AdaptiveController

export step_size_control

"""
Proportional-Integral (PI) adaptive step-size controller.
"""
function step_size_control(err, tol, dt; safety=0.9, min_scale=0.2, max_scale=5.0)
    if err == 0
        return dt * max_scale
    end
    scale = safety * (tol / err)^0.2
    scale = clamp(scale, min_scale, max_scale)
    return dt * scale
end

end # module
