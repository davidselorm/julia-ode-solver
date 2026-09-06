function adaptive_tolerance_check(err, tol, dt)
    if err < tol
        return dt * 1.25
    else
        return dt * 0.5
    end
end
