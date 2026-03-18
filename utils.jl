include("constants.jl")

function gradient_simple(x, y)
    # simple nearest neighbour gradient for 
    # arguments: x: independent variable data points, y = f(x) with df/dx being the gradient
    # output: df/dx at the middle points between x
    gradient = diff(f) ./ diff(x)
    x_middle = x[1:end-1] + diff(x)/2
    return x_middle, gradient
end


function E_ev(v_mps, m_ev)
    E = m_ev / 2 * v_mps .* abs.(v_mps) / c.c^2
    return E
end

function plasma_freq(ne, m_kg)
    return sqrt(ne * c.qe^2 / (m_kg * c.eps0 ))
end

function thermal_velocity(T, m)
    return sqrt(2*c.kb*T/m)
end
