include("constants.jl")

using NCDatasets

function load_psd(ncfile)
    # Loads one PSD file produced by AURORA.jl. AURORA now writes its output as
    # netCDF (psd.nc) instead of the old MATLAB .mat files. We read the reduced
    # parallel-velocity distribution F and its coordinate grids, renaming the
    # netCDF variables to the field names the rest of the code expects:
    #   altitude -> h_atm, time -> t_run.
    # Returns a Dict{String,Any} (like the old `matread`), so it stays compatible
    # with the NamedTuple conversion and the .mat output written further down.
    NCDataset(ncfile, "r") do ds
        Dict{String, Any}(
            "h_atm"        => Array(ds["altitude"]),
            "t_run"        => Array(ds["time"]),
            "vpar_centers" => Array(ds["vpar_centers"]),
            "vpar_edges"   => Array(ds["vpar_edges"]),
            "F"            => Array(ds["F"]),
        )
    end
end

function gradient_simple(x, y)
    # simple nearest neighbour gradient for 
    # arguments: x: independent variable data points, y = f(x) with df/dx being the gradient
    # output: df/dx at the middle points between x
    gradient = diff(y) ./ diff(x)
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
