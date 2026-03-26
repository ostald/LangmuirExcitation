using MAT
using DataInterpolations
using Peaks
using GLMakie #, MakieExtra
#using Makie
using PCHIPInterpolation

include("load_guisdap.jl")
include("utils.jl")
include("constants.jl")
c = physical_constants

showplots = false

#load guisdap data
guisdap_dir = "analyzed_parameters/2018-12-07_folke_6.4@42mb"
files = filter(x -> x[end-3:end] == ".mat", readdir(guisdap_dir))

# 3 different time intervals to try out:
interval1 = 143:161
interval2 = 292:310
interval3 = 246:264

#for interval in [interval1, interval2, interval3]

interval = interval2
load_files = files[interval]

if showplots
    file = files[290]
    ex_time, h_param, ne, Te = load_guisdap(joinpath(guisdap_dir, file))
    println(ex_time)
    #t_start = ex_time[1, :]
    #t_end = ex_time[2, :]

    if showplots
        lines(ne, h_param/1e3)
    end
end

t_start, t_end, h_median, ne_median, Te_median, h_mad, ne_mad, Te_mad = 
    average_guisdap(joinpath.(guisdap_dir, load_files))
gd = load_guisdap.(joinpath.(guisdap_dir, load_files))

#interpolate guisdap data to psd height resolution
ne_itp = Interpolator(h_median, ne_median).(psd.h_atm)
Te_itp = Interpolator(h_median, Te_median).(psd.h_atm)

if showplots
    ne_mad_lower = copy(ne_mad)
    ne_mad_lower[(ne_median - ne_mad) .< 0] .= ne_median[(ne_median - ne_mad) .< 0] ./10
    fig, ax, lin = scatterlines(ne_median, h_median/1e3, 
        marker = 'x', 
        label = "median",
        axis = (xscale = log10, limits = ((1e6, 2e12), nothing)
    ),
    )
    errorbars!(ne_median, h_median/1e3, ne_mad_lower, ne_mad, direction = :x; color = :red) # same low and high error
    for m in gd
        scatter!(m.ne, m.h_param/1e3, color = "black", alpha = 0.3
        )
    end
    lines!(ne_itp, psd.h_atm/1e3)
    axislegend(ax)
    display(fig)
end

if showplots
    Te_mad_lower = copy(Te_mad)
    Te_mad_lower[(Te_median - Te_mad) .< 0] .= Te_median[(Te_median - Te_mad) .< 0] ./10
    fig, ax, lin = scatterlines(Te_median, h_median/1e3, 
        marker = 'x', 
        label = "median",
        #axis = (xscale = log10, #limits = ((1e6, 2e12), nothing)
    #),
    )
    errorbars!(Te_median, h_median/1e3, Te_mad_lower, Te_mad, direction = :x; color = :red) # same low and high error
    for m in gd
        scatter!(m.Te, m.h_param/1e3, color = "black", alpha = 0.3
        )
    end
    lines!(Te_itp, psd.h_atm/1e3)
    axislegend(ax)
    display(fig)
end

# calculate plasma frequency and thermal velocity 
# necessary for resonance condition ang growth rate
wp = plasma_freq.(ne_itp, c.me)
vth = thermal_velocity.(Te_itp, c.me)


psd_dir = "Afven Train PSB/psd/"
psd_files = filter(x-> contains(x, ".mat"), readdir(psd_dir))
# loop over all files
for psd_f in psd_files
    println("Processing ", psd_f)
    psd_file = joinpath(psd_dir, psd_f)
    #psd_file = "Afven Train PSB/psd/psd-30.mat"
    # Load the .mat file
    psd_data = matread(psd_file)

    # Convert psd_data into a named tuple
    psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data)
        )
    #println(keys(psd))

    showplots = false
    if showplots
        fig, ax, hm = heatmap(psd.vpar_centres, psd.h_atm/1e3, psd.F[:, :, 50], 
            colorscale = log10,
            colorrange = (1e20, 1e28),
            axis = (xlabel = "Velocity [m/s]",
                ylabel = "Height [km]",
                #xscale = Makie.Symlog10(0.01),
                ),
            )
        #xlims!(-100, 100)
        Colorbar(fig[1, 2], hm)
        display(fig)
    end

    dfpardvpar = stack(
        [gradient_simple(psd.vpar_centres, psd.F[:, ih, it])[2]
            for ih in eachindex(psd.h_atm),
                it in eachindex(psd.t_run)
            ]
        )
    v_middle = psd.vpar_centres[1:end-1] .+ diff(psd.vpar_centres)

    Ekin = E_ev(psd.vpar_centres, c.me_ev)
    E    = E_ev(v_middle, c.me_ev)

    grad = copy(dfpardvpar[abs.(E) .< 200, :, :])
    v = v_middle[abs.(E) .< 200]
    E    = E_ev(v, c.me_ev)
    grad[abs.(E) .< 3, :, :] .= NaN

    dfdvmax = Matrix{Any}(undef, size(grad, 2), size(grad, 3))
    i_dfdvmax = Matrix{Any}(undef, size(grad, 2), size(grad, 3))
    v_dfdvmax = Matrix{Any}(undef, size(grad, 2), size(grad, 3))
    showplots = false
    for it in axes(psd.t_run, 1)
        for ih in axes(psd.h_atm, 1)
            #ih = 395
            #it = 36
            imax, hmax = findmaxima(grad[:, ih, it]) |> peakheights(; min = 0)

            if showplots    
                lines(E, grad[:, ih, it])
            end
            if showplots
                fig, ax, hm = heatmap(E, psd.h_atm/1e3, grad[:, :, it], colorrange = (-1e20, 1e20),)
                xlims!(-200, 200)
                Colorbar(fig[1, 2], hm)
                display(fig)
            end
            if showplots
                fig, ax, lin = lines(E, grad[:, ih, it])
                scatter!(E[imax], hmax)
                display(fig)
            end

            dfdvmax[ih, it] = hmax
            i_dfdvmax[ih, it] = imax
            v_dfdvmax[ih, it] = v[imax]

            #im = 3
            #vmax = v[imax[im]]
            #dfdvmax = hmax[im]
            #E_ev(vmax, c.me_ev)
        end
    end

    if false
        it = 36
        showplots = true
        if showplots
            fig, ax, hm = heatmap(E, psd.h_atm/1e3, grad[:, :, it], colorrange = (-1e20, 1e20),)
            xlims!(-200, 200)
            Colorbar(fig[1, 2], hm, label = "df∥/d|v|∥")
            E_profiles = E_ev.(v_dfdvmax[:, it], c.me_ev)
            nmax = maximum(length.(E_profiles[:]))
            E_mat = stack(rpad_array.(E_profiles, nmax, NaN))
            [scatter!(Ep, psd.h_atm/1e3, marker = 'x') for Ep in eachrow(E_mat)]
            display(fig)
        end

        if showplots
            fig, ax, hm = heatmap(Ekin, psd.h_atm/1e3, psd.F[:, :, it], 
               colorrange = (1e20, 1e28),
                #colorrange = (1e24, 1e26), 
                colorscale = log10,
            )
            xlims!(-200, 200)
            Colorbar(fig[1, 2], hm, label = "f∥")
            E_profiles = E_ev.(v_dfdvmax[:, it], c.me_ev)
            nmax = maximum(length.(E_profiles[:]))
            E_mat = stack(rpad_array.(E_profiles, nmax, NaN))
            [scatter!(Ep, psd.h_atm/1e3, marker = 'x') for Ep in eachrow(E_mat)]
            display(fig)
        end

        if showplots
            ih = 370
            vmax = only(v_dfdvmax[ih, it])
            vth_ = vth[ih]
            wp_ = wp[ih]
            k = -2000:2000
            fig, ax, lin = lines(k, k*vmax)
            lines!(k, sqrt.(wp_.^2 .+ 3/2*vth_.^2 .*k.^2))
            display(fig)
        end

        #solve for k and omega
        #k^2 vmax^2 .- 3/2*vth.^2*k.^2 3/2 = wp^2
        #k^2 = wp^2 / (vmax^2 .- 3/2*vth.^2)
        k_growth = wp[ih] / sqrt(vmax^2 - 3/2*vth[ih]^2)
        w_growth = k_growth * vmax

        #growth rate!
        #double check!!
        gamma = pi * wp[ih]^2 / k_growth^2 * dfdvmax[ih, it]

        #check with eiscat wavenumber
    end

    k_growth = Vector{Any}(undef, size(psd.t_run, 1))
    v_max = Vector{Any}(undef, size(psd.t_run, 1))
    gamma = Vector{Any}(undef, size(psd.t_run, 1))

    for it in axes(psd.t_run, 1)
        v_profile = v_dfdvmax[:, it]
        v_max[it] = v_profile
        dfdvm_profile = dfdvmax[:, it]
        #E_profiles = E_ev.(v_profile, c.me_ev)
        nmax = maximum(length.(v_profile[:]))
        v_mat = stack(rpad_array.(v_profile, nmax, NaN))'
        dfdvm_mat = stack(rpad_array.(dfdvm_profile, nmax, NaN))'
        #E_mat = stack(rpad_array.(E_profiles, nmax, NaN))
        k_g = wp ./ sqrt.(v_mat.^2 .- 3/2*vth .^2)
        k_growth[it] = k_g
        w_g = k_g .* v_mat
        gamma[it] = pi * wp .^2 ./ k_g .^2 .* dfdvm_mat
    end

    #save output for each file
    psd_data["k_growth"] = k_growth
    psd_data["v_phase"] = v_max
    psd_data["growth_rate"] = gamma
    
    if !isdir(joinpath(psd_dir, "LE"))
        mkdir(joinpath(psd_dir, "LE"))
    end

    if !isdir(joinpath(psd_dir, "LE", string(interval)))
        mkdir(joinpath(psd_dir, "LE", string(interval)))
    end

    matwrite(joinpath(psd_dir, "LE", string(interval), psd_f[1:end-4]*"_LE.mat"), psd_data)
end
end
"""
it = 1:60 per file:     done
ih= 1:406               done

ifiles                  done
"""