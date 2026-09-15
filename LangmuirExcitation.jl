using MAT
using DataInterpolations
using Peaks
using WGLMakie
#`ssh -L 9384:localhost:9384 oliver@avaruus.uit.no`
#scatter(0, 0)
using PCHIPInterpolation
using Formatting

include("load_guisdap.jl")
include("utils.jl")
include("constants.jl")
c = physical_constants

showplots = false

# output directory:
save_dir = "/nfs/revontuli/data/oliver/LangmuirExcitation/results"

# intput directory:
psd_dir = "/nfs/revontuli/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_train_330-345s_pchip/analysis"
psd_file = "psd.nc"

# load psd data
psd_data = load_psd(joinpath(psd_dir, psd_file))
psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data))
h_atm = psd.h_atm

#define guisdap data
guisdap_dir = "analyzed_parameters/2018-12-07_folke_6.4@42mb"
files = filter(x -> x[end-3:end] == ".mat", readdir(guisdap_dir))

# 3 different time intervals of guisdap to try out:
interval1 = 143:161
interval2 = 292:310
interval3 = 246:264

#for interval in [interval1, interval2, interval3]
    interval = interval1
    load_files = files[interval]

    if showplots
        file = files[290]
        gd = load_guisdap(joinpath(guisdap_dir, file))
        ex_time, h_param, Te, ne = gd
        println(ex_time)
        #t_start = ex_time[1, :]
        #t_end = ex_time[2, :]

        if showplots
            lines(ne, h_param/1e3,
                axis = (xlabel = "Electron Density [m⁻³]",
                    ylabel = "Height [km]",
                    xscale = log10,),
                )
        end
    end

    ## guisdap data part
    #load guisdap data
    t_start, t_end, h_median, ne_median, Te_median, h_mad, ne_mad, Te_mad = 
        average_guisdap(joinpath.(guisdap_dir, load_files))
    gd = load_guisdap.(joinpath.(guisdap_dir, load_files))


    #interpolate guisdap data to psd height resolution
    ne_itp = Interpolator(h_median, ne_median).(h_atm)
    Te_itp = Interpolator(h_median, Te_median).(h_atm)

    if showplots
        ne_mad_lower = copy(ne_mad)
        ne_mad_lower[(ne_median - ne_mad) .< 0] .= ne_median[(ne_median - ne_mad) .< 0] ./10
        fig, ax, lin = scatterlines(ne_median, h_median/1e3, 
            marker = 'x', 
            label = "median",
            axis = (xscale = log10, 
                limits = ((1e6, 2e12), nothing),
                xlabel = "Electron Density [m⁻³]",
                ylabel = "Height [km]",
            ),
        )
        errorbars!(ne_median, h_median/1e3, ne_mad_lower, ne_mad, direction = :x; color = :red) # same low and high error
        for m in gd
            scatter!(m.ne, m.h_param/1e3, color = "black", alpha = 0.3
            )
        end
        lines!(ne_itp, h_atm/1e3)
        axislegend(ax)
        display(fig)
    end

    if showplots
        Te_mad_lower = copy(Te_mad)
        Te_mad_lower[(Te_median - Te_mad) .< 0] .= Te_median[(Te_median - Te_mad) .< 0] ./10
        fig, ax, lin = scatterlines(Te_median, h_median/1e3, 
            marker = 'x', 
            label = "median",
            axis = (xlabel = "Electron Temperature [K]",
                ylabel = "Height [km]",
            ),
        )
        errorbars!(Te_median, h_median/1e3, Te_mad_lower, Te_mad, direction = :x; color = :red) # same low and high error
        for m in gd
            scatter!(m.Te, m.h_param/1e3, color = "black", alpha = 0.3
            )
        end
        lines!(Te_itp, h_atm/1e3)
        axislegend(ax)
        display(fig)
    end

    # calculate plasma frequency and thermal velocity 
    # necessary for resonance condition and growth rate
    wp = plasma_freq.(ne_itp, c.me)
    vth = thermal_velocity.(Te_itp, c.me)

   
    #showplots = true
    if showplots
            it = 1030
            ih = 389
            fig, ax, hm = heatmap(psd.vpar_centers, 
                psd.h_atm/1e3, 
                psd.F[:, :, it], 
                colorscale = log10,
                colorrange = (1e-2, 1e2),
                axis = (xlabel = "Velocity [m/s]",
                    ylabel = "Height [km]",
                    #xscale = Makie.Symlog10(0.01),
                    title = format("t = {:.3f} s", psd.t_run[it]),
                ),
            )
            #xlims!(-500, 500)
            sleep(2)
            Colorbar(fig[1, 2], hm, label = "red. par. velocity distribution [?]")
            #display(fig)
        

            it_iterator = range(1, length(psd.t_run))
            framerate = 30*2
            F = Observable(psd.F[:, :, 1])

            fig, ax, hm = heatmap(E_ev(psd.vpar_centers, c.me_ev),
                psd.h_atm/1e3, 
                F, 
                colorscale = log10,
                colorrange = (1e-5, 1e2),
                axis = (xlabel = "Energy [eV]",
                    ylabel = "Height [km]",
                    #xscale = Makie.Symlog10(0.01),
                    title = format("t = {:.3f} s", psd.t_run[it]),
                ),
            )
            xlims!(-400, 400)
            sleep(2)
            Colorbar(fig[1, 2], hm)
            #display(fig)

            record(fig, "fpar_E_v2_500ev.mp4", it_iterator;
                framerate = framerate) do it
                F[] = psd.F[:, 140:160, it]
                ax.title = format("t = {:.3f} s", psd.t_run[it])
            end
    end

    ## psd data part
    # calculate gradient in parellel velocity 
    @time dfpardvpar = stack(
        [gradient_simple(psd.vpar_centers, psd.F[:, ih, it])[2]
            for ih in eachindex(psd.h_atm),
                it in eachindex(psd.t_run)
            ]
        )
    
    v_middle = psd.vpar_centers[1:end-1] .+ diff(psd.vpar_centers)
    dfpardvpar = dfpardvpar .* sign.(v_middle)

    Ekin = E_ev(psd.vpar_centers, c.me_ev)
    E    = E_ev(v_middle, c.me_ev)

    if showplots
            fig, ax, hm = heatmap(Ekin, 
                psd.h_atm/1e3, 
                dfpardvpar[:, :, it], 
                #colorscale = log10,
                colorrange = (-1e-7, 1e-7),
                axis = (xlabel = "kinetic Energy [eV]",
                    ylabel = "Height [km]",
                    #xscale = Makie.Symlog10(0.01),
                    ),
                )
            sleep(2)
            xlims!(-200, 200)
            Colorbar(fig[1, 2], hm)
            #display(fig)
    end

    # select subset (between 3 and 200 eV)
    grad = copy(dfpardvpar[abs.(E) .< 200, :, :])
    v = v_middle[abs.(E) .< 200]
    E    = E_ev(v, c.me_ev)
    grad[abs.(E) .< 3, :, :] .= NaN

    #prepare for saving
    dfdvmax = 0 ./zeros(size(grad));
    iv_dfdvmax = copy(dfdvmax);
    v_dfdvmax = copy(dfdvmax);
    k_growth = copy(dfdvmax);
    gamma_dfdvmax = copy(dfdvmax);

    ## mixed psd and guisdap evaluation
    @time for it in axes(psd.t_run, 1)
        for ih in axes(psd.h_atm, 1)
            showplots = false
            #ih = 395
            #it = 2
            iv_max, h_dfdv_max = findmaxima(grad[:, ih, it])

            #filter for h_dfdv_max larger than 0 (i.e. positive gradient)
            iv_max = iv_max[h_dfdv_max .> 0]
            h_dfdv_max = h_dfdv_max[h_dfdv_max .> 0]

            v_max = v[iv_max]

            if showplots    
                    lines(E, grad[:, ih, it])
            end
            if showplots
                    fig, ax, hm = heatmap(E, psd.h_atm/1e3, grad[:, :, it],
                        colorrange = (-1e-7, 1e-7),)
                    xlims!(-200, 200)
                    Colorbar(fig[1, 2], hm)
                    display(fig)
            end
            if showplots
                fig, ax, lin = lines(E, grad[:, ih, it])
                scatter!(E[iv_max], h_dfdv_max)
                ylims!(-1e-6, 1e-6)
                display(fig)
            end

            dfdvmax[iv_max, ih, it] = h_dfdv_max
            iv_dfdvmax[iv_max, ih, it] = iv_max
            v_dfdvmax[iv_max, ih, it] = v_max

            k_g = wp[ih] ./ sqrt.(v_max.^2 .- 3/2*vth[ih] .^2)
            w_g = k_g .* v_max
            gamma = pi * wp[ih] .^2 ./ k_g .^2 .* h_dfdv_max

            k_growth[iv_max, ih, it] = k_g
            gamma_dfdvmax[iv_max, ih, it] = gamma
            # what about w_g, what is it used for?

            #im = 3
            #vmax = v[iv_max[im]]
            #dfdvmax = h_dfdv_max[im]
            #E_ev(vmax, c.me_ev)
        end
    end

    if false
            showplots = true
            if showplots
                fig, ax, hm = heatmap(E, 
                    psd.h_atm/1e3, 
                    grad[:, :, it], 
                    colorrange = (-1e-7, 1e-7),
                    )
                xlims!(-200, 200)
                Colorbar(fig[1, 2], hm, label = "df∥/d|v|∥")
                E_profiles = E_ev.(v_dfdvmax[:, it], c.me_ev)
                nmax = maximum(length.(E_profiles[:]))
                E_mat = stack(rpad_array.(E_profiles, nmax, NaN))
                [scatter!(Ep, psd.h_atm/1e3, marker = 'x', color = "black") for Ep in eachrow(E_mat)]
                display(fig)
            end

            if showplots
                fig, ax, hm = heatmap(Ekin, 
                    psd.h_atm/1e3, 
                    psd.F[:, :, it], 
                    colorrange = (1e-2, 1e0),
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
                vmax = v_dfdvmax[ih, it][1]
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

    # save results in struct for easy saving
    # (not used with nc files)
    psd_data["k_growth"] = k_growth
    psd_data["v_dfdvmax"] = v_dfdvmax
    psd_data["growth_rate"] = gamma_dfdvmax
    #psd_data


    ## save output
    # 1. copy nc file for backup
    run(`cp -f $(joinpath(psd_dir, psd_file)) $(joinpath(save_dir, psd_file))`)

    # 2. copy psd to psd_le
    if !isdir(joinpath(save_dir, string(interval)))
        mkdir(joinpath(save_dir, string(interval)))
    end
    run(`cp -f $(joinpath(save_dir, psd_file)) $(joinpath(save_dir, string(interval), "psd_le.nc"))`)

    # 3. modify contents
    ds = NCDataset(joinpath(save_dir, string(interval), "psd_le.nc"), "a")

    defDim(ds, "vpar_grad_red", length(v))
    vp_gc_r = defVar(ds, "vpar_gradient_center_reduced", Float64, ("vpar_grad_red",);
                attrib=["units" => "m s-1", "long_name" => "v_parallel bin centers of gradient reduced"])
    vp_gc_r[:] = v
    
    v_k_growth = defVar(ds, "k_growth", Float64, ("vpar_grad_red", "altitude", "time"),
        attrib=["units" => "m-1", "long_name" => "magnitude of k vector"])
    v_k_growth[:] = k_growth[:]
    
    v_gamma_dfdvmax = defVar(ds, "gamma_dfdvmax", Float64, ("vpar_grad_red", "altitude", "time"),
        attrib=["units" => "m-2 s-1 (?)", "long_name" => "growth rate"])
    v_gamma_dfdvmax[:] = gamma_dfdvmax[:]

    defVar(ds, "iv_dfdvmax_", iv_dfdvmax, ("vpar_grad_red", "altitude", "time"),)
    v_iv_dfdvmax_ = defVar(ds, "i_dfdvmax_", Float64, ("vpar_grad_red", "altitude", "time"),
        attrib=["long_name" => "index of maximum dF_par/dv_par"])
    v_iv_dfdvmax_[:] = iv_dfdvmax[:]

    close(ds);

end