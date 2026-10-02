using MAT
using DataInterpolations
using Peaks
using WGLMakie
using PCHIPInterpolation
using Formatting

include("load_guisdap.jl")
include("utils.jl")
include("plotting_util.jl")
include("constants.jl")
c = physical_constants

## define evaluation parameters
showplots = false
it = 1030
ih = 389

# output directory:
save_dir = "/nfs/revontuli/data/oliver/LangmuirExcitation/results"

# intput directory:
psd_dir = "/nfs/revontuli/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_train_330-345s_pchip/analysis"
psd_file = "psd.nc"
# copy nc file for backup
#run(`cp -f $(joinpath(psd_dir, psd_file)) $(joinpath(save_dir, psd_file))`)

#define guisdap data
guisdap_dir = "analyzed_parameters/2018-12-07_folke_6.4@42mb"
files = filter(x -> x[end-3:end] == ".mat", readdir(guisdap_dir))

# 3 different time intervals of guisdap to try out:
interval1 = 143:161
interval2 = 292:310
interval3 = 246:264

#for interval in [interval1, interval2, interval3]
    interval = interval3
    save_dir_int = joinpath(save_dir, string(interval))
    load_files = files[interval]

    if showplots
        file = files[290]
        gd = load_guisdap(joinpath(guisdap_dir, file))
        ex_time, h_param, Te, ne = gd
        println(ex_time)
        fig, ax, ln = plot_ne_profile(ne, h_param/1e3, 
            title_add = parse_ex_time(ex_time)
        )
        save(joinpath(save_dir_int, "ne_profile.jpg"), fig)
    end

    ## guisdap data part
    #load guisdap data
    t_start, t_end, h_median, ne_median, Te_median, h_mad, ne_mad, Te_mad = 
        average_guisdap(joinpath.(guisdap_dir, load_files))
    gd = load_guisdap.(joinpath.(guisdap_dir, load_files))

    # load psd data
    psd_data = load_psd(joinpath(psd_dir, psd_file))
    psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data))
    h_atm = psd.h_atm

    #interpolate guisdap data to psd height resolution
    ne_itp = Interpolator(h_median, ne_median).(h_atm)
    Te_itp = Interpolator(h_median, Te_median).(h_atm)

    if showplots
        ne_gd = vcat([m.ne for m in gd]...)
        h_gd = vcat([m.h_param/1e3 for m in gd]...)

        fig, ax, lin = plot_median_ne_profile_errorbars(
            ne_median, h_median/1e3, ne_mad,
            title_add = "\n"*parse_ex_time(gd[1].ex_time) 
                * " | " * parse_ex_time(gd[end].ex_time)
            )
        save(joinpath(save_dir_int, "median_ne_density.jpg"), fig)
                        
        scatter!(ne_gd, h_gd, color = "black", alpha = 0.2, markersize = 3,
            label = "Measured")

        lines!(ax, ne_itp, h_atm/1e3, linestyle = :dash, color = "orange",
            label = "Interp.")
        axislegend(ax)
        #display(fig)
    end

    if showplots
        Te_gd = vcat([m.Te for m in gd]...)
        h_gd = vcat([m.h_param/1e3 for m in gd]...)

        fig, ax, lin = plot_median_Te_profile_errorbars(
            Te_median, h_median/1e3, Te_mad,
            title_add = "\n"*parse_ex_time(gd[1].ex_time) 
                * " | " * parse_ex_time(gd[end].ex_time)
            )
        save(joinpath(save_dir_int, "median_Te_density.jpg"), fig)
        
        scatter!(Te_gd, h_gd, color = "black", alpha = 0.2, markersize = 3,
            label = "Measured")

        lines!(ax, Te_itp, h_atm/1e3, linestyle = :dash, color = "orange",
            label = "Interp.")
        axislegend(ax)
        #display(fig)
    end

    # calculate plasma frequency and thermal velocity 
    # necessary for resonance condition and growth rate
    wp = plasma_freq.(ne_itp, c.me)
    vth = thermal_velocity.(Te_itp, c.me)


    ## psd data part
    # calculate gradient in parellel velocity 
    @time dfpardvpar = stack(
        [diff(psd.F[:, ih, it]) ./ diff(psd.vpar_centers)
            for ih in eachindex(psd.h_atm),
                it in eachindex(psd.t_run)
        ]
    );    
    v_middle = psd.vpar_centers[1:end-1] .+ diff(psd.vpar_centers)/2
    sign_v = sign.(v_middle)
    sign_v[sign_v .== 0] .= 1
    dfpardvpar = dfpardvpar .* sign_v

    psd_grad = (psd..., v_middle = v_middle, dfpardvpar = dfpardvpar)


    #showplots = true
    if showplots
        it = 1030
        ih = 389

        fig, ax, hm = plot_Fpar_frame(psd, it)
        save(joinpath(save_dir_int, "example_F.jpg"), fig)


        plot_Fpar_frame_sliders(psd, it)

        t_range = : ;
        t_range = 1100:1200 ;
        make_movie(psd_grad, save_dir_int; t_range = t_range)

    end

    #=
    @time dfpardvpar = stack(
        [gradient_simple(psd.vpar_centers, psd.F[:, ih, it])[2]
            for ih in eachindex(psd.h_atm),
                it in eachindex(psd.t_run[1:1000])
            ]
        );
    =#

    if showplots
        plot_suprathermal(v_middle)
    end

    if showplots
        it = 1030
        ih = 389
        plot_thermal_suparthermal_gradient(psd_grad, ih, it)
    end

    
    Ekin = E_ev(psd_grad.vpar_centers, c.me_ev)
    E    = E_ev(psd_grad.v_middle, c.me_ev)


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

    if false #old code, incompatible now
            showplots = true
            #plotting gradiend dfdv heatmap, superimposed with scatter of where dfdv has a maximum
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

            #plotting phase space F heatmap, superimposed with scatter of where dfdv has a maximum
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

            # intersection of langmuir dispersion relation and ??
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
    psd_data["v"] = v
    psd_data["k_growth"] = k_growth
    psd_data["v_dfdvmax"] = v_dfdvmax
    psd_data["gamma_dfdvmax"] = gamma_dfdvmax
    psd_data["iv_dfdvmax"] = iv_dfdvmax

    #psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data))

    #psd_data
    psd_grad = (psd...,
        v_middle = v_middle, 
        dfpardvpar = dfpardvpar,
        grad = grad,
        v = v,
        E = E,
        dfdvmax = dfdvmax,  
        iv_dfdvmax = iv_dfdvmax,
        v_dfdvmax = v_dfdvmax,
        k_growth = k_growth,
        gamma_dfdvmax = gamma_dfdvmax,
        )


    ## save output
    # 2. copy psd to psd_les
    if !isdir(joinpath(save_dir_int))
        mkdir(joinpath(save_dir_int))
    end
    run(`cp -f $(joinpath(save_dir, psd_file)) $(joinpath(save_dir_int, "psd_le.nc"))`)

    # 3. modify contents
    ds = NCDataset(joinpath(save_dir_int, "psd_le.nc"), "a")

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

    v_iv_dfdvmax_ = defVar(ds, "iv_dfdvmax", Float64, ("vpar_grad_red", "altitude", "time"),
        attrib=["long_name" => "index of maximum dF_par/dv_par"])
    v_iv_dfdvmax_[:] = iv_dfdvmax[:]

    close(ds);


    #if showplots
        plot_Fpar_frame_sliders(psd_grad, it)

        fig, ax, hm = plot_dfdv(psd_grad, it)
        


    end

end


# to do:
# 1. also save guisdap data in nc file
# 2. restrict velocity/energy interval from the beginning?
#   => related: only save dfpardvpar for small velocities? ie only grad
#   => related: rename grad
# 3. add timestamp to all plots that need it
# 4. check maximum detection. looks iffy
#   => realted: put gradient on a colorscale that is diverging from 0 (e.g. :bam), so 0 is easily found
#       maybe the maximum detetction is ok, it's just that some maxima are below 0 and therefore not contributing to wave growth? double check
# 5. add "if key exists in psd then draw that too"
