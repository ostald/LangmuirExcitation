using MAT
using DataInterpolations
using Peaks
using CairoMakie
using PCHIPInterpolation
using Formatting

#using WGLMakie
#WGLMakie.activate!()

include("load_guisdap.jl")
include("utils.jl")
include("plotting_util.jl")
include("constants.jl")
c = physical_constants

## define evaluation parameters
showplots = true
evaluate = true

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
interval2 = 246:264
interval3 = 292:310


if evaluate == true
for interval in [interval1, interval2, interval3]
#    interval = interval3
    println("processing interval ", string(interval))
    showplots = false

    save_dir_int = joinpath(save_dir, string(interval))
    if !isdir(joinpath(save_dir_int))
        println("Making directory ", save_dir_int)
        mkdir(joinpath(save_dir_int))
    end

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
    println("Loading guisdap files")
    t_start, t_end, h_median, ne_median, Te_median, h_mad, ne_mad, Te_mad = 
        average_guisdap(joinpath.(guisdap_dir, load_files))
    gd = load_guisdap.(joinpath.(guisdap_dir, load_files))

    # load psd data
    println("Loading psd data")
    psd = netCDF_to_namedTuple(joinpath(psd_dir, psd_file))

    #interpolate guisdap data to psd height resolution
    ne_itp = Interpolator(h_median, ne_median).(psd.altitude)
    Te_itp = Interpolator(h_median, Te_median).(psd.altitude)

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

        lines!(ax, ne_itp, psd.altitude/1e3, linestyle = :dash, color = "orange",
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

        lines!(ax, Te_itp, psd.altitude/1e3, linestyle = :dash, color = "orange",
            label = "Interp.")
        axislegend(ax)
        #display(fig)
    end

    # calculate plasma frequency and thermal velocity 
    # necessary for resonance condition and growth rate
    wp = plasma_freq.(ne_itp, c.me)
    vth = thermal_velocity.(Te_itp, c.me)


    ## psd data part

    #cut data where Energy is in thermal range and replace F'
    E = E_ev(psd.vpar_centers, c.me_ev)
    F_st = psd.F
    F_st[abs.(E) .< 2, :, :] .= NaN 
    psd = (psd..., F_st = F_st)

    # calculate gradient in parellel velocity 
    println("Calculating gradient")
    dfpardvpar = stack(
        [diff(psd.F_st[:, ih, it]) ./ diff(psd.vpar_centers)
            for ih in eachindex(psd.altitude),
                it in eachindex(psd.time)
        ]
    );   

    v_dfdv      = psd.vpar_centers[1:end-1] .+ diff(psd.vpar_centers)/2
    sign_v_dfdv = sign.(v_dfdv)
    sign_v_dfdv[sign_v_dfdv .== 0] .= 1
    dfpardvpar  = dfpardvpar .* sign_v_dfdv
    E_kin       = E_ev(v_dfdv, c.me_ev)


    # select subset (between 3 and 200 eV)
    dfdv_s  = copy(dfpardvpar[abs.(E_kin) .< 200, :, :])
    v_s     = v_dfdv[abs.(E_kin) .< 200]
    E_kin_s = E_ev(v_s, c.me_ev)
    dfdv_s[abs.(E_kin_s) .< 3, :, :] .= NaN

    #prepare for saving
    dfdvmax = 0 ./zeros(size(dfdv_s));
    iv_growth = copy(dfdvmax);
    v_growth = copy(dfdvmax);
    k_growth = copy(dfdvmax);
    gamma_growth = copy(dfdvmax);

    ## mixed psd and guisdap evaluation
    println("Evaluate excitation of Langmuir waves")
    for it in axes(psd.time, 1)
        for ih in axes(psd.altitude, 1)
            showplots = false
            #ih = 395
            #it = 2
            iv_max, h_dfdv_max = findmaxima(dfdv_s[:, ih, it])

            #filter for h_dfdv_max larger than 0 (i.e. positive gradient)
            iv_max = iv_max[h_dfdv_max .> 0]
            h_dfdv_max = h_dfdv_max[h_dfdv_max .> 0]

            v_max = v_s[iv_max]

            if showplots    
                    lines(E_kin_s, dfdv_s[:, ih, it])
            end
            if showplots
                    fig, ax, hm = heatmap(E_kin_s, psd.altitude/1e3, dfdv_s[:, :, it],
                        colorrange = (-1e-7, 1e-7),)
                    xlims!(-200, 200)
                    Colorbar(fig[1, 2], hm)
                    display(fig)
            end
            if showplots
                fig, ax, lin = lines(E_kin_s, dfdv_s[:, ih, it])
                scatter!(E_kin_s[iv_max], h_dfdv_max)
                ylims!(-1e-6, 1e-6)
                display(fig)
            end

            dfdvmax[iv_max, ih, it] = h_dfdv_max
            iv_growth[iv_max, ih, it] = iv_max
            v_growth[iv_max, ih, it] = v_max

            k_g = wp[ih] ./ sqrt.(v_max.^2 .- 3/2*vth[ih] .^2)
            w_g = k_g .* v_max
            gamma = pi * wp[ih] .^2 ./ k_g .^2 .* h_dfdv_max
            #growth rate formula see eg.
            # Bump-on-tail instability in space plasmas
            # Susmita Sarkar ; Samit Paul; Raicharan Denra

            k_growth[iv_max, ih, it] = k_g
            gamma_growth[iv_max, ih, it] = gamma
            # what about w_g, what is it used for?

            #im = 3
            #vmax = v[iv_max[im]]
            #dfdvmax = h_dfdv_max[im]
            #E_ev(vmax, c.me_ev)
        end
    end


    # save results in named tuple for easy access
    psd_le = (psd...,
        v_dfdv = v_dfdv, 
        dfpardvpar = dfpardvpar,
        dfdv_s = dfdv_s,
        v_s = v_s,
        E_kin_s = E_kin_s,
        dfdvmax = dfdvmax,  
        iv_growth = iv_growth,
        v_growth = v_growth,
        k_growth = k_growth,
        gamma_growth = gamma_growth,
        )


    
    ## save output
    println("saving")
    # 2. copy psd to psd_les
    run(`cp -f $(joinpath(save_dir, psd_file)) $(joinpath(save_dir_int, "psd_le.nc"))`)

    # 3. modify contents
    ds = NCDataset(joinpath(save_dir_int, "psd_le.nc"), "a")

    defDim(ds, "vpar_dfdv_s", length(v_s))

    vp_gc_r = defVar(ds, "v_s", Float64, ("vpar_dfdv_s",);
                attrib=["units" => "m s-1", "long_name" => "v_parallel bin centers of gradient reduced"])
    vp_gc_r[:] = v_s

    ncvar_dfdv_s = defVar(ds, "dfdv_s", Float64, ("vpar_dfdv_s", "altitude", "time"),
                attrib=["units" => "s2 m-5", "long_name" => "gradient dF_par/dv_par"])
    ncvar_dfdv_s[:] = dfdv_s
    
    ncvar_k_growth = defVar(ds, "k_growth", Float64, ("vpar_dfdv_s", "altitude", "time"),
        attrib=["units" => "m-1", "long_name" => "magnitude of k vector"])
    ncvar_k_growth[:] = k_growth[:]
    
    ncvar_gamma_growth = defVar(ds, "gamma_growth", Float64, ("vpar_dfdv_s", "altitude", "time"),
        attrib=["units" => "m-3", "long_name" => "growth rate"])
    ncvar_gamma_growth[:] = gamma_growth[:]

    ncvar_iv_growth = defVar(ds, "iv_growth", Float64, ("vpar_dfdv_s", "altitude", "time"),
        attrib=["long_name" => "index of maximum dF_par/dv_par"])
    ncvar_iv_growth[:] = iv_growth[:]

    close(ds);

end
end



#if showplots == true
println("Evaluation finished, starting plotting")


# output directory:
save_dir = "/nfs/revontuli/data/oliver/LangmuirExcitation/results"

# 3 different time intervals of guisdap to try out:
interval1 = 143:161
interval2 = 246:264
interval3 = 292:310

radar_freq = 500e6 #Hz
klim = (radar_freq * 2 *pi)/c.c


for interval in [interval1, interval2, interval3]
#if true
    #interval = interval1
    println("processing interval ", string(interval))
    save_dir_int = joinpath(save_dir, string(interval))

    #psd = netCDF_to_namedTuple(joinpath(save_dir, "psd.nc"))
    #psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data))

    println("loading netcdf file")
    psd = netCDF_to_namedTuple(joinpath(save_dir_int, "psd_le.nc"))
    #psd_le = NamedTuple{Tuple(Symbol(k) for k in keys(ds))}(Array(ds[k]) for k in keys(ds))
    #close(ds);

    E = E_ev(psd.vpar_centers, c.me_ev)
    F_st = psd.F
    F_st[abs.(E) .< 2, :, :] .= NaN 
    psd = (psd..., F_st = F_st)

    println("Calculating df/dv")
    #reconstruct some variables
    dfpardvpar = stack(
        [diff(psd.F[:, ih, it]) ./ diff(psd.vpar_centers)
            for ih in eachindex(psd.altitude),
                it in eachindex(psd.time)
        ]
    );

    v_dfdv      = psd.vpar_centers[1:end-1] .+ diff(psd.vpar_centers)/2
    sign_v_dfdv = sign.(v_dfdv)
    sign_v_dfdv[sign_v_dfdv .== 0] .= 1
    dfpardvpar  = dfpardvpar .* sign_v_dfdv
    E_kin       = E_ev(v_dfdv, c.me_ev)
    E_kin_s     = E_ev(psd.v_s, c.me_ev)


    #psd2 = merge(psd, psd_le)
    psd_le = (psd...,
        v_dfdv = v_dfdv,
        dfpardvpar = dfpardvpar, 
        E_kin_s = E_kin_s,
        )

    println("plotting")
    showplots = true
#    if showplots
        it = 1030
        ih = 389


        fig, ax, hm = plot_Fpar_frame(psd_le, it)
        save(joinpath(save_dir_int, "example_F.png"), fig)
        
        plot_Fpar_frame_sliders(psd_le, it)

        fig, ax, hm = plot_dfdv(psd_le, it)
        save(joinpath(save_dir_int, "dfdv.png"), fig)

        t_range = : ;
        #t_range = 1100:1200 ;
        make_movie(psd_le, save_dir_int; t_range = t_range)

        v_dfdv      = psd.vpar_centers[1:end-1] .+ diff(psd.vpar_centers)/2
        #plot_suprathermal(v_dfdv)

        fig, ax, ln = plot_thermal_suparthermal_gradient(psd_le, ih, it)
        save(joinpath(save_dir_int, "thermal_suparthermal_gradient.png"), fig)

        fig, ax, hm = plot_k_hm(psd_le, it, klim)
        save(joinpath(save_dir_int, "k_vs_v_and_h_hm.png"), fig)

        fig, ax, sc = plot_k_scatter(psd_le, it, klim)
        save(joinpath(save_dir_int, "k_h.png"), fig)

        fig, ax, hm = k_t_h(psd_le, klim)
        save(joinpath(save_dir_int, "k_t_h.png"), fig)
        #make it interactive!!

        fig, ax, hm = g_t_h(psd_le)
        save(joinpath(save_dir_int, "g_t_h.png"), fig)

  #  end

end
end




# to do:
# 1. also save guisdap data in nc file
# 2. restrict velocity/energy interval from the beginning?
#   => related: only save dfpardvpar for small velocities? ie only grad? yes!!
# 3. add timestamp to all plots that need it
# 4. check maximum detection. looks iffy
#   => related: put gradient on a colorscale that is diverging from 0 (e.g. :bam), so 0 is easily found
#       maybe the maximum detetction is ok, it's just that some maxima are below 0 and therefore not contributing to wave growth? double check
# 5. add "if key exists in psd then draw that too"
# 6. make nice x axis all the way through all plots (Energy and nsqrt)
# 7. check thermal distribtution and cutoff energy
# 8. make k number plot
# 9. plot_k_scatter: plot all available k, not just maximum
#10. mask out 3eV threshold for low energy