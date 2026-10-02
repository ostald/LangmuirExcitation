
function nsqrt(x)
    sign_x = sign(x)
    if sign_x == 0 sign_x = 1 end
    return sqrt(abs(x)) * sign_x
end

function square(x)
    sign_x = sign(x)
    if sign_x == 0 sign_x = 1 end
    return x^2 * sign_x
end

Makie.defaultlimits(::typeof(nsqrt)) = (-Inf, Inf)
Makie.defined_interval(::typeof(nsqrt)) = Makie.OpenInterval(-Inf, Inf)
Makie.inverse_transform(::typeof(nsqrt)) = square


function plot_ne_profile(ne, h; title_add = "")
    #plots electron density ne against height h
    # ne: electron density [m-3]
    # h: height [km], same length as ne
    # title_add: string to add to the title
    fig, ax, ln = lines(ne, h,
        axis = (xlabel = "Electron Density [m⁻³]",
            ylabel = "Height [km]",
            xscale = log10,
            title = "Electron Density Profile " * title_add 
            ),
        )
    display(fig)
    return fig, ax, ln
end

function parse_ex_time(ex_time)
    et2 = lpad.(string.(round.(Int, ex_time)), 2, "0")
    et3 = join(et2[1, 1:3], "-") * " " * join(et2[1, 4:6], ":")
    return et3
end

function plot_median_ne_profile_errorbars(ne_median, h, ne_error; title_add = "")
    fig, ax, lin = scatterlines(ne_median, h, 
        marker = 'x', 
        label = "Median",
        axis = (xscale = log10, 
            limits = ((1e8, 2e12), nothing),
            xlabel = "Electron Density [m⁻³]",
            ylabel = "Height [km]",
            title = "Median Electron Density Profile " * title_add
        ),
    )

    ne_err_lower = copy(ne_error)
    ne_err_lower[(ne_median - ne_error) .< 0] .= ne_median[(ne_median - ne_error) .< 0] .* (1-1e-12)        
    errorbars!(ne_median, h, ne_err_lower, ne_err_lower, 
        direction = :x; color = (:red, 0.8), linewidth = 0.4,
        whiskerwidth = 3, label = "Error") 

    axislegend(ax)
    display(fig)
    return fig, ax, lin
end


function plot_median_Te_profile_errorbars(Te_median, h, Te_error; title_add = "")
    fig, ax, lin = scatterlines(Te_median, h, 
        marker = 'x', 
        label = "Median",
        axis = (xlabel = "Electron Temperature [K]",
            ylabel = "Height [km]",
            title = "Median Electron Temperature Profile " * title_add
        ),
    )

    errorbars!(Te_median, h, Te_error, Te_error, 
        direction = :x; color = (:red, 0.8), linewidth = 0.4,
        whiskerwidth = 3, label = "Error") 

    axislegend(ax)
    display(fig)
    return fig, ax, lin
end

function plot_Fpar_frame(psd, it)
    #replace velocity with energy and scal nsqrt?
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
    Colorbar(fig[1, 2], hm, label = "red. par. velocity distribution [?]")

    if :iv_dfdvmax in keys(psd_grad)
        scatter!(ax,
            psd_grad.v, 
            psd_grad.h_atm/1e3, 
            psd_grad.iv_dfdvmax[:, :, it], 
            marker = 'x',
            color = "black",)
    end

    display(fig)
    return fig, ax, hm        
end

function plot_Fpar_frame_sliders(psd, it)
    #replace velocity with energy and scal nsqrt?
    fig = Figure()
    ax = Axis(fig[1, 1], 
        title = format("t = {:.3f} s", psd.t_run[it]),
        xlabel = "Velocity [m/s]",
        ylabel = "Height [km]",
        #xscale = Makie.Symlog10(0.01),
    )

    # Choose sensible bounds for the sliders
    data_min, data_max = extrema(psd.F[:, :, it])
    data_min = 1e-10

    slider_range = 10 .^ range(log10(data_min), log10(data_max), length = 500)

    rs_v = IntervalSlider(fig[2, 1], range = slider_range,
        startvalues = (data_min, data_max))

    color_range =  rs_v.interval

    hm = heatmap!(ax, 
        psd.vpar_centers, 
        psd.h_atm/1e3, 
        psd.F[:, :, it],
        #colorscale = log10,
        colormap = :viridis, 
        colorrange = color_range,
    )

    Colorbar(fig[1, 2], hm, label = "red. par. velocity distribution [?]")

    if :iv_dfdvmax in keys(psd_grad)
        scatter!(ax,
            psd_grad.v, 
            psd_grad.h_atm/1e3, 
            psd_grad.iv_dfdvmax[:, :, it], 
            marker = 'x',
            color = "black",)
    end

    display(fig)
    return fig, ax, hm
end


function make_movie(psd_grad, save_dir_int; t_range = :)
    tt = psd_grad.t_run[t_range]
    it_iterator = range(1, length(psd_grad.t_run))
    if tt == psd_grad.t_run
        nothing
    else 
        it_iterator = t_range
    end
    framerate = 30*2
    F = Observable(psd_grad.F[:, :, 1])
    dfdv = Observable(psd_grad.dfpardvpar[:, :, 1])
    maxdfdv = Observable(psd_grad.iv_dfdvmax[:, :, 1])

    colorrange_F = (1e-5, 1e2)

    fig, ax, hm = heatmap(E_ev(psd_grad.vpar_centers, c.me_ev),
        psd_grad.h_atm/1e3, 
        F, 
        colorscale = log10,
        colorrange = colorrange_F,
        axis = (#xlabel = "v∥ [eV]",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-5000, 5000), nothing),    # x limit must be specified!
            title = format("t = {:.3f} s", psd_grad.t_run[1]),
        ),
        figure = (size = (600, 600),),
    )
    #xlims!(-400, 400)
    #Colorbar(fig[1, 2], hm)
    #display(fig)

    if :iv_dfdvmax in keys(psd_grad)
        scatter!(ax,
            psd_grad.E, 
            psd_grad.h_atm/1e3, 
            maxdfdv,
            marker = 'x',
            color = "black",)
    end


    ax2, hm2 = heatmap(fig[2,1], E_ev(psd_grad.vpar_centers, c.me_ev),
        psd_grad.h_atm/1e3, 
        F, 
        colorscale = log10,
        colorrange = colorrange_F,
        axis = (#xlabel = "v∥ [eV]",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-400, 400), nothing),    # x limit must be specified!
            #title = format("t = {:.3f} s", psd_grad.t_run[1]),
        ),
    )
    xlims!(-400, 400)
    Colorbar(fig[1:2, 2], hm)

    if :iv_dfdvmax in keys(psd_grad)
        scatter!(ax2,
            psd_grad.E, 
            psd_grad.h_atm/1e3, 
            maxdfdv,
            marker = 'x',
            color = "black",)
    end


    ax3, hm3 = heatmap(fig[3,1],
        E_ev(psd_grad.vpar_centers, c.me_ev),
        psd_grad.h_atm/1e3, 
        dfdv,
        #colorscale = log10,
        colorrange = (-1e-7, 1e-7),
        axis = (xlabel = "v∥ [eV]",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-400, 400), nothing),    # x limit must be specified!
            #title = format("t = {:.3f} s", psd_grad.t_run[1]),
        ),
    )
    Colorbar(fig[3, 2], hm)

    if :iv_dfdvmax in keys(psd_grad)
        scatter!(ax3,
            psd_grad.E, 
            psd_grad.h_atm/1e3, 
            maxdfdv,
            marker = 'x',
            color = "black",)
    end


    record(fig, joinpath(save_dir_int, "fpar_E_v2_500ev.mp4"), it_iterator;
        framerate = framerate) do it
        F[] = psd_grad.F[:, :, it]
        dfdv[] = psd_grad.dfpardvpar[:, :, it]
        maxdfdv[] = psd_grad.iv_dfdvmax[:, :, it]
        ax.title = format("t = {:.3f} s", psd_grad.t_run[it])
    end
end


function plot_thermal_suparthermal_gradient(psd_grad, ih, it)
    #s plot maxwellian distribution for field aligned together with suprathermmal
    # second panel for gradients
    fig, ax, ln = lines(
        E_ev(psd_grad.vpar_centers, c.me_ev),
        psd_grad.F[:, ih, it],
        label = "subrathermal";
        axis = (xscale = nsqrt,
            yscale = log10,
            limits = ((-100.0, 100.0), (1e-1, 1e2)),
            ylabel = "fₑ(v∥) [s m⁻⁴]",
            #xtickformat = values -> ["" for value in values],
        ),
    )
    lines!(ax, 
        E_ev(psd_grad.v_middle, c.me_ev),
        maxwellian_fa(c.me, 2500, v_middle)*2.5e11,
        label = "thermal"
    )
    axislegend(ax)
    
    ax2, ln = lines(fig[2, 1], E_ev(psd_grad.v_middle, c.me_ev),
        gradient_simple(psd_grad.vpar_centers, psd_grad.F[:, ih, it])[2],
        label = "simple nearest neighbour";
        axis = (
            xscale = nsqrt,
            xlabel = "vₑ [eV]",
            ylabel = "df∥/dv∥ [m⁻³]",
            limits = ((-100.0, 100.0), nothing),
        ),
    )
    lines!(ax2, E_ev(psd_grad.vpar_centers, c.me_ev),
        smoothing_gradient(psd_grad.F[:, ih, it], mean(diff(psd_grad.vpar_centers))),
        label = "smooth gradient";
    )
    axislegend(ax2)
    linkxaxes!(ax, ax2)
    display(fig)
    save(joinpath(save_dir_int, "thermal_suparthermal_gradient.jpg"), fig)
end


function plot_dfdv(psd_grad, it)
    Ekin = E_ev(psd_grad.vpar_centers, c.me_ev)
    fig, ax, hm = heatmap(Ekin, 
        psd_grad.h_atm/1e3, 
        psd_grad.dfpardvpar[:, :, it], 
        #colorscale = log10,
        colorrange = (-1e-7, 1e-7),
        axis = (xlabel = "v∥ [eV]",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-200, 100), nothing),    # x limit must be specified!
            title = format("t = {:.3f} s", psd.t_run[it]),
            ),
        )
    Colorbar(fig[1, 2], hm)
    if :iv_dfdvmax in keys(psd_grad)
        scatter!(ax,
            psd_grad.E, 
            psd_grad.h_atm/1e3, 
            psd_grad.iv_dfdvmax[:, :, it], 
            marker = 'x',
            color = "black",)
    end
    display(fig)
    save(joinpath(save_dir_int, "dfdv.jpg"), fig)
    return fig, ax, hm
end


function plot_suprathermal(v_middle)
    scatterlines(v_middle, maxwellian_fa(c.me, 2500, v_middle)*2e11,
        axis = (xticks = [-v_abs.([1, 3, 10, 30, 100, 300, 1000], c.me_ev);
                v_abs.([0, 1, 3, 10, 30, 100, 300, 1000], c.me_ev)
                ],
            xtickformat = values -> ["$(round(Int, E_ev(value, c.me_ev)))" for value in values],
            yscale = log10,
            limits = ((-v_abs.(100.0, c.me_ev), v_abs.(100.0, c.me_ev)), (1e-1, 1e2)),
            xlabel = "vₑ [eV]",
            ylabel = "fₑ(v∥) [s m⁻⁴]",
        )
    )
end