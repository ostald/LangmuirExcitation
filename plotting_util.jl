using MathTeXEngine
set_texfont_family!(FontFamily("TeXGyreHeros"))

using NaNStatistics

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


function add_line_colorbar(value, fig, cb, colorrange)
    rect = lift(cb.blockscene.plots[2][1], cb.blockscene.plots[2][2]) do (xmin, xmax), (ymin, ymax)
        if any(!isfinite, (xmin, xmax, ymin, ymax))
            return BBox(0.0, 1.0, 0.0, 1.0)
        end
        BBox(xmin, xmax, ymin, ymax)
    end

    inside_ax = Axis(fig,
        bbox = rect,
        limits = (0, 1, colorrange...),
        yscale = cb.scale
    )

    hidedecorations!(inside_ax)
    hidespines!(inside_ax)
    hlines!(inside_ax, value, color = :red, linewidth = 3)
    #return fig, ax, hm
end


function rescale_colorbar_ticks!(fig, cb, power::Int; digits = 2)
    # Get the current tick values from the colorbar
    tickvals = cb.axis.tickvalues.val

    # If there are no ticks, stop
    isempty(tickvals) && return cb

    # Rescale the tick labels
    scaled = tickvals ./ 10^float(power)
    labels = [string(round(x; digits = digits)) for x in scaled]

    # Apply the relabeled ticks while keeping the same positions
    cb.ticks = (tickvals, labels)

    # Put the multiplier above the colorbar
    Label(fig[1, 2, Top()],
        rich("×10", superscript(string(power))), 
        halign = :left,
        padding = (0.0f0, 0.0f0, 7.0f0, 0.0f0),
    )

    return cb
end


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
    # parses guisdap timestamps in a string
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
    E_kin = E_ev(psd.vpar_centers, c.me_ev)
    fig, ax, hm = heatmap(E_kin,
        psd.altitude/1e3, 
        psd.F_st[:, :, it], 
        colorscale = log10,
        colorrange = (1e-2, 1e2),
        axis = (xlabel = L"\mathrm{v_{\parallel}\;[eV]}",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = (extrema(E_kin), nothing),
            title = format("t = {:.3f} s", psd.time[it]),
            xticks = [-4000, -2000, -500, 0, 500, 2000, 4000],
            xminorticks = -5000:500:5000,
            xminorticksvisible = true,
        ),
    )

    Colorbar(fig[1, 2], hm, label = L"\mathrm{F_{\parallel}\;[s\,m^{-4}]}")

    if :iv_growth in keys(psd)
        scatter!(ax,
            E_ev(psd.v_s, c.me_ev), 
            psd.altitude/1e3, 
            psd.iv_growth[:, :, it], 
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
        title = format("t = {:.3f} s", psd.time[it]),
        xlabel = L"\mathrm{v_\parallel \,  [eV]}",
        ylabel = "Height [km]",
        xscale = nsqrt,
        limits = (extrema(E_ev(psd.vpar_centers, c.me_ev)), nothing),
        #xticks = [-4000, -2000, -500, 0, 500, 2000, 4000],
        #xminorticks = -5000:500:5000,
        xticks = WilkinsonTicks(5, k_min = 5),
        xminorticksvisible = true,
    )

    # Choose sensible bounds for the sliders
    data_min, data_max = nanextrema(psd.F_st[:, :, it])
    data_min = 1e-10

    slider_range = 10 .^ range(log10(data_min), log10(data_max), length = 500)

    rs_v = IntervalSlider(fig[2, 1], range = slider_range,
        startvalues = (data_min, data_max))

    color_range =  rs_v.interval

    hm = heatmap!(ax, 
        E_ev(psd.vpar_centers, c.me_ev),
        psd.altitude/1e3, 
        psd.F_st[:, :, it],
        colorscale = log10,
        colormap = :viridis, 
        colorrange = color_range,
    )

    Colorbar(fig[1, 2], hm, label = L"\mathrm{F_{\parallel}\;[s\,m^{-4}]}")

    if :iv_growth in keys(psd)
        scatter!(ax,
            E_ev(psd.v_s, c.me_ev),
            psd.altitude/1e3, 
            psd.iv_growth[:, :, it], 
            marker = 'x',
            color = "black",)
    end

    display(fig)
    return fig, ax, hm
end


function make_movie(psd_grad, save_dir_int; t_range = :)
    tt = psd_grad.time[t_range]
    it_iterator = range(1, length(psd_grad.time))
    if tt == psd_grad.time
        nothing
    else 
        it_iterator = t_range
    end
    framerate = 30*2
    F = Observable(psd_grad.F_st[:, :, it_iterator[1]])
    dfdv = Observable(psd_grad.dfpardvpar[:, :, it_iterator[1]])
    maxdfdv = Observable(psd_grad.iv_growth[:, :, it_iterator[1]])
    k1 = Observable(psd_grad.k_growth[:, :, it_iterator[1]])
    k2 = Observable(vec(psd_grad.k_growth[:, :, it_iterator[1]]))
    gam = Observable(vec(psd_grad.gamma_growth[:, :, it_iterator[1]]))

    colorrange_F = (1e-5, 1e2)

    fig, ax, hm = heatmap(E_ev(psd_grad.vpar_centers, c.me_ev),
        psd_grad.altitude/1e3, 
        F, 
        colorscale = log10,
        colorrange = colorrange_F,
        axis = (#xlabel = "v∥ [eV]",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-5000, 5000), nothing),    # x limit must be specified!
            title = format("t = {:.3f} s", psd_grad.time[1]),
            xticks = [-5000, -2500, -1000, -200, 0, 200, 1000, 2500, 5000],
            xminorticksvisible = true,
            xminorticks = vcat(vcat(0:50:190, 200:200:1000, 1000:500:5000),
                vcat(0:50:190, 200:200:1000, 1000:500:5000) .*-1),
        ),
        figure = (size = (600, 800),),
    )
    vlines!(ax, [-400, 400], linestyle = :dash, color = "white")
    #xlims!(-400, 400)
    #Colorbar(fig[1, 2], hm)
    #display(fig)

    if :iv_growth in keys(psd_grad)
        scatter!(ax,
            psd_grad.E_kin_s, 
            psd_grad.altitude/1e3, 
            maxdfdv,
            marker = 'x',
            color = "black",)
    end

    ax2, hm2 = heatmap(fig[2,1], E_ev(psd_grad.vpar_centers, c.me_ev),
        psd_grad.altitude/1e3, 
        F, 
        colorscale = log10,
        colorrange = colorrange_F,
        axis = (#xlabel = "v∥ [eV]",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-400, 400), nothing),    # x limit must be specified!
            xticks = [-400, -200, -50, -10, 0, 10, 50, 200, 400],
            #xticklabelsvisible = false,
            xminorticksvisible = true,
            xminorticks = vcat(vcat(0:2:9, 10:10:49, 50:50:400),
                vcat(0:2:9, 10:10:49, 50:50:400) .*-1),
        ),
    )
    #xlims!(-400, 400)
    Colorbar(fig[1:2, 2], hm, label = L"\mathrm{F_{\parallel} \; [s \, m^{-4}]}")

    if :iv_growth in keys(psd_grad)
        scatter!(ax2,
            psd_grad.E_kin_s, 
            psd_grad.altitude/1e3, 
            maxdfdv,
            marker = 'x',
            color = "black",)
    end

    cm = Makie.to_colormap(:vik);

    ax3, hm3 = heatmap(fig[3,1],
        E_ev(psd_grad.vpar_centers, c.me_ev),
        psd_grad.altitude/1e3, 
        dfdv,
        colormap = cm[60:end-60],
        #colorscale = log10,
        colorrange = (-1.001e-7, 1.001e-7),
        axis = (#xlabel = L"\mathrm{v_{\parallel} \; [eV]}",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-400, 400), nothing),    # x limit must be specified!
            xticks = [-400, -200, -50, -10, 0, 10, 50, 200, 400],
            #xticklabelsvisible = false,
            xminorticksvisible = true,
            xminorticks = vcat(vcat(0:2:9, 10:10:49, 50:50:400), vcat(0:2:9, 10:10:49, 50:50:400) .*-1),
        ),
    )
    cb = Colorbar(fig[3, 2], hm3, label = L"\mathrm{dF_{\parallel}/dv_{\parallel} \; [s^2 \, m^{-5}]}")

    #rescale_colorbar_ticks!(fig[3, 2], cb, -7; digits = 2)

    if :iv_growth in keys(psd_grad)
        scatter!(ax3,
            psd_grad.E_kin_s, 
            psd_grad.altitude/1e3, 
            maxdfdv,
            marker = 'x',
            color = "black",)
    end

    ax4, hm4 = heatmap(fig[4,1],
        E_ev(psd_grad.v_s, c.me_ev),
        psd_grad.altitude/1e3, 
        k1,
        colorrange = (0, nanmaximum(psd_grad.k_growth)/2),
        axis = (xlabel = L"\mathrm{v_{\parallel} \; [eV]}",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-400, 400), nothing),    # x limit must be specified!
            xticks = [-400, -200, -50, -10, 0, 10, 50, 200, 400],
            xminorticksvisible = true,
            xminorticks = vcat(vcat(0:2:9, 10:10:49, 50:50:400),
                vcat(0:2:9, 10:10:49, 50:50:400) .*-1),
        ),
    )
    cb = Colorbar(fig[4, 2], hm4, label = L"\mathrm{k \; [m^{-1}]}")
    add_line_colorbar(klim, fig[4, 2], cb, (0, nanmaximum(psd_grad.k_growth)/2))


    nk = size(psd_grad.k_growth, 1)
    maxk = nanmaximum(psd_grad.k_growth)
    ax5, sc5 = scatter(fig[5,1],
        k2,
        repeat(psd_grad.altitude'/1e3, outer = (nk, 1))[:],
        color = gam,
        colorscale = log10,
        colorrange = (1e4, 1e8),
        axis = (xlabel = L"\mathrm{k \; [m^{-1}]}",
            ylabel = "Height [km]",
            #title = format("t = {:.3f} s", psd_grad.time[it]),
            limits = ((maxk*(-0.04), maxk*1.1), 
                (psd_grad.altitude[1]/1e3, psd_grad.altitude[end]/1e3)),
            xticks = WilkinsonTicks(5, k_min = 5),
            xminorticksvisible = true,
        ),
    )
    #scatter!(ax5, 0, psd_grad.altitude[1]/1e3, alpha = 0)
    #scatter!(ax5, maximum(psd_grad.k_growth)*1.1, psd_grad.altitude[end]/1e3, alpha = 0)
    vlines!(ax5, klim, color = :red, linewidth = 3)
    cb5 = Colorbar(fig[5, 2], sc5, label = L"\mathrm{Growth Rate \; [m^{-3}]}")


    record(fig, joinpath(save_dir_int, "fpar_E_v2_500ev.mp4"), it_iterator;
        framerate = framerate) do it
        F[] = psd_grad.F_st[:, :, it]
        dfdv[] = psd_grad.dfpardvpar[:, :, it]
        maxdfdv[] = psd_grad.iv_growth[:, :, it]
        k1[] = psd_grad.k_growth[:, :, it]
        k2[] = vec(psd_grad.k_growth[:, :, it])
        gam[] = vec(psd_grad.gamma_growth[:, :, it])

        ax.title = format("t = {:.3f} s", psd_grad.time[it])
    end
end


function plot_k_hm(psd_grad, it, klim)
    fig, ax, hm = heatmap(
        E_ev(psd_grad.v_s, c.me_ev),
        psd_grad.altitude/1e3, 
        psd_grad.k_growth[:, :, it],
        axis = (xlabel = L"\mathrm{v_{\parallel} \; [eV]}",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-400, 400), nothing),    # x limit must be specified!
            title = format("t = {:.3f} s", psd_grad.time[it]),
        ),
    )
    cb = Colorbar(fig[1, 2], hm, label = L"\mathrm{k \; [m^{-1}]}")
    add_line_colorbar(klim, fig, cb, (0, nanmaximum(psd_grad.k_growth[:, :, it])))
    display(fig)
    return fig, ax, hm
end

function plot_k_scatter(psd_grad, it, klim)
    nk = size(psd_grad.k_growth, 1)
    fig, ax, sc = scatter(
        psd_grad.k_growth[:, :, it][:],
        repeat(psd_grad.altitude'/1e3, outer = (nk, 1))[:],
        color = psd_grad.gamma_growth[:, :, it][:],
        colorscale = log10,
        colorrange = (1e4, 1e8),
        axis = (xlabel = L"\mathrm{k \; [m^{-1}]}",
            ylabel = "Height [km]",
            title = format("t = {:.3f} s", psd_grad.time[it]),
            #limits = ((nothing, maximum(k)*1.1), (psd_grad.altitude[1]/1e3, psd_grad.altitude[end]/1e3))
        ),
    )
    scatter!(ax, 0, psd_grad.altitude[1]/1e3, alpha = 0)
    scatter!(ax, maximum(psd_grad.k_growth)*1.1, psd_grad.altitude[end]/1e3, alpha = 0)
    vlines!(ax, klim, color = :red, linewidth = 3)

    cb = Colorbar(fig[1, 2],
        sc, 
        label = L"\mathrm{Growth \, Rate \; [m^{-3}]}",
        )

    display(fig)
    return fig, ax, sc
end

function plot_histogram_growth_log(psd_le)
    g_no_nan = psd_le.gamma_growth[.!isnan.(psd_le.gamma_growth)]
    
    #=
    fig, ax, hist = WGLMakie.hist(g_no_nan, bins = 20,
        axis = (
            ylabel = "Counts",
            xlabel = "Growth Rate ",
            #xscale = log10,
            #yscale = log10,
            #limits = ((1e6, 1e9), (1e2, 1e6)),
        ),
    )

    ax2, hist2 = WGLMakie.hist(fig[1, 2], g_no_nan, bins = 20,
        label = "the first bar should extend to 0",
        axis = (
            ylabel = "Counts",
            xlabel = "Growth Rate [m⁻³]",
            xscale = log10,
            yscale = log10,
            limits = ((1e6, 1e9), (1e2, 1e6)),
            title = "this is a shitty plot"

        ),
    )
    axislegend(ax2)
    =#

    fig, ax, hist = WGLMakie.hist(log10.(g_no_nan), bins = 100,
        axis = (
            ylabel = "Counts",
            xlabel = "Growth Rate [m⁻³]",
            xtickformat = values -> [rich("10", superscript("$(round(Int, value))")) for value in values],
            #xscale = log10,
            yscale = log10,
            #limits = ((1e6, 1e9), (1e2, 1e6)),
        ),
    )

    ax4, hist4 = WGLMakie.hist(fig[2, 1], log10.(g_no_nan[g_no_nan .> 1e-5]), bins = 100,
        axis = (
            ylabel = "Counts",
            xlabel = "Growth Rate [m⁻³]",
            xtickformat = values -> [rich("10", superscript("$(round(Int, value))")) for value in values],
            #xscale = log10,
            yscale = log10,
            #limits = ((1e-1, nothing), nothing),
        ),
    )

    #=
    WGLMakie.hist(log10.(g_no_nan[g_no_nan .> 1e-5]), bins = 80,
        axis = (
            ylabel = "Counts",
            xlabel = "Growth Rate [m⁻³]",
            xtickformat = values -> [rich("10", superscript("$(round(Int, value))")) for value in values],
            #xscale = log10,
            yscale = log10,
            #limits = ((1e6, 1e9), (1e2, 1e6)),
        ),
    )
    =#

    display(fig)
    return fig
end

function plot_thermal_suparthermal_gradient(psd_grad, ih, it)
    #s plot maxwellian distribution for field aligned together with suprathermmal
    # second panel for gradients
    fig, ax, ln = lines(
        E_ev(psd_grad.vpar_centers, c.me_ev),
        psd_grad.F_st[:, ih, it],
        label = "suprathermal";
        axis = (xscale = nsqrt,
            yscale = log10,
            limits = ((-100.0, 100.0), (1e-1, 1e2)),
            ylabel = L"\mathrm{F_{\parallel} \; [s \, m^{-4}]}",
            #xtickformat = values -> ["" for value in values],
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(5),
            title = "$(psd_grad.time[it]) s, $(round(psd_grad.altitude[ih]/1e3)) km"
        ),
    )
    lines!(ax, 
        E_ev(psd_grad.v_s, c.me_ev),
        maxwellian_fa(c.me, 2500, psd_grad.v_s)*2.5e11,
        label = "thermal"
    )
    axislegend(ax)

    ax2, ln = lines(fig[2, 1], E_ev(psd_grad.v_dfdv, c.me_ev),
        psd_grad.dfpardvpar[:, ih, it],
        label = "simple nearest neighbour";
        axis = (
            xscale = nsqrt,
            xlabel = L"\mathrm{v_{\parallel} \; [eV]}",
            ylabel = L"\mathrm{dF_{\parallel}/dv_{\parallel} \; [s^2 \, m^{-5}]}",
            limits = ((-100.0, 100.0), nothing),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(5)
        ),
    )

    smooth_gr =  smoothing_gradient(psd_grad.F_st[:, ih, it], mean(diff(psd_grad.vpar_centers)))
    sig = sign.(psd_grad.vpar_centers)
    sig[sig .== 0] .= 1
    smooth_gr = smooth_gr .* sig

    lines!(ax2, E_ev(psd_grad.vpar_centers, c.me_ev),
        smooth_gr,
        label = "smooth gradient";
    )
    axislegend(ax2)
    linkxaxes!(ax, ax2)
    display(fig)
    return fig, ax, ln
end


function plot_dfdv(psd_grad, it, cb_power = -7)
    Ekin = E_ev(psd_grad.vpar_centers, c.me_ev)
    cm = Makie.to_colormap(:vik)
    fig, ax, hm = heatmap(Ekin, 
        psd_grad.altitude/1e3, 
        psd_grad.dfpardvpar[:, :, it], 
        #colorscale = log10,
        colormap = cm[60:end-60],
        colorrange = (-1.001e-7, 1.001e-7),
        axis = (xlabel = L"\mathrm{v_\parallel \, [eV]}",
            ylabel = "Height [km]",
            xscale = nsqrt,
            limits = ((-400, 400), nothing),    # x limit must be specified!
            title = format("t = {:.3f} s", psd_grad.time[it]),
            xticks = [-400, -200, -50, -10, 0, 10, 50, 200, 400],
            xminorticksvisible = true,
            xminorticks = vcat(vcat(0:2:9, 10:10:49, 50:50:400), vcat(0:2:9, 10:10:49, 50:50:400) .*-1),
            ),
        )
    cb = Colorbar(fig[1, 2], hm, 
        label = L"\mathrm{dF_{\parallel}/dv_{\parallel} \; [s^2 \, m^{-5}]}")

    rescale_colorbar_ticks!(fig, cb, cb_power; digits = 2)

#    Label(fig[1, 2, Top()], 
#        rich("    ×10", superscript(string(cb_power))), 
#        halign = :left,
#        #valign = :bottom,
#        padding = (0.0f0, 0.0f0, 7.0f0, 0.0f0),
#    )
    #text!(fig, x, y+10, text = "x10^-7")

    if :iv_growth in keys(psd_grad)
        scatter!(ax,
            psd_grad.E_kin_s, 
            psd_grad.altitude/1e3, 
            psd_grad.iv_growth[:, :, it], 
            marker = 'x',
            color = "black",)
    end
    display(fig)
    return fig, ax, hm
end


function plot_suprathermal(vel_mps)
    scatterlines(vel_mps, maxwellian_fa(c.me, 2500, vel_mps)*2e11,
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



function k_t_h(psd_le, klim)
    k_max = dropdims(nanmaximum(psd_le.k_growth, dims = 1), dims = 1)
    fig, ax, hm = heatmap(
        psd_le.time,
        psd_le.altitude/1e3, 
        k_max',
        colorrange = (klim/2, nanmaximum(k_max[:])),
        #lowclip = "white",
        axis = (xlabel = "Time [s]",
            ylabel = "Height [km]",
            #limits = ((9, 9.5), nothing),
        ),
    )
    cb = Colorbar(fig[1, 2], hm, label = "k [m⁻¹]", )
    add_line_colorbar(klim, fig, cb, (0, nanmaximum(psd_le.k_growth)))
    display(fig)
    return fig, ax, hm 
end



function g_t_h(psd_le)
    g_max = dropdims(nanmaximum(psd_le.gamma_growth, dims = 1), dims = 1)
    fig, ax, hm = heatmap(
        psd_le.time,
        psd_le.altitude/1e3, 
        g_max',
        colorrange = (0, nanmaximum(g_max[:])),
        #lowclip = "white",
        axis = (xlabel = "Time [s]",
            ylabel = "Height [km]",
            #limits = ((9, 9.5), nothing),
        ),
    )
    cb = Colorbar(fig[1, 2], hm, label = "Growth Rate [m⁻³]", )
    display(fig)
    return fig, ax, hm 
end