using MAT
using DataInterpolations
using Peaks
using GLMakie #, MakieExtra
#using Makie

include("utils.jl")
include("constants.jl")
c = physical_constants

showplots = true

psd_dir = "Afven Train PSB/psd/"
#for psd_file in readdir(psd_dir)
# loop over all files

# Load the .mat file
psd_file = "Afven Train PSB/psd/psd-30.mat"
psd_data = matread(psd_file)

# Convert psd_data into a named tuple
psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data))
println(keys(psd))

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

if showplots
    ih = 395
    it = 36
    lines(E, grad[:, ih, it])
end

if showplots
    fig, ax, hm = heatmap(E, psd.h_atm/1e3, grad[:, :, it], colorrange = (-1e20, 1e20),)
    xlims!(-200, 200)
    Colorbar(fig[1, 2], hm)
    display(fig)
end

imax, hmax = findmaxima(grad[:, ih, it])
lines(E, grad[:, ih, it])
scatter!(E[imax], hmax)


#for (imax, hmax) in findmaxima(grad[:, ih, it])
#   if hmax <= 0 continue end

im = 3
vmax = v[imax[im]]
dfdvmax = hmax[1]
E_ev(vmax, c.me_ev)


#load guisdap data
guisdap_dir = "analyzed_parameters/2018-12-07_folke_6.4@42mb"
files = filter(x -> x[end-3:end] == ".mat", readdir(guisdap_dir))
file = files[290]

ex_time, h_param, ne, Te = load_guisdap(joinpath(guisdap_dir, file))
println(ex_time)
#t_start = ex_time[1, :]
#t_end = ex_time[2, :]

if showplots
    lines(ne, h_ne/1e3)
end

h_median, ne_median, Te_median, h_mad, ne_mad, Te_mad = average_guisdap(joinpath.(guisdap_dir, files[290:292]))

gd = load_guisdap.(joinpath.(guisdap_dir, files[290:292]))
n_files = size(gd, 1)
t_start = gd[1].ex_time[1, :]
t_end = gd[end].ex_time[end, :]
interval = t_end - t_start

nh = maximum([size(m.h_param, 1) for m in gd])

h_mat = stack(rpad_array.([m.h_param for m in gd], nh, NaN))
h_mean  = [mean(row[.!isnan.(row)]) for row in eachrow(h_mat)]
h_median = [median(row[.!isnan.(row)]) for row in eachrow(h_mat)]

ne_mat = stack(rpad_array.([m.ne for m in gd], nh, NaN))
ne_mean  = [mean(row[.!isnan.(row)]) for row in eachrow(ne_mat)]
ne_median = [median(row[.!isnan.(row)]) for row in eachrow(ne_mat)]

Te_mat = stack(rpad_array.([m.Te for m in gd], nh, NaN))
Te_mean  = [mean(row[.!isnan.(row)]) for row in eachrow(Te_mat)]
Te_median = [median(row[.!isnan.(row)]) for row in eachrow(Te_mat)]

h_mad  = [mad(row[.!isnan.(row)]) for row in eachrow(h_mat)]
ne_mad = [mad(row[.!isnan.(row)]) for row in eachrow(ne_mat)]
Te_mad = [mad(row[.!isnan.(row)]) for row in eachrow(Te_mat)]

h_mad2  = [mad(row[.!isnan.(row)], center = mean(row[.!isnan.(row)])) for row in eachrow(h_mat)]
ne_mad2 = [mad(row[.!isnan.(row)], center = mean(row[.!isnan.(row)])) for row in eachrow(ne_mat)]
Te_mad2 = [mad(row[.!isnan.(row)], center = mean(row[.!isnan.(row)])) for row in eachrow(Te_mat)]

ne_mad[(ne_median .- ne_mad) .< 0] .= ne_median[(ne_median .- ne_mad) .< 0] ./ 10
#ne_mad2[(ne_mean .- ne_mad) .< 0] .= ne_mean[(ne_mean .- ne_mad) .< 0] ./ 10


if showplots
    fig, ax, lin = scatterlines(ne_mean, h_mean/1e3, 
        marker = 'x', 
        label = "mean",
        axis = (xscale = log10, limits = ((1e6, 2e12), nothing)
    ),
    )
    #errorbars!(ne_mean, h_mean/1e3, ne_mad2, direction = :x; color = :green) # same low and high error
    scatterlines!(ne_median, h_median/1e3, marker = 'x')
    errorbars!(ne_median, h_median/1e3, ne_mad, direction = :x; color = :red) # same low and high error
    for m in gd
        scatter!(m.ne, m.h_param/1e3, #color = "black", alpha = 0.2
        )
    end
    scatter!(gd[2].ne, gd[2].h_param/1e3, color = "black")
    #xlims!(1e6, 1e12)
    #xscale(log10)
    axislegend(ax)
    display(fig)
end

ih = 25

wp = plasma_freq(ne[ih], c.me)
vth = thermal_velocity.(Te, c.me)

vth = vth[ih]

k = -2000:2000
lines(k, k*vmax)
lines!(k, sqrt.(wp^2 .+ 3/2*vth.^2*k.^2))

#solve for k and omega
#k^2 vmax^2 .- 3/2*vth.^2*k.^2 3/2 = wp^2
#k^2 = wp^2 / (vmax^2 .- 3/2*vth.^2)
k_growth = wp / sqrt(vmax^2 - 3/2*vth^2)
w_growth = k_growth * vmax

#growth rate!
#double check!!
gamma = pi * wp^2 / k_growth^2 * dfdvmax

#check with iescat wavenumber
