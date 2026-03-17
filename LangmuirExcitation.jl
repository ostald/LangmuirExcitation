using MAT
using GLMakie
using DataInterpolations
using Peaks



# Load the .mat file
file_path = "Afven Train PSB/psd/psd-30.mat"
mat_data = matread(file_path)

# Inspect the contents
println(keys(mat_data))  # List all variables in the file

# Iterate through all keys in the dictionary and assign them as variables
for (key, value) in mat_data
    eval(Meta.parse("const $key = mat_data[\"$key\"]"))
end

c = 3e8
me = 511e3

using GLMakie, MakieExtra
using Makie
fig, ax, hm = heatmap(vpar_centres, h_atm/1e3, F[:, :, 50], 
    colorscale = log10,
    colorrange = (1e20, 1e28),
    axis = (xlabel = "Velocity [m/s]",
        ylabel = "Height [km]",
        #xscale = Makie.Symlog10(0.01),
        ),
    )
#xlims!(-100, 100)
Colorbar(fig[1, 2], hm)

it = 36
ih = 395

f = F[:, ih, it]
v = vpar_centres


function gradient_fpar_var(v_par, f_par)
    # for equidistand independent variable array, calculate gradient from neighbours
    grad = diff(f_par) ./ diff(abs.(v_par))
    return grad
end

me = 511e3 #eV
function E_ev(v_mps, m_ev)
    c = 299792458 #m/s
    E = m_ev / 2 * v_mps .* abs.(v_mps) / c^2
    return E
end

for (it, _) in enumerate(t_run)
    for (ih, _) in enumerate(h_atm)
        gradient_fpar_var(vpar_centres, F[:, ih, it])[2]
    end
end

dfpardvpar = stack([gradient_fpar_var(vpar_centres, F[:, ih, it]) for ih in eachindex(h_atm), it in eachindex(t_run)])
v_middle = vpar_centres[1:end-1] .+ diff(vpar_centres)

Ekin = E_ev(vpar_centres, me)
E    = E_ev(v_middle, me)

#lines!(E, grad[:, ih, it])

grad = copy(dfpardvpar[abs.(E) .< 200, :, :])
v = v_middle[abs.(E) .< 200]
E    = E_ev(v, me)

grad[abs.(E) .< 3, :, :] .= 0
fig, ax, hm = heatmap(E, h_atm/1e3, grad[:, :, it], colorrange = (-1e20, 1e20),)
xlims!(-200, 200)
Colorbar(fig[1, 2], hm)

#grad[abs.(E) .< 3] .= 0
#grad[abs.(E) .> 200].= 0




imax, hmax = findmaxima(grad[:, ih, it])
lines(E, grad[:, ih, it])
scatter!(E[imax], hmax)


#for (imax, hmax) in findmaxima(grad[:, ih, it])
#   if hmax <= 0 continue end

im = 3
vmax = v[imax[im]]
dfdvmax = hmax[1]
E_ev(vmax, me)

qe = 1.60217663e-19 #C
me_kg = 9.1093837e-31 #kg
function plasma_freq(ne, qe, m_kg)
    eps0 = 8.8541878188e-12
    return sqrt(ne*qe^2 / ( m_kg * eps0 ))
end

guisdap_dir = "analyzed_parameters/2018-12-07_folke_6.4@42mb"
files = filter(x -> x[end-3:end] == ".mat", readdir(guisdap_dir))
file = files[1]
guisdap_data = matread(joinpath(guisdap_dir, file))
ne = guisdap_data["r_pp"]
h_ne = guisdap_data["r_pprange"]
ex_time = guisdap_data["r_time"]
t_start = ex_time[1, :]
t_end = ex_time[2, :]

wp = plasma_freq(ne[1], qe, me_kg)
vth = 1200 #m/s

k = -200e3:200e3
lines(k, k*vmax)
lines!(k, sqrt.(wp^2 .+ 3*vth^2*k.^2))

