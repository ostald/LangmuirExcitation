using WGLMakie
using MAT
using Serialization
include("utils.jl")

klim = (224e6 * 2 *pi)/3e8 *2

dir = "results/143:161"
res_file = joinpath(dir, "psd_le.nc")

@time psd_data = NCDataset(res_file, "r") do ds
        Dict{String, Any}(
            "h_atm"        => Array(ds["altitude"]),
            "t_run"        => Array(ds["time"]),
            #"vpar_centers" => Array(ds["vpar_centers"]),
            #"vpar_edges"   => Array(ds["vpar_edges"]),
            #"F"            => Array(ds["F"]),
            "v"            => Array(ds["vpar_gradient_center_reduced"]),
            "k_growth"     => Array(ds["k_growth"]),
            "gamma_dfdvmax"=> Array(ds["gamma_dfdvmax"]),
            "iv_dfdvmax"   => Array(ds["iv_dfdvmax"]),
        )
    end
psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data));


"""
res_file2 = "/nfs/revontuli/data/oliver/LangmuirExcitation/AlfvenTrainPSD/PSD/psd/LE/143:161/kmat.bin"
io = open(res_file2, "r")
tplot, h_atm, kmat = deserialize(io)
close(io)
"""

kmat = permutedims(psd.k_growth, (3, 2, 1))
tplot = psd.t_run
h_atm = psd.h_atm

kmat[kmat .< klim/2] .= NaN
tmat = ones(size(kmat)) .* tplot;
hmat = permutedims(permutedims(ones(size(kmat)), (2, 1, 3)) .* h_atm, (2, 1, 3));


##
#using CairoMakie
#CairoMakie.activate!()
fig = Figure()
sleep(2)
ax = Axis3(fig[1, 1], 
    xlabel = "Time [s]", 
    ylabel = "wavenumber [m-1]",
    zlabel = "Height [km]"
    )
sleep(1)

scatter!(ax, tmat[.!isnan.(kmat)], kmat[.!isnan.(kmat)], hmat[.!isnan.(kmat)]/1e3, color =kmat[.!isnan.(kmat)])
ax.azimuth = pi*1.1
ax.elevation = pi*0.05
#display(fig)
##

"""
io = open(res_file, "r")
tplot, h_atm, kmat = deserialize(io)
close(io)
"""

kmat[isnan.(kmat)] .= 0

fig, ax, hm = heatmap(tplot, 
    h_atm/1e3, 
    dropdims(maximum(abs.(kmat), dims = 3), dims = 3),
    axis = (xlabel = "Time [s]",
    ylabel = "Height [km]"),
    )
Colorbar(fig[1, 2], hm, label = "k [m⁻¹]")

fig, ax, hm = heatmap(tplot, 
    h_atm/1e3, 
    dropdims(maximum(abs.(kmat), dims = 3), dims = 3),
    colorrange = (klim/2, maximum(kmat[:])),
    lowclip = "white",
    axis = (xlabel = "Time [s]",
    ylabel = "Height [km]"),
    )
Colorbar(fig[1, 2], hm, label = "k [m⁻¹]", )

##
    for i in 1:100
        ax.azimuth = pi*0.98 + i/1000
        sleep(0.01)
    end


    for _ in 1:1
        for i in 1:10
            ax.azimuth = pi*0.98 + i/100
            sleep(0.05)
        end
        for i in 1:10
            ax.azimuth = pi*0.98 + 0.1 - i/100
            sleep(0.05)
        end
    end


    for _ in 1:1
        for i in 1:100
            ax.elevation = pi*0.02 + i/100
            sleep(0.01)
        end
        for i in 1:100
            ax.elevation = pi*0.02 + 1 - i/100
            sleep(0.01)
        end
    end

