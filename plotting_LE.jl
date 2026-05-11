using WGLMakie
using MAT
using Serialization

klim = (224e6 * 2 *pi)/3e8 *2

dir = "AlfvenTrainPSD/PSD/psd/LE/143:161"
res_file = joinpath(dir, "kmat.bin")
psd_files = readdir(dir)
psd_files = filter(x-> contains(x, ".mat"), psd_files)
# Sort files numerically by extracting the number from filenames
psd_files = sort(psd_files, by = x -> parse(Int, match(r"psd-(\d+)", x).captures[1]))

#load sample
psd_data = matread(joinpath(dir, psd_files[1]))
psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data)
    )


psd_files = psd_files

nfiles = length(psd_files)
nt = length(psd.t_run)
nh = length(psd.h_atm)
kmat = zeros(nfiles*nt, nh, 10)*NaN;
tplot = zeros(nfiles*nt)*NaN;
h_atm = psd.h_atm


if !isfile(res_file)
    for (ifile, psd_file) in enumerate(psd_files)
        #psd_file = "psd-60_LE.mat"
        #ifile = 6
        println("Processing ", psd_file)
        #psd_file = joinpath(psd_dir, psd_f)

        # Load the .mat file
        psd_data = matread(joinpath(dir, psd_file))

        # Convert psd_data into a named tuple
        psd = NamedTuple{Tuple(Symbol(k) for k in keys(psd_data))}(values(psd_data)
            )
        #println(keys(psd))

        psd.t_run
        psd.h_atm
        psd.k_growth

        tplot[(ifile-1)*nt + 1:(ifile)*nt] = psd.t_run

        for it in axes(psd.t_run, 1)
            #it = 60
            kmat[(ifile-1)*nt + it, axes(psd.k_growth[it])...] = psd.k_growth[it]
        end
    end

    open(res_file, "w") do io
        serialize(io, [tplot, h_atm, kmat])
    end
end

io = open(res_file, "r")
tplot, h_atm, kmat = deserialize(io)
close(io)

kmat[kmat .< klim/2] .= NaN
tmat = ones(size(kmat)) .* tplot;
hmat = permutedims(permutedims(ones(size(kmat)), (2, 1, 3)) .* h_atm, (2, 1, 3));


##
fig = Figure()
sleep(1)
ax = Axis3(fig[1, 1], 
    xlabel = "Time [s]", 
    ylabel = "wavenumber [m-1]",
    zlabel = "Height [km]"
    )
sleep(1)

scatter!(ax, tmat[.!isnan.(kmat)], kmat[.!isnan.(kmat)], hmat[.!isnan.(kmat)]/1e3, color =kmat[.!isnan.(kmat)])
ax.azimuth = pi*1.1
ax.elevation = pi*0.05

##

io = open(res_file, "r")
tplot, h_atm, kmat = deserialize(io)
close(io)

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

