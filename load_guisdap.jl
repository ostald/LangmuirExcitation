using StatsBase

function load_guisdap(matfile)
    #loads one singe matfile
    guisdap_data = matread(matfile)
    ex_time = guisdap_data["r_time"]
    println(ex_time)
    h_param = dropdims(guisdap_data["r_h"], dims = 2)*1000
    ne = guisdap_data["r_param"][:, 1]
    Ti = guisdap_data["r_param"][:, 2]
    TeTi = guisdap_data["r_param"][:, 3]
    Te = TeTi .* Ti
    gd = (ex_time = ex_time,
        h_param = h_param,
        Te = Te, 
        ne = ne)
    return gd
end


function average_guisdap(files, mode = "median")
    # load several matfiles, averages over the contents
    gd = load_guisdap.(files)
    n_files = size(gd, 1)
    t_start = gd[1].ex_time[1, :]
    t_end = gd[end].ex_time[end, :]

    nh = maximum([size(m.h_param, 1) for m in gd])


    """
    gd_new = NamedTuple()

    for key in keys(gd[1])
    # make loop for all keys in gd, as all data is processed the same way
        if key == :Symbol(:ex_time) continue end
        stacked_k = stack(rpad_array.([m.key for m in gd], nh, NaN))
        median_k  = [median(filter_nan(row)) for row in eachrow(stacked)]
        #mean_k   = [mean(filter_nan(row))   for row in eachrow(stacked)]
        mad_k     = [mad(filter_nan(row))    for row in eachrow(stacked)]
        #mad_k     = [mad(filter_nan(row), center = mean(filter_nan(row))) for row in eachrow(stacked)]
        gd_new(key) = median_k

    end
    """


    h_mat = stack(rpad_array.([m.h_param for m in gd], nh, NaN))
    ne_mat = stack(rpad_array.([m.ne for m in gd], nh, NaN))
    Te_mat = stack(rpad_array.([m.Te for m in gd], nh, NaN))

    if mode == "median"
        h_median  = [median(row[.!isnan.(row)]) for row in eachrow(h_mat)]
        ne_median = [median(row[.!isnan.(row)]) for row in eachrow(ne_mat)]
        Te_median = [median(row[.!isnan.(row)]) for row in eachrow(Te_mat)]

        h_mad  = [mad(row[.!isnan.(row)]) for row in eachrow(h_mat)]
        ne_mad = [mad(row[.!isnan.(row)]) for row in eachrow(ne_mat)]
        Te_mad = [mad(row[.!isnan.(row)]) for row in eachrow(Te_mat)]
        return t_start, t_end, h_median, ne_median, Te_median, h_mad, ne_mad, Te_mad
    elseif mode == "mean"
        h_mean  = [mean(row[.!isnan.(row)]) for row in eachrow(h_mat)]
        ne_mean = [mean(row[.!isnan.(row)]) for row in eachrow(ne_mat)]
        Te_mean = [mean(row[.!isnan.(row)]) for row in eachrow(Te_mat)]

        h_mad  = [mad(row[.!isnan.(row)], center = mean(row[.!isnan.(row)])) for row in eachrow(h_mat)]
        ne_mad = [mad(row[.!isnan.(row)], center = mean(row[.!isnan.(row)])) for row in eachrow(ne_mat)]
        Te_mad = [mad(row[.!isnan.(row)], center = mean(row[.!isnan.(row)])) for row in eachrow(Te_mat)]

        return t_start, t_end, h_mean, ne_mean, Te_mean, h_mad, ne_mad, Te_mad
    end
end


function rpad_array(arr, length, value)
    # pads an array arr on the rigth side with value until it that the required length length
    padded_arr = fill(value, length)
    padded_arr[axes(arr, 1)] .= arr
    return padded_arr
end

function filter_nan(arr)
    arr_filtered = arr[.!isnan.(arr)]
    return arr_filtered
end