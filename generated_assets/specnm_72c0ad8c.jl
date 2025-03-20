### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a9019660-1f03-4132-bccf-09bdb1421ad9
begin
	using CondaPkg
	CondaPkg.add("pytest")
	CondaPkg.add("sympy")
	CondaPkg.add("flake8")
	CondaPkg.add("Cython")
	CondaPkg.add("mpi4py")
	CondaPkg.add("petsc")
	CondaPkg.add("petsc4py")
	CondaPkg.add("slepc4py")
	CondaPkg.add_pip("specnm", version="@https://gitlab.com/JohKem1/specnm/-/archive/main/specnm-main.zip")
end

# ╔═╡ 447b4e82-fe73-11ef-30b1-69824c8e3d24
using PythonCall, PlutoPlotly

# ╔═╡ 9f419394-528a-4bde-98a5-d62787c17fa8
using PlutoUI, FFTW, StatsBase

# ╔═╡ cd3c5ff7-7992-41f7-9a5e-29934a1355ca
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ bdae265d-3a96-4ecc-a6dd-c166357e801c
md"""
# 🌍 Exploring Free Oscillations
In this notebook, we will explore **free oscillations of spherically symmetric bodies**, which describe how a planet vibrates after a large disturbance (e.g., an earthquake).  
The numerical tool **`specnm`** allows us to compute **gravito-elastic normal modes** by solving the **radial ordinary differential equations (ODEs)** using a **spectral element discretization**.

`specnm`: 
A spectral element approach to computing normal modes, 
J Kemper, M van Driel, F Munch, A Khan, D Giardini,
Geophysical Journal International, Volume 229, Issue 2, May 2022, Pages 915–932.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 82cc9219-7633-41c6-91e1-17968904b2b6
md"## Real Data"

# ╔═╡ 5ec5e8d2-681b-4ce4-a581-478f00b91dc9
md"## Specnm"

# ╔═╡ 9b0a9a80-e65f-4385-a377-372e408b19ad
md"## Appendix"

# ╔═╡ 53495d28-cf5e-4108-9921-0ca016f1b24d
CondaPkg.add_pip("obspy")

# ╔═╡ 79406ee0-6026-4da8-a29d-245048c27e47
obspy = pyimport("obspy")

# ╔═╡ acd7956a-855c-43f0-a353-7d5533f6aaf1
begin
	# Import ObsPy modules correctly
	# obspy = pyimport("obspy.core")
	fdsn = pyimport("obspy.clients.fdsn")
end

# ╔═╡ 335162fc-52c3-4c3d-bb88-d2e2f4fde37f
client = fdsn.Client("IRIS")

# ╔═╡ 58baaf3c-a5b3-4995-9f4e-c54946f7e798
specnm = pyimport("specnm")

# ╔═╡ 5852bfac-a8ff-405e-8ed1-6438c6827091
md"### Select Earth Model"

# ╔═╡ d99b1935-4f90-4b82-809c-b6a801c37e0d
model_fname = "../specnm_models/prem_ani"

# ╔═╡ 1b22f312-9197-4044-ba9b-1f12789b88ff
ray = specnm.rayleigh(model_fname, fmax=0.1)

# ╔═╡ d9c71796-e79f-404d-89ef-54adfbd3335c
ray_out = ray.rayleigh_problem(attenuation_mode="elastic", fmax=0.005)

# ╔═╡ beddae25-17bc-48e9-8eef-41f21a08fb10
begin
	ray_angular_orders = pyconvert(Array, ray_out["angular orders"])
	ray_frequencies = pyconvert(Array, ray_out["frequencies"] * 1000.0) # in mHz
end;

# ╔═╡ 69699338-d461-4e05-bbb4-a874ad4ad970
ray_overtones = let 
	angular_orders = pyconvert(Vector{Int}, ray_out["angular orders"])

# Get unique values and their counts
l_countmap = countmap(angular_orders)  # Dictionary of {l => count}
l_uniq = collect(keys(l_countmap))  # Unique angular orders
l_count = collect(values(l_countmap))  # Their counts

# Compute overtones
# Compute overtones
overtones = vcat([lu != 1 ? collect(0:lc-1) : collect(2:lc+1) for (lu, lc) in zip(l_uniq, l_count)])

# Convert to a single vector
ray_overtones = collect(reduce(vcat, overtones))
end;

# ╔═╡ 55b7b9a8-caea-4ec6-b47d-937231eaaee8
lov = specnm.love(model_fname, fmax=0.1)

# ╔═╡ 3751319d-d1b8-484a-8a70-ad46a6a65634
lov_out = lov.love_problem(attenuation_mode="elastic", fmax=0.005)

# ╔═╡ 21d48160-8f70-4b15-9534-a7b0f5131335
@bind ray_lov_selected Select([[ray, ray_out] => "Spheroidal", [lov, lov_out] => "Toroidal"])

# ╔═╡ c692b5d4-e136-47f4-9cc6-432d4e5ef3fe
begin
	lov_angular_orders = pyconvert(Array, lov_out["angular orders"])
	lov_frequencies = pyconvert(Array, lov_out["frequencies"] * 1000.0) # in mHz
end;

# ╔═╡ c1238b16-770c-4a31-9e04-c8c9fb45e1f5
lov_overtones = let 
	angular_orders = pyconvert(Vector{Int}, lov_out["angular orders"])

# Get unique values and their counts
l_countmap = countmap(angular_orders)  # Dictionary of {l => count}
l_uniq = collect(keys(l_countmap))  # Unique angular orders
l_count = collect(values(l_countmap))  # Their counts

# Compute overtones
overtones = vcat([lu != 1 ? collect(0:lc-1) : collect(1:lc) for (lu, lc) in zip(l_uniq, l_count)])

# Convert to a single vector
lov_overtones = collect(reduce(vcat, overtones))
end;

# ╔═╡ 526c0fd9-4935-409e-82e0-059f5c084b57
function read_model(filename::String)
    # Read all lines from the file
    lines = readlines(filename)

    # Extract numeric data, skipping header and comment lines
    data = []
    for line in lines
        # Skip lines that contain non-numeric text
        words = split(line)
        if length(words) > 0 && all(w -> tryparse(Float64, w) !== nothing, words)
            push!(data, parse.(Float64, words))
        end
    end

    # Convert list to matrix
    data_matrix = hcat(data...)

    # Extract relevant columns
    radius = data_matrix[1, :]  # Column 1: Radius (m)
    rho = data_matrix[2, :]     # Column 2: Density (kg/m³)
    vp = data_matrix[3, :]      # Column 3: P-wave velocity (m/s)
    vs = data_matrix[4, :]      # Column 4: S-wave velocity (m/s)

    # Convert radius to depth (Earth's radius is ~6371 km)
    depth = (6371000 .- radius) ./ 1000  # Convert meters to kilometers

    return depth, rho, vp, vs
end

# ╔═╡ b5bfa4ac-f20a-4d1a-b870-fe4a3fa52da3
md"### Get Earthquakes"

# ╔═╡ a5491b1e-ce27-4fb5-82de-853899f75006
begin
	starttime = obspy.UTCDateTime("2000-01-01")
	endtime = obspy.UTCDateTime("2024-12-31")
	min_magnitude = 7.5
end

# ╔═╡ 0ecb0dea-e87e-4602-956c-9fc4a9375cb8
# Fetch events from IRIS
catalog = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=min_magnitude)

# ╔═╡ 7db650a6-5155-494f-9fc3-c5d664f98b32
# Store earthquake details in a structured format
earthquake_list = [
    (
        "Mw $(event.magnitudes[0].mag) at $(event.origins[0].time)",
        event.origins[0].time, event.origins[0].latitude, 
        event.origins[0].longitude, event.origins[0].depth / 1000, 
        event.magnitudes[0].mag
    ) for event in catalog
]

# ╔═╡ 5a632ee0-e405-432f-a35e-e37defe34552
earthquakes = [
    (
        event.origins[0].time, event.origins[0].latitude,
        event.origins[0].longitude, event.magnitudes[0].mag
    ) for event in catalog
]

# ╔═╡ 849d1a59-2c2f-4e44-b45b-266382255b1a
@bind clicked_eq let
	p = PlutoPlot(Plot(scattergeo(
    lon=pyconvert.(Float64, [e[3] for e in earthquakes]), lat=pyconvert.(Float64, [e[2] for e in earthquakes]),
    text=pyconvert.(String, [string("M", e[4], " at ", e[1]) for e in earthquakes]),
    mode="markers", marker=attr(size=6, color="red")
), Layout(title="Earthquakes (click to select)")))
	add_plotly_listener!(p,"plotly_click", "
	(e) => {

	console.log(e)
    let dt = e.points[0]
	PLOT.value = [dt.lat, dt.lon]
	PLOT.dispatchEvent(new CustomEvent('input'))
}
	")
	p
end

# ╔═╡ 902f56d9-041d-4b4a-af72-1adcf20a4db1
# Extract details of the selected earthquake (filtered using lat/lon)
selected_event_index = findall(e -> pyconvert(Float32, e[2]) ≈ clicked_eq[1] && pyconvert(Float32, e[3]) ≈ clicked_eq[2], earthquakes)[1]

# ╔═╡ ea6fe2a5-42f3-4c6f-b9e9-f4e4172463da
begin
	# Find details of the selected earthquake
	selected_details = earthquake_list[selected_event_index]
	
	eq_time = selected_details[2]
	eq_lat = selected_details[3]
	eq_lon = selected_details[4]
end

# ╔═╡ 2afaeeec-6fe2-48db-81de-1fa4ccb0fc1b
md"### Get Stations"

# ╔═╡ 73042ddf-790f-4a11-98b2-f9b6b9d29fe0
begin
	# Fetch GSN stations that recorded this earthquake
	network = "IU"  # IU = Global Seismographic Network (GSN)
	station_list = client.get_stations(network=network, latitude=eq_lat, longitude=eq_lon, maxradius=180)
	
end

# ╔═╡ 521d6111-e4fb-4532-9acb-561a35aa5607
# Extract station metadata
stations = [
    (
        s.code, s.latitude, s.longitude, s.elevation
    ) for net in station_list for s in net.stations
]

# ╔═╡ 74e0c240-dda8-44bd-bf0b-fb862a6ea52d
# Create a dropdown for selecting a station (if stations exist)
station_names = length(stations) > 0 ? [s[1] for s in stations] : ["No stations found"]

# ╔═╡ 2d7e4963-1cee-413c-8541-e8e77962c1fe
@bind clicked_station let
	p = PlutoPlot(Plot(scattergeo(
    lon=pyconvert.(Float32, [e[3] for e in stations]), lat=pyconvert.(Float32, [e[2] for e in stations]),
    text=pyconvert.(String, [string(e) for e in station_names]),
    mode="markers", marker=attr(size=6, symbol="triangle-down", color="blue")
), Layout(title="GSN Stations (click to select)")))
	add_plotly_listener!(p,"plotly_click", "
	(e) => {

	console.log(e)
    let dt = e.points[0]
	PLOT.value = [dt.lat, dt.lon]
	PLOT.dispatchEvent(new CustomEvent('input'))
}
	")
	p
end

# ╔═╡ 5360665a-81ba-4cdb-8fdc-55829c6f4255
# Extract details of the selected earthquake (filtered using lat/lon)
selected_station_index = findall(e -> pyconvert(Float32, e[2]) ≈ clicked_station[1] && pyconvert(Float32, e[3]) ≈ clicked_station[2], stations)[1]

# ╔═╡ aa28ca2f-662b-4789-a6fe-579f85d27ada
selected_station = stations[selected_station_index]

# ╔═╡ ba1041f2-8173-4436-9413-ba8c8c85af20
 inventory = client.get_stations(
        network="IU",
        station=selected_station[1],
        level="channel"
    )

# ╔═╡ 55968681-0b8c-4bb0-8ee9-1a9d81ee34e6
   available_channels = [c.code for net in inventory for sta in net.stations for c in sta.channels]


# ╔═╡ 580d562e-4c3f-4adb-94cf-65048ee7ff35
begin
	freqmin = 0.2 * inv(1000)  # Lower cutoff frequency (Hz)
	freqmax = 5.0 * inv(1000)    # Upper cutoff frequency (Hz)
end

# ╔═╡ 7493f1d5-2017-416b-bc0c-183767bad68a
begin
	starttime_data = eq_time + 0.5 * 3600  # Roughly start half a day after earthquake time
	endtime_data = starttime_data + 48 * 3600  # 6 hours of data
	st = client.get_waveforms(
	        "IU",  # Network (IU, IC, II)
	        selected_station[1],  # Station code
	        "*",  # Any location
	        "LHZ",  # Vertical component, 
		attach_response=true,
	        starttime_data, endtime_data
	    )
	trace = st[0]
	trace.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=true)
	data = pyconvert(Array, trace.normalize().detrend().data)
	sampling_rate = pyconvert(Float32, trace.stats.sampling_rate)
	
end

# ╔═╡ b57d88b9-bdd5-447b-ba24-30a7b1ea1e0b
plot(data)

# ╔═╡ c168a471-f21f-4175-b7a0-5ae44d5d30ee
begin
	# Compute FFT
	n = length(data)
	freqs = rfftfreq(n, 1/sampling_rate)
	
	# Convert frequency to mHz (as normal modes are in mHz)
	freqs_mHz = freqs .* 1000

	
	ifmax = findfirst(x->x>freqmax*1000.0, freqs_mHz)
	ifmin = findfirst(x->x>freqmin*1000.0, freqs_mHz)
	freqs_mHz = freqs_mHz[ifmin:ifmax]
	freqs = freqs[1:ifmax]
	
	spectrum = abs.(rfft(data))[ifmin:ifmax]
	
end

# ╔═╡ 1a67d073-bf07-4b01-beea-6c11013965be
let
	# Create spectrum plot
	spectrum_plot = [scatter(
	    x=freqs_mHz[freqs_mHz .> 0], y=spectrum[freqs_mHz .> 0],
	    mode="lines", line=attr(color="blue"),
		name="",
	)]

	# Add vertical lines at each mode frequency
for (f, l, n) in zip(ray_frequencies, ray_angular_orders, ray_overtones)
    push!(spectrum_plot, scatter(
        x=[f, f],  # Vertical line at frequency f
        y=[0, maximum(spectrum)],  # Span full range of l
        mode="lines",
		name="$((l,n))",
        line=attr(color="gray", width=0.5),
        showlegend=false
    ))
end

	
	
	# Interactive Plotly figure
	plot(spectrum_plot, Layout(
	    title="Normal Mode Spectrum",
	    xaxis_title="Frequency (mHz)", yaxis_title="Amplitude",
	))
end

# ╔═╡ 89c81d3c-67c1-4cc7-b892-d64208d845f2
md"### Plots"

# ╔═╡ ef32df8e-5053-4c00-8ab8-2087f7771051
pyconvert(Array, ray_out["frequencies"] * 1000.0)

# ╔═╡ 39710dda-1054-491e-b9d5-7d8c4cb8e3a6
function spectrum_plot(cls, cls_out;
                      plot_title=true, plot_efs=true, exclude_slichter=false,
                      plot_wavenumber=false, plot_phasevelocity=false,
                      init_mode=nothing)
        oplotcolors = ["#e31a1c", "#33a02c", "#ff7f00", "#1f78b4"]

        # addt = 1 + Int(exclude_slichter) * Int(cls.sph_type > 1)

        ls = cls_out["angular orders"]
        if plot_wavenumber
            xlabel = "wave number k in km"
            xs = pyconvert(Array, cls_out["wave numbers"] .* cls_out["radius"][end] / 1000.0)
        else
            xlabel = "angular degree l"
            xs = pyconvert(Array, cls_out["angular orders"])
        end

        if plot_phasevelocity
            ylabel = "phase velocity in km/s"
            fcp = pyconvert(Array, cls_out["phase velocities"] / 1000.0)
        else
            ylabel = "frequency in mHz"
            fcp = pyconvert(Array, cls_out["frequencies"] * 1000.0)
        end

        return Plot(scatter(x=xs, y=fcp, mode="markers", marker_size=4, marker_color="black"), Layout(xaxis=attr(title=xlabel), yaxis=attr(title=ylabel), title="Type: $(cls.mode), Modelname: $(cls.modelname)"))

    end

# ╔═╡ 58041dd4-fd1c-4223-a00d-d4e98a7ef412
@bind clicked_mode let
	p = PlutoPlot(spectrum_plot(ray_lov_selected...))
	add_plotly_listener!(p,"plotly_click", "
	(e) => {

	console.log(e)
    let dt = e.points[0]
	PLOT.value = [dt.x, dt.y]
	PLOT.dispatchEvent(new CustomEvent('input'))
}
	")
	p
end

# ╔═╡ ecee2903-a103-498a-9aaf-b5d920f531f4
"""
Plot eigenfunctions
"""
function eigenfunction_plot(cls, cls_out, clicked_mode)
	eigenfunctions = pyconvert(Matrix, cls_out["eigenfunctions"])
	xs = pyconvert(Array, cls_out["angular orders"])
	fcp = pyconvert(Array, cls_out["frequencies"] * 1000.0) # in mHz
	 if pyconvert(String, cls.mode) == "spheroidal"
				
                if pyconvert(Int, cls.sph_type) == 3
                    plabels = ["U", "V", "P"]
                    efs = [eigenfunctions[:, 1:3:end],
                           eigenfunctions[:, 2:3:end],
                           eigenfunctions[:, 3:3:end]]
                else
                    plabels = ["U", "V"]
                    efs = [eigenfunctions[:, 1:2:end],
                           eigenfunctions[:, 2:2:end]]
                end
            elseif pyconvert(String, cls.mode) == "radial"
                plabels = ["U"]
                efs = [eigenfunctions]
            elseif pyconvert(String, cls.mode) == "toroidal"
                plabels = ["W"]
                efs = [eigenfunctions]
            else
                error("Not implemented for type $(cls.mode)")
            end
	rad_grid = pyconvert(Array, cls.r / 1000.00)
	# depth_grid = pyconvert(Float64, cls.radius / 1000.00) .- 
	depth_grid = -1.0 .* (rad_grid  .- pyconvert(Float64, cls.radius / 1000.00))
	dataid = findfirst(x -> isapprox([x...], clicked_mode, atol=1e-3), collect(zip(pyconvert(Array, cls_out["angular orders"]), pyconvert(Array, cls_out["frequencies"] .* 1000.0))))
	traces = map(efs, plabels) do e, label
		scatter(x=e[dataid, :], y=depth_grid, name=label)
	end
	plot(traces, Layout(title="Eigen function; f (mHz): $(round(fcp[dataid], digits=3))<br> 
	period (min): $(round(inv(fcp[dataid])*1000.0/60.0, digits=3))<br>l=$(xs[dataid])", width=300, height=550, yaxis=attr(range=[7000, -100], title="depth in km")))
	
	# plot(scatter(x=f[dataid, :], y=cls.r / 1000.0, mode="lines", line_color=oplotcolors[e], name=plabels[e]), )
end

# ╔═╡ 71af98fe-12cb-42bf-b447-5afde0482e47
"""
Plot vp, vs and rho as a function of depth
"""
function plot_Earth_model(model_fname)
	depth, rho, vp, vs = read_model(model_fname)
    # Create traces for each variable
    trace_rho = scatter(x=rho, y=depth, mode="lines", name="Density (kg/m³)", line=attr(color="red"))
    trace_vp = scatter(x=vp, y=depth, mode="lines", name="P-Wave Velocity (m/s)", line=attr(color="blue"))
    trace_vs = scatter(x=vs, y=depth, mode="lines", name="S-Wave Velocity (m/s)", line=attr(color="green"))

    # Create figure with reversed depth axis
    fig = plot([trace_rho, trace_vp, trace_vs],
        Layout(title="Earth Model",
               xaxis_title="Property Value",
               yaxis_title="depth in km",
               yaxis=attr(range=[7000, -100]),  # Depth increases downward
               legend=attr(x=0.9, y=1), width=400, height=550,
        )
    )
    return fig
end

# ╔═╡ a9bbb105-3138-481d-9a73-0ec2b8303c36
PlutoUI.ExperimentalLayout.hbox([plot_Earth_model(model_fname), eigenfunction_plot(ray_lov_selected..., clicked_mode)])

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CondaPkg = "~0.2.26"
FFTW = "~1.8.1"
PlutoPlotly = "~0.6.2"
PlutoUI = "~0.7.61"
PythonCall = "~0.9.24"
StatsBase = "~0.34.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "18456e1aa235134cfb42b1b67a74a9aa4043e0c3"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CondaPkg]]
deps = ["JSON3", "Markdown", "MicroMamba", "Pidfile", "Pkg", "Preferences", "Scratch", "TOML", "pixi_jll"]
git-tree-sha1 = "44d759495ed1711e3c0ca469f8d609429318b332"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.26"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "7de7c78d681078f027389e067864a8d53bd7c3c9"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "1d322381ef7b087548321d3f878cb4c9bd8f8f9b"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.1"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MicroMamba]]
deps = ["Pkg", "Scratch", "micromamba_jll"]
git-tree-sha1 = "011cab361eae7bcd7d278f0a7a00ff9c69000c51"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.14"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "90af5c9238c1b3b25421f1fdfffd1e8fca7a7133"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.20"

    [deps.PlotlyBase.extensions]
    DataFramesExt = "DataFrames"
    DistributionsExt = "Distributions"
    IJuliaExt = "IJulia"
    JSON3Ext = "JSON3"

    [deps.PlotlyBase.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "ColorSchemes", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "PrecompileTools", "Reexport", "ScopedValues", "Scratch", "TOML"]
git-tree-sha1 = "9ebe25fc4703d4112cc418834d5e4c9a4b29087d"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.2"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "Requires", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "feab249add2d40873acbd6b286b450bd30b083dd"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.24"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "1147f140b4c8ddab224c94efa9569fc23d63ab44"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnsafePointers]]
git-tree-sha1 = "c81331b3b2e60a982be57c046ec91f599ede674a"
uuid = "e17b2a0c-0bdf-430a-bd0c-3a23cae4ff39"
version = "1.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.micromamba_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "b4a5a3943078f9fd11ae0b5ab1bdbf7718617945"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "1.5.8+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.pixi_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "f349584316617063160a947a82638f7611a8ef0f"
uuid = "4d7b5844-a134-5dcd-ac86-c8f19cd51bed"
version = "0.41.3+0"
"""

# ╔═╡ Cell order:
# ╠═cd3c5ff7-7992-41f7-9a5e-29934a1355ca
# ╟─bdae265d-3a96-4ecc-a6dd-c166357e801c
# ╟─21d48160-8f70-4b15-9534-a7b0f5131335
# ╟─58041dd4-fd1c-4223-a00d-d4e98a7ef412
# ╟─a9bbb105-3138-481d-9a73-0ec2b8303c36
# ╟─82cc9219-7633-41c6-91e1-17968904b2b6
# ╟─849d1a59-2c2f-4e44-b45b-266382255b1a
# ╟─2d7e4963-1cee-413c-8541-e8e77962c1fe
# ╠═b57d88b9-bdd5-447b-ba24-30a7b1ea1e0b
# ╟─1a67d073-bf07-4b01-beea-6c11013965be
# ╟─5ec5e8d2-681b-4ce4-a581-478f00b91dc9
# ╠═1b22f312-9197-4044-ba9b-1f12789b88ff
# ╠═55b7b9a8-caea-4ec6-b47d-937231eaaee8
# ╠═3751319d-d1b8-484a-8a70-ad46a6a65634
# ╠═d9c71796-e79f-404d-89ef-54adfbd3335c
# ╠═beddae25-17bc-48e9-8eef-41f21a08fb10
# ╠═c692b5d4-e136-47f4-9cc6-432d4e5ef3fe
# ╠═69699338-d461-4e05-bbb4-a874ad4ad970
# ╠═c1238b16-770c-4a31-9e04-c8c9fb45e1f5
# ╟─9b0a9a80-e65f-4385-a377-372e408b19ad
# ╠═a9019660-1f03-4132-bccf-09bdb1421ad9
# ╠═53495d28-cf5e-4108-9921-0ca016f1b24d
# ╠═79406ee0-6026-4da8-a29d-245048c27e47
# ╠═acd7956a-855c-43f0-a353-7d5533f6aaf1
# ╠═335162fc-52c3-4c3d-bb88-d2e2f4fde37f
# ╠═447b4e82-fe73-11ef-30b1-69824c8e3d24
# ╠═9f419394-528a-4bde-98a5-d62787c17fa8
# ╠═58baaf3c-a5b3-4995-9f4e-c54946f7e798
# ╟─5852bfac-a8ff-405e-8ed1-6438c6827091
# ╠═d99b1935-4f90-4b82-809c-b6a801c37e0d
# ╠═526c0fd9-4935-409e-82e0-059f5c084b57
# ╟─b5bfa4ac-f20a-4d1a-b870-fe4a3fa52da3
# ╠═a5491b1e-ce27-4fb5-82de-853899f75006
# ╠═0ecb0dea-e87e-4602-956c-9fc4a9375cb8
# ╠═7db650a6-5155-494f-9fc3-c5d664f98b32
# ╠═5a632ee0-e405-432f-a35e-e37defe34552
# ╠═902f56d9-041d-4b4a-af72-1adcf20a4db1
# ╠═ea6fe2a5-42f3-4c6f-b9e9-f4e4172463da
# ╟─2afaeeec-6fe2-48db-81de-1fa4ccb0fc1b
# ╠═73042ddf-790f-4a11-98b2-f9b6b9d29fe0
# ╠═521d6111-e4fb-4532-9acb-561a35aa5607
# ╠═74e0c240-dda8-44bd-bf0b-fb862a6ea52d
# ╠═5360665a-81ba-4cdb-8fdc-55829c6f4255
# ╠═aa28ca2f-662b-4789-a6fe-579f85d27ada
# ╠═ba1041f2-8173-4436-9413-ba8c8c85af20
# ╠═55968681-0b8c-4bb0-8ee9-1a9d81ee34e6
# ╠═7493f1d5-2017-416b-bc0c-183767bad68a
# ╠═580d562e-4c3f-4adb-94cf-65048ee7ff35
# ╠═c168a471-f21f-4175-b7a0-5ae44d5d30ee
# ╟─89c81d3c-67c1-4cc7-b892-d64208d845f2
# ╠═ef32df8e-5053-4c00-8ab8-2087f7771051
# ╠═39710dda-1054-491e-b9d5-7d8c4cb8e3a6
# ╠═ecee2903-a103-498a-9aaf-b5d920f531f4
# ╠═71af98fe-12cb-42bf-b447-5afde0482e47
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
