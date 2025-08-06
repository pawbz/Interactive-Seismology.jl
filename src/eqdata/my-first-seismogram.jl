### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 2cc62eeb-edc2-4fcb-876b-010d8a180162
begin
    using CondaPkg
    CondaPkg.add_pip("obspy")
end

# ╔═╡ 53f5276e-7954-43d5-bcfb-6898bc6e85ab
using PythonCall, PlutoPlotly

# ╔═╡ efb48ae5-c01c-4e58-ba08-22f0c885ae4a
using PlutoUI, FFTW, StatsBase

# ╔═╡ 277155cf-77d8-4cd6-a122-346be2ff1516
using Colors

# ╔═╡ ca9f95ec-71d2-11f0-2e94-7da041469cdf
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ 5a72389a-1277-4e18-a532-ce81a4b18c74
md"""
# My First Seismogram
This notebook lets users interactively select a global earthquake and a seismic station, then downloads and visualizes the corresponding seismogram data. It also marks theoretical seismic phase arrivals on the trace using the TauP toolkit for educational exploration.



##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 18eb286d-263e-42ee-94fc-26f3223ce108
md"""
### Parameters

#### Select Receivers (distance to the earthquake in degrees)
- Min Distance: $(@bind minradius Slider(0:0.5:180, show_value=true, default=10))
- Max Distance: $(@bind maxradius Slider(0:0.5:180, show_value=true, default=90))

#### Frequency Filter (Hz)
- Minimum Frequency: $(@bind freqmin Slider(0.01:0.01:0.1, show_value=true, default=0.05))
- Maximum Frequency: $(@bind freqmax Slider(0.2:0.01:1.0, show_value=true, default=0.5))

#### Total Duration (seconds)
- $(@bind total_duration Slider(1000:10:5000, show_value=true, default=4000))
"""

# ╔═╡ d7023542-1db3-42fa-a4c6-f716f8a6fcfb
md"## Random Earthquake and Station Selection"

# ╔═╡ eddf7e8c-21e3-4d0d-9cfb-dacb24047b44
md"## Selected Indices"

# ╔═╡ 0a78aecc-932f-4aa9-9172-9d79736b5301
md"## Station List"

# ╔═╡ 9ef37a61-295d-45ca-b686-9592cb417e5a
md"## Data Download"

# ╔═╡ 98b22766-a9b2-419e-8ec2-6287f65e2722
md"## TauP"

# ╔═╡ d77d897a-ca34-4a00-a5c2-fdc988c6b75c
md"## Earthquake List"

# ╔═╡ d2aba787-dae6-4da5-8f16-0cadfcc38371
md"## Appendix"

# ╔═╡ 4afe5b07-8c6c-463f-aa77-c866d378ea61
obspy = pyimport("obspy")

# ╔═╡ 68087281-2e24-40a5-872e-6d549bdd2eaa
begin
    starttime = obspy.UTCDateTime("2000-01-01")
    endtime = obspy.UTCDateTime("2024-12-31")
    min_magnitude = 7.5
end

# ╔═╡ 1d74bb4f-eb8b-4e6a-bb12-9f8287bfa05b
taup = pyimport("obspy.taup")

# ╔═╡ 4b0309a1-e369-4954-8b40-e58fbe41754f
model = taup.TauPyModel(model="iasp91")

# ╔═╡ 20bf25f7-a56d-4e01-a4cd-ec7506028755
UTCDateTime = obspy.UTCDateTime

# ╔═╡ 4ee50f18-9d36-43c3-aced-6b33a661a959
fdsn = pyimport("obspy.clients.fdsn")

# ╔═╡ 08001b34-7cd7-4b72-9fc8-b4a8c869c460
iris = pyimport("obspy.clients.iris")

# ╔═╡ b4199996-5f1b-41b7-8732-6971a8d1824c
client = fdsn.Client("IRIS")

# ╔═╡ 00670b47-a0a1-43aa-b304-3a63bd1615c2
# Fetch events from IRIS
catalog = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=min_magnitude)

# ╔═╡ ce54afb1-9a66-4d88-aced-e8b94b162776
# Store earthquake details in a structured format
earthquake_list = [
    (
        "Mw $(event.magnitudes[0].mag) at $(event.origins[0].time)",
        event.origins[0].time, event.origins[0].latitude,
        event.origins[0].longitude, event.origins[0].depth / 1000,
        event.magnitudes[0].mag
    ) for event in catalog
]

# ╔═╡ dd849802-4614-4f30-bbce-d87660df3840
earthquakes = [
    (
        event.origins[0].time, event.origins[0].latitude,
        event.origins[0].longitude, event.magnitudes[0].mag
    ) for event in catalog
]

# ╔═╡ 1939e9df-8c71-45b9-981f-e410d0a9965b
@bind clicked_eq let
    p = PlutoPlot(Plot(scattergeo(
            lon=pyconvert.(Float64, [e[3] for e in earthquakes]), lat=pyconvert.(Float64, [e[2] for e in earthquakes]),
            text=pyconvert.(String, [string("M", e[4], " at ", e[1]) for e in earthquakes]),
            mode="markers", marker=attr(size=6, color="red")
        ), Layout(title="Earthquakes (click on the map to select)")))
    add_plotly_listener!(
        p,
        "plotly_click",
        "
	(e) => {

	console.log(e)
    let dt = e.points[0]
	PLOT.value = [dt.lat, dt.lon]
	PLOT.dispatchEvent(new CustomEvent('input'))
}
	"
    )
    p
end

# ╔═╡ 3e20d6bd-a212-44a2-b7ff-0abba9a36a90
default_eq = let
    eq = rand(earthquakes)
    pyconvert.(Float64, [eq[2], eq[3]])
end

# ╔═╡ 834f3dae-81dd-4ff4-a891-2e62a1e10ce8
# Extract details of the selected earthquake (filtered using lat/lon)
begin
    selected_eq = (clicked_eq === nothing) ? default_eq : clicked_eq
    selected_earthquake_index = findall(e -> pyconvert(Float32, e[2]) ≈ selected_eq[1] && pyconvert(Float32, e[3]) ≈ selected_eq[2], earthquakes)[1]
end

# ╔═╡ c0ba9e33-0365-44ae-8ed9-f3d67d0e20c4
begin
    # Find details of the selected earthquake
    selected_earthquake_details = earthquake_list[selected_earthquake_index]

    eq_time = selected_earthquake_details[2]
    eq_lat = selected_earthquake_details[3]
    eq_dep = selected_earthquake_details[5]
    eq_lon = selected_earthquake_details[4]
end

# ╔═╡ ad11d31e-2aed-4fc3-9e77-17c40c907e83
station_list = let
    # Fetch GSN stations that recorded this earthquake
    network = "IU"  # IU = Global Seismographic Network (GSN)
    station_list = client.get_stations(network=network, latitude=eq_lat, longitude=eq_lon, minradius=minradius, maxradius=maxradius)
end

# ╔═╡ 5c05410d-d364-4f32-a584-5601874887d0
# Extract station metadata
stations = [
    (
        s.code, s.latitude, s.longitude, s.elevation
    ) for net in station_list for s in net.stations
]

# ╔═╡ 97de4846-f1ba-4269-bb25-68068dca274c
default_station = let
    st = rand(stations)
    pyconvert.(Float64, [st[2], st[3]])
end

# ╔═╡ caf0a5ef-a65b-4e11-bedb-c6a5463f6387
# Create a dropdown for selecting a station (if stations exist)
station_names = length(stations) > 0 ? [s[1] for s in stations] : ["No stations found"]

# ╔═╡ ad4d1627-f3ec-4a34-b109-bf07b1ce5bc0
@bind clicked_station let
    p = PlutoPlot(Plot(scattergeo(
            lon=pyconvert.(Float32, [e[3] for e in stations]), lat=pyconvert.(Float32, [e[2] for e in stations]),
            text=pyconvert.(String, [string(e) for e in station_names]),
            mode="markers", marker=attr(size=6, symbol="triangle-down", color="blue")
        ), Layout(title="GSN Stations (click on the map to select)")))
    add_plotly_listener!(
        p,
        "plotly_click",
        "
	(e) => {

	console.log(e)
    let dt = e.points[0]
	PLOT.value = [dt.lat, dt.lon]
	PLOT.dispatchEvent(new CustomEvent('input'))
}
	"
    )
    p
end

# ╔═╡ 1617897a-dd18-4c80-b12a-c2d063965e09
begin
    selected_station = (clicked_station === nothing) ? default_station : clicked_station
    # Extract details of the selected receiver
    selected_station_index = findall(e -> pyconvert(Float32, e[2]) ≈ selected_station[1] && pyconvert(Float32, e[3]) ≈ selected_station[2], stations)[1]
end

# ╔═╡ 44944cf9-b2c7-4264-b766-6163e3946d86
selected_station_details = stations[selected_station_index]

# ╔═╡ 0ec9801f-515e-4e42-9bc0-8c96c1e8937f
md"""
### Current Selections

- **Selected Earthquake**: $(selected_earthquake_details)

- **Selected Receiver**: $(selected_station_details)
"""

# ╔═╡ 4be1db35-014e-4a56-8006-a7fe2ed358f1
println(station_list)

# ╔═╡ fd936886-101c-4244-898d-dfcf1333555b
begin
    starttime_data = eq_time  # Roughly start half a day after earthquake time
    endtime_data = starttime_data + total_duration
    traces = client.get_waveforms(
        "IU",  # Network (IU, IC, II)
        selected_station_details[1],  # Station code
        "*",  # Any location
        "BH?",  # Vertical component, 
        attach_response=true,
        starttime_data, endtime_data
    )
    traces = traces.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=true)
    traces = traces.normalize().detrend()
    # sampling_rate = pyconvert(Float32, first(traces.stats.sampling_rate)
end

# ╔═╡ 04687d41-3a81-4413-acb8-8c2293e39139
channels = pyconvert.(String, [s.stats.channel for s in traces])

# ╔═╡ 78e9243e-4981-4489-b120-6fa38352ded3
@bind selected_channel_index Select([i => channels[i] for i in 1:length(channels)])

# ╔═╡ 84c0bc78-7f66-4a6d-88ee-af6ea5154e49
begin
    # Compute distance and arrivals
    distance_m, az, baz = obspy.geodetics.gps2dist_azimuth(eq_lat, eq_lon, selected_station_details[2], selected_station_details[3])
    epi_dist_deg = obspy.geodetics.kilometer2degrees(distance_m / 1000.0)

    arrivals = model.get_ray_paths(source_depth_in_km=eq_dep,
        distance_in_degree=epi_dist_deg)
end

# ╔═╡ fc1f9bcd-2a4e-4468-a066-cf3540578e27
let

    trace = traces[selected_channel_index-1]
    times = pyconvert(Vector, trace.times("relative"))
    data = pyconvert(Vector, trace.data)

    # Base trace
    seis_trace = scatter(x=times, y=data, mode="lines", name="Seismogram")
    phase_colors = distinguishable_colors(length(arrivals))
    # Mark arrivals
    arrival_traces = map(1:length(arrivals)) do i
        arr_time = pyconvert(Float64, arrivals[i-1].time)
        arr_name = arrivals[i-1].name
        if 0.0 < arr_time < total_duration
            shape = scatter(x=[arr_time, arr_time],
                y=[minimum(data), maximum(data)],
                mode="lines",
                name=pyconvert(String, arr_name),
                line=attr(dash="solid", color=phase_colors[i]))
            return shape
        else
            return nothing
        end
    end

    # Filter out any `nothing` values
    arrival_traces = filter(!isnothing, arrival_traces)

    # Combine traces
    plot([seis_trace; arrival_traces...], Layout(title="Seismogram at Station $(selected_station_details[1])",
        xaxis_title="Time Relative to Earthquake Origin (s)",
        yaxis_title="Amplitude"))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Colors = "~0.12.11"
CondaPkg = "~0.2.29"
FFTW = "~1.9.0"
PlutoPlotly = "~0.6.4"
PlutoUI = "~0.7.68"
PythonCall = "~0.9.26"
StatsBase = "~0.34.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "0034747ae103e5bc15fd622859acb9755389869e"

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
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

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
git-tree-sha1 = "3a3dfb30697e96a440e4149c8c51bf32f818c0f3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.17.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CondaPkg]]
deps = ["JSON3", "Markdown", "MicroMamba", "Pidfile", "Pkg", "Preferences", "Scratch", "TOML", "pixi_jll"]
git-tree-sha1 = "93e81a68a84dba7e652e61425d982cd71a1a0835"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.29"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

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
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "797762812ed063b9b94f6cc7742bc8883bb5e69e"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.9.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

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
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "411eccfe8aba0814ffa0fdf4860913ed09c34975"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.3"

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
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

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
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

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
git-tree-sha1 = "28278bb0053da0fd73537be94afd1682cc5a0a83"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.21"

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
git-tree-sha1 = "232630fee92e588c11c2b260741b4fa70784b4c5"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.4"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ec9e63bd098c50e4ad28e7cb95ca7a4860603298"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.68"

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
git-tree-sha1 = "f03464b21983fb5af2f8cea99106b8d8f48ac69d"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.26"

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
git-tree-sha1 = "7f44eef6b1d284465fafc66baf4d9bdcc239a15b"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.4.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

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
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "b81c5035922cc89c2d9523afc6c54be512411466"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.5"

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
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

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
git-tree-sha1 = "0fc001395447da85495b7fef1dfae9789fdd6e31"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.11"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

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
git-tree-sha1 = "2ca2ac0b23a8e6b76752453e08428b3b4de28095"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "1.5.12+0"

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
# ╠═ca9f95ec-71d2-11f0-2e94-7da041469cdf
# ╟─5a72389a-1277-4e18-a532-ce81a4b18c74
# ╟─18eb286d-263e-42ee-94fc-26f3223ce108
# ╟─0ec9801f-515e-4e42-9bc0-8c96c1e8937f
# ╟─1939e9df-8c71-45b9-981f-e410d0a9965b
# ╟─ad4d1627-f3ec-4a34-b109-bf07b1ce5bc0
# ╟─78e9243e-4981-4489-b120-6fa38352ded3
# ╟─fc1f9bcd-2a4e-4468-a066-cf3540578e27
# ╟─d7023542-1db3-42fa-a4c6-f716f8a6fcfb
# ╠═3e20d6bd-a212-44a2-b7ff-0abba9a36a90
# ╠═97de4846-f1ba-4269-bb25-68068dca274c
# ╟─eddf7e8c-21e3-4d0d-9cfb-dacb24047b44
# ╠═834f3dae-81dd-4ff4-a891-2e62a1e10ce8
# ╠═c0ba9e33-0365-44ae-8ed9-f3d67d0e20c4
# ╠═1617897a-dd18-4c80-b12a-c2d063965e09
# ╠═44944cf9-b2c7-4264-b766-6163e3946d86
# ╟─0a78aecc-932f-4aa9-9172-9d79736b5301
# ╠═ad11d31e-2aed-4fc3-9e77-17c40c907e83
# ╠═5c05410d-d364-4f32-a584-5601874887d0
# ╠═caf0a5ef-a65b-4e11-bedb-c6a5463f6387
# ╠═4be1db35-014e-4a56-8006-a7fe2ed358f1
# ╟─9ef37a61-295d-45ca-b686-9592cb417e5a
# ╠═fd936886-101c-4244-898d-dfcf1333555b
# ╠═04687d41-3a81-4413-acb8-8c2293e39139
# ╟─98b22766-a9b2-419e-8ec2-6287f65e2722
# ╠═84c0bc78-7f66-4a6d-88ee-af6ea5154e49
# ╠═4b0309a1-e369-4954-8b40-e58fbe41754f
# ╟─d77d897a-ca34-4a00-a5c2-fdc988c6b75c
# ╠═68087281-2e24-40a5-872e-6d549bdd2eaa
# ╠═00670b47-a0a1-43aa-b304-3a63bd1615c2
# ╠═ce54afb1-9a66-4d88-aced-e8b94b162776
# ╠═dd849802-4614-4f30-bbce-d87660df3840
# ╟─d2aba787-dae6-4da5-8f16-0cadfcc38371
# ╠═2cc62eeb-edc2-4fcb-876b-010d8a180162
# ╠═4afe5b07-8c6c-463f-aa77-c866d378ea61
# ╠═1d74bb4f-eb8b-4e6a-bb12-9f8287bfa05b
# ╠═20bf25f7-a56d-4e01-a4cd-ec7506028755
# ╠═4ee50f18-9d36-43c3-aced-6b33a661a959
# ╠═08001b34-7cd7-4b72-9fc8-b4a8c869c460
# ╠═b4199996-5f1b-41b7-8732-6971a8d1824c
# ╠═53f5276e-7954-43d5-bcfb-6898bc6e85ab
# ╠═efb48ae5-c01c-4e58-ba08-22f0c885ae4a
# ╠═277155cf-77d8-4cd6-a122-346be2ff1516
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
