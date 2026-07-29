### A Pluto.jl notebook ###
# v0.20.19

#> [frontmatter]
#> title = "Love Wave Dispersion in Real Seismogram Data"
#> tags = ["fundamentals"]
#> layout = "layout.jlhtml"
#> description = "Interactive exploration of Love wave dispersion using real earthquake data from US seismic networks"

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

# ╔═╡ bd18bcc4-a581-11f0-a887-f7159f142a24
TableOfContents(include_definitions=true)

# ╔═╡ bd18bd30-a581-11f0-8c32-5df4ea3e1a17
md"""
# Love Wave Dispersion in Real Seismogram Data

This interactive notebook demonstrates **Love wave dispersion** using real earthquake recordings from the US seismic network. You'll explore how surface waves propagate at different velocities depending on their frequency, creating the classic dispersive wave trains observed in seismology.

1. **Observe real Love waves** in earthquake recordings
2. **Understand dispersion** through frequency filtering
3. **Analyze wave propagation** across different epicentral distances
4. **Connect theory to observation** using actual seismic data

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,  
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 78677a5a-a586-11f0-b3b3-7fe00db28da7
md"""
**Select Earthquake by Index:** $(@bind selected_earthquake_index Slider(1:length(earthquake_catalog), default=1, show_value=true))

**Selected:** $(earthquake_catalog[selected_earthquake_index].name)
"""

# ╔═╡ 8553843a-89df-4800-b967-02caabb12512
md"""
##### Frequency Band Controls

Adjust the frequency band to observe dispersion effects:

**Low Frequency (Hz):** $(@bind freq_low Slider(0.005:0.005:0.1, default=0.02, show_value=true))

**High Frequency (Hz):** $(@bind freq_high Slider(0.05:0.01:0.5, default=0.1, show_value=true))

**Show Raw Data:** $(@bind show_raw CheckBox(default=false))

**Normalize Amplitudes:** $(@bind normalize_amplitudes CheckBox(default=true))
"""

# ╔═╡ 40c1ebd2-a584-11f0-a0ec-5d4ca09a95cf
begin
    if @isdefined(seismic_data) && !isempty(seismic_data) && !isnothing(earthquake_info)
        
        # Create the main dispersion plot
        fig = plot(Layout(
            title="Love Wave Dispersion: $(earthquake_info.name) (M $(earthquake_info.mag))<br>Frequency Band: $(freq_low)-$(freq_high) Hz",
            xaxis_title="Time after Origin (seconds)",
            yaxis_title="Distance (km) + Normalized Amplitude",
            height=700,
            showlegend=false,
            font=attr(size=10)
        ))
        
        # Plot each station's Love wave data
        for (i, station_data) in enumerate(seismic_data)
            time = station_data["time"]
            distance = station_data["distance_km"]
            
            # Select data to display
            amplitude = show_raw ? station_data["transverse_raw"] : station_data["transverse"]
            
            # Additional normalization for raw data
            if show_raw && normalize_amplitudes && maximum(abs.(amplitude)) > 0
                amplitude ./= maximum(abs.(amplitude))
            end
            
            # Scale amplitude for plotting clarity
            scale_factor = 40  # Adjust trace separation
            y_plot = distance .+ amplitude .* scale_factor
            
            # Color coding based on distance
            color_intensity = i / length(seismic_data)
            color = "rgb($(Int(floor(255*color_intensity))), $(Int(floor(255*(1-color_intensity)))), 100)"
            
            add_trace!(fig, scatter(
                x=time,
                y=pyconvert(Array{Float64}, y_plot),
                mode="lines",
                name="$(station_data["station"])_$(round(pyconvert(Float64, distance), digits=0))km",
                line=attr(color=color, width=1.2),
                hovertemplate="Station: $(station_data["station"])<br>" *
                            "Network: $(station_data["network"])<br>" *
                            "Distance: $(round(pyconvert(Float64, distance), digits=1)) km<br>" *
                            "Azimuth: $(round(pyconvert(Float64, station_data["azimuth"]), digits=1))°<br>" *
                            "Time: %{x:.1f} s<br>" *
                            "<extra></extra>"
            ))
        end
        
        # Add theoretical Love wave velocity curves
        if length(seismic_data) >= 2
            min_dist = pyconvert(Float64, minimum([s["distance_km"] for s in seismic_data]))
            max_dist = pyconvert(Float64, maximum([s["distance_km"] for s in seismic_data]))
            min_time = 0  # Start from origin time
            
            # Theoretical Love wave velocities (km/s)
            velocities = [3.0, 3.5, 4.0, 4.5, 5.0]
            velocity_colors = ["red", "blue", "green", "orange", "purple"]
            
            for (v, color) in zip(velocities, velocity_colors)
                t_min = min_dist / v
                t_max = max_dist / v
                
                if pyconvert(Float64, t_min) >= min_time  # Only show if visible in time window
                    add_trace!(fig, scatter(
                        x=[t_min, t_max],
                        y=[min_dist, max_dist],
                        mode="lines",
                        name="$(v) km/s",
                        line=attr(color=color, width=2, dash="dash"),
                        opacity=0.7,
                        showlegend=true
                    ))
                end
            end
        end
	end
        
        fig
end

# ╔═╡ 7031673d-23af-41d9-a254-fc57b27dd417
md"# Appendix"

# ╔═╡ bd18b862-a581-11f0-9c3a-d75bd0f4186b
begin
	using CondaPkg
	CondaPkg.add_pip("obspy")
	CondaPkg.add_pip("matplotlib")
	
	using PlutoUI, PlutoPlotly
	using PythonCall
	using Dates, LinearAlgebra
end

# ╔═╡ bd18d1a6-a581-11f0-a393-05b8ff32c765
"""
Process seismic data to extract Love waves
"""
function process_love_waves(waveforms, freq_low, freq_high, normalize_amp=true)
    processed_data = []
    
    for wf_data in waveforms
        try
            st = wf_data["stream"]
            sta_info = wf_data["station_info"]
            origin_time = wf_data["origin_time"]
            
            # Find horizontal components
            north_tr = nothing
            east_tr = nothing
            
            for tr in st
                if occursin("N", tr.stats.channel) || occursin("1", tr.stats.channel)
                    north_tr = tr
                elseif occursin("E", tr.stats.channel) || occursin("2", tr.stats.channel)
                    east_tr = tr
                end
            end
            
            if north_tr === nothing || east_tr === nothing
                continue
            end
            
            # Ensure same sampling rate and length
            if north_tr.stats.sampling_rate != east_tr.stats.sampling_rate
                continue
            end
            
            min_len = min(length(north_tr.data), length(east_tr.data))
            north_tr.data = north_tr.data[1:min_len]
            east_tr.data = east_tr.data[1:min_len]
            
            # Rotate to radial-transverse
            radial, transverse = rotate_to_rt(north_tr, east_tr, sta_info["back_azimuth"])
            
            if radial === nothing || transverse === nothing
                continue
            end
            
            # Create time array
            dt = 1.0 / north_tr.stats.sampling_rate
            t_start = (north_tr.stats.starttime - origin_time)
            time_array = collect(t_start:dt:(t_start + (length(transverse)-1)*dt))
            
            # Filter for Love waves (transverse component)
            nyquist = north_tr.stats.sampling_rate / 2
            if freq_high >= nyquist
                freq_high = nyquist * 0.9
            end
            
            # Design Butterworth filter
            sos = scipy_signal.butter(4, [freq_low, freq_high], btype="band", 
                                    fs=north_tr.stats.sampling_rate, output="sos")
            transverse_filtered = scipy_signal.sosfilt(sos, transverse)
            
            # Normalize amplitude if requested
            if normalize_amp && maximum(abs.(transverse_filtered)) > 0
                transverse_filtered ./= maximum(abs.(transverse_filtered))
            end
            
            push!(processed_data, Dict(
                "transverse" => transverse_filtered,
                "transverse_raw" => transverse,
                "time" => time_array,
                "distance_km" => sta_info["distance_km"],
                "station" => "$(sta_info["network"]).$(sta_info["station"])",
                "sampling_rate" => north_tr.stats.sampling_rate
            ))
            
        catch e
            println("Error processing station: $e")
            continue
        end
    end
    
    # Sort by distance
    sort!(processed_data, by = x -> x["distance_km"])
    
    return processed_data
end

# ╔═╡ bd18e4ea-a581-11f0-9752-d7101c5581f0
md"""
### **Data Processing Steps:**

1. **Component Rotation**: N-E → Radial-Transverse coordinates
2. **Frequency Filtering**: Isolate specific frequency bands
3. **Amplitude Normalization**: Compare relative wave shapes
4. **Time Alignment**: Relative to earthquake origin time
"""

# ╔═╡ bfd89e80-a583-11f0-8c1e-5d357ed7ab93
begin
	# Earthquake catalog for Love wave analysis
	earthquake_catalog = [
		(lat=35.7695, lon=-117.5993, mag=7.1, time="2019-07-06T03:19:53", depth=8.0, name="2019 Ridgecrest, CA"),
		(lat=61.346, lon=-149.955, mag=7.1, time="2018-11-30T17:29:29", depth=46.7, name="2018 Anchorage, AK"),
		(lat=38.215, lon=-122.312, mag=6.0, time="2014-08-24T10:20:44", depth=11.3, name="2014 Napa, CA"),
		(lat=37.936, lon=-77.933, mag=5.8, time="2011-08-23T17:51:04", depth=6.0, name="2011 Virginia"),
		(lat=40.640, lon=-125.13, mag=7.0, time="2014-03-10T05:18:13", depth=7.0, name="2014 Ferndale, CA"),
		(lat=55.308, lon=-134.79, mag=7.5, time="2013-01-05T08:58:20", depth=9.8, name="2013 Alaska"),
	]
	
	# Create earthquake plot
	# for (i, eq) in 
	eq_traces = map(enumerate(earthquake_catalog)) do (i, eq)
		marker_color = eq.mag >= 7.0 ? "red" : eq.mag >= 6.5 ? "orange" : "yellow"
		marker_size = 8 + 2 * eq.mag
		
		trace = scatter(
			x=[eq.lon], y=[eq.lat],
			mode="markers",
			marker=attr(
				size=marker_size,
				color=marker_color,
				line=attr(width=2, color="black"),
				symbol="circle"
			),
			name="M$(eq.mag) $(eq.name)",
			text=["M$(eq.mag) $(eq.name)<br>$(eq.time)<br>Depth: $(eq.depth) km"],
			hovertemplate="%{text}<extra></extra>",
			customdata=[i]
		)
		trace
	end
	
	earthquake_map = plot(
		eq_traces,
		Layout(
			title="Select Earthquake for Love Wave Analysis",
			geo=attr(
				projection_type="natural earth",
				showland=true, landcolor="lightgray",
				showocean=true, oceancolor="lightblue",
				showcountries=true, countrycolor="white",
				showlakes=true, lakecolor="lightblue",
				scope="usa",
				center=attr(lat=40, lon=-100),
				projection_scale=1
			),
			height=400,
			showlegend=false
		)
	)
end

# ╔═╡ bfd8a3ee-a583-11f0-9977-c9f4d0d5cfa5
begin
	# Display the earthquake map
	earthquake_map
end

# ╔═╡ bfd8a45c-a583-11f0-857b-d7918378d9e2
let
	if @isdefined(selected_earthquake_index)
		selected_eq = earthquake_catalog[selected_earthquake_index]
		
		Markdown.parse("""
		### 📍 Selected Earthquake: $(selected_eq.name)
		- **Magnitude:** M$(selected_eq.mag)
		- **Location:** $(round(selected_eq.lat, digits=3))°N, $(round(abs(selected_eq.lon), digits=3))°W  
		- **Time:** $(selected_eq.time)
		- **Depth:** $(selected_eq.depth) km
		
		This earthquake will be used for Love wave dispersion analysis.
		""")
	else
		md"👆 **Select an earthquake above to begin analysis**"
	end
end

# ╔═╡ 3a9ecf88-9562-4185-ad78-9c12d0eab09b
md"""
## 📊 Data Collection Parameters

**Number of Stations:** $(@bind n_stations Slider(6:2:20, default=10, show_value=true))

**Maximum Distance (km):** $(@bind max_distance Slider(300:100:2000, default=2000, show_value=true))

**Minimum Distance (km):** $(@bind min_distance Slider(50:50:500, default=200, show_value=true))
"""

# ╔═╡ 8e440781-d910-4cce-b318-95d8316af8be
begin
	# Initialize Python environment and import required libraries
	obspy = pyimport("obspy")
	obspy_clients = pyimport("obspy.clients.fdsn")
	geodetics = pyimport("obspy.geodetics")
	signal = pyimport("obspy.signal")
	numpy = pyimport("numpy")
	scipy_signal = pyimport("scipy.signal")
	
	# Initialize FDSN client
	iris_client = obspy_clients.Client("IRIS")
	
	md"✅ Python libraries loaded successfully"
end

# ╔═╡ 6c89c2c1-e2ea-4550-9ccf-84b1e77d13cf
# waveforms_result = download_love_wave_data(selected_eq, n_stations, min_distance, max_distance)

# ╔═╡ 52e85d7f-8504-45b9-b1a4-7526ecbdf9d4
"""
Download seismic data for Love wave analysis following specnm.jl patterns
"""
function download_love_wave_data(eq_params, n_stations, min_dist, max_dist)
    try
        # Parse origin time
        origin_time = obspy.UTCDateTime(eq_params.time)
        eq_lat = eq_params.lat
        eq_lon = eq_params.lon
        
        # Define time window for surface waves
        start_time = origin_time - 120  # 2 minutes before
        end_time = origin_time + max_dist/3.0 + 600  # Surface wave travel time + buffer
        
        # Get available stations within geographic bounds
        lat_range = 25.0  # degrees
        lon_range = 35.0  # degrees
        
        inventory = iris_client.get_stations(
            network="US,TA,IU,II,N4",  # Focus on high-quality networks
            station="*", 
            location="*",
            channel="BH*,HH*",
            starttime=start_time,
            endtime=end_time,
            minlatitude=max(eq_lat - lat_range, -90),
            maxlatitude=min(eq_lat + lat_range, 90),
            minlongitude=max(eq_lon - lon_range, -180),
            maxlongitude=min(eq_lon + lon_range, 180),
            level="station"
        )
        
        # Calculate distances and collect station information
        stations_data = []
        for network in inventory
            for station in network
				@show station
                try
                    distance_m, azimuth, back_azimuth = geodetics.gps2dist_azimuth(
                        eq_lat, eq_lon, station.latitude, station.longitude
                    )
                    distance_km = distance_m / 1000.0
                    if min_dist <= pyconvert(Float64, distance_km) <= max_dist
                        push!(stations_data, Dict(
                            "network" => network.code,
                            "station" => station.code,
                            "latitude" => station.latitude,
                            "longitude" => station.longitude,
                            "distance_km" => distance_km,
                            "azimuth" => azimuth,
                            "back_azimuth" => back_azimuth
                        ))
                    end
                catch e
                    continue  # Skip stations with geodetic calculation issues
                end
            end
        end
        
        # Sort by distance and select best stations
        sort!(stations_data, by = x -> x["distance_km"])
        selected_stations = stations_data[1:min(n_stations, length(stations_data))]
        
        # Download waveform data
        successful_downloads = []
        for sta_info in selected_stations
            try
                # Request 3-component data
                st = iris_client.get_waveforms(
                    sta_info["network"], sta_info["station"], "*", "BH*,HH*",
                    start_time, end_time
                )
                
                # Basic preprocessing similar to specnm.jl
                st.detrend("demean")
                st.detrend("linear")
                st.taper(0.05)
                
                # Pre-filter to focus on surface wave band
                st.filter("bandpass", freqmin=0.005, freqmax=0.5, corners=2)
                
                # Store with metadata
                push!(successful_downloads, Dict(
                    "stream" => st,
                    "station_info" => sta_info,
                    "origin_time" => origin_time
                ))
                
                if length(successful_downloads) >= n_stations
                    break
                end
                
            catch e
                println("Download failed for $(sta_info["station"]): $e")
                continue
            end
        end
        
        return successful_downloads, length(stations_data)
        
    catch e
        println("Data download error: $e")
        return [], 0
    end
end

# ╔═╡ ce884ca0-a583-11f0-bafc-6d33854697d4
# begin
# 	# Only attempt data download if earthquake is selected
# 	if !isnothing(earthquake_clicked) && haskey(earthquake_clicked, "points") && !isempty(earthquake_clicked["points"])
# 		eq_index = earthquake_clicked["points"][1]["customdata"] + 1
# 		selected_eq = earthquake_catalog[eq_index]
		
# 		try
# 			md"🔄 **Downloading seismic data...** (this may take a moment)"
# 		catch e
# 			md"❌ **Error:** Please try again or select a different earthquake."
# 		end
# 	else
# 		md"⏳ **Waiting for earthquake selection...**"
# 	end
# end

# ╔═╡ ea24374e-a583-11f0-b189-b906b9fd5b91
waveforms_result = let
	selected_eq = earthquake_catalog[selected_earthquake_index]
	
	
		waveforms_result = download_love_wave_data(selected_eq, n_stations, min_distance, max_distance)
		
		if length(waveforms_result[1]) > 0
			Markdown.parse("""
			✅ **Data Download Complete**
			- **Downloaded:** $(length(waveforms_result[1])) stations
			- **Available stations:** $(waveforms_result[2]) total
			- **Distance range:** $(min_distance) - $(max_distance) km
			""")
		else
			md"❌ **No data available** - Try adjusting distance parameters or selecting a different earthquake"
		end
	
	waveforms_result
end

# ╔═╡ ea243a14-a583-11f0-8564-3d89bcf5f405
# ╠═╡ disabled = true
#=╠═╡
# let
# 	selected_eq = earthquake_catalog[selected_earthquake_index]
# 	waveforms_data = try
# 		download_love_wave_data(selected_eq, n_stations, min_distance, max_distance)[1]
# 	catch e
# 		println("Download error: $e")
# 		[]
# 	end
	
# 	if !isempty(waveforms_data)
# 		md"✅ **Downloaded $(length(waveforms_data)) stations successfully**"
# 	else
# 		md"❌ **No data downloaded - check parameters**"
# 	end
# end
  ╠═╡ =#

# ╔═╡ 13e2b1b5-c187-4929-ac97-1c5bdfd1c1e0
"""
Rotate horizontal components to radial and transverse following specnm.jl patterns
"""
function rotate_to_rt(north_data, east_data, back_azimuth)
    try
        # Convert back azimuth to radians
        baz_rad = deg2rad(pyconvert(Float64, back_azimuth))
        
        # Rotation matrix for N-E to R-T
        # R = N*cos(baz) + E*sin(baz) (positive away from source)
        # T = -N*sin(baz) + E*cos(baz) (positive for left-lateral motion)
        radial = north_data .* cos(baz_rad) .+ east_data .* sin(baz_rad)
        transverse = -north_data .* sin(baz_rad) .+ east_data .* cos(baz_rad)
        
        return radial, transverse
        
    catch e
        return nothing, nothing
    end
end

# ╔═╡ a50fcf5d-ef3f-4d75-a20e-6d36ccaa08af
seismic_data = process_love_wave_signals(waveforms_data, freq_low, freq_high, normalize_amplitudes)

# ╔═╡ 6aaafe83-5ddf-42b0-a401-59cd93ba56fa
"""
Process seismic data for Love wave analysis with enhanced error handling
"""
function process_love_wave_signals(waveforms_data, freq_low, freq_high, normalize_amp=true)
    if isempty(waveforms_data)
        return []
    end
    
    processed_data = []
    
    for wf_data in waveforms_data
        try
            st = wf_data["stream"]
			st = st.resample(inv(4*freq_high))
            sta_info = wf_data["station_info"]
            origin_time = wf_data["origin_time"]
            
            # Find horizontal components (N and E)
            north_tr = nothing
            east_tr = nothing
            
            for tr in st
                channel = pyconvert(String, tr.stats.channel)
                if endswith(channel, 'N') || endswith(channel, '1')
                    north_tr = tr
                elseif endswith(channel, 'E') || endswith(channel, '2')
                    east_tr = tr
                end
            end
        
            if isnothing(north_tr) || isnothing(east_tr)
                continue  # Skip if missing horizontal components
            end
            
            # Extract data arrays
            north_data = pyconvert(Array, north_tr.data)
            east_data = pyconvert(Array, east_tr.data)

            # Ensure same length
            min_length = min(length(north_data), length(east_data))
            north_data = north_data[1:min_length]
            east_data = east_data[1:min_length]
            
            # Rotate to radial and transverse components
            radial, transverse = rotate_to_rt(north_data, east_data, sta_info["back_azimuth"])
            
            if isnothing(transverse)
                continue
            end

            # Store raw transverse component
            transverse_raw = copy(transverse)
            
            # Apply frequency filtering to isolate Love waves
            sampling_rate = pyconvert(Float64, north_tr.stats.sampling_rate)
            nyquist = sampling_rate / 2.0
            
            freq_high_adj = min(freq_high, nyquist * 0.95)  # Avoid Nyquist frequency issues
            
            # Design Butterworth bandpass filter
            filter_order = 4
            low_norm = freq_low / nyquist
            high_norm = freq_high_adj / nyquist
            
            # Use scipy for robust filtering
            sos = scipy_signal.butter(filter_order, [low_norm, high_norm], btype="band", output="sos")
            transverse_filtered = pyconvert(Array, scipy_signal.sosfilt(sos, transverse))
            
            # Create time array relative to earthquake origin
            start_time_utc = pyconvert(Float64, north_tr.stats.starttime.timestamp)
            origin_time_utc = pyconvert(Float64, origin_time.timestamp)
            
            dt = 1.0 / sampling_rate
            time_relative = collect(0:dt:(length(transverse_filtered)-1)*dt) .+ (start_time_utc - origin_time_utc)
            
            # Normalize amplitudes if requested
            if normalize_amp && maximum(abs.(transverse_filtered)) > 0
                transverse_filtered ./= maximum(abs.(transverse_filtered))
            end
            
            # Store processed data
            push!(processed_data, Dict(
                "station" => sta_info["station"],
                "network" => sta_info["network"],
                "distance_km" => sta_info["distance_km"],
                "azimuth" => sta_info["azimuth"],
                "back_azimuth" => sta_info["back_azimuth"],
                "time" => time_relative,
                "transverse" => transverse_filtered,
                "transverse_raw" => transverse_raw,
                "sampling_rate" => sampling_rate
            ))
            
        catch e
            println("Processing failed for station: $e")
            continue
        end
    end
    
    # Sort by distance for better plotting
    sort!(processed_data, by = x -> x["distance_km"])
    
    return processed_data
end

# ╔═╡ cfe0821b-c768-4824-b8b6-f1b9808641af
waveforms_data =waveforms_result[1]

# ╔═╡ a5620ac8-5a4b-4c55-9f5e-980c0a6d0995
earthquake_info = earthquake_catalog[selected_earthquake_index]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"

[compat]
CondaPkg = "~0.2.33"
PlutoPlotly = "~0.6.5"
PlutoUI = "~0.7.72"
PythonCall = "~0.9.28"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "505cfdf1a79ea0ac44d582f31df7d2e3265451c1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

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
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CondaPkg]]
deps = ["JSON3", "Markdown", "MicroMamba", "Pidfile", "Pkg", "Preferences", "Scratch", "TOML", "pixi_jll"]
git-tree-sha1 = "bd491d55b97a036caae1d78729bdb70bf7dababc"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.33"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

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

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

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
git-tree-sha1 = "8acd04abc9a636ef57004f4c2e6f3f6ed4611099"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.5"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "34510e11cabd7964291f32f14d28b367e9960e6e"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.28"

    [deps.PythonCall.extensions]
    CategoricalArraysExt = "CategoricalArrays"
    PyCallExt = "PyCall"

    [deps.PythonCall.weakdeps]
    CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"

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
git-tree-sha1 = "c3b2323466378a2ba15bea4b2f73b081e022f473"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.5.0"

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

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

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
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

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
# ╠═bd18bcc4-a581-11f0-a887-f7159f142a24
# ╟─bd18bd30-a581-11f0-8c32-5df4ea3e1a17
# ╟─78677a5a-a586-11f0-b3b3-7fe00db28da7
# ╟─8553843a-89df-4800-b967-02caabb12512
# ╟─40c1ebd2-a584-11f0-a0ec-5d4ca09a95cf
# ╟─7031673d-23af-41d9-a254-fc57b27dd417
# ╠═bd18b862-a581-11f0-9c3a-d75bd0f4186b
# ╠═bd18d1a6-a581-11f0-a393-05b8ff32c765
# ╟─bd18e4ea-a581-11f0-9752-d7101c5581f0
# ╠═bfd89e80-a583-11f0-8c1e-5d357ed7ab93
# ╠═bfd8a3ee-a583-11f0-9977-c9f4d0d5cfa5
# ╠═bfd8a45c-a583-11f0-857b-d7918378d9e2
# ╟─3a9ecf88-9562-4185-ad78-9c12d0eab09b
# ╠═8e440781-d910-4cce-b318-95d8316af8be
# ╠═6c89c2c1-e2ea-4550-9ccf-84b1e77d13cf
# ╠═52e85d7f-8504-45b9-b1a4-7526ecbdf9d4
# ╠═ce884ca0-a583-11f0-bafc-6d33854697d4
# ╠═ea24374e-a583-11f0-b189-b906b9fd5b91
# ╠═ea243a14-a583-11f0-8564-3d89bcf5f405
# ╠═13e2b1b5-c187-4929-ac97-1c5bdfd1c1e0
# ╠═a50fcf5d-ef3f-4d75-a20e-6d36ccaa08af
# ╠═6aaafe83-5ddf-42b0-a401-59cd93ba56fa
# ╠═cfe0821b-c768-4824-b8b6-f1b9808641af
# ╠═a5620ac8-5a4b-4c55-9f5e-980c0a6d0995
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
