### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "PREM Earth Free Oscillations"
#> tags = ["normalmodes"]
#> layout = "layout.jlhtml"
#> description = "Free oscillations of spherically symmetric bodies"

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

# ╔═╡ a751c6dc-414b-4ce8-b1b1-5aef8e23ab63
begin
	# `specnm`'s native deps (scipy/petsc4py/mpi4py/Cython) don't resolve cleanly
	# through CondaPkg's managed micromamba env, so this notebook deliberately points
	# PythonCall at a system Python with `specnm` installed manually instead.
	ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
	local python_exe = "/opt/miniconda3/bin/python3"
	isfile(python_exe) || error("specnm.jl expects a system Python at $python_exe with the `specnm` package installed manually (see the Environment & Package Setup section for why CondaPkg is disabled here) -- edit `python_exe` above if this machine's Python lives elsewhere.")
	ENV["JULIA_PYTHONCALL_EXE"] = python_exe
	using PythonCall
end

# ╔═╡ a9019660-1f03-4132-bccf-09bdb1421ad9
using CondaPkg

# ╔═╡ 447b4e82-fe73-11ef-30b1-69824c8e3d24
using PlutoPlotly

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

# ╔═╡ f2a1c3d4-7e6b-4a5c-9d8e-1b2c3d4e5f6a
md"""
##### What are spheroidal and toroidal modes?

A planet's free oscillations split into two independent families, both indexed by
angular degree `l` (spherical-harmonic degree, more nodal lines as `l` grows) and
overtone number `n` (radial order -- `n=0` is the fundamental, higher `n` adds radial
nodes):

- **Spheroidal (`ₙSₗ`)**: motion has both radial and horizontal components, coupled
  to gravity and compression -- shown as **U** (radial displacement), **V**
  (horizontal displacement), and, for a self-gravitating solve, **P** (perturbation
  to the gravitational potential).
- **Toroidal (`ₙTₗ`)**: purely horizontal, rotational shearing motion with no radial
  component and no coupling to gravity -- shown as **W** (horizontal displacement).
  Toroidal modes don't exist in a fluid (they need shear rigidity), so they vanish
  wherever `vs = 0` in the Earth model.

The eigenfunction panel plots each component versus depth for whichever mode you
click in the spectrum -- their zero-crossings are exactly the mode's radial nodes,
and the amplitude near the surface versus near the core tells you how deeply that
mode samples the planet.
"""

# ╔═╡ 82cc9219-7633-41c6-91e1-17968904b2b6
md"## Real Data"

# ╔═╡ 921c2136-563a-4219-98f5-21cafa67e516
md"""
Start time after earthquake (min) $(@bind starttime_cut Slider(range(1, 10 * 24 * 60, length=1000), show_value=true, default=12*60))
"""

# ╔═╡ fe54d35f-8017-4bc0-ac6c-7dee815d1477
md"""
End time after earthquake (min) $(@bind endtime_cut Slider(range(starttime_cut, 10 * 24 * 60, length=1000), show_value=true, default=5*24*60))
"""

# ╔═╡ 89c81d3c-67c1-4cc7-b892-d64208d845f2
md"### Plots"

# ╔═╡ 9b0a9a80-e65f-4385-a377-372e408b19ad
md"## Appendix"

# ╔═╡ 2afaeeec-6fe2-48db-81de-1fa4ccb0fc1b
md"### Get Stations"

# ╔═╡ 79406ee0-6026-4da8-a29d-245048c27e47
obspy = pyimport("obspy")

# ╔═╡ acd7956a-855c-43f0-a353-7d5533f6aaf1
fdsn = pyimport("obspy.clients.fdsn")

# ╔═╡ 335162fc-52c3-4c3d-bb88-d2e2f4fde37f
client = fdsn.Client("IRIS")

# ╔═╡ 58baaf3c-a5b3-4995-9f4e-c54946f7e798
specnm = pyimport("specnm")

# ╔═╡ 5852bfac-a8ff-405e-8ed1-6438c6827091
md"### Select Earth Model & Solve"

# ╔═╡ d99b1935-4f90-4b82-809c-b6a801c37e0d
"""
All 13 model files under `src/specnm_models/`, spanning specnm's three
auto-detected formats (deck, layered, poly).
"""
const SPECNM_MODEL_NAMES = ["europa", "homo-full-sphere", "homo-full-sphere-att",
    "homo-full-sphere-dbl", "homo-full-sphere2", "mars", "mars_lvl", "prem_ani",
    "prem_iso_one_crust", "prem_noat", "prem_noocean", "prem_noocean_noat", "vpremoon"]

# ╔═╡ 1b22f312-9197-4044-ba9b-1f12789b88ff
"""
    specnm_overtones(angular_orders, l1_start)

Overtone number `n` for each mode in `angular_orders` (sorted by degree `l` as
`specnm` returns them): within each `l`, overtones run `0, 1, 2, ...`, except at
`l == 1` where the lowest `l1_start` overtones don't exist (translation/rotation
modes for spheroidal, `l1_start = 2`; the missing `0T1` for toroidal, `l1_start = 1`)
and numbering starts from `l1_start` instead.
"""
function specnm_overtones(angular_orders::AbstractVector{<:Integer}, l1_start::Integer)
    l_countmap = countmap(angular_orders)
    l_uniq = collect(keys(l_countmap))
    l_count = collect(values(l_countmap))
    overtones = vcat([lu != 1 ? collect(0:lc-1) : collect(l1_start:lc+l1_start-1) for (lu, lc) in zip(l_uniq, l_count)])
    return collect(reduce(vcat, overtones))
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
# Same @bind-nested-argument ordering bug as cls_selected/clicked_mode above --
# `earthquakes` is referenced only inside a nested list comprehension, so Pluto's
# static dependency analysis can miss the edge on a fresh restart.
begin
	earthquakes
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
begin
	stations; station_names
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
end;

# ╔═╡ 7493f1d5-2017-416b-bc0c-183767bad68a
begin
	starttime_data = eq_time  # Roughly start half a day after earthquake time
	endtime_data = starttime_data + 10 * 24 * 3600  # 10 hours of data
	st = client.get_waveforms(
	        "IU",  # Network (IU, IC, II)
	        selected_station[1],  # Station code
	        "*",  # Any location
	        "LHE",  # Vertical component, 
		attach_response=true,
	        starttime_data, endtime_data
	    )
	trace = st[0]
	trace.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=true)
	data = pyconvert(Array, trace.normalize().detrend().data)
	sampling_rate = pyconvert(Float32, trace.stats.sampling_rate)
	
end

# ╔═╡ 3a325fab-9cf2-4cb8-a645-f0731e40b1f5
istarttime_cut = round(Int, starttime_cut * 60 / sampling_rate)

# ╔═╡ a7d64a11-17c4-4aa0-a10f-e47d39b1eebb
iendtime_cut = round(Int, endtime_cut * 60 / sampling_rate)

# ╔═╡ c168a471-f21f-4175-b7a0-5ae44d5d30ee
begin
	data_cut = data[istarttime_cut:iendtime_cut]
	
	# Compute FFT
	n = length(data_cut)
	freqs = rfftfreq(n, 1/sampling_rate)
	
	# Convert frequency to mHz (as normal modes are in mHz)
	freqs_mHz = freqs .* 1000

	
	ifmax = findfirst(x->x>freqmax*1000.0, freqs_mHz)
	ifmin = findfirst(x->x>freqmin*1000.0, freqs_mHz)
	freqs_mHz = freqs_mHz[ifmin:ifmax]
	freqs = freqs[1:ifmax]
	
	spectrum = abs.(rfft(data_cut))[ifmin:ifmax]
	
end;

# ╔═╡ b57d88b9-bdd5-447b-ba24-30a7b1ea1e0b
plot(data_cut)

# ╔═╡ 39710dda-1054-491e-b9d5-7d8c4cb8e3a6
"""
    specnm_type_payload(cls, cls_out, angular_orders, frequencies, overtones)

Serialize one solved class's (spheroidal or toroidal) full mode catalogue and
eigenfunctions to a JSON object string for [`SpecnmBrowseView`](@ref): `l`/`f`/`n`
per mode, component `labels` (`["U","V"]`, `["U","V","P"]`, or `["W"]`), the shared
`depth` grid, and `ef` -- one flat **mode-major** array per component
(`ef[c][mode*nradial + i]`, so the widget can slice out one mode's curve by index
without any further Julia round-trip).
"""
function specnm_type_payload(cls, cls_out, angular_orders, frequencies, overtones)
	num(x) = isfinite(x) ? string(round(Float64(x), digits=5)) : "0"
	jsonarr(v) = "[" * join(num.(v), ",") * "]"
	jsonintarr(v) = "[" * join(round.(Int, v), ",") * "]"

	mode = pyconvert(String, cls.mode)
	eigenfunctions = pyconvert(Matrix, cls_out["eigenfunctions"])
	rad_grid = pyconvert(Array, cls.r ./ 1000.0)
	radius_planet = pyconvert(Float64, cls.radius / 1000.0)
	depth_grid = radius_planet .- rad_grid

	if mode == "spheroidal"
		sph_type = pyconvert(Int, cls.sph_type)
		if sph_type == 3
			labels = ["U", "V", "P"]
			efs = [eigenfunctions[:, 1:3:end], eigenfunctions[:, 2:3:end], eigenfunctions[:, 3:3:end]]
		else
			labels = ["U", "V"]
			efs = [eigenfunctions[:, 1:2:end], eigenfunctions[:, 2:2:end]]
		end
	else
		labels = ["W"]
		efs = [eigenfunctions]
	end

	ef_flat = [jsonarr(vec(permutedims(ef))) for ef in efs]
	string(
		"{\"l\":", jsonintarr(angular_orders),
		",\"f\":", jsonarr(frequencies),
		",\"n\":", jsonintarr(overtones),
		",\"labels\":[", join(["\"$l\"" for l in labels], ","), "]",
		",\"depth\":", jsonarr(depth_grid),
		",\"nradial\":", length(depth_grid),
		",\"ef\":[", join(ef_flat, ","), "]",
		"}",
	)
end

# ╔═╡ ecee2903-a103-498a-9aaf-b5d920f531f4
md"### The Interactive Widgets"

# ╔═╡ 71af98fe-12cb-42bf-b447-5afde0482e47
begin
    """
        SpecnmSolveInput(; model="prem_ani", fmax=0.01)

    Tier A -- the expensive, gated widget. Picks an Earth model + max frequency; the
    actual `specnm.rayleigh`/`specnm.love` solve only fires on the **Compute** button
    click (never on dropdown/slider drag), since mesh construction alone takes 40+
    seconds and the full solve several minutes. Bound value is a `Dict{String,Any}`
    with keys `"model"`, `"fmax"` (Pluto's default JS-object bond transport).
    """
    struct SpecnmSolveInput
        model::String
        fmax::Float64
    end

    SpecnmSolveInput(; model="prem_ani", fmax=0.01) = SpecnmSolveInput(model, Float64(fmax))

    Base.get(w::SpecnmSolveInput) = Dict{String,Any}("model" => w.model, "fmax" => w.fmax)

    function Base.show(io::IO, ::MIME"text/html", w::SpecnmSolveInput)
        model_options = join(["<option value=\"$m\"" * (m == w.model ? " selected" : "") * ">$m</option>" for m in SPECNM_MODEL_NAMES], "\n          ")
        write(io, """
<div id="ssiwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#ssiwidget){width:min(80vw,1000px)!important;margin-left:calc((100% - min(80vw,1000px))/2)!important}
    #ssiwidget{width:100%;box-sizing:border-box;color:#d1d5db;font:14px sans-serif}
    #ssiwidget .ssi-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #ssiwidget .ssi-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #ssiwidget .ssi-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #ssiwidget .ssi-controls{width:100%;display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:8px;font:14px sans-serif}
    #ssiwidget .ssi-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
    #ssiwidget .ssi-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #ssiwidget .ssi-control-row{display:grid;grid-template-columns:minmax(70px,90px) minmax(70px,1fr) minmax(50px,70px);gap:6px;align-items:center;margin:6px 0}
    #ssiwidget .ssi-control-row input[type=range]{width:100%;min-width:0}
    #ssiwidget .ssi-value{color:#d1d5db;text-align:right;font-variant-numeric:tabular-nums}
    #ssiwidget select{width:100%;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:5px 6px}
    #ssiwidget .ssi-actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    #ssiwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 14px;font-size:14px;cursor:pointer}
    #ssiwidget button:disabled{opacity:0.5;cursor:default}
  </style>
  <div class="ssi-title">
    <div class="ssi-title-desc">Pick an Earth model and frequency cutoff, then solve for its normal modes.</div>
    <div class="ssi-title-hint">the spectral-element solve takes anywhere from tens of seconds to several minutes -- it only runs on Compute, never while you're still adjusting the controls</div>
  </div>
  <div class="ssi-controls">
    <div class="ssi-control-group">
      <div class="ssi-control-title">Earth Model</div>
      <select id="ssi-model">
          $(model_options)
      </select>
    </div>
    <div class="ssi-control-group">
      <div class="ssi-control-title">Max Frequency</div>
      <label class="ssi-control-row"><span>fmax (Hz)</span><input type="range" id="ssi-fmax" min="0.004" max="0.02" step="0.001" value="$(w.fmax)"><span id="ssi-fmax-v" class="ssi-value">$(w.fmax)</span></label>
    </div>
    <div class="ssi-control-group">
      <div class="ssi-control-title">Compute</div>
      <div class="ssi-actions">
        <button id="ssi-compute" type="button">Compute</button>
        <span id="ssi-status" style="font-size:13px">idle</span>
      </div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  let model = "$(w.model)", fmax = $(w.fmax)

  par.querySelector('#ssi-fmax').addEventListener('input', e=>{
    fmax = parseFloat(e.target.value)
    par.querySelector('#ssi-fmax-v').textContent = fmax.toFixed(3)
  })
  par.querySelector('#ssi-model').addEventListener('change', e=>{ model = e.target.value })

  const computeBtn = par.querySelector('#ssi-compute')
  const statusEl = par.querySelector('#ssi-status')
  computeBtn.addEventListener('click', ()=>{
    computeBtn.disabled = true
    statusEl.style.color = '#9ca3af'
    statusEl.textContent = 'computing... (can take a few minutes)'
    par.value = {model, fmax}
    par.dispatchEvent(new CustomEvent('input'))
  })

  window.addEventListener('specnm-solve-done', e=>{
    const d = e.detail ? JSON.parse(e.detail) : null
    computeBtn.disabled = false
    if(!d){ return }
    if(d.ok){ statusEl.style.color = '#4ade80'; statusEl.textContent = 'done' }
    else { statusEl.style.color = '#f87171'; statusEl.textContent = 'error: ' + (d.error || 'solve failed') }
  })
</script>
""")
    end

    # Forces the correct execution order on a fresh restart (Pluto's static
    # dependency analysis doesn't reliably detect that the bind cell below depends
    # on this cell defining SpecnmSolveInput) -- same pattern as `_epi_ready` in
    # earth-internal-structure.jl.
    const _ssi_ready = true
end

# ╔═╡ d9c71796-e79f-404d-89ef-54adfbd3335c
begin
	_ssi_ready
	@bind specnm_solve SpecnmSolveInput()
end

# ╔═╡ beddae25-17bc-48e9-8eef-41f21a08fb10
begin
	model_fname = specnm_solve isa AbstractDict ? "../specnm_models/" * specnm_solve["model"] : "../specnm_models/prem_ani"
	fmax_value = specnm_solve isa AbstractDict ? specnm_solve["fmax"] : 0.01
end;

# ╔═╡ 69699338-d461-4e05-bbb4-a874ad4ad970
"""
The actual (expensive) solve, gated behind `SpecnmSolveInput`'s Compute button --
`model_fname`/`fmax_value` only change on a button click, never on drag/hover, so
this cell (and the several-minutes-long spectral-element solve inside it) only
reruns when the user explicitly asks for it. Wrapped in try/catch so a bad
model/config (several of `SPECNM_MODEL_NAMES` are synthetic test fixtures with
unverified solver behavior) surfaces a friendly in-widget error instead of an
unreadable Python stacktrace cascading through a dozen downstream cells.
"""
specnm_result = let
	try
		ray = specnm.rayleigh(model_fname, fmax=fmax_value)
		ray_out = ray.rayleigh_problem(attenuation_mode="elastic", fmax=fmax_value / 2)
		ray_angular_orders = pyconvert(Array, ray_out["angular orders"])
		ray_frequencies = pyconvert(Array, ray_out["frequencies"] * 1000.0) # in mHz
		ray_overtones = specnm_overtones(pyconvert(Vector{Int}, ray_out["angular orders"]), 2)

		lov = specnm.love(model_fname, fmax=fmax_value)
		# NOTE: attenuation_mode differs from Rayleigh's ("elastic" vs "full") in the
		# original notebook -- kept as-is (not silently "fixed") since it may be an
		# intentional choice by the notebook's author; worth confirming with them.
		lov_out = lov.love_problem(attenuation_mode="full", fmax=fmax_value / 2)
		lov_angular_orders = pyconvert(Array, lov_out["angular orders"])
		lov_frequencies = pyconvert(Array, lov_out["frequencies"] * 1000.0) # in mHz
		lov_overtones = specnm_overtones(pyconvert(Vector{Int}, lov_out["angular orders"]), 1)

		(; ok=true, error=nothing, ray, ray_out, ray_angular_orders, ray_frequencies, ray_overtones,
			lov, lov_out, lov_angular_orders, lov_frequencies, lov_overtones)
	catch e
		(; ok=false, error=sprint(showerror, e), ray=nothing, ray_out=nothing,
			ray_angular_orders=Int[], ray_frequencies=Float64[], ray_overtones=Int[],
			lov=nothing, lov_out=nothing, lov_angular_orders=Int[], lov_frequencies=Float64[], lov_overtones=Int[])
	end
end;

# ╔═╡ 55b7b9a8-caea-4ec6-b47d-937231eaaee8
# Aliases restoring the plain top-level names (`ray`, `ray_out`, ...) the rest of the
# notebook already uses -- keeps `cls_selected` and friends unchanged below.
begin
	ray = specnm_result.ray
	ray_out = specnm_result.ray_out
	ray_angular_orders = specnm_result.ray_angular_orders
	ray_frequencies = specnm_result.ray_frequencies
	ray_overtones = specnm_result.ray_overtones
	lov = specnm_result.lov
	lov_out = specnm_result.lov_out
	lov_angular_orders = specnm_result.lov_angular_orders
	lov_frequencies = specnm_result.lov_frequencies
	lov_overtones = specnm_result.lov_overtones
end;

# ╔═╡ 1a67d073-bf07-4b01-beea-6c11013965be
let
	# Create spectrum plot
	spectrum_plot = [scatter(
	    x=freqs_mHz[freqs_mHz .> 0], y=spectrum[freqs_mHz .> 0],
	    mode="lines", line=attr(color="blue"),
		name="",
	)]

	# Add vertical lines at each mode frequency (spheroidal, matching the dropdown's
	# former default selection -- the Spheroidal/Toroidal picker moved into
	# SpecnmBrowseView above and no longer produces a `cls_selected` variable here)
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

# ╔═╡ 3751319d-d1b8-484a-8a70-ad46a6a65634
specnm_result.error === nothing ? nothing :
	Markdown.MD(Markdown.Admonition("danger", "Solve failed", [Markdown.Paragraph(specnm_result.error)]))

# ╔═╡ c692b5d4-e136-47f4-9cc6-432d4e5ef3fe
# Push solve status back to the SpecnmSolveInput widget so its "computing..." label
# resolves to "done"/"error: ..." once this cell (re)runs.
let
	payload = string("{\"ok\":", specnm_result.ok, ",\"error\":",
		specnm_result.error === nothing ? "null" : repr(replace(specnm_result.error, "\"" => "'")),
		"}")
	HTML("""<script>window.dispatchEvent(new CustomEvent('specnm-solve-done', {detail: $(repr(payload))}));</script>""")
end

# ╔═╡ 58041dd4-fd1c-4223-a00d-d4e98a7ef412
# Push the full mode catalogue (both families) + eigenfunctions + Earth-model profile
# to SpecnmBrowseView in one shot after a (re)solve -- ~85k floats for prem_ani,
# small enough that all of the widget's browsing below (family toggle, click-select,
# hover) is pure client-side array indexing, zero further Julia round-trips.
let
	num(x) = isfinite(x) ? string(round(Float64(x), digits=5)) : "0"
	jsonarr(v) = "[" * join(num.(v), ",") * "]"
	if specnm_result.ok
		em_radius = pyconvert(Float64, specnm_result.ray.radius / 1000.0)
		em_depth = em_radius .- pyconvert(Array, specnm_result.ray.r ./ 1000.0)
		em_rho = pyconvert(Array, specnm_result.ray.rho)
		em_vp = pyconvert(Array, specnm_result.ray.vp)
		em_vs = pyconvert(Array, specnm_result.ray.vs)
		payload = string(
			"{\"spheroidal\":", specnm_type_payload(specnm_result.ray, specnm_result.ray_out, specnm_result.ray_angular_orders, specnm_result.ray_frequencies, specnm_result.ray_overtones),
			",\"toroidal\":", specnm_type_payload(specnm_result.lov, specnm_result.lov_out, specnm_result.lov_angular_orders, specnm_result.lov_frequencies, specnm_result.lov_overtones),
			",\"earth_model\":{\"depth\":", jsonarr(em_depth), ",\"rho\":", jsonarr(em_rho), ",\"vp\":", jsonarr(em_vp), ",\"vs\":", jsonarr(em_vs), "}",
			",\"model_name\":\"", basename(model_fname), "\"",
			"}",
		)
		HTML("""<script>window.dispatchEvent(new CustomEvent('specnm-browse-data', {detail: $(repr(payload))}));</script>""")
	end
end

# ╔═╡ a9bbb105-3138-481d-9a73-0ec2b8303c36
begin
    """
        SpecnmBrowseView()

    Tier B -- the cheap, live-browsing widget. Entirely fed by a `CustomEvent` pushed
    once after [`SpecnmSolveInput`](@ref)'s solve completes (see the push cell above);
    every interaction below (spheroidal/toroidal toggle, hover, click-to-select a
    mode) is then pure client-side array indexing with zero further Julia calls.
    Three canvases: the mode spectrum (angular order `l` vs frequency, colored by
    overtone `n`), the selected mode's eigenfunction(s) vs depth, and the solved
    Earth model's density/velocity profile vs depth.
    """
    struct SpecnmBrowseView end

    function Base.show(io::IO, ::MIME"text/html", w::SpecnmBrowseView)
        write(io, """
<div id="sbvwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#sbvwidget){width:min(90vw,1400px)!important;margin-left:calc((100% - min(90vw,1400px))/2)!important}
    #sbvwidget{width:100%;box-sizing:border-box;color:#d1d5db;font:14px sans-serif}
    #sbvwidget .sbv-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #sbvwidget .sbv-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #sbvwidget .sbv-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #sbvwidget .sbv-workspace{display:flex;gap:12px;align-items:flex-start;justify-content:center;width:100%;flex-wrap:wrap}
    #sbvwidget canvas{background:#000;border:1px solid #374151;border-radius:6px;display:block;cursor:crosshair}
    #sbvwidget .sbv-controls{width:100%;margin-top:8px;display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:8px;font:14px sans-serif}
    #sbvwidget .sbv-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
    #sbvwidget .sbv-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #sbvwidget .sbv-actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    #sbvwidget label{color:#d1d5db}
  </style>
  <div class="sbv-title">
    <div class="sbv-title-desc">Every mode from the last Compute, browsed for free -- no re-solve needed.</div>
    <div class="sbv-title-hint">click a point in the spectrum to see its eigenfunction &middot; hover to identify a mode &middot; toggle spheroidal/toroidal</div>
  </div>
  <div class="sbv-workspace">
    <canvas id="sbvspectrum" width="440" height="440"></canvas>
    <canvas id="sbveigen" width="260" height="440"></canvas>
    <canvas id="sbvmodel" width="260" height="440"></canvas>
  </div>
  <div class="sbv-controls">
    <div class="sbv-control-group">
      <div class="sbv-control-title">Mode Family</div>
      <div class="sbv-actions">
        <label><input type="radio" name="sbv-type" id="sbv-type-s" checked> Spheroidal</label>
        <label><input type="radio" name="sbv-type" id="sbv-type-t"> Toroidal</label>
      </div>
    </div>
    <div class="sbv-control-group">
      <div class="sbv-control-title">Status</div>
      <div id="sbv-status" style="font-size:13px;color:#9ca3af">waiting for a Compute run</div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.9, par.clientWidth || window.innerWidth*0.9, 1400)
  const totalW = Math.max(760, availW)
  const SPW = Math.round(Math.min(totalW*0.45, 480)), SIDEW = Math.round(Math.min(totalW*0.25, 280))
  const H = Math.round(Math.min(window.innerHeight - 380, 480, 440))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)

  let data = null       // {spheroidal:{...}, toroidal:{...}, earth_model:{...}, model_name}
  let curType = 'spheroidal'
  let selIdx = -1, hoverIdx = -1

  const spCvs = par.querySelector('#sbvspectrum'), sctx = spCvs.getContext('2d')
  const eiCvs = par.querySelector('#sbveigen'), ectx = eiCvs.getContext('2d')
  const moCvs = par.querySelector('#sbvmodel'), mctx = moCvs.getContext('2d')
  function hidpi(cv, cx, w, h){ cv.width=Math.round(w*DPR); cv.height=Math.round(h*DPR); cv.style.width=w+'px'; cv.style.height=h+'px'; cx.setTransform(DPR,0,0,DPR,0,0) }
  hidpi(spCvs, sctx, SPW, H); hidpi(eiCvs, ectx, SIDEW, H); hidpi(moCvs, mctx, SIDEW, H)

  const PAD_L = 42, PAD_R = 12, PAD_T = 10, PAD_B = 26

  // Categorical palette cycling by overtone number, readable on black.
  const PALETTE = ['#38bdf8','#f97316','#4ade80','#facc15','#f472b6','#a78bfa','#fb7185','#2dd4bf','#fbbf24','#818cf8']
  function branchColor(n){ return PALETTE[((n%PALETTE.length)+PALETTE.length)%PALETTE.length] }

  function cur(){ return data ? data[curType] : null }

  function drawSpectrum(){
    sctx.clearRect(0,0,SPW,H)
    const d = cur()
    if(!d){
      sctx.fillStyle = '#6b7280'; sctx.font = '13px sans-serif'
      sctx.fillText('press Compute above to solve a model', 12, 20)
      return
    }
    const x0=PAD_L, x1=SPW-PAD_R, y0=PAD_T, y1=H-PAD_B
    const lmax = Math.max(...d.l), fmax = Math.max(...d.f)*1.05
    const X = l => x0 + (l/lmax)*(x1-x0)
    const Y = f => y1 - (f/fmax)*(y1-y0)

    sctx.strokeStyle = '#1f2937'; sctx.fillStyle = '#9ca3af'; sctx.font = '10px sans-serif'; sctx.textAlign = 'right'
    for(let i=0;i<=5;i++){ const f = fmax*i/5, y = Y(f); sctx.beginPath(); sctx.moveTo(x0,y); sctx.lineTo(x1,y); sctx.stroke(); sctx.fillText(f.toFixed(1), x0-6, y+3) }
    sctx.textAlign = 'center'
    for(let i=0;i<=5;i++){ const l = lmax*i/5, x = X(l); sctx.fillText(Math.round(l)+'', x, y1+14) }
    sctx.fillStyle = '#e5e7eb'; sctx.font = '12px sans-serif'
    sctx.fillText((curType==='spheroidal'?'Spheroidal':'Toroidal') + ' spectrum -- ' + (data.model_name||''), (x0+x1)/2, 16)
    sctx.font = '11px sans-serif'; sctx.fillStyle = '#9ca3af'
    sctx.fillText('angular degree l', (x0+x1)/2, H-4)
    sctx.save(); sctx.translate(12, (y0+y1)/2); sctx.rotate(-Math.PI/2); sctx.fillText('frequency (mHz)', 0, 0); sctx.restore()

    for(let i=0;i<d.l.length;i++){
      const isSel = i===selIdx, isHov = i===hoverIdx
      sctx.beginPath()
      sctx.arc(X(d.l[i]), Y(d.f[i]), isSel?5:(isHov?4:2.5), 0, 2*Math.PI)
      sctx.fillStyle = branchColor(d.n[i])
      sctx.globalAlpha = (selIdx===-1 || isSel || isHov) ? 0.9 : 0.35
      sctx.fill()
      sctx.globalAlpha = 1
      if(isSel){ sctx.strokeStyle = '#f5f3ef'; sctx.lineWidth = 1.4; sctx.stroke() }
    }

    if(hoverIdx>=0){
      const label = (curType==='spheroidal'?'S':'T') + ': l='+d.l[hoverIdx]+' n='+d.n[hoverIdx]+' f='+d.f[hoverIdx].toFixed(3)+'mHz'
      sctx.font='12px sans-serif'
      const tw = sctx.measureText(label).width
      const tx = Math.min(X(d.l[hoverIdx])+8, SPW-tw-14), ty = 34
      sctx.fillStyle = 'rgba(11,18,32,0.9)'; sctx.fillRect(tx-5, ty-13, tw+10, 18)
      sctx.strokeStyle = '#374151'; sctx.strokeRect(tx-5, ty-13, tw+10, 18)
      sctx.fillStyle = '#e5e7eb'; sctx.textAlign='left'; sctx.fillText(label, tx, ty)
    }
  }

  function nearestIndex(mx, my){
    const d = cur(); if(!d) return -1
    const x0=PAD_L, x1=SPW-PAD_R, y0=PAD_T, y1=H-PAD_B
    const lmax = Math.max(...d.l), fmax = Math.max(...d.f)*1.05
    const X = l => x0 + (l/lmax)*(x1-x0), Y = f => y1 - (f/fmax)*(y1-y0)
    let best=-1, bestD=10
    for(let i=0;i<d.l.length;i++){
      const dd = Math.hypot(mx-X(d.l[i]), my-Y(d.f[i]))
      if(dd<bestD){ bestD=dd; best=i }
    }
    return best
  }

  function drawEigen(){
    ectx.clearRect(0,0,SIDEW,H)
    const d = cur()
    const x0=PAD_L, x1=SIDEW-PAD_R, y0=PAD_T+14, y1=H-PAD_B
    if(!d || selIdx<0){
      ectx.fillStyle = '#6b7280'; ectx.font = '12px sans-serif'
      ectx.fillText('click a mode', 10, 18)
      return
    }
    const nr = d.nradial, depth = d.depth
    const depthMax = depth[depth.length-1]
    const Y = dep => y0 + (dep/depthMax)*(y1-y0)
    let vmax = 1e-9
    for(const comp of d.ef){ for(let i=0;i<nr;i++) vmax = Math.max(vmax, Math.abs(comp[selIdx*nr+i])) }
    const X = v => (x0+x1)/2 + (v/vmax)*(x1-x0)*0.48

    ectx.strokeStyle = '#1f2937'; ectx.fillStyle = '#9ca3af'; ectx.font='10px sans-serif'; ectx.textAlign='right'
    for(let i=0;i<=5;i++){ const dep=depthMax*i/5, y=Y(dep); ectx.beginPath(); ectx.moveTo(x0,y); ectx.lineTo(x1,y); ectx.stroke(); ectx.fillText(Math.round(dep)+'', x0-6, y+3) }
    ectx.strokeStyle = '#374151'; ectx.beginPath(); ectx.moveTo(X(0),y0); ectx.lineTo(X(0),y1); ectx.stroke()

    const colors = ['#38bdf8','#f97316','#4ade80']
    d.ef.forEach((comp, c) => {
      ectx.strokeStyle = colors[c%colors.length]; ectx.lineWidth = 1.8
      ectx.beginPath()
      for(let i=0;i<nr;i++){ const px=X(comp[selIdx*nr+i]), py=Y(depth[i]); i===0?ectx.moveTo(px,py):ectx.lineTo(px,py) }
      ectx.stroke()
    })
    ectx.textAlign='center'; ectx.fillStyle='#e5e7eb'; ectx.font='12px sans-serif'
    const nS = (curType==='spheroidal') ? 'S' : 'T'
    ectx.fillText(d.n[selIdx]+nS+d.l[selIdx]+'  f='+d.f[selIdx].toFixed(3)+'mHz', (x0+x1)/2, 16)
    ectx.font='11px sans-serif'
    d.labels.forEach((lab,c)=>{ ectx.fillStyle = colors[c%colors.length]; ectx.fillText(lab, x0+16+c*24, H-6) })
    ectx.save(); ectx.translate(12, (y0+y1)/2); ectx.rotate(-Math.PI/2); ectx.fillStyle='#9ca3af'; ectx.fillText('depth (km)', 0, 0); ectx.restore()
  }

  function drawModel(){
    mctx.clearRect(0,0,SIDEW,H)
    const em = data ? data.earth_model : null
    const x0=PAD_L, x1=SIDEW-PAD_R, y0=PAD_T+14, y1=H-PAD_B
    if(!em){
      mctx.fillStyle = '#6b7280'; mctx.font = '12px sans-serif'
      mctx.fillText('no model solved yet', 10, 18)
      return
    }
    const depthMax = em.depth[em.depth.length-1]
    const Y = dep => y0 + (dep/depthMax)*(y1-y0)
    const series = [['rho','#f97316','density (kg/m3)'], ['vp','#38bdf8','vp (m/s)'], ['vs','#4ade80','vs (m/s)']]

    mctx.strokeStyle = '#1f2937'; mctx.fillStyle = '#9ca3af'; mctx.font='10px sans-serif'; mctx.textAlign='right'
    for(let i=0;i<=5;i++){ const dep=depthMax*i/5, y=Y(dep); mctx.beginPath(); mctx.moveTo(x0,y); mctx.lineTo(x1,y); mctx.stroke(); mctx.fillText(Math.round(dep)+'', x0-6, y+3) }

    series.forEach(([key,color,label], si)=>{
      const vals = em[key]
      const vmax = Math.max(...vals)*1.05
      const X = v => x0 + (v/vmax)*(x1-x0)
      mctx.strokeStyle = color; mctx.lineWidth = 1.8
      mctx.beginPath()
      for(let i=0;i<vals.length;i++){ const px=X(vals[i]), py=Y(em.depth[i]); i===0?mctx.moveTo(px,py):mctx.lineTo(px,py) }
      mctx.stroke()
      mctx.fillStyle = color; mctx.font='11px sans-serif'; mctx.textAlign='left'
      mctx.fillText(label, x0+4, 16+si*12)
    })
    mctx.textAlign='center'; mctx.fillStyle='#e5e7eb'; mctx.font='12px sans-serif'
    mctx.fillText('Earth model', (x0+x1)/2, y0-2)
    mctx.save(); mctx.translate(12, (y0+y1)/2); mctx.rotate(-Math.PI/2); mctx.fillStyle='#9ca3af'; mctx.font='11px sans-serif'; mctx.fillText('depth (km)', 0, 0); mctx.restore()
  }

  function redraw(){ drawSpectrum(); drawEigen(); drawModel() }
  redraw()

  spCvs.addEventListener('mousemove', e=>{
    const i = nearestIndex(e.offsetX, e.offsetY)
    if(i !== hoverIdx){ hoverIdx = i; drawSpectrum() }
  })
  spCvs.addEventListener('mouseleave', ()=>{ if(hoverIdx!==-1){ hoverIdx=-1; drawSpectrum() } })
  spCvs.addEventListener('click', e=>{
    const i = nearestIndex(e.offsetX, e.offsetY)
    if(i>=0){ selIdx = i; drawSpectrum(); drawEigen() }
  })

  par.querySelector('#sbv-type-s').addEventListener('change', ()=>{ curType='spheroidal'; selIdx=-1; hoverIdx=-1; redraw() })
  par.querySelector('#sbv-type-t').addEventListener('change', ()=>{ curType='toroidal'; selIdx=-1; hoverIdx=-1; redraw() })

  window.addEventListener('specnm-browse-data', e=>{
    data = JSON.parse(e.detail)
    par.querySelector('#sbv-status').textContent = 'model: ' + (data.model_name||'') + '  (' + data.spheroidal.l.length + ' spheroidal, ' + data.toroidal.l.length + ' toroidal modes)'
    selIdx = -1; hoverIdx = -1
    redraw()
  })
</script>
""")
    end
end

# ╔═╡ bf3349da-ad56-4c4b-9b5c-c2e138abd1c0
SpecnmBrowseView()

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
CondaPkg = "~0.2.33"
FFTW = "~1.10.0"
PlutoPlotly = "~0.6.6"
PlutoUI = "~0.7.83"
PythonCall = "~0.9.31"
StatsBase = "~0.34.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "9d3d67d856891af255e148e59dc0d57e31573b21"

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
version = "1.3.0+1"

[[deps.CondaPkg]]
deps = ["JSON3", "Markdown", "MicroMamba", "Pidfile", "Pkg", "Preferences", "Scratch", "TOML", "pixi_jll"]
git-tree-sha1 = "bd491d55b97a036caae1d78729bdb70bf7dababc"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.33"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

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
version = "1.7.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

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
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

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
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "411eccfe8aba0814ffa0fdf4860913ed09c34975"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.3"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

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
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

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
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

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
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

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
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "6256ab3ee24ef079b3afa310593817e069925eeb"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.23"

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
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "6122f9423393a2294e26a4efdf44960c5f8acb70"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.78"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "982f3f017f08d31202574ef6bdcf8b3466430bea"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.31"

    [deps.PythonCall.extensions]
    CategoricalArraysExt = "CategoricalArrays"
    PyCallExt = "PyCall"

    [deps.PythonCall.weakdeps]
    CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
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

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

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
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "9297459be9e338e546f5c4bedb59b3b5674da7f1"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

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
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

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
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.micromamba_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "2ca2ac0b23a8e6b76752453e08428b3b4de28095"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "1.5.12+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.pixi_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "f349584316617063160a947a82638f7611a8ef0f"
uuid = "4d7b5844-a134-5dcd-ac86-c8f19cd51bed"
version = "0.41.3+0"
"""

# ╔═╡ Cell order:
# ╠═cd3c5ff7-7992-41f7-9a5e-29934a1355ca
# ╟─bdae265d-3a96-4ecc-a6dd-c166357e801c
# ╟─d9c71796-e79f-404d-89ef-54adfbd3335c
# ╠═bf3349da-ad56-4c4b-9b5c-c2e138abd1c0
# ╠═f2a1c3d4-7e6b-4a5c-9d8e-1b2c3d4e5f6a
# ╟─82cc9219-7633-41c6-91e1-17968904b2b6
# ╟─849d1a59-2c2f-4e44-b45b-266382255b1a
# ╟─2d7e4963-1cee-413c-8541-e8e77962c1fe
# ╠═b57d88b9-bdd5-447b-ba24-30a7b1ea1e0b
# ╟─921c2136-563a-4219-98f5-21cafa67e516
# ╟─fe54d35f-8017-4bc0-ac6c-7dee815d1477
# ╠═1a67d073-bf07-4b01-beea-6c11013965be
# ╟─89c81d3c-67c1-4cc7-b892-d64208d845f2
# ╟─9b0a9a80-e65f-4385-a377-372e408b19ad
# ╟─2afaeeec-6fe2-48db-81de-1fa4ccb0fc1b
# ╠═a751c6dc-414b-4ce8-b1b1-5aef8e23ab63
# ╠═a9019660-1f03-4132-bccf-09bdb1421ad9
# ╠═447b4e82-fe73-11ef-30b1-69824c8e3d24
# ╠═9f419394-528a-4bde-98a5-d62787c17fa8
# ╠═79406ee0-6026-4da8-a29d-245048c27e47
# ╠═acd7956a-855c-43f0-a353-7d5533f6aaf1
# ╠═335162fc-52c3-4c3d-bb88-d2e2f4fde37f
# ╠═58baaf3c-a5b3-4995-9f4e-c54946f7e798
# ╟─5852bfac-a8ff-405e-8ed1-6438c6827091
# ╠═d99b1935-4f90-4b82-809c-b6a801c37e0d
# ╠═1b22f312-9197-4044-ba9b-1f12789b88ff
# ╠═beddae25-17bc-48e9-8eef-41f21a08fb10
# ╠═69699338-d461-4e05-bbb4-a874ad4ad970
# ╠═55b7b9a8-caea-4ec6-b47d-937231eaaee8
# ╠═3751319d-d1b8-484a-8a70-ad46a6a65634
# ╟─c692b5d4-e136-47f4-9cc6-432d4e5ef3fe
# ╟─b5bfa4ac-f20a-4d1a-b870-fe4a3fa52da3
# ╠═a5491b1e-ce27-4fb5-82de-853899f75006
# ╠═0ecb0dea-e87e-4602-956c-9fc4a9375cb8
# ╠═7db650a6-5155-494f-9fc3-c5d664f98b32
# ╠═5a632ee0-e405-432f-a35e-e37defe34552
# ╠═902f56d9-041d-4b4a-af72-1adcf20a4db1
# ╠═ea6fe2a5-42f3-4c6f-b9e9-f4e4172463da
# ╠═73042ddf-790f-4a11-98b2-f9b6b9d29fe0
# ╠═521d6111-e4fb-4532-9acb-561a35aa5607
# ╠═74e0c240-dda8-44bd-bf0b-fb862a6ea52d
# ╠═5360665a-81ba-4cdb-8fdc-55829c6f4255
# ╠═aa28ca2f-662b-4789-a6fe-579f85d27ada
# ╠═ba1041f2-8173-4436-9413-ba8c8c85af20
# ╠═55968681-0b8c-4bb0-8ee9-1a9d81ee34e6
# ╠═580d562e-4c3f-4adb-94cf-65048ee7ff35
# ╠═7493f1d5-2017-416b-bc0c-183767bad68a
# ╠═3a325fab-9cf2-4cb8-a645-f0731e40b1f5
# ╠═a7d64a11-17c4-4aa0-a10f-e47d39b1eebb
# ╠═c168a471-f21f-4175-b7a0-5ae44d5d30ee
# ╠═39710dda-1054-491e-b9d5-7d8c4cb8e3a6
# ╟─58041dd4-fd1c-4223-a00d-d4e98a7ef412
# ╟─ecee2903-a103-498a-9aaf-b5d920f531f4
# ╠═71af98fe-12cb-42bf-b447-5afde0482e47
# ╟─a9bbb105-3138-481d-9a73-0ec2b8303c36
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
