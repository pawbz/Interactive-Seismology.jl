### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "My First Seismogram"
#> date = "2025-08-06"
#> layout = "layout.jlhtml"
#> description = "Pick a real earthquake and station, then download and filter the seismogram"
#> tags = ["eqdata"]

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
begin
    using PythonCall, PlutoPlotly
    using PlutoUI, FFTW, StatsBase
    using Colors
    using DelimitedFiles

    # Forces TableOfContents() below to depend on this cell -- a bare call to a
    # PlutoUI-exported function carries no argument for Pluto's static analysis to
    # detect, so on a fresh restart it can run before `using PlutoUI` has (same class
    # of ordering bug documented for the two widgets' @bind cells below).
    const _pkgs_ready = true
end

# ╔═╡ ca9f95ec-71d2-11f0-2e94-7da041469cdf
begin
    _pkgs_ready
    PlutoUI.TableOfContents(include_definitions=true)
end

# ╔═╡ 5a72389a-1277-4e18-a532-ce81a4b18c74
md"""
# My First Seismogram

A seismogram is what one station recorded during one earthquake -- which means *getting* one
is really three decisions: **which earthquake**, **which station**, and **how to look at the
recording** (how long a window, which frequencies). This notebook walks through exactly those
three steps, downloading real data from live earthquake and seismic-network catalogs at each
one.

!!! note "Data sources"
	Earthquake catalogs and station/waveform data live in different places these days:
	the earthquake catalog below comes from the **USGS** FDSN event service, and the
	station metadata + waveforms come from **EarthScope** (the network formerly known as
	the IRIS Data Management Center) -- see the Appendix's `event_client`/`waveform_client`
	split for why a single shared client no longer works for both.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000001
md"""
## Step 1: Pick an Earthquake

Not every earthquake makes a useful teaching example -- too small and the signal is buried in
noise at a global station, too close and there's no interesting wave propagation to see. Drag
a box on the map below to a region you're curious about, set a minimum magnitude, and hit
**Search**. Click any dot to select that earthquake -- dot size scales with magnitude.
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000003
md"""
**Selected earthquake**: $(selected_earthquake_details)
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000004
md"""
## Step 2: Pick a Station

A station only gives you a seismogram for this earthquake if it was actually **operating**
during it -- picking blind can land on a station with no data for that time window at all.
Drag a box over a region (GSN/`IU` stations only), hit **Search**: this queries for stations
that were online across the whole event window, not just ones that happen to be nearby.
Click one to select it. The star marks the epicenter, and once you pick a station the line
between them is the great-circle path the waves actually traveled -- the epicentral distance
is exactly what set `minradius`/`maxradius` in this notebook's earlier version, now read off
directly instead of being a search parameter.
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000006
md"""
**Selected station**: $(selected_station_details)
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000007
md"""
## Step 3: Download and Filter

The **total duration** below sets how much of the recording to pull -- long enough to catch
the phases you care about, without downloading more than needed. Once you have a station
picked, hit **Download Waveform**. The **frequency band** sliders don't re-download anything:
they just re-filter whatever's already on your machine, so drag them freely.
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000008
md"""
Total Duration (seconds): $(@bind total_duration Slider(500:10:5000, show_value=true, default=4000))

$(@bind download_click CounterButton("Download Waveform"))
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000009
md"""
Minimum Frequency (Hz): $(@bind freqmin Slider(0.01:0.01:0.1, show_value=true, default=0.05))

Maximum Frequency (Hz): $(@bind freqmax Slider(0.2:0.01:1.0, show_value=true, default=0.5))
"""

# ╔═╡ d2aba787-dae6-4da5-8f16-0cadfcc38371
md"""
## Appendix
"""

# ╔═╡ 6a1e0001-0000-4000-8000-00000000000c
md"### Package Imports"

# ╔═╡ 6a1e0001-0000-4000-8000-00000000000d
md"### FDSN Clients"

# ╔═╡ 6a1e0001-0000-4000-8000-00000000000e
obspy = pyimport("obspy")

# ╔═╡ 6a1e0001-0000-4000-8000-00000000000f
taup = pyimport("obspy.taup")

# ╔═╡ 6a1e0001-0000-4000-8000-000000000010
UTCDateTime = obspy.UTCDateTime

# ╔═╡ 6a1e0001-0000-4000-8000-000000000011
fdsn = pyimport("obspy.clients.fdsn")

# ╔═╡ 6a1e0001-0000-4000-8000-000000000012
begin
    # PythonCall's embedded Python does not inherit `SSL_CERT_FILE` the way a normal
    # shell-activated `python3` invocation does, so any `requests`-based network call
    # (which is what `fdsn.Client(...)`'s service-discovery step makes) fails
    # deterministically with `CERTIFICATE_VERIFY_FAILED: unable to get local issuer
    # certificate` -- confirmed directly: the identical `Client("USGS")` call succeeds
    # every time from a plain `python3` in this same CondaPkg environment, and fails
    # every time through PythonCall, until this is set. Pointing both at `certifi`'s own
    # bundle (already a transitive obspy dependency, so always present) fixes it.
    certifi = pyimport("certifi")
    cafile = pyconvert(String, certifi.where())
    ENV["SSL_CERT_FILE"] = cafile
    ENV["REQUESTS_CA_BUNDLE"] = cafile
end

# ╔═╡ 0b12e2c0-9183-11f1-898e-dbccf4ea72db
begin
    # EarthScope (the network formerly known as the IRIS Data Management Center) retired
    # its FDSN `event` web service on 2026-06-01 -- a single shared `Client("IRIS")` can no
    # longer serve earthquake catalogs at all, only station/waveform data. USGS still runs
    # a full `event` service, so the fix is two clients instead of one: `event_client` for
    # the catalog, `waveform_client` for everything else. Confirmed directly against both
    # live services.
    cafile   # bare reference: forces this cell to depend on the SSL_CERT_FILE fix above
    event_client = fdsn.Client("USGS")
    waveform_client = fdsn.Client("IRIS")
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000013
taup_model = taup.TauPyModel(model="iasp91")

# ╔═╡ 6a1e0001-0000-4000-8000-000000000014
md"### Layer 1: Flat-Map Projection & Coastlines"

# ╔═╡ 6a1e0001-0000-4000-8000-000000000015
begin
    # Natural Earth 1:110m coastline (public domain), the same vendored dataset used by
    # `earth-internal-structure.jl`'s globe widget -- drawn here as a flat equirectangular
    # reference layer so a viewer can place a dragged search box geographically.
    coast_raw, _coast_header = readdlm(joinpath(@__DIR__, "data", "coastlines_110m.csv"), ',';
        comments=true, comment_char='#', header=true)
    coast_line_id = Int.(coast_raw[:, 1])
    coast_lon = Float64.(coast_raw[:, 2])
    coast_lat = Float64.(coast_raw[:, 3])
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000016
"""
    coast_js_literal(line_id, lon, lat)

Groups the flat `(line_id, lon, lat)` coastline arrays into one `[[lon,lat],...]` array per
polyline and renders that as a `[[[lon,lat],...],...]` JS literal, so a widget's JS side can
iterate "one polyline at a time" without redoing the grouping. Shared by both map widgets
below (identical to `earth-internal-structure.jl`'s helper of the same name).
"""
function coast_js_literal(line_id, lon, lat)
    segs = Dict{Int,Vector{Tuple{Float64,Float64}}}()
    for i in eachindex(line_id)
        push!(get!(() -> Tuple{Float64,Float64}[], segs, line_id[i]), (lon[i], lat[i]))
    end
    ids = sort(collect(keys(segs)))
    "[" * join(["[" * join(["[$(lo),$(la)]" for (lo, la) in segs[id]], ",") * "]" for id in ids], ",") * "]"
en

# ╔═╡ 6a1e0001-0000-4000-8000-000000000017
md"### Layer 2: Earthquake Search"

# ╔═╡ 6a1e0001-0000-4000-8000-000000000018
begin
    """
        EarthquakeExplorerInput(; latmin=-60.0, latmax=60.0, lonmin=-180.0, lonmax=180.0, minmag=6.5)

    Drag a rectangle directly on the map to set the lat/lon search box (mousedown a corner,
    drag to the opposite corner) -- the box itself is the only location control, there are no
    separate lat/lon sliders. Set a minimum magnitude, press **Search** to query the USGS
    earthquake catalog for that box/magnitude (network call, gated behind the button -- never
    fires just from dragging or adjusting the slider), then click the nearest dot to select
    that earthquake. Dot size scales with magnitude. Self-contained: the fetched catalog is
    rendered and hit-tested entirely client-side once pushed in via the `eqx-results` event;
    Julia only ever sees the two explicit actions (`search`, `select`) dispatched through the
    bound value.
    """
    struct EarthquakeExplorerInput
        latmin::Float64
        latmax::Float64
        lonmin::Float64
        lonmax::Float64
        minmag::Float64
    end

    EarthquakeExplorerInput(; latmin=-60.0, latmax=60.0, lonmin=-180.0, lonmax=180.0, minmag=6.5) =
        EarthquakeExplorerInput(Float64(latmin), Float64(latmax), Float64(lonmin), Float64(lonmax), Float64(minmag))

    Base.get(w::EarthquakeExplorerInput) = Dict{String,Any}(
        "action" => "init", "latmin" => w.latmin, "latmax" => w.latmax,
        "lonmin" => w.lonmin, "lonmax" => w.lonmax, "minmag" => w.minmag,
    )

    function Base.show(io::IO, ::MIME"text/html", w::EarthquakeExplorerInput)
        write(io, """
<div id="eqxwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#eqxwidget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #eqxwidget .eqx-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #eqxwidget .eqx-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #eqxwidget .eqx-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #eqxwidget canvas{cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block;max-width:100%}
    #eqxwidget .eqx-controls{display:flex;gap:10px;flex-wrap:wrap;width:100%;margin-top:12px}
    #eqxwidget .eqx-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px;flex:1 1 240px;min-width:240px}
    #eqxwidget .eqx-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #eqxwidget .eqx-control-row{display:grid;grid-template-columns:minmax(90px,120px) minmax(70px,1fr) minmax(40px,56px);gap:6px;align-items:center;margin:5px 0}
    #eqxwidget .eqx-control-row input[type=range]{width:100%;min-width:0;vertical-align:middle}
    #eqxwidget .eqx-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums;font-size:13px}
    #eqxwidget .eqx-readout{font-size:13px;line-height:1.7}
    #eqxwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px;cursor:pointer}
    #eqxwidget button:disabled{opacity:0.5;cursor:default}
  </style>
  <div class="eqx-title">
    <div class="eqx-title-desc">Every earthquake is a source location and an origin time -- pick one to see how its waves travel.</div>
    <div class="eqx-title-hint">drag a box on the map &middot; set min magnitude &middot; Search, then click a dot to select</div>
  </div>
  <canvas id="eqxMap" width="900" height="450"></canvas>
  <div class="eqx-controls">
    <div class="eqx-control-group">
      <div class="eqx-control-title">Search</div>
      <label class="eqx-control-row"><span>min mag</span><input type="range" id="eqx-minmag" min="4" max="8.5" step="0.1" value="$(w.minmag)"><span id="eqx-minmag-v" class="eqx-value">$(w.minmag)</span></label>
      <div style="display:flex;gap:10px;align-items:center;margin-top:6px">
        <button id="eqx-search" type="button">Search</button>
        <span id="eqx-status" style="font-size:13px">drag a box above, then Search</span>
      </div>
    </div>
    <div class="eqx-control-group" style="flex:1 1 280px">
      <div class="eqx-control-title">Selected Earthquake</div>
      <div class="eqx-readout" id="eqx-readout">none yet</div>
    </div>
    <div class="eqx-control-group" style="flex:1 1 240px">
      <div class="eqx-control-title">Legend</div>
      <div class="eqx-readout"><span style="color:#ef4444">&#9679;</span> earthquake (size = magnitude) &middot; <span style="color:#facc15">&#9679;</span> selected &middot; <span style="color:#facc15">box</span> = search region</div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1400)
  const W = Math.round(Math.max(500, availW))
  const H = Math.round(W/2)
  const DPR = Math.min(window.devicePixelRatio || 1, 2)
  const cvs = par.querySelector('#eqxMap'), ctx = cvs.getContext('2d')
  cvs.width = Math.round(W*DPR); cvs.height = Math.round(H*DPR)
  cvs.style.width = W+'px'; cvs.style.height = H+'px'
  ctx.setTransform(DPR,0,0,DPR,0,0)

  const COAST = $(coast_js_literal(coast_line_id, coast_lon, coast_lat))

  let latmin=$(w.latmin), latmax=$(w.latmax), lonmin=$(w.lonmin), lonmax=$(w.lonmax), minmag=$(w.minmag)
  let eqData = null   // filled in by the 'eqx-results' push from Julia
  let selectedIdx = -1
  let dragStart = null, dragMoved = false

  function lonlatToXY(lon, lat){ return [(lon+180)/360*W, (90-lat)/180*H] }
  function xyToLonLat(x, y){ return [x/W*360-180, 90-y/H*180] }

  function drawCoast(){
    ctx.strokeStyle = 'rgba(107,114,128,0.55)'; ctx.lineWidth = 1
    for(const line of COAST){
      ctx.beginPath()
      let prevLon = null, started = false
      for(const [lon,lat] of line){
        const [x,y] = lonlatToXY(lon,lat)
        if(prevLon !== null && Math.abs(lon-prevLon) > 180){ ctx.moveTo(x,y) }
        else if(!started){ ctx.moveTo(x,y) } else { ctx.lineTo(x,y) }
        started = true; prevLon = lon
      }
      ctx.stroke()
    }
  }

  function drawBox(){
    const [x0,y0] = lonlatToXY(lonmin, latmax)
    const [x1,y1] = lonlatToXY(lonmax, latmin)
    ctx.strokeStyle = '#facc15'; ctx.lineWidth = 2
    ctx.strokeRect(Math.min(x0,x1), Math.min(y0,y1), Math.abs(x1-x0), Math.abs(y1-y0))
    ctx.fillStyle = 'rgba(250,204,21,0.08)'
    ctx.fillRect(Math.min(x0,x1), Math.min(y0,y1), Math.abs(x1-x0), Math.abs(y1-y0))
  }

  function drawEvents(){
    if(!eqData) return
    for(let i=0;i<eqData.length;i++){
      const [lon,lat,mag] = eqData[i]
      const [x,y] = lonlatToXY(lon,lat)
      const r = Math.max(3, Math.min(13, 2+mag*1.4))
      ctx.beginPath(); ctx.arc(x,y,r,0,2*Math.PI)
      ctx.fillStyle = i===selectedIdx ? '#facc15' : 'rgba(239,68,68,0.75)'
      ctx.fill()
      if(i===selectedIdx){ ctx.strokeStyle='#0a0f18'; ctx.lineWidth=1.5; ctx.stroke() }
    }
  }

  function redraw(){
    ctx.clearRect(0,0,W,H)
    ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
    drawCoast()
    drawEvents()
    drawBox()
  }
  redraw()

  function fmt(v){ return v.toFixed(2) }
  function updateBoxReadout(){
    // status line already shows search feedback; box bounds shown inline near legend instead
  }

  const searchBtn = par.querySelector('#eqx-search')
  const statusEl = par.querySelector('#eqx-status')
  searchBtn.addEventListener('click', ()=>{
    searchBtn.disabled = true
    statusEl.textContent = 'searching…'
    par.value = {action:'search', latmin, latmax, lonmin, lonmax, minmag}
    par.dispatchEvent(new CustomEvent('input'))
  })

  window.addEventListener('eqx-results', e=>{
    const d = JSON.parse(e.detail)
    searchBtn.disabled = false
    if(!d.ok){ statusEl.textContent = 'error: ' + d.error; eqData = []; redraw(); return }
    eqData = d.events
    selectedIdx = -1
    statusEl.textContent = d.count ? (d.count+' earthquakes found — click one') : 'no earthquakes found, try widening the box or lowering magnitude'
    redraw()
  })

  function trySelect(mx, my){
    if(!eqData || eqData.length===0) return
    let best=-1, bestD=Infinity
    for(let i=0;i<eqData.length;i++){
      const [lon,lat] = eqData[i]
      const [x,y] = lonlatToXY(lon,lat)
      const d = Math.hypot(x-mx,y-my)
      if(d<bestD){ bestD=d; best=i }
    }
    if(bestD < 16){
      selectedIdx = best
      const [lon,lat,mag,depth,time] = eqData[best]
      par.querySelector('#eqx-readout').innerHTML =
        '<b>M'+mag.toFixed(1)+'</b> at '+time+'<br>('+fmt(lat)+'°, '+fmt(lon)+'°), depth '+depth.toFixed(0)+' km'
      par.value = {action:'select', lat, lon, depth, mag, time}
      par.dispatchEvent(new CustomEvent('input'))
      redraw()
    }
  }

  cvs.addEventListener('mousedown', e=>{ dragStart=[e.offsetX,e.offsetY]; dragMoved=false })
  cvs.addEventListener('mousemove', e=>{
    if(!dragStart) return
    const dx=e.offsetX-dragStart[0], dy=e.offsetY-dragStart[1]
    if(Math.hypot(dx,dy) > 6) dragMoved = true
    if(dragMoved){
      const [lonA,latA] = xyToLonLat(dragStart[0], dragStart[1])
      const [lonB,latB] = xyToLonLat(e.offsetX, e.offsetY)
      latmin=Math.min(latA,latB); latmax=Math.max(latA,latB)
      lonmin=Math.min(lonA,lonB); lonmax=Math.max(lonA,lonB)
      redraw()
    }
  })
  cvs.addEventListener('mouseup', e=>{
    if(!dragStart) return
    if(!dragMoved) trySelect(e.offsetX, e.offsetY)
    dragStart=null; dragMoved=false
  })
  window.addEventListener('mouseup', ()=>{ dragStart=null; dragMoved=false })

  par.querySelector('#eqx-minmag').addEventListener('input', e=>{
    minmag = parseFloat(e.target.value)
    par.querySelector('#eqx-minmag-v').textContent = minmag.toFixed(1)
  })
</script>
""")
    end

    const _eqx_ready = true
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000002
begin
    _eqx_ready
    @bind eqx EarthquakeExplorerInput()
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000019
let
    if eqx isa AbstractDict && eqx["action"] == "search"
        payload = try
            cat = event_client.get_events(
                starttime=UTCDateTime("1990-01-01"), endtime=UTCDateTime("2030-01-01"),
                minmagnitude=eqx["minmag"],
                minlatitude=eqx["latmin"], maxlatitude=eqx["latmax"],
                minlongitude=eqx["lonmin"], maxlongitude=eqx["lonmax"],
            )
            events = [
                (
                    lon=pyconvert(Float64, ev.origins[0].longitude),
                    lat=pyconvert(Float64, ev.origins[0].latitude),
                    mag=pyconvert(Float64, ev.magnitudes[0].mag),
                    depth=pyconvert(Float64, ev.origins[0].depth) / 1000,
                    time=pyconvert(String, string(ev.origins[0].time)),
                ) for ev in cat
            ]
            eventsjson = "[" * join(["[$(e.lon),$(e.lat),$(e.mag),$(e.depth),\"$(e.time)\"]" for e in events], ",") * "]"
            "{\"ok\":true,\"count\":$(length(events)),\"events\":$eventsjson}"
        catch e
            errmsg = replace(sprint(showerror, e), "\"" => "'", "\n" => " ")
            "{\"ok\":false,\"error\":\"$errmsg\"}"
        end
        HTML("""<script>window.dispatchEvent(new CustomEvent('eqx-results', {detail: $(repr(payload))}));</script>""")
    end
end

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001a
selected_earthquake_details = if eqx isa AbstractDict && eqx["action"] == "select"
    (time=eqx["time"], lat=eqx["lat"], lon=eqx["lon"], depth=eqx["depth"], mag=eqx["mag"])
else
    "none yet -- search and click an earthquake above"
en

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001b
md"### Layer 3: Station Search"

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001c
begin
    """
        StationExplorerInput(; latmin=-60.0, latmax=60.0, lonmin=-180.0, lonmax=180.0)

    Same box-drag-then-Search shape as `EarthquakeExplorerInput`, for GSN (`IU` network)
    stations instead of earthquakes -- no magnitude control, since stations don't have one.
    The epicenter of the currently-selected earthquake (read directly from the notebook-level
    `selected_earthquake_details`, baked into this widget's HTML at construction time the same
    way `coast_js_literal(...)` bakes in the coastline data) is drawn as a star, so picking a
    *new* earthquake naturally produces a *fresh* widget instance -- any station search/
    selection from a previous earthquake is correctly discarded rather than silently stale.
    Once a station is clicked, the great-circle path to the epicenter is drawn and its length
    reported live (client-side spherical trig, for instant feedback -- the authoritative
    `epi_dist_deg` used for the TauP call is computed separately in Julia via
    `obspy.geodetics`).
    """
    struct StationExplorerInput
        latmin::Float64
        latmax::Float64
        lonmin::Float64
        lonmax::Float64
    end

    StationExplorerInput(; latmin=-60.0, latmax=60.0, lonmin=-180.0, lonmax=180.0) =
        StationExplorerInput(Float64(latmin), Float64(latmax), Float64(lonmin), Float64(lonmax))

    Base.get(w::StationExplorerInput) = Dict{String,Any}(
        "action" => "init", "latmin" => w.latmin, "latmax" => w.latmax,
        "lonmin" => w.lonmin, "lonmax" => w.lonmax,
    )

    function Base.show(io::IO, ::MIME"text/html", w::StationExplorerInput)
        epi_lat_js = selected_earthquake_details isa NamedTuple ? string(selected_earthquake_details.lat) : "null"
        epi_lon_js = selected_earthquake_details isa NamedTuple ? string(selected_earthquake_details.lon) : "null"
        write(io, """
<div id="stxwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#stxwidget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #stxwidget .stx-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #stxwidget .stx-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #stxwidget .stx-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #stxwidget canvas{cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block;max-width:100%}
    #stxwidget .stx-controls{display:flex;gap:10px;flex-wrap:wrap;width:100%;margin-top:12px}
    #stxwidget .stx-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px;flex:1 1 240px;min-width:240px}
    #stxwidget .stx-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #stxwidget .stx-readout{font-size:13px;line-height:1.7}
    #stxwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px;cursor:pointer}
    #stxwidget button:disabled{opacity:0.5;cursor:default}
  </style>
  <div class="stx-title">
    <div class="stx-title-desc">A station only gives you data for this earthquake if it was actually recording during it.</div>
    <div class="stx-title-hint">drag a box on the map &middot; Search (only returns stations online for this event) &middot; click one</div>
  </div>
  <canvas id="stxMap" width="900" height="450"></canvas>
  <div class="stx-controls">
    <div class="stx-control-group">
      <div class="stx-control-title">Search</div>
      <div style="display:flex;gap:10px;align-items:center">
        <button id="stx-search" type="button">Search</button>
        <span id="stx-status" style="font-size:13px">drag a box above, then Search</span>
      </div>
    </div>
    <div class="stx-control-group" style="flex:1 1 280px">
      <div class="stx-control-title">Selected Station</div>
      <div class="stx-readout" id="stx-readout">none yet</div>
    </div>
    <div class="stx-control-group" style="flex:1 1 240px">
      <div class="stx-control-title">Legend</div>
      <div class="stx-readout"><span style="color:#facc15">&#9733;</span> epicenter &middot; <span style="color:#3b82f6">&#9650;</span> station &middot; <span style="color:#22d3ee">line</span> = great-circle path</div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1400)
  const W = Math.round(Math.max(500, availW))
  const H = Math.round(W/2)
  const DPR = Math.min(window.devicePixelRatio || 1, 2)
  const cvs = par.querySelector('#stxMap'), ctx = cvs.getContext('2d')
  cvs.width = Math.round(W*DPR); cvs.height = Math.round(H*DPR)
  cvs.style.width = W+'px'; cvs.style.height = H+'px'
  ctx.setTransform(DPR,0,0,DPR,0,0)

  const COAST = $(coast_js_literal(coast_line_id, coast_lon, coast_lat))
  const epiLat = $(epi_lat_js), epiLon = $(epi_lon_js)

  let latmin=$(w.latmin), latmax=$(w.latmax), lonmin=$(w.lonmin), lonmax=$(w.lonmax)
  let stData = null
  let selectedIdx = -1
  let dragStart = null, dragMoved = false

  function lonlatToXY(lon, lat){ return [(lon+180)/360*W, (90-lat)/180*H] }
  function xyToLonLat(x, y){ return [x/W*360-180, 90-y/H*180] }
  function llToXYZ(lon, lat){
    const th=(90-lat)*Math.PI/180, ph=lon*Math.PI/180
    return [Math.sin(th)*Math.cos(ph), Math.sin(th)*Math.sin(ph), Math.cos(th)]
  }
  function xyzToLL(v){
    const n=Math.hypot(v[0],v[1],v[2])
    const th=Math.acos(Math.max(-1,Math.min(1,v[2]/n)))
    return [Math.atan2(v[1],v[0])*180/Math.PI, 90-th*180/Math.PI]
  }
  function slerp(a,b,f){
    const dot=a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    const ang=Math.acos(Math.max(-1,Math.min(1,dot)))
    if(ang<1e-9) return a
    const s=Math.sin(ang), w1=Math.sin((1-f)*ang)/s, w2=Math.sin(f*ang)/s
    return [w1*a[0]+w2*b[0], w1*a[1]+w2*b[1], w1*a[2]+w2*b[2]]
  }

  function drawCoast(){
    ctx.strokeStyle = 'rgba(107,114,128,0.55)'; ctx.lineWidth = 1
    for(const line of COAST){
      ctx.beginPath()
      let prevLon=null, started=false
      for(const [lon,lat] of line){
        const [x,y] = lonlatToXY(lon,lat)
        if(prevLon!==null && Math.abs(lon-prevLon)>180){ ctx.moveTo(x,y) }
        else if(!started){ ctx.moveTo(x,y) } else { ctx.lineTo(x,y) }
        started=true; prevLon=lon
      }
      ctx.stroke()
    }
  }

  function drawBox(){
    const [x0,y0] = lonlatToXY(lonmin, latmax)
    const [x1,y1] = lonlatToXY(lonmax, latmin)
    ctx.strokeStyle = '#facc15'; ctx.lineWidth = 2
    ctx.strokeRect(Math.min(x0,x1), Math.min(y0,y1), Math.abs(x1-x0), Math.abs(y1-y0))
    ctx.fillStyle = 'rgba(250,204,21,0.08)'
    ctx.fillRect(Math.min(x0,x1), Math.min(y0,y1), Math.abs(x1-x0), Math.abs(y1-y0))
  }

  function drawEpicenter(){
    if(epiLat===null) return
    const [x,y] = lonlatToXY(epiLon, epiLat)
    ctx.fillStyle = '#facc15'
    ctx.beginPath()
    for(let i=0;i<5;i++){
      const a = -Math.PI/2 + i*2*Math.PI/5
      const a2 = a + Math.PI/5
      const p1 = [x+Math.cos(a)*9, y+Math.sin(a)*9]
      const p2 = [x+Math.cos(a2)*3.5, y+Math.sin(a2)*3.5]
      i===0 ? ctx.moveTo(p1[0],p1[1]) : ctx.lineTo(p1[0],p1[1])
      ctx.lineTo(p2[0],p2[1])
    }
    ctx.closePath(); ctx.fill()
  }

  function drawStations(){
    if(!stData) return
    for(let i=0;i<stData.length;i++){
      const [code,lat,lon] = stData[i]
      const [x,y] = lonlatToXY(lon,lat)
      ctx.beginPath()
      ctx.moveTo(x, y-6); ctx.lineTo(x-6, y+5); ctx.lineTo(x+6, y+5); ctx.closePath()
      ctx.fillStyle = i===selectedIdx ? '#facc15' : 'rgba(59,130,246,0.8)'
      ctx.fill()
      if(i===selectedIdx){ ctx.strokeStyle='#0a0f18'; ctx.lineWidth=1.5; ctx.stroke() }
    }
  }

  function drawGreatCircle(lon2,lat2){
    if(epiLat===null) return
    const A = llToXYZ(epiLon,epiLat), B = llToXYZ(lon2,lat2)
    ctx.strokeStyle = '#22d3ee'; ctx.lineWidth = 2
    ctx.beginPath()
    let prevLon=null, started=false
    const N=64
    for(let i=0;i<=N;i++){
      const [lon,lat] = xyzToLL(slerp(A,B,i/N))
      const [x,y] = lonlatToXY(lon,lat)
      if(prevLon!==null && Math.abs(lon-prevLon)>180){ ctx.moveTo(x,y) }
      else if(!started){ ctx.moveTo(x,y) } else { ctx.lineTo(x,y) }
      started=true; prevLon=lon
    }
    ctx.stroke()
  }

  function redraw(){
    ctx.clearRect(0,0,W,H)
    ctx.fillStyle = '#000'; ctx.fillRect(0,0,W,H)
    drawCoast()
    if(selectedIdx>=0){ const [,lat,lon] = stData[selectedIdx]; drawGreatCircle(lon,lat) }
    drawEpicenter()
    drawStations()
    drawBox()
  }
  redraw()

  const searchBtn = par.querySelector('#stx-search')
  const statusEl = par.querySelector('#stx-status')
  searchBtn.addEventListener('click', ()=>{
    searchBtn.disabled = true
    statusEl.textContent = 'searching…'
    par.value = {action:'search', latmin, latmax, lonmin, lonmax}
    par.dispatchEvent(new CustomEvent('input'))
  })

  window.addEventListener('stx-results', e=>{
    const d = JSON.parse(e.detail)
    searchBtn.disabled = false
    if(!d.ok){ statusEl.textContent = 'error: ' + d.error; stData = []; redraw(); return }
    stData = d.stations
    selectedIdx = -1
    statusEl.textContent = d.count ? (d.count+' stations found — click one') : 'no stations were online for this event in that box'
    redraw()
  })

  function trySelect(mx, my){
    if(!stData || stData.length===0) return
    let best=-1, bestD=Infinity
    for(let i=0;i<stData.length;i++){
      const [,lat,lon] = stData[i]
      const [x,y] = lonlatToXY(lon,lat)
      const d = Math.hypot(x-mx,y-my)
      if(d<bestD){ bestD=d; best=i }
    }
    if(bestD < 16){
      selectedIdx = best
      const [code,lat,lon,elevation] = stData[best]
      let distTxt = ''
      if(epiLat!==null){
        const A=llToXYZ(epiLon,epiLat), B=llToXYZ(lon,lat)
        const dot=A[0]*B[0]+A[1]*B[1]+A[2]*B[2]
        const distDeg = Math.acos(Math.max(-1,Math.min(1,dot)))*180/Math.PI
        distTxt = '<br>'+distDeg.toFixed(1)+'&deg; from epicenter'
      }
      par.querySelector('#stx-readout').innerHTML =
        '<b>'+code+'</b> ('+lat.toFixed(2)+'°, '+lon.toFixed(2)+'°), '+elevation.toFixed(0)+' m elev'+distTxt
      par.value = {action:'select', code, lat, lon, elevation}
      par.dispatchEvent(new CustomEvent('input'))
      redraw()
    }
  }

  cvs.addEventListener('mousedown', e=>{ dragStart=[e.offsetX,e.offsetY]; dragMoved=false })
  cvs.addEventListener('mousemove', e=>{
    if(!dragStart) return
    const dx=e.offsetX-dragStart[0], dy=e.offsetY-dragStart[1]
    if(Math.hypot(dx,dy) > 6) dragMoved = true
    if(dragMoved){
      const [lonA,latA] = xyToLonLat(dragStart[0], dragStart[1])
      const [lonB,latB] = xyToLonLat(e.offsetX, e.offsetY)
      latmin=Math.min(latA,latB); latmax=Math.max(latA,latB)
      lonmin=Math.min(lonA,lonB); lonmax=Math.max(lonA,lonB)
      redraw()
    }
  })
  cvs.addEventListener('mouseup', e=>{
    if(!dragStart) return
    if(!dragMoved) trySelect(e.offsetX, e.offsetY)
    dragStart=null; dragMoved=false
  })
  window.addEventListener('mouseup', ()=>{ dragStart=null; dragMoved=false })
</script>
""")
    end

    const _stx_ready = true
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000005
begin
    _stx_ready
    @bind stx StationExplorerInput()
end

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001d
let
    if stx isa AbstractDict && stx["action"] == "search"
        payload = if !(selected_earthquake_details isa NamedTuple)
            "{\"ok\":false,\"error\":\"Select an earthquake first (Step 1).\"}"
        else
            try
                eq_time_for_search = UTCDateTime(selected_earthquake_details.time)
                sta = waveform_client.get_stations(network="IU",
                    minlatitude=stx["latmin"], maxlatitude=stx["latmax"],
                    minlongitude=stx["lonmin"], maxlongitude=stx["lonmax"],
                    starttime=eq_time_for_search, endtime=eq_time_for_search + 3600)
                stations = [
                    (code=pyconvert(String, s.code), lat=pyconvert(Float64, s.latitude),
                        lon=pyconvert(Float64, s.longitude), elevation=pyconvert(Float64, s.elevation))
                    for net in sta for s in net.stations
                ]
                stationsjson = "[" * join(["[\"$(s.code)\",$(s.lat),$(s.lon),$(s.elevation)]" for s in stations], ",") * "]"
                "{\"ok\":true,\"count\":$(length(stations)),\"stations\":$stationsjson}"
            catch e
                errmsg = replace(sprint(showerror, e), "\"" => "'", "\n" => " ")
                "{\"ok\":false,\"error\":\"$errmsg\"}"
            end
        end
        HTML("""<script>window.dispatchEvent(new CustomEvent('stx-results', {detail: $(repr(payload))}));</script>""")
    end
end

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001e
selected_station_details = if stx isa AbstractDict && stx["action"] == "select"
    (code=stx["code"], lat=stx["lat"], lon=stx["lon"], elevation=stx["elevation"])
else
    "none yet -- search and click a station above"
en

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001f
md"### Layer 4: Waveform Download & Phase Arrivals"

# ╔═╡ 6a1e0001-0000-4000-8000-000000000020
"""
Gated behind the **Download Waveform** button (via the bare `download_click` reference, the
same "reference a `CounterButton`'s value to gate a `let` block" trick the notebook's own
`nsample` cell used previously) -- dragging the duration slider or picking a different channel
never re-triggers a download, only pressing the button does. Wrapped in try/catch because a
station having no data for a given window is a *normal*, common outcome (confirmed live: e.g.
`IU.AFI` returns `FDSNNoDataException`/HTTP 204 for plenty of real event windows), not a bug --
it should surface as a friendly message, not crash every cell downstream.
"""
download_result = let
    download_click
    if !(selected_earthquake_details isa NamedTuple) || !(selected_station_details isa NamedTuple)
        (data=nothing, error=nothing)
    else
        try
            eq_time = UTCDateTime(selected_earthquake_details.time)
            tr = waveform_client.get_waveforms("IU", selected_station_details.code, "*", "BH?",
                attach_response=true, starttime=eq_time, endtime=eq_time + total_duration)
            tr = tr.normalize().detrend()
            (data=tr, error=nothing)
        catch e
            (data=nothing, error=sprint(showerror, e))
        end
    end
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000021
raw_traces = download_result.data

# ╔═╡ 6a1e0001-0000-4000-8000-000000000022
download_result.error === nothing ? md"" : md"""
!!! danger "Download failed"
	$(download_result.error)

	Try a different station -- not every station in the box was necessarily recording on
	every channel for the full event window.
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000023
"""
Cheap, local, and *live*: unlike the download above, re-filtering an already-downloaded trace
needs no network access, so `freqmin`/`freqmax` are plain reactive sliders -- dragging them
never re-triggers `download_result`.
"""
filtered_traces = raw_traces === nothing ? nothing :
    raw_traces.copy().filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=true)

# ╔═╡ 6a1e0001-0000-4000-8000-000000000024
channel_options = filtered_traces === nothing ? [1 => "—"] :
    [i => pyconvert(String, filtered_traces[i-1].stats.channel) for i in 1:length(filtered_traces)]

# ╔═╡ 6a1e0001-0000-4000-8000-00000000000b
begin
    channel_options   # bare reference: forces detection of this cell's dependency on
    # `channel_options`, which otherwise sits only inside `@bind`'s nested Select(...)
    # argument -- the same "@bind can't see nested-argument dependencies" gotcha as
    # the widgets below, just via a plain PlutoUI element instead of a custom struct.
    md"""
    Channel: $(@bind selected_channel_index Select(channel_options, default=1))
    """
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000025
begin
    if selected_earthquake_details isa NamedTuple && selected_station_details isa NamedTuple
        distance_m, az, baz = obspy.geodetics.gps2dist_azimuth(
            selected_earthquake_details.lat, selected_earthquake_details.lon,
            selected_station_details.lat, selected_station_details.lon)
        epi_dist_deg = pyconvert(Float64, obspy.geodetics.kilometer2degrees(distance_m / 1000.0))
        arrivals = taup_model.get_ray_paths(source_depth_in_km=selected_earthquake_details.depth,
            distance_in_degree=epi_dist_deg)
    else
        epi_dist_deg = NaN
        arrivals = []
    end
end

# ╔═╡ 6a1e0001-0000-4000-8000-00000000000a
let
    if filtered_traces === nothing
        md"Pick an earthquake and a station above, then press **Download Waveform**."
    else
        trace = filtered_traces[selected_channel_index]
        times = pyconvert(Vector, trace.times("relative"))
        data = pyconvert(Vector, trace.data)

        seis_trace = scatter(x=times, y=data, mode="lines", name="Seismogram", line=attr(color="#38bdf8"))
        phase_colors = distinguishable_colors(length(arrivals))
        arrival_traces = map(1:length(arrivals)) do i
            arr_time = pyconvert(Float64, arrivals[i-1].time)
            arr_name = pyconvert(String, arrivals[i-1].name)
            if 0.0 < arr_time < total_duration
                scatter(x=[arr_time, arr_time], y=[minimum(data), maximum(data)],
                    mode="lines", name=arr_name,
                    line=attr(dash="solid", color="#" * hex(phase_colors[i])))
            else
                nothing
            end
        end
        arrival_traces = filter(!isnothing, arrival_traces)

        plot([seis_trace; arrival_traces...], Layout(
            title="Seismogram at $(selected_station_details.code) ($(round(epi_dist_deg, digits=1))° away)",
            xaxis_title="Time Relative to Earthquake Origin (s)",
            yaxis_title="Amplitude",
            paper_bgcolor="#000", plot_bgcolor="#000",
            font=attr(color="#d1d5db"),
            xaxis=attr(gridcolor="#1f2937"), yaxis=attr(gridcolor="#1f2937"),
        ))
    end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Colors = "~0.13.1"
CondaPkg = "~0.2.36"
FFTW = "~1.10.0"
PlutoPlotly = "~0.6.6"
PlutoUI = "~0.7.83"
PythonCall = "~0.9.35"
StatsBase = "~0.34.12"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "e8271f6fcad4fb71b36b9d685d3fe53cbcf99239"

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
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

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
deps = ["JSON", "Markdown", "MicroMamba", "Pidfile", "Pkg", "Preferences", "Scratch", "TOML", "pixi_jll"]
git-tree-sha1 = "2b1afb8ae65a0758795b00adafb37f97e67ef0e9"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.36"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

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
git-tree-sha1 = "6866aec60ef98e3164cd8d6855225684207e9dff"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.12+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

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
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

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
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c89d196f5ffb64bfbf80985b699ea913b0d2c211"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.6.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
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
git-tree-sha1 = "bba2d9aa057d8f126415de240573e86a8f39d2a1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.1"

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
git-tree-sha1 = "535656ce55266bfed0575cd051acc4f36dc869a0"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.15"

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
git-tree-sha1 = "94ba93778373a53bfd5a0caaf7d809c445292ff4"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.2"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "32a4e09c5f29402573d673901778a0e03b0807b9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.6"

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
git-tree-sha1 = "2b9e3d771adfe535a4fdda855f4741fdaacd3f7f"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.6"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Preferences", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "2b67e030054dd9438a00e3d7f59927e839b00569"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.35"

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
git-tree-sha1 = "67a144433c4ce877ee6d1ada69a124d6b1ecf7be"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.6.2"

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
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

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
git-tree-sha1 = "e4d7a1a0edc20af42689ea6f4f3587a2175d50ee"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.12"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "82bee338d650aa515f31866c460cb7e3bcef90b8"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
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
git-tree-sha1 = "0f38a06c83f0007bbab3cf911262841c9a0f07e0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.13.0"

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
git-tree-sha1 = "3b0738bd7c5645641845da25cbd99800b8718689"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.2"

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
git-tree-sha1 = "717df6f6892af4ee13279a73aa58474e58a88667"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "2.3.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "da8c1f6eee04831f14edcfa5dae611d309807e57"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.3.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.pixi_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "3667b0931a7fe50f0a5554c61af00e5640019e21"
uuid = "4d7b5844-a134-5dcd-ac86-c8f19cd51bed"
version = "0.63.2+0"
"""

# ╔═╡ Cell order:
# ╠═ca9f95ec-71d2-11f0-2e94-7da041469cdf
# ╟─5a72389a-1277-4e18-a532-ce81a4b18c74
# ╟─6a1e0001-0000-4000-8000-000000000001
# ╟─6a1e0001-0000-4000-8000-000000000002
# ╟─6a1e0001-0000-4000-8000-000000000003
# ╟─6a1e0001-0000-4000-8000-000000000004
# ╟─6a1e0001-0000-4000-8000-000000000005
# ╟─6a1e0001-0000-4000-8000-000000000006
# ╟─6a1e0001-0000-4000-8000-000000000007
# ╟─6a1e0001-0000-4000-8000-000000000008
# ╟─6a1e0001-0000-4000-8000-00000000000b
# ╟─6a1e0001-0000-4000-8000-000000000009
# ╟─6a1e0001-0000-4000-8000-00000000000a
# ╟─d2aba787-dae6-4da5-8f16-0cadfcc38371
# ╟─6a1e0001-0000-4000-8000-00000000000c
# ╠═2cc62eeb-edc2-4fcb-876b-010d8a180162
# ╠═53f5276e-7954-43d5-bcfb-6898bc6e85ab
# ╟─6a1e0001-0000-4000-8000-00000000000d
# ╠═6a1e0001-0000-4000-8000-00000000000e
# ╠═6a1e0001-0000-4000-8000-00000000000f
# ╠═6a1e0001-0000-4000-8000-000000000010
# ╠═6a1e0001-0000-4000-8000-000000000011
# ╠═6a1e0001-0000-4000-8000-000000000012
# ╠═0b12e2c0-9183-11f1-898e-dbccf4ea72db
# ╠═6a1e0001-0000-4000-8000-000000000013
# ╟─6a1e0001-0000-4000-8000-000000000014
# ╠═6a1e0001-0000-4000-8000-000000000015
# ╠═6a1e0001-0000-4000-8000-000000000016
# ╟─6a1e0001-0000-4000-8000-000000000017
# ╠═6a1e0001-0000-4000-8000-000000000018
# ╠═6a1e0001-0000-4000-8000-000000000019
# ╠═6a1e0001-0000-4000-8000-00000000001a
# ╟─6a1e0001-0000-4000-8000-00000000001b
# ╠═6a1e0001-0000-4000-8000-00000000001c
# ╠═6a1e0001-0000-4000-8000-00000000001d
# ╠═6a1e0001-0000-4000-8000-00000000001e
# ╟─6a1e0001-0000-4000-8000-00000000001f
# ╠═6a1e0001-0000-4000-8000-000000000020
# ╠═6a1e0001-0000-4000-8000-000000000021
# ╟─6a1e0001-0000-4000-8000-000000000022
# ╠═6a1e0001-0000-4000-8000-000000000023
# ╠═6a1e0001-0000-4000-8000-000000000024
# ╠═6a1e0001-0000-4000-8000-000000000025
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
