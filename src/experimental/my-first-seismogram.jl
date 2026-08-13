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
    # Pin below Python 3.14: PythonCall + this Malt-based worker process (the exact
    # mechanism Pluto uses to run cells) segfaults inside CPython itself
    # (`slot_tp_init`/`type_call`) when constructing an `obspy` FDSN `Client` under
    # Python 3.14.6 -- reproduced independently outside Pluto with a bare `Malt.Worker`,
    # not present under Python 3.12. Very likely a native-extension/C-API compatibility
    # gap in a 3.14 release that's only a couple of months old.
    CondaPkg.add("python"; version="3.12.*")
    CondaPkg.add_pip("obspy")
end

# ╔═╡ 53f5276e-7954-43d5-bcfb-6898bc6e85ab
begin
    using PythonCall
    using PlutoUI
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
## Pick an Earthquake, Pick a Station, See the Seismogram

Not every earthquake makes a useful teaching example -- too small and the signal is buried in
noise at a global station, too close and there's no interesting wave propagation to see. On
the **left** map, drag a box to a region you're curious about, set a minimum magnitude, and
hit **Search**; click any dot to select that earthquake (dot size scales with magnitude).

A station only gives you a seismogram for this earthquake if it was actually **operating**
during it -- picking blind can land on a station with no data for that time window at all.
Once an earthquake is selected, the **right** map's Search queries for GSN (`IU`) stations
that were online across the whole event window, not just ones that happen to be nearby. The
star marks the epicenter; once you pick a station the line between them is the great-circle
path the waves actually traveled -- the epicentral distance is exactly what set
`minradius`/`maxradius` in this notebook's earlier version, now read off directly instead of
being a search parameter.

Once you have a station, hit **Download** -- that's the only step that touches the network
again. The **channel**, **duration**, and **frequency band** controls below the seismogram
panel only ever re-slice or re-filter what's already on your machine (duration can trim the
downloaded window shorter, but not extend it past what you asked for at download time), so
drag them as freely as you like.
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

# ╔═╡ 1294ec0f-34e3-4f9e-9e1a-0d86c6bb7674
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
    coast_raw, _coast_header = readdlm(joinpath(@__DIR__, "..", "assets", "data", "coastlines_110m.csv"), ',';
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
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000017
md"### Layer 2: Earthquake & Station Search"

# ╔═╡ 6a1e0001-0000-4000-8000-000000000018
begin
    """
        SeismogramExplorerInput(; eqlatmin=-60.0, eqlatmax=60.0, eqlonmin=-180.0, eqlonmax=180.0,
                                   minmag=6.5, stlatmin=-60.0, stlatmax=60.0, stlonmin=-180.0, stlonmax=180.0,
                                   duration=2000.0, periodmin=2.0, periodmax=20.0)

    One widget, start to finish. **Left**/**right** maps: drag a box, press **Search** (USGS
    for earthquakes, EarthScope/GSN `IU` for stations -- station search only returns stations
    that were actually operating across the event window), click a dot/triangle to select.
    The station panel's epicenter star, great-circle path, and Mohr... *(no relation --
    just the traction path)* are all local state, so picking a new earthquake naturally
    clears any stale station search. **Bottom row**: once both are selected, press
    **Download** -- the *only* step that touches the network again, fetching every channel
    for the current duration setting (pinned to one location code, since a station can have
    several co-located instruments) plus, whenever both horizontals are available, ObsPy-rotated
    **radial/transverse** components (`->ZNE` first if the horizontals aren't already true-north
    aligned, then `NE->RT` using the source-to-station back-azimuth) appended to the channel list
    alongside the originals. After that, the **channel** dropdown, **duration**
    slider, and **period band** (0.1-200 s, the natural seismological unit -- log-scaled
    dual-handle slider, since a linear one would give almost no resolution below ~10 s over
    that range) are pure client-side redraws: duration re-slices the already-downloaded
    samples (it can only shorten the window, not extend past what was actually fetched), and
    the period band runs a real zero-phase Butterworth-style bandpass (cascaded 2-pole
    high/low sections via `filtfilt`, the standard bilinear-transform biquad design -- not
    bit-identical to ObsPy's SciPy-based filter, but the same qualitative behavior) entirely
    in JavaScript, converting period back to frequency only at the point the filter needs it,
    so none of the three ever re-triggers a download. Phase arrivals (TauP, depth + distance
    only) are pushed in once both ends are picked and drawn as vertical markers over the
    trace.
    """
    struct SeismogramExplorerInput
        eqlatmin::Float64
        eqlatmax::Float64
        eqlonmin::Float64
        eqlonmax::Float64
        minmag::Float64
        stlatmin::Float64
        stlatmax::Float64
        stlonmin::Float64
        stlonmax::Float64
        duration::Float64
        periodmin::Float64
        periodmax::Float64
    end

    SeismogramExplorerInput(; eqlatmin=-60.0, eqlatmax=60.0, eqlonmin=-180.0, eqlonmax=180.0, minmag=6.5,
        stlatmin=-60.0, stlatmax=60.0, stlonmin=-180.0, stlonmax=180.0,
        duration=2000.0, periodmin=2.0, periodmax=20.0) =
        SeismogramExplorerInput(Float64(eqlatmin), Float64(eqlatmax), Float64(eqlonmin), Float64(eqlonmax), Float64(minmag),
            Float64(stlatmin), Float64(stlatmax), Float64(stlonmin), Float64(stlonmax),
            Float64(duration), Float64(periodmin), Float64(periodmax))

    Base.get(w::SeismogramExplorerInput) = Dict{String,Any}(
        "action" => "init",
        "eqlatmin" => w.eqlatmin, "eqlatmax" => w.eqlatmax, "eqlonmin" => w.eqlonmin, "eqlonmax" => w.eqlonmax,
        "minmag" => w.minmag,
        "stlatmin" => w.stlatmin, "stlatmax" => w.stlatmax, "stlonmin" => w.stlonmin, "stlonmax" => w.stlonmax,
        "duration" => w.duration, "periodmin" => w.periodmin, "periodmax" => w.periodmax,
    )

    function Base.show(io::IO, ::MIME"text/html", w::SeismogramExplorerInput)
        write(io, """
<div id="sexwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#sexwidget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #sexwidget .sex-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #sexwidget .sex-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #sexwidget .sex-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #sexwidget .sex-panels{display:flex;gap:14px;align-items:flex-start;justify-content:center;width:100%}
    #sexwidget .sex-col{flex:1 1 0;min-width:0;display:flex;flex-direction:column;align-items:center}
    #sexwidget .sex-panel-label{font-size:12px;color:#6b7280;margin-bottom:4px;text-transform:uppercase;letter-spacing:.04em}
    #sexwidget canvas{cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block;max-width:100%}
    #sexwidget .sex-seis-row{width:100%;margin-top:14px;display:flex;flex-direction:column;align-items:center}
    #sexwidget .sex-controls{display:flex;gap:10px;flex-wrap:wrap;width:100%;margin-top:12px}
    #sexwidget .sex-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px;flex:1 1 220px;min-width:220px}
    #sexwidget .sex-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #sexwidget .sex-control-row{display:grid;grid-template-columns:minmax(60px,90px) minmax(70px,1fr) minmax(40px,56px);gap:6px;align-items:center;margin:5px 0}
    #sexwidget .sex-control-row input[type=range]{width:100%;min-width:0;vertical-align:middle}
    #sexwidget .sex-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums;font-size:13px}
    #sexwidget .sex-readout{font-size:13px;line-height:1.6}
    #sexwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px;cursor:pointer}
    #sexwidget button:disabled{opacity:0.5;cursor:default}
    #sexwidget select{background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:4px 6px;font-size:13px}
    #sexwidget .sex-dualrange{position:relative;height:22px}
    #sexwidget .sex-dualrange-track{position:absolute;left:0;right:0;top:9px;height:4px;background:#374151;border-radius:2px}
    #sexwidget .sex-dualrange-fill{position:absolute;top:9px;height:4px;background:#38bdf8;border-radius:2px;pointer-events:none}
    #sexwidget .sex-dualrange input[type=range]{position:absolute;left:0;top:0;width:100%;height:22px;margin:0;background:transparent;-webkit-appearance:none;appearance:none;pointer-events:none}
    #sexwidget .sex-dualrange input[type=range]::-webkit-slider-runnable-track{background:transparent;height:4px}
    #sexwidget .sex-dualrange input[type=range]::-webkit-slider-thumb{-webkit-appearance:none;appearance:none;pointer-events:auto;width:14px;height:14px;border-radius:50%;background:#38bdf8;border:1px solid #0a0f18;margin-top:3px;cursor:pointer}
    #sexwidget .sex-dualrange input[type=range]::-moz-range-track{background:transparent;height:4px;border:none}
    #sexwidget .sex-dualrange input[type=range]::-moz-range-thumb{pointer-events:auto;width:14px;height:14px;border-radius:50%;background:#38bdf8;border:1px solid #0a0f18;cursor:pointer}
    @media (max-width: 900px){
      #sexwidget .sex-panels{flex-direction:column;align-items:center}
      #sexwidget .sex-col{width:100%;max-width:560px}
    }
  </style>
  <div class="sex-title">
    <div class="sex-title-desc">My First Seismogram</div>
    <div class="sex-title-hint">drag a box &middot; Search &middot; click a dot &middot; Download &middot; channel/duration/band are free to explore</div>
  </div>
  <div class="sex-panels">
    <div class="sex-col">
      <div class="sex-panel-label">Earthquakes (drag a box, then Search)</div>
      <canvas id="sexEqMap" width="600" height="300"></canvas>
    </div>
    <div class="sex-col">
      <div class="sex-panel-label">Stations (drag a box, then Search)</div>
      <canvas id="sexStMap" width="600" height="300"></canvas>
    </div>
  </div>
  <div class="sex-seis-row">
    <div class="sex-panel-label">Seismogram</div>
    <canvas id="sexSeis" width="1214" height="220"></canvas>
  </div>
  <div class="sex-controls">
    <div class="sex-control-group">
      <div class="sex-control-title">Earthquake Search</div>
      <label class="sex-control-row"><span>min mag</span><input type="range" id="sex-minmag" min="4" max="8.5" step="0.1" value="$(w.minmag)"><span id="sex-minmag-v" class="sex-value">$(w.minmag)</span></label>
      <div style="display:flex;gap:10px;align-items:center;margin-top:6px">
        <button id="sex-eq-search" type="button">Search</button>
        <span id="sex-eq-status" style="font-size:13px">drag a box above, then Search</span>
      </div>
    </div>
    <div class="sex-control-group">
      <div class="sex-control-title">Station Search</div>
      <div style="display:flex;gap:10px;align-items:center">
        <button id="sex-st-search" type="button" disabled>Search</button>
        <span id="sex-st-status" style="font-size:13px">select an earthquake first</span>
      </div>
    </div>
    <div class="sex-control-group" style="flex:1 1 320px">
      <div class="sex-control-title">Waveform</div>
      <div class="sex-control-row" style="grid-template-columns:60px 1fr">
        <span>channel 1</span><select id="sex-channel" disabled><option>-- download first --</option></select>
      </div>
      <div class="sex-control-row" style="grid-template-columns:60px 1fr">
        <span>channel 2</span><select id="sex-channel2" disabled><option>-- download first --</option></select>
      </div>
      <label class="sex-control-row"><span>duration</span><input type="range" id="sex-duration" min="100" max="5000" step="50" value="$(w.duration)"><span id="sex-duration-v" class="sex-value">$(w.duration) s</span></label>
      <div class="sex-control-row" style="grid-template-columns:60px 1fr">
        <span>period</span>
        <div class="sex-dualrange" id="sex-period-wrap">
          <div class="sex-dualrange-track"></div>
          <div class="sex-dualrange-fill" id="sex-period-fill"></div>
          <input type="range" id="sex-period-min" min="-1" max="2.30103" step="0.005" value="$(log10(w.periodmin))">
          <input type="range" id="sex-period-max" min="-1" max="2.30103" step="0.005" value="$(log10(w.periodmax))">
        </div>
      </div>
      <div class="sex-control-row" style="grid-template-columns:60px 1fr">
        <span></span><span id="sex-period-v" class="sex-value">$(round(w.periodmin, digits=2))–$(round(w.periodmax, digits=1)) s</span>
      </div>
      <div style="display:flex;gap:10px;align-items:center;margin-top:6px">
        <button id="sex-download" type="button" disabled>Download</button>
        <span id="sex-dl-status" style="font-size:13px">select a station first</span>
      </div>
    </div>
    <div class="sex-control-group" style="flex:1 1 260px">
      <div class="sex-control-title">Selection</div>
      <div class="sex-readout" id="sex-readout">no earthquake or station selected yet</div>
    </div>
    <div class="sex-control-group" style="flex:1 1 220px">
      <div class="sex-control-title">Legend</div>
      <div class="sex-readout"><span style="color:#ef4444">&#9679;</span> earthquake (size=magnitude) &middot; <span style="color:#3b82f6">&#9650;</span> station &middot; <span style="color:#facc15">&#9679;/&#9733;</span> selected/epicenter &middot; <span style="color:#22d3ee">line</span> = great-circle path &middot; <span style="color:#facc15">dashed</span> = phase arrival</div>
      <div class="sex-readout"><span style="color:#fb923c">&#8594;</span> R (radial, away from source) &middot; <span style="color:#a78bfa">&#8594;</span> T (transverse) &middot; <span style="color:#38bdf8">&#9678;</span> Z (vertical, out of page) at the selected station</div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  // Pluto's own @bind wiring listens for 'input'/'change' on `par` (the same node our own code
  // calls `par.dispatchEvent(new CustomEvent('input'))` on for eq-search/st-search/download/
  // select actions) -- and it listens in the CAPTURE phase, so a plain `e.stopPropagation()`
  // inside a descendant control's own bubble-phase listener runs too late to stop it: Pluto's
  // capture-phase handler on `par` already fired on the way DOWN to the control, before the
  // control's own handler even runs. Any purely-local control (channel selects, duration/period/
  // minmag sliders) that fires a native 'input'/'change' event otherwise bubbles that event
  // straight into Pluto re-sending the *stale* last-dispatched `sex` value (still the download
  // action) to Julia -- re-running the download gate cell (a real, redundant network re-fetch)
  // and re-firing `sex-waveform`, which silently resets the channel selects/zoom/etc back to
  // their post-download defaults. Registering our OWN capture-phase listener on `par`, as the
  // very first thing here (before Pluto's own bond-scanning script has had a chance to run,
  // since that only happens after this whole cell-output script finishes executing), wins the
  // registration-order race: same node, same phase, first-registered runs first, so
  // stopImmediatePropagation() here reliably blocks Pluto's listener from ever seeing these.
  //
  // IMPORTANT: capture-phase listeners on an ANCESTOR (`par`) always run before ANY listener on
  // a descendant target, no matter what phase that descendant listener itself is registered
  // with -- there's no way to let the event continue on to a normal `el.addEventListener(...)`
  // on the control itself once we've decided to intercept it here. So the actual per-control
  // logic below is registered into `SEX_LOCAL_HANDLERS` (id -> function) and invoked *from
  // inside* this one delegated capture handler, instead of via separate listeners on each
  // control -- that's what makes the controls responsive again while still keeping Pluto from
  // ever seeing the event.
  const SEX_LOCAL_HANDLERS = {}
  ;['input','change'].forEach(evt => par.addEventListener(evt, e => {
    const id = e.target && e.target.id
    if(id && SEX_LOCAL_HANDLERS[id]){ e.stopImmediatePropagation(); SEX_LOCAL_HANDLERS[id](e) }
  }, true))
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1400)
  const W = Math.round(Math.max(280, (availW-14)/2))
  const H = Math.round(W/2)
  const FULLW = Math.round(2*W+14)
  // Tall enough to comfortably stack two traces (channel 1/2) when both are selected -- a
  // single trace just uses the same height with more headroom around it, no separate case.
  const SEISH = Math.round(Math.min(340, Math.max(180, FULLW*0.28)))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)

  function hidpi(cv, cx, w, h){ cv.width=Math.round(w*DPR); cv.height=Math.round(h*DPR); cv.style.width=w+'px'; cv.style.height=h+'px'; cx.setTransform(DPR,0,0,DPR,0,0) }
  const eqCvs = par.querySelector('#sexEqMap'), eqCtx = eqCvs.getContext('2d')
  const stCvs = par.querySelector('#sexStMap'), stCtx = stCvs.getContext('2d')
  const seisCvs = par.querySelector('#sexSeis'), seisCtx = seisCvs.getContext('2d')
  hidpi(eqCvs, eqCtx, W, H); hidpi(stCvs, stCtx, W, H); hidpi(seisCvs, seisCtx, FULLW, SEISH)

  const COAST = $(coast_js_literal(coast_line_id, coast_lon, coast_lat))

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
  function drawCoast(ctx){
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
  function drawBox(ctx, latmin,latmax,lonmin,lonmax){
    const [x0,y0] = lonlatToXY(lonmin, latmax)
    const [x1,y1] = lonlatToXY(lonmax, latmin)
    ctx.strokeStyle = '#facc15'; ctx.lineWidth = 2
    ctx.strokeRect(Math.min(x0,x1), Math.min(y0,y1), Math.abs(x1-x0), Math.abs(y1-y0))
    ctx.fillStyle = 'rgba(250,204,21,0.08)'
    ctx.fillRect(Math.min(x0,x1), Math.min(y0,y1), Math.abs(x1-x0), Math.abs(y1-y0))
  }
  function drawGreatCircle(ctx, lon1,lat1,lon2,lat2){
    const A = llToXYZ(lon1,lat1), B = llToXYZ(lon2,lat2)
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
  // The radial/transverse split (see the ObsPy `rotate('NE->RT', ...)` step in the Appendix) is
  // a geometric fact about the source-station geometry, not just a filter -- R is the horizontal
  // motion pointing away from the epicenter along the great-circle path (continuing past the
  // station), T is 90° from it in the horizontal plane, and Z is straight up, perpendicular to
  // both (drawn as the standard "out of the page" dot-in-circle glyph since a flat map can't show
  // a genuinely vertical arrow). Approximating the tangent at the station with the flat-map
  // epicenter->station direction is a small, deliberate simplification -- it's exact along the
  // great-circle's own path and only drifts from the true spherical tangent far off that path,
  // which is not where these arrows are drawn.
  function drawComponentAxes(ctx, stX, stY, eqLon, eqLat, lon, lat){
    const [ex, ey] = lonlatToXY(eqLon, eqLat)
    let dx = stX-ex, dy = stY-ey
    const len = Math.hypot(dx,dy) || 1
    dx/=len; dy/=len                // R: epicenter -> station, continuing outward
    const tx = -dy, ty = dx          // T: perpendicular to R, in-plane
    const ARR = 26
    function arrow(ux,uy,color,label){
      const x2 = stX+ux*ARR, y2 = stY+uy*ARR
      ctx.strokeStyle=color; ctx.fillStyle=color; ctx.lineWidth=2
      ctx.beginPath(); ctx.moveTo(stX,stY); ctx.lineTo(x2,y2); ctx.stroke()
      const ang = Math.atan2(uy,ux)
      ctx.beginPath(); ctx.moveTo(x2,y2)
      ctx.lineTo(x2-7*Math.cos(ang-0.4), y2-7*Math.sin(ang-0.4))
      ctx.lineTo(x2-7*Math.cos(ang+0.4), y2-7*Math.sin(ang+0.4))
      ctx.closePath(); ctx.fill()
      ctx.font='11px sans-serif'
      ctx.fillText(label, x2+ux*11, y2+uy*11+4)
    }
    arrow(dx,dy,'#fb923c','R')
    arrow(tx,ty,'#a78bfa','T')
    const zx = stX-dx*16+tx*16, zy = stY-dy*16+ty*16   // offset so it doesn't sit on the marker
    ctx.strokeStyle='#38bdf8'; ctx.fillStyle='#38bdf8'; ctx.lineWidth=1.5
    ctx.beginPath(); ctx.arc(zx,zy,7,0,2*Math.PI); ctx.stroke()
    ctx.beginPath(); ctx.arc(zx,zy,2,0,2*Math.PI); ctx.fill()
    ctx.font='11px sans-serif'; ctx.fillText('Z', zx+10, zy+4)
  }

  // ---------- Earthquake panel state ----------
  let eqLatmin=$(w.eqlatmin), eqLatmax=$(w.eqlatmax), eqLonmin=$(w.eqlonmin), eqLonmax=$(w.eqlonmax), minmag=$(w.minmag)
  let eqData = null, eqSelectedIdx = -1
  let eqDragStart = null, eqDragMoved = false

  function eqRadius(mag){ return Math.max(3, Math.min(18, 3 + (mag-4)*2.8)) }

  function drawEqEvents(){
    if(!eqData) return
    for(let i=0;i<eqData.length;i++){
      const [lon,lat,mag] = eqData[i]
      const [x,y] = lonlatToXY(lon,lat)
      const r = eqRadius(mag)
      eqCtx.beginPath(); eqCtx.arc(x,y,r,0,2*Math.PI)
      eqCtx.fillStyle = i===eqSelectedIdx ? '#facc15' : 'rgba(239,68,68,0.7)'
      eqCtx.fill()
      if(i===eqSelectedIdx){ eqCtx.strokeStyle='#0a0f18'; eqCtx.lineWidth=1.5; eqCtx.stroke() }
    }
  }
  function redrawEq(){
    eqCtx.clearRect(0,0,W,H); eqCtx.fillStyle='#000'; eqCtx.fillRect(0,0,W,H)
    drawCoast(eqCtx); drawEqEvents(); drawBox(eqCtx, eqLatmin,eqLatmax,eqLonmin,eqLonmax)
  }
  redrawEq()

  const eqSearchBtn = par.querySelector('#sex-eq-search')
  const eqStatusEl = par.querySelector('#sex-eq-status')
  eqSearchBtn.addEventListener('click', ()=>{
    eqSearchBtn.disabled = true
    eqStatusEl.textContent = 'searching…'
    par.value = {action:'eq-search', latmin:eqLatmin, latmax:eqLatmax, lonmin:eqLonmin, lonmax:eqLonmax, minmag}
    par.dispatchEvent(new CustomEvent('input'))
  })
  window.addEventListener('sex-eq-results', e=>{
    const d = JSON.parse(e.detail)
    eqSearchBtn.disabled = false
    if(!d.ok){ eqStatusEl.textContent = 'error: ' + d.error; eqData = []; redrawEq(); return }
    eqData = d.events
    eqSelectedIdx = -1
    eqStatusEl.textContent = d.count ? (d.count+' found — click one') : 'none found, widen the box or lower magnitude'
    redrawEq()
  })

  function fmt(v){ return v.toFixed(2) }
  function updateReadout(){
    let html = ''
    if(eqSelectedIdx>=0){
      const [lon,lat,mag,depth,time] = eqData[eqSelectedIdx]
      html += '<b>Earthquake</b>: M'+mag.toFixed(1)+' at '+time+'<br>('+fmt(lat)+'°, '+fmt(lon)+'°), depth '+depth.toFixed(0)+' km<br>'
    } else {
      html += 'no earthquake selected yet<br>'
    }
    if(stSelectedIdx>=0){
      const [code,lat,lon,elevation] = stData[stSelectedIdx]
      let distTxt = ''
      if(eqSelectedIdx>=0){
        const [eqLon,eqLat] = eqData[eqSelectedIdx]
        const A=llToXYZ(eqLon,eqLat), B=llToXYZ(lon,lat)
        const dot=A[0]*B[0]+A[1]*B[1]+A[2]*B[2]
        const distDeg = Math.acos(Math.max(-1,Math.min(1,dot)))*180/Math.PI
        distTxt = ', '+distDeg.toFixed(1)+'&deg; away'
      }
      html += '<b>Station</b>: '+code+' ('+lat.toFixed(2)+'°, '+lon.toFixed(2)+'°)'+distTxt
    } else {
      html += 'no station selected yet'
    }
    par.querySelector('#sex-readout').innerHTML = html
  }

  function tryEqSelect(mx, my){
    if(!eqData || eqData.length===0) return
    let best=-1, bestD=Infinity
    for(let i=0;i<eqData.length;i++){
      const [lon,lat,mag] = eqData[i]
      const [x,y] = lonlatToXY(lon,lat)
      const d = Math.hypot(x-mx,y-my)
      if(d<bestD){ bestD=d; best=i }
    }
    const [,,bestMag] = eqData[best]
    if(bestD < eqRadius(bestMag)+6){
      eqSelectedIdx = best
      const [lon,lat,mag,depth,time] = eqData[best]
      updateReadout()
      par.value = {action:'eq-select', lat, lon, depth, mag, time}
      par.dispatchEvent(new CustomEvent('input'))
      stSearchBtn.disabled = false
      stStatusEl.textContent = 'drag a box above, then Search'
      stData = null; stSelectedIdx = -1
      dlBtn.disabled = true
      dlStatusEl.textContent = 'select a station first'
      redrawEq(); redrawSt()
    }
  }

  eqCvs.addEventListener('mousedown', e=>{ eqDragStart=[e.offsetX,e.offsetY]; eqDragMoved=false })
  eqCvs.addEventListener('mousemove', e=>{
    if(!eqDragStart) return
    const dx=e.offsetX-eqDragStart[0], dy=e.offsetY-eqDragStart[1]
    if(Math.hypot(dx,dy) > 6) eqDragMoved = true
    if(eqDragMoved){
      const [lonA,latA] = xyToLonLat(eqDragStart[0], eqDragStart[1])
      const [lonB,latB] = xyToLonLat(e.offsetX, e.offsetY)
      eqLatmin=Math.min(latA,latB); eqLatmax=Math.max(latA,latB)
      eqLonmin=Math.min(lonA,lonB); eqLonmax=Math.max(lonA,lonB)
      redrawEq()
    }
  })
  eqCvs.addEventListener('mouseup', e=>{
    if(!eqDragStart) return
    if(!eqDragMoved) tryEqSelect(e.offsetX, e.offsetY)
    eqDragStart=null; eqDragMoved=false
  })
  window.addEventListener('mouseup', ()=>{ eqDragStart=null; eqDragMoved=false; stDragStart=null; stDragMoved=false })

  SEX_LOCAL_HANDLERS['sex-minmag'] = e=>{
    minmag = parseFloat(e.target.value)
    par.querySelector('#sex-minmag-v').textContent = minmag.toFixed(1)
  }

  // ---------- Station panel state ----------
  let stLatmin=$(w.stlatmin), stLatmax=$(w.stlatmax), stLonmin=$(w.stlonmin), stLonmax=$(w.stlonmax)
  let stData = null, stSelectedIdx = -1
  let stDragStart = null, stDragMoved = false

  function drawEpicenter(){
    if(eqSelectedIdx<0) return
    const [lon,lat] = eqData[eqSelectedIdx]
    const [x,y] = lonlatToXY(lon,lat)
    stCtx.fillStyle = '#facc15'
    stCtx.beginPath()
    for(let i=0;i<5;i++){
      const a=-Math.PI/2+i*2*Math.PI/5, a2=a+Math.PI/5
      const p1=[x+Math.cos(a)*8, y+Math.sin(a)*8], p2=[x+Math.cos(a2)*3, y+Math.sin(a2)*3]
      i===0 ? stCtx.moveTo(p1[0],p1[1]) : stCtx.lineTo(p1[0],p1[1])
      stCtx.lineTo(p2[0],p2[1])
    }
    stCtx.closePath(); stCtx.fill()
  }
  function drawStations(){
    if(!stData) return
    for(let i=0;i<stData.length;i++){
      const [code,lat,lon] = stData[i]
      const [x,y] = lonlatToXY(lon,lat)
      stCtx.beginPath()
      stCtx.moveTo(x, y-6); stCtx.lineTo(x-6, y+5); stCtx.lineTo(x+6, y+5); stCtx.closePath()
      stCtx.fillStyle = i===stSelectedIdx ? '#facc15' : 'rgba(59,130,246,0.8)'
      stCtx.fill()
      if(i===stSelectedIdx){ stCtx.strokeStyle='#0a0f18'; stCtx.lineWidth=1.5; stCtx.stroke() }
    }
  }
  function redrawSt(){
    stCtx.clearRect(0,0,W,H); stCtx.fillStyle='#000'; stCtx.fillRect(0,0,W,H)
    drawCoast(stCtx)
    if(stSelectedIdx>=0 && eqSelectedIdx>=0){
      const [,lat,lon] = stData[stSelectedIdx]
      const [eqLon,eqLat] = eqData[eqSelectedIdx]
      drawGreatCircle(stCtx, eqLon,eqLat, lon,lat)
    }
    drawEpicenter(); drawStations(); drawBox(stCtx, stLatmin,stLatmax,stLonmin,stLonmax)
    if(stSelectedIdx>=0 && eqSelectedIdx>=0){
      const [,lat,lon] = stData[stSelectedIdx]
      const [eqLon,eqLat] = eqData[eqSelectedIdx]
      const [stX,stY] = lonlatToXY(lon,lat)
      drawComponentAxes(stCtx, stX,stY, eqLon,eqLat, lon,lat)
    }
  }
  redrawSt()

  const stSearchBtn = par.querySelector('#sex-st-search')
  const stStatusEl = par.querySelector('#sex-st-status')
  const dlBtn = par.querySelector('#sex-download')
  const dlStatusEl = par.querySelector('#sex-dl-status')
  function eqPayloadFields(){
    if(eqSelectedIdx<0) return {}
    const [eqLon,eqLat,eqMag,eqDepth,eqTime] = eqData[eqSelectedIdx]
    return {eqlat:eqLat, eqlon:eqLon, eqdepth:eqDepth, eqmag:eqMag, eqtime:eqTime}
  }
  function stPayloadFields(){
    if(stSelectedIdx<0) return {}
    const [stCode,stLat,stLon,stElevation] = stData[stSelectedIdx]
    return {stcode:stCode, stlat:stLat, stlon:stLon, stelevation:stElevation}
  }
  stSearchBtn.addEventListener('click', ()=>{
    if(eqSelectedIdx<0) return
    stSearchBtn.disabled = true
    stStatusEl.textContent = 'searching…'
    par.value = Object.assign({action:'st-search', latmin:stLatmin, latmax:stLatmax, lonmin:stLonmin, lonmax:stLonmax}, eqPayloadFields())
    par.dispatchEvent(new CustomEvent('input'))
  })
  window.addEventListener('sex-st-results', e=>{
    const d = JSON.parse(e.detail)
    stSearchBtn.disabled = (eqSelectedIdx<0)
    if(!d.ok){ stStatusEl.textContent = 'error: ' + d.error; stData = []; redrawSt(); return }
    stData = d.stations
    stSelectedIdx = -1
    stStatusEl.textContent = d.count ? (d.count+' found — click one') : 'none were online for this event in that box'
    redrawSt()
  })

  function trySelectSt(mx, my){
    if(!stData || stData.length===0) return
    let best=-1, bestD=Infinity
    for(let i=0;i<stData.length;i++){
      const [,lat,lon] = stData[i]
      const [x,y] = lonlatToXY(lon,lat)
      const d = Math.hypot(x-mx,y-my)
      if(d<bestD){ bestD=d; best=i }
    }
    if(bestD < 14){
      stSelectedIdx = best
      const [code,lat,lon,elevation] = stData[best]
      updateReadout()
      par.value = Object.assign({action:'st-select', code, lat, lon, elevation}, eqPayloadFields())
      par.dispatchEvent(new CustomEvent('input'))
      dlBtn.disabled = false
      dlStatusEl.textContent = 'ready to download'
      redrawSt()
    }
  }

  stCvs.addEventListener('mousedown', e=>{ stDragStart=[e.offsetX,e.offsetY]; stDragMoved=false })
  stCvs.addEventListener('mousemove', e=>{
    if(!stDragStart) return
    const dx=e.offsetX-stDragStart[0], dy=e.offsetY-stDragStart[1]
    if(Math.hypot(dx,dy) > 6) stDragMoved = true
    if(stDragMoved){
      const [lonA,latA] = xyToLonLat(stDragStart[0], stDragStart[1])
      const [lonB,latB] = xyToLonLat(e.offsetX, e.offsetY)
      stLatmin=Math.min(latA,latB); stLatmax=Math.max(latA,latB)
      stLonmin=Math.min(lonA,lonB); stLonmax=Math.max(lonA,lonB)
      redrawSt()
    }
  })
  stCvs.addEventListener('mouseup', e=>{
    if(!stDragStart) return
    if(!stDragMoved) trySelectSt(e.offsetX, e.offsetY)
    stDragStart=null; stDragMoved=false
  })

  // ---------- Waveform / seismogram panel state ----------
  let waveform = null           // [{name, sr, data:[...]}, ...] pushed after Download
  let selectedChannelIdx = 0
  let selectedChannelIdx2 = -1  // -1 = "none" -- second trace is optional, off by default
  let duration = $(w.duration)
  // Stored as period (seconds, the natural seismological unit) -- converted to Hz only at
  // the point the Butterworth filter actually needs it. periodMin is the *shorter* period
  // (higher frequency) edge of the band, periodMax the *longer* period (lower frequency) edge.
  let periodMin = $(w.periodmin), periodMax = $(w.periodmax)
  const PERIOD_LOG_MIN = -1, PERIOD_LOG_MAX = Math.log10(200)   // 0.1 s .. 200 s
  let arrivalsData = null       // [[name, time], ...]
  // Plotly-style horizontal zoom: drag a range on the seismogram to zoom into it, double-click
  // to reset. `viewEnd===null` means "full (duration-trimmed) window" -- resetting on a fresh
  // download or a duration-slider change (both change what "full" even means) but *not* on a
  // channel switch, so zooming in on one component and flipping to another compares the same
  // time window.
  let viewStart = 0, viewEnd = null
  // Updated by drawSeismogram() every call; the drag handlers use these to invert a pixel
  // position back to a time, since the pixel<->time mapping depends on the current zoom.
  let seisX0 = 0, seisX1 = 0, seisVStart = 0, seisVEnd = 0

  window.addEventListener('sex-arrivals', e=>{
    const d = JSON.parse(e.detail)
    arrivalsData = d.arrivals
    drawSeismogram()
  })

  dlBtn.addEventListener('click', ()=>{
    dlBtn.disabled = true
    dlStatusEl.textContent = 'downloading…'
    par.value = Object.assign({action:'download', duration}, eqPayloadFields(), stPayloadFields())
    par.dispatchEvent(new CustomEvent('input'))
  })
  window.addEventListener('sex-waveform', e=>{
    const d = JSON.parse(e.detail)
    dlBtn.disabled = false
    if(!d.ok){ dlStatusEl.textContent = 'error: ' + d.error; return }
    waveform = d.channels
    selectedChannelIdx = 0
    selectedChannelIdx2 = -1
    viewStart = 0; viewEnd = null
    const sel = par.querySelector('#sex-channel')
    sel.disabled = false
    sel.innerHTML = ''
    waveform.forEach((ch,i)=>{
      const opt = document.createElement('option')
      opt.value = i; opt.textContent = ch.name
      sel.appendChild(opt)
    })
    sel.value = 0
    const sel2 = par.querySelector('#sex-channel2')
    sel2.disabled = false
    sel2.innerHTML = '<option value="-1">-- none --</option>'
    waveform.forEach((ch,i)=>{
      const opt = document.createElement('option')
      opt.value = i; opt.textContent = ch.name
      sel2.appendChild(opt)
    })
    sel2.value = -1
    dlStatusEl.textContent = 'downloaded ' + waveform.length + ' channel(s)'
    drawSeismogram()
  })
  // Every control below only needs to update local JS state and redraw -- it must NEVER reach
  // Julia. Registered into `SEX_LOCAL_HANDLERS` (see the delegated capture listener near the
  // top of this script) rather than via `addEventListener` on the control itself, because a
  // capture-phase listener on `par` -- an ancestor -- always runs before any listener on a
  // descendant control, so once we intercept there, a separate listener on the control would
  // just never fire.
  SEX_LOCAL_HANDLERS['sex-channel'] = e=>{
    selectedChannelIdx = parseInt(e.target.value)
    drawSeismogram()
  }
  SEX_LOCAL_HANDLERS['sex-channel2'] = e=>{
    selectedChannelIdx2 = parseInt(e.target.value)
    drawSeismogram()
  }
  SEX_LOCAL_HANDLERS['sex-duration'] = e=>{
    duration = parseFloat(e.target.value)
    par.querySelector('#sex-duration-v').textContent = duration.toFixed(0) + ' s'
    viewStart = 0; viewEnd = null   // the base window changed, any previous zoom no longer applies
    drawSeismogram()
  }
  // Two overlapping native <input type=range> elements (the standard no-library way to build
  // a dual-handle range slider) sharing one log10(period) scale -- a linear scale would leave
  // almost no usable resolution below ~10 s across a 0.1-200 s span. A thin fill bar between
  // the two thumbs (position/width set from the same log-percent math) makes the selected band
  // visible at a glance, same as any standard dual-range control.
  function periodPct(v){ return (Math.log10(v)-PERIOD_LOG_MIN)/(PERIOD_LOG_MAX-PERIOD_LOG_MIN)*100 }
  function updatePeriodUI(){
    const fill = par.querySelector('#sex-period-fill')
    const a = periodPct(periodMin), b = periodPct(periodMax)
    fill.style.left = a+'%'; fill.style.width = Math.max(0,b-a)+'%'
    par.querySelector('#sex-period-v').textContent = periodMin.toFixed(2)+'–'+periodMax.toFixed(1)+' s'
  }
  updatePeriodUI()
  SEX_LOCAL_HANDLERS['sex-period-min'] = e=>{
    let v = Math.pow(10, parseFloat(e.target.value))
    if(v > periodMax) v = periodMax   // keep the two handles from crossing
    periodMin = v
    e.target.value = Math.log10(v)
    updatePeriodUI()
    drawSeismogram()
  }
  SEX_LOCAL_HANDLERS['sex-period-max'] = e=>{
    let v = Math.pow(10, parseFloat(e.target.value))
    if(v < periodMin) v = periodMin
    periodMax = v
    e.target.value = Math.log10(v)
    updatePeriodUI()
    drawSeismogram()
  }

  // 2-pole Butterworth high/low-pass biquads (RBJ Audio-EQ-Cookbook design, Q=1/sqrt(2) for
  // the maximally-flat Butterworth response) cascaded into a bandpass and run forward+backward
  // (filtfilt) for zero phase -- the same shape of filter ObsPy applies, implemented here in
  // pure JS so the frequency sliders never need to call back into Julia/Python.
  function butter2(freq, sr, type){
    const Q = Math.SQRT1_2
    const w0 = 2*Math.PI*Math.max(1e-6, Math.min(freq, sr*0.49))/sr
    const alpha = Math.sin(w0)/(2*Q), cosw0 = Math.cos(w0)
    let b0,b1,b2,a0,a1,a2
    if(type==='lowpass'){ b0=(1-cosw0)/2; b1=1-cosw0; b2=(1-cosw0)/2 }
    else { b0=(1+cosw0)/2; b1=-(1+cosw0); b2=(1+cosw0)/2 }
    a0=1+alpha; a1=-2*cosw0; a2=1-alpha
    return {b:[b0/a0,b1/a0,b2/a0], a1:a1/a0, a2:a2/a0}
  }
  function applyBiquad(x, c){
    const y = new Float64Array(x.length)
    for(let n=0;n<x.length;n++){
      let v = c.b[0]*x[n]
      if(n>=1) v += c.b[1]*x[n-1] - c.a1*y[n-1]
      if(n>=2) v += c.b[2]*x[n-2] - c.a2*y[n-2]
      y[n] = v
    }
    return y
  }
  function filtfiltPass(x, c){
    const fwd = applyBiquad(x, c)
    const rev = applyBiquad(fwd.slice().reverse(), c)
    return rev.reverse()
  }
  function bandpass(x, fmin, fmax, sr){
    if(!(fmax>fmin) || x.length<12) return x
    let y = filtfiltPass(x, butter2(fmin, sr, 'highpass'))
    y = filtfiltPass(y, butter2(fmax, sr, 'lowpass'))
    return y
  }

  // Simple min/max-per-pixel-bucket decimation so a redraw stays fast regardless of how many
  // samples are in the current window (a raw per-sample lineTo would be fine up to tens of
  // thousands of points, but this keeps every slider drag equally smooth at any duration).
  function decimateMinMax(arr, targetPx){
    if(arr.length <= targetPx*2){
      const out = new Array(arr.length)
      for(let i=0;i<arr.length;i++) out[i] = [i, arr[i]]
      return out
    }
    const bucket = arr.length/targetPx
    const out = []
    for(let p=0;p<targetPx;p++){
      const s=Math.floor(p*bucket), e=Math.max(s+1, Math.floor((p+1)*bucket))
      let mn=Infinity, mx=-Infinity, mnI=s, mxI=s
      for(let i=s;i<e && i<arr.length;i++){ if(arr[i]<mn){mn=arr[i];mnI=i} if(arr[i]>mx){mx=arr[i];mxI=i} }
      if(mnI<=mxI){ out.push([mnI,mn]); out.push([mxI,mx]) } else { out.push([mxI,mx]); out.push([mnI,mn]) }
    }
    return out
  }

  // Classic "nice number" tick step (1/2/5 x 10^n) so the seismogram's time axis always lands
  // on round values no matter the zoom window, instead of just labeling the two window edges.
  function niceTicks(lo, hi, targetCount){
    const range = hi-lo
    if(!(range>0)) return [lo]
    const roughStep = range/targetCount
    const mag = Math.pow(10, Math.floor(Math.log10(roughStep)))
    const norm = roughStep/mag
    const step = (norm<1.5?1:norm<3?2:norm<7?5:10)*mag
    const out = []
    for(let t=Math.ceil(lo/step)*step; t<=hi+step*1e-6; t+=step) out.push(t)
    return out
  }

  function drawSeismogram(){
    seisCtx.clearRect(0,0,FULLW,SEISH)
    seisCtx.fillStyle='#000'; seisCtx.fillRect(0,0,FULLW,SEISH)
    if(!waveform){
      seisCtx.fillStyle='#6b7280'; seisCtx.font='13px sans-serif'
      seisCtx.fillText('pick an earthquake and station above, then Download', 12, 20)
      return
    }
    // Channel 2 is optional (-1 = "none") -- when active, both traces share one time axis but
    // get their own vertically-offset band (stacked, not overlaid) and their own color/autoscale,
    // the same "record section" idiom seismologists use to compare two components at a glance.
    const ch1 = waveform[selectedChannelIdx]
    const hasCh2 = selectedChannelIdx2 >= 0 && !!waveform[selectedChannelIdx2]
    const ch2 = hasCh2 ? waveform[selectedChannelIdx2] : null
    function prepChannel(ch){
      const sr = ch.sr
      const nSamp = Math.max(2, Math.min(ch.data.length, Math.round(duration*sr)))
      const filtered = bandpass(ch.data.slice(0, nSamp), 1/periodMax, 1/periodMin, sr)
      return {sr, filtered, fullTMax: nSamp/sr}
    }
    const p1 = prepChannel(ch1)
    const p2 = hasCh2 ? prepChannel(ch2) : null
    const fullTMax = hasCh2 ? Math.min(p1.fullTMax, p2.fullTMax) : p1.fullTMax

    // The current zoom window, clamped to the (possibly just-changed) full trimmed trace.
    let vEnd = (viewEnd===null || viewEnd>fullTMax) ? fullTMax : viewEnd
    let vStart = Math.max(0, Math.min(viewStart, vEnd-0.01))

    const PAD_L=46, PAD_R=12, PAD_T=16, PAD_B=26
    const x0=PAD_L, x1=FULLW-PAD_R, y0=PAD_T, y1=SEISH-PAD_B
    seisX0=x0; seisX1=x1; seisVStart=vStart; seisVEnd=vEnd
    const X = t => x0 + ((t-vStart)/(vEnd-vStart))*(x1-x0)

    // One below the other with a small gap between -- not overlaid -- when channel 2 is active;
    // a single full-height band otherwise (identical to the pre-channel-2 layout).
    const bandGap = 14
    const bandH = hasCh2 ? (y1-y0-bandGap)/2 : (y1-y0)
    const bands = hasCh2 ? [{y0, y1:y0+bandH}, {y0:y0+bandH+bandGap, y1}] : [{y0, y1}]

    if(arrivalsData){
      for(const [name,t] of arrivalsData){
        if(t>vStart && t<vEnd){
          const px = X(t)
          seisCtx.setLineDash([3,3])
          seisCtx.strokeStyle='rgba(250,204,21,0.6)'; seisCtx.lineWidth=1
          seisCtx.beginPath(); seisCtx.moveTo(px,y0); seisCtx.lineTo(px,y1); seisCtx.stroke()
          seisCtx.setLineDash([])
          seisCtx.save(); seisCtx.translate(px+2,y0+10); seisCtx.fillStyle='#facc15'; seisCtx.font='11px sans-serif'
          seisCtx.fillText(name, 0, 0); seisCtx.restore()
        }
      }
    }

    // Different color per trace (blue/orange, the repo's diverging accent pair) so two
    // components stay visually distinguishable even where their amplitudes overlap in time.
    function drawTrace(prep, band, color, label){
      const iStart = Math.max(0, Math.floor(vStart*prep.sr))
      const iEnd = Math.min(prep.filtered.length, Math.ceil(vEnd*prep.sr))
      const visible = prep.filtered.slice(iStart, iEnd)
      let amax = 1e-9
      for(const v of visible) amax = Math.max(amax, Math.abs(v))
      const mid = (band.y0+band.y1)/2, halfH = (band.y1-band.y0)/2
      const Y = v => mid - (v/amax)*halfH*0.92
      seisCtx.strokeStyle='#1f2937'; seisCtx.lineWidth=1
      seisCtx.beginPath(); seisCtx.moveTo(x0,mid); seisCtx.lineTo(x1,mid); seisCtx.stroke()
      seisCtx.strokeStyle=color; seisCtx.lineWidth=1.2
      seisCtx.beginPath()
      const pts = decimateMinMax(visible, Math.round(x1-x0))
      pts.forEach((p,k)=>{ const px=X((iStart+p[0])/prep.sr), py=Y(p[1]); k===0?seisCtx.moveTo(px,py):seisCtx.lineTo(px,py) })
      seisCtx.stroke()
      seisCtx.fillStyle=color; seisCtx.font='12px sans-serif'
      seisCtx.fillText(label, 10, band.y0+14)
    }
    drawTrace(p1, bands[0], '#38bdf8', ch1.name)
    if(hasCh2) drawTrace(p2, bands[1], '#fb923c', ch2.name)

    // Continuous time axis: a "nice" tick step (1/2/5 x 10^n) chosen from the visible span so
    // ticks land on round numbers at any zoom level, not just the two window endpoints.
    for(const t of niceTicks(vStart, vEnd, Math.max(2, Math.round((x1-x0)/70)))){
      const px = X(t)
      seisCtx.strokeStyle='#1f2937'; seisCtx.lineWidth=1
      seisCtx.beginPath(); seisCtx.moveTo(px,y0); seisCtx.lineTo(px,y1); seisCtx.stroke()
      seisCtx.strokeStyle='#4b5563'
      seisCtx.beginPath(); seisCtx.moveTo(px,y1); seisCtx.lineTo(px,y1+4); seisCtx.stroke()
      seisCtx.fillStyle='#9ca3af'; seisCtx.font='11px sans-serif'; seisCtx.textAlign='center'
      seisCtx.fillText(t.toFixed(0), px, y1+16)
      seisCtx.textAlign='left'
    }
    seisCtx.strokeStyle='#4b5563'; seisCtx.lineWidth=1
    seisCtx.beginPath(); seisCtx.moveTo(x0,y1); seisCtx.lineTo(x1,y1); seisCtx.stroke()
    seisCtx.fillStyle='#6b7280'; seisCtx.font='11px sans-serif'; seisCtx.textAlign='right'
    seisCtx.fillText('s', x1, y1+16)
    seisCtx.textAlign='left'

    if(vStart>0 || vEnd<fullTMax){
      seisCtx.fillStyle='#6b7280'; seisCtx.font='11px sans-serif'; seisCtx.textAlign='right'
      seisCtx.fillText('drag to zoom · double-click to reset', x1, 14)
      seisCtx.textAlign='left'
    }
  }
  drawSeismogram()

  // Plotly-style click-drag horizontal zoom, double-click to reset -- mirrors the map
  // panels' drag-vs-click disambiguation (a few px of movement means "drag", not "click").
  let seisDragStartPx = null, seisDragMoved = false
  seisCvs.addEventListener('mousedown', e=>{ seisDragStartPx = e.offsetX; seisDragMoved = false })
  seisCvs.addEventListener('mousemove', e=>{
    if(seisDragStartPx===null || !waveform) return
    if(Math.abs(e.offsetX-seisDragStartPx) > 3) seisDragMoved = true
    if(!seisDragMoved) return
    drawSeismogram()
    const xA = Math.min(seisDragStartPx, e.offsetX), xB = Math.max(seisDragStartPx, e.offsetX)
    seisCtx.fillStyle='rgba(250,204,21,0.15)'; seisCtx.fillRect(xA,0,xB-xA,SEISH)
    seisCtx.strokeStyle='rgba(250,204,21,0.7)'; seisCtx.lineWidth=1; seisCtx.strokeRect(xA,0,xB-xA,SEISH)
  })
  seisCvs.addEventListener('mouseup', e=>{
    if(seisDragStartPx===null) return
    const xA = Math.min(seisDragStartPx, e.offsetX), xB = Math.max(seisDragStartPx, e.offsetX)
    seisDragStartPx = null
    if(!seisDragMoved || xB-xA < 5 || seisX1<=seisX0){ seisDragMoved=false; drawSeismogram(); return }
    seisDragMoved = false
    const tA = seisVStart + (xA-seisX0)/(seisX1-seisX0)*(seisVEnd-seisVStart)
    const tB = seisVStart + (xB-seisX0)/(seisX1-seisX0)*(seisVEnd-seisVStart)
    viewStart = Math.max(0, tA); viewEnd = tB
    drawSeismogram()
  })
  seisCvs.addEventListener('dblclick', ()=>{
    viewStart = 0; viewEnd = null
    drawSeismogram()
  })
  window.addEventListener('mouseup', ()=>{ seisDragStartPx = null; seisDragMoved = false })

  updateReadout()
</script>
""")
    end

    const _sex_ready = true
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000002
begin
    _sex_ready
    @bind sex SeismogramExplorerInput()
end

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001a
"""
Every dispatch from the widget that follows an earthquake selection (`st-search`,
`st-select`, `download`) also echoes back that earthquake's own fields (`eqlat`, `eqlon`,
...) alongside its own -- checking `haskey(sex, "eqtime")` directly (rather than
enumerating every action name that carries it) means this stays correct as new actions get
added, instead of quietly reverting to "none yet" the way a forgotten action name once did.
Pluto cells can't be self-referential (there's no way to read this variable's *previous*
value from inside its own definition), so the widget re-sending the earthquake fields on
every action is what lets this stay a plain, stateless, input-only computation instead of
needing a mutable `Ref`.
"""
selected_earthquake_details = if sex isa AbstractDict && sex["action"] == "eq-select"
    (time=sex["time"], lat=sex["lat"], lon=sex["lon"], depth=sex["depth"], mag=sex["mag"])
elseif sex isa AbstractDict && haskey(sex, "eqtime")
    (time=sex["eqtime"], lat=sex["eqlat"], lon=sex["eqlon"], depth=sex["eqdepth"], mag=sex["eqmag"])
else
    "none yet -- search and click an earthquake above"
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000003
md"""
**Selected earthquake**: $(selected_earthquake_details)
"""

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001e
# Same generic-`haskey` shape as `selected_earthquake_details` above, and for the same
# reason: `download` echoes the station back too (as `stcode`/`stlat`/...), so this must not
# be keyed to the `st-select` action name alone.
selected_station_details = if sex isa AbstractDict && sex["action"] == "st-select"
    (code=sex["code"], lat=sex["lat"], lon=sex["lon"], elevation=sex["elevation"])
elseif sex isa AbstractDict && haskey(sex, "stcode")
    (code=sex["stcode"], lat=sex["stlat"], lon=sex["stlon"], elevation=sex["stelevation"])
else
    "none yet -- search and click a station above"
end

# ╔═╡ 6a1e0001-0000-4000-8000-000000000006
md"""
**Selected station**: $(selected_station_details)
"""

# ╔═╡ 6a1e0001-0000-4000-8000-000000000019
let
    if sex isa AbstractDict
        act = sex["action"]
        if act == "eq-search"
            payload = try
                cat = event_client.get_events(
                    starttime=UTCDateTime("1990-01-01"), endtime=UTCDateTime("2030-01-01"),
                    minmagnitude=sex["minmag"],
                    minlatitude=sex["latmin"], maxlatitude=sex["latmax"],
                    minlongitude=sex["lonmin"], maxlongitude=sex["lonmax"],
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
            HTML("""<script>window.dispatchEvent(new CustomEvent('sex-eq-results', {detail: $(repr(payload))}));</script>""")
        elseif act == "st-search"
            payload = if !(selected_earthquake_details isa NamedTuple)
                "{\"ok\":false,\"error\":\"Select an earthquake first.\"}"
            else
                try
                    eq_time_for_search = UTCDateTime(selected_earthquake_details.time)
                    sta = waveform_client.get_stations(network="IU",
                        minlatitude=sex["latmin"], maxlatitude=sex["latmax"],
                        minlongitude=sex["lonmin"], maxlongitude=sex["lonmax"],
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
            HTML("""<script>window.dispatchEvent(new CustomEvent('sex-st-results', {detail: $(repr(payload))}));</script>""")
        elseif act == "download"
            payload = if !(selected_earthquake_details isa NamedTuple) || !(selected_station_details isa NamedTuple)
                "{\"ok\":false,\"error\":\"Select an earthquake and a station first.\"}"
            else
                try
                    eq_time = UTCDateTime(selected_earthquake_details.time)
                    dur = sex["duration"]
                    tr = waveform_client.get_waveforms("IU", selected_station_details.code, "*", "BH?",
                        attach_response=true, starttime=eq_time, endtime=eq_time + dur)
                    # A "*" location wildcard can return several co-located instruments at the same
                    # station (e.g. both "00" and "10"), each contributing its own Z/N/E -- pin down
                    # to one location code so the channel list isn't full of duplicates and rotation
                    # below has one unambiguous horizontal pair to work with.
                    locs = String[]
                    for i in 0:pyconvert(Int, tr.__len__())-1
                        loc = pyconvert(String, tr[i].stats.location)
                        loc in locs || push!(locs, loc)
                    end
                    if length(locs) > 1
                        chosen_loc = "" in locs ? "" : ("00" in locs ? "00" : first(locs))
                        tr = tr.select(location=chosen_loc)
                    end
                    tr = tr.normalize().detrend()
                    n = pyconvert(Int, tr.__len__())
                    parts = String[]
                    chan_lastletters = String[]
                    for i in 0:n-1
                        trace = tr[i]
                        chan = pyconvert(String, trace.stats.channel)
                        sr = pyconvert(Float64, trace.stats.sampling_rate)
                        data = pyconvert(Vector{Float64}, trace.data)
                        dstr = join(round.(data, sigdigits=6), ",")
                        push!(parts, "{\"name\":\"$chan\",\"sr\":$sr,\"data\":[$dstr]}")
                        push!(chan_lastletters, chan[end:end])
                    end
                    # Radial/transverse are geometric rotations of the two horizontal components
                    # toward/across the source-to-station great circle -- only meaningful (and only
                    # computable) once both horizontals are present, so this is a best-effort add-on
                    # to the requested BH? channels, never a reason to fail the download itself.
                    try
                        _, _, baz_py = obspy.geodetics.gps2dist_azimuth(
                            selected_earthquake_details.lat, selected_earthquake_details.lon,
                            selected_station_details.lat, selected_station_details.lon)
                        baz = pyconvert(Float64, baz_py)
                        rt = tr.copy()
                        if !("N" in chan_lastletters && "E" in chan_lastletters)
                            # Horizontals aren't already true-north-aligned (e.g. borehole `--1`/`--2`
                            # orientations) -- rotate to ZNE first using channel orientation metadata.
                            inv = waveform_client.get_stations(network="IU", station=selected_station_details.code,
                                channel="BH?", level="channel", starttime=eq_time, endtime=eq_time + dur)
                            rt.rotate(method="->ZNE", inventory=inv)
                        end
                        rt.rotate(method="NE->RT", back_azimuth=baz)
                        for comp in ("R", "T")
                            sel = rt.select(component=comp)
                            if pyconvert(Int, sel.__len__()) > 0
                                rtrace = sel[0]
                                rchan = pyconvert(String, rtrace.stats.channel)
                                rsr = pyconvert(Float64, rtrace.stats.sampling_rate)
                                rdata = pyconvert(Vector{Float64}, rtrace.data)
                                rdstr = join(round.(rdata, sigdigits=6), ",")
                                push!(parts, "{\"name\":\"$rchan\",\"sr\":$rsr,\"data\":[$rdstr]}")
                            end
                        end
                    catch
                        # No usable horizontal pair (single-component station, missing orientation
                        # metadata, etc.) -- R/T just won't be offered for this station.
                    end
                    "{\"ok\":true,\"channels\":[$(join(parts, ","))]}"
                catch e
                    errmsg = replace(sprint(showerror, e), "\"" => "'", "\n" => " ")
                    "{\"ok\":false,\"error\":\"$errmsg\"}"
                end
            end
            HTML("""<script>window.dispatchEvent(new CustomEvent('sex-waveform', {detail: $(repr(payload))}));</script>""")
        end
    end
end

# ╔═╡ 6a1e0001-0000-4000-8000-00000000001f
md"### Layer 4: Phase Arrivals"

# ╔═╡ 95cc2dde-244a-407b-b4af-2fa5e9658926
# Reactively depends on both `selected_earthquake_details` and `selected_station_details` --
# fires the moment *either* changes, so the widget's phase-arrival overlay is ready before
# the waveform is even downloaded. Only needs depth + epicentral distance (via
# `obspy.geodetics`/TauP, no JS equivalent exists), so unlike the waveform itself this never
# needs to be re-fetched for a channel/duration/frequency change.
let
    if selected_earthquake_details isa NamedTuple && selected_station_details isa NamedTuple
        distance_m, az, baz = obspy.geodetics.gps2dist_azimuth(
            selected_earthquake_details.lat, selected_earthquake_details.lon,
            selected_station_details.lat, selected_station_details.lon)
        dist_deg = pyconvert(Float64, obspy.geodetics.kilometer2degrees(distance_m / 1000.0))
        arr = taup_model.get_ray_paths(source_depth_in_km=selected_earthquake_details.depth,
            distance_in_degree=dist_deg)
        n = pyconvert(Int, arr.__len__())
        parts = String[]
        for i in 0:n-1
            a = arr[i]
            nm = pyconvert(String, a.name)
            t = pyconvert(Float64, a.time)
            push!(parts, "[\"$nm\",$(round(t, digits=2))]")
        end
        payload = "{\"epiDistDeg\":$(round(dist_deg, digits=3)),\"arrivals\":[$(join(parts, ","))]}"
        HTML("""<script>window.dispatchEvent(new CustomEvent('sex-arrivals', {detail: $(repr(payload))}));</script>""")
    end
end


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"

[compat]
CondaPkg = "~0.2.36"
PlutoUI = "~0.7.83"
PythonCall = "~0.9.35"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "e0199d6bee1fd7c839bdff0edfab18a187098f19"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

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

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

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
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "65979512c25a0727f050e6e4be40f0fd9ec893f7"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.7.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

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
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MicroMamba]]
deps = ["Pkg", "Scratch", "micromamba_jll"]
git-tree-sha1 = "535656ce55266bfed0575cd051acc4f36dc869a0"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.15"

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

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "3de8f5e6e90ebfa8d6d1f86997d6cdcd6a912ff3"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.7"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

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

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
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

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "a99e557661bfcee04af1ba688ab6c211b25327f9"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.3"

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
# ╟─6a1e0001-0000-4000-8000-000000000006
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
# ╠═1294ec0f-34e3-4f9e-9e1a-0d86c6bb7674
# ╠═6a1e0001-0000-4000-8000-000000000013
# ╟─6a1e0001-0000-4000-8000-000000000014
# ╠═6a1e0001-0000-4000-8000-000000000015
# ╠═6a1e0001-0000-4000-8000-000000000016
# ╟─6a1e0001-0000-4000-8000-000000000017
# ╠═6a1e0001-0000-4000-8000-000000000018
# ╠═6a1e0001-0000-4000-8000-000000000019
# ╠═6a1e0001-0000-4000-8000-00000000001a
# ╠═6a1e0001-0000-4000-8000-00000000001e
# ╟─6a1e0001-0000-4000-8000-00000000001f
# ╠═95cc2dde-244a-407b-b4af-2fa5e9658926
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
