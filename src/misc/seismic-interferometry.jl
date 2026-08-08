### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Seismic Interferometry"
#> layout = "layout.jlhtml"
#> tags = ["misc"]

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

# ╔═╡ 290b50eb-0faf-4a5f-86bf-36628e61ff06
begin
    using PlutoUI, Statistics, PlutoPlotly, FourierTools, FFTW, StatsBase, MLUtils, ColorSchemes
    import PlutoUI: combine
end

# ╔═╡ ec238280-1edc-442f-837e-d3c807ea5fc6
using LinearAlgebra

# ╔═╡ 7ede27ab-f524-4f20-8bcb-5be70600c892
using LaTeXStrings

# ╔═╡ f986cb42-7cbe-11ee-22de-cbc3714ed81d
TableOfContents(include_definitions=true)

# ╔═╡ bd04a8af-7666-47b7-87cc-a212ec8adcbd
md"""
## Seismic Interferometry Demo
The notebook is an interactive sandbox for seismic interferometry with ambient noise sources. The main source/receiver canvas lets you control the geometry: two receivers sit on the horizontal axis, and you can place, clear, drag, or preset distributions of sources around them. The controls change receiver spacing d, wave period T, bandwidth, source distance scale, spray width, and a stationary-zone threshold. As these change, the notebook regenerates synthetic seismograms, cross-correlates receiver records, and averages selected source contributions.
The top comparison panel shows the recovered interstation response against the ideal “true” virtual-source response. The orange trace is the averaged cross-correlation, and the blue trace is the true reference response. The legend labels are clickable, so you can inspect either trace alone or both together. Markers show zero lag and expected positive/negative arrivals, making it easy to see whether the causal and acausal branches line up with the physical receiver travel time.
The source colors show which sources contribute strongly to the recovered response. With the stationary-zone threshold, weak or geometrically poor contributors are faded out, while stronger contributors are colored by contribution strength. This helps demonstrate the central idea of interferometry: not every noise source helps equally, and source geometry controls whether the cross-correlation reconstructs a clean Green’s-function-like response.
The metrics panel quantifies the two branches separately. For both causal and acausal arrivals, it estimates group velocity U, phase velocity c, phase mismatch Δφ, three-station residuals R3U and R3c, and a mode-consistency residual Rmode. Invalid or too-weak branches are shown as --, using the relative amplitude gate we added so tiny tails or numerical leakage do not produce misleading numbers.
There is also an optional phase correction mode. When enabled, the notebook estimates the phase mismatch between the averaged trace and the reference trace near the expected arrival windows, then rotates valid branches accordingly before reporting measurements. If a branch is too weak, it is not phase-corrected and its phase-related output stays invalid, which keeps the UI honest rather than inventing metrics from noise.
Overall, the notebook is meant to make the full interferometry workflow visible: choose source geometry, generate receiver records, cross-correlate, average, compare to the ideal response, and evaluate whether the recovered causal/acausal branches are physically meaningful.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India

"""

# ╔═╡ 4a592f7e-fa42-4724-8281-1534333816d0
md"---"

# ╔═╡ b8671b46-d45b-4ec7-99f2-b1a2cb0aee12
function average_selected_cross_correlations(cross, contributions, threshold)
    nlag, nsrc = size(cross)
    avg = zeros(Float32, nlag, 1)
    nselected = 0
    use_threshold = threshold > 0 && length(contributions) == nsrc

    @inbounds for j in 1:nsrc
        if !use_threshold || contributions[j] > threshold
            nselected += 1
            for i in 1:nlag
                avg[i, 1] += cross[i, j]
            end
        end
    end

    if nselected > 0
        invn = 1 / nselected
        @inbounds for i in 1:nlag
            avg[i, 1] *= invn
        end
    end

    return avg
end

# ╔═╡ 1bc0ce25-4de0-4628-a308-8106544fc0ec
md"""
## Cross-Correlation
"""

# ╔═╡ 5427df02-9278-41c9-9b3b-b5d1ab5b4c24
stationary_phase_source_distance = 1000.0

# ╔═╡ e0ae191c-b674-46c7-bf05-bfecade313a5
md"""
## Generate Wavelets
"""

# ╔═╡ 9a74489f-8e6d-45ce-ab57-464ee36d3316
function delaywav!(out, w, t, dt)
    fill!(out, 0)
    sample_center = t / dt + 1
    wavelet_center = (length(w) + 1) / 2
    @inbounds for i in eachindex(out)
        j = i - sample_center + wavelet_center
        j0 = floor(Int, j)
        a = j - j0
        if 1 <= j0 < length(w)
            out[i] = (1 - a) * w[j0] + a * w[j0 + 1]
        elseif j0 == length(w)
            out[i] = w[j0]
        end
    end
    return out
end

# ╔═╡ f0e27c51-233b-45ef-be1e-34b624b5a4b2
function delaywav(w, t, nt, dt)
    y = zeros(Float32, nt)
    delaywav!(y, w, t, dt)
    return y
end

# ╔═╡ 577d411c-88b5-495d-a49a-d2379336f533
md"""
## Generate Random Source Locations
"""

# ╔═╡ 8d45d812-fb70-45c5-a123-4f93844fef6d
md"""
## Receivers
"""

# ╔═╡ 3b144853-0d82-46e2-a4ee-301516f3ce10
vel = 4

# ╔═╡ 4e20077c-3687-47a1-8ab0-048b99fb0a46
"""
    get_traveltime(rec, src)

Calculate the travel time of a seismic wave between a receiver and a source.

# Arguments
- `rec::Tuple{Float64, Float64}`: Coordinates of the receiver (x, y).
- `src::Tuple{Float64, Float64}`: Coordinates of the source (x, y).

# Returns
- `tm::Float64`: Travel time of the seismic wave.

# Notes
- The function assumes a constant velocity `vel` which should be defined in the scope where this function is used.
"""
function get_traveltime(rec, src, velocity)
    dis = sqrt((rec[1] - src[1])^2 + (rec[2] - src[2])^2)
    tm = dis / velocity
    return tm
end

# ╔═╡ e04a82d2-9eda-419f-a6b3-0af28431e554
get_traveltime(rec, src) = get_traveltime(rec, src, vel)

# ╔═╡ 1ce924ce-5a35-4a84-835a-3bdf16b876dc
"""
    generate_seismograms(rec1, rec2, srcloc)

Generate seismograms for two receivers from multiple source locations.

# Arguments
- `rec1`: The first receiver location.
- `rec2`: The second receiver location.
- `srcloc`: An array of source locations.

# Returns
- `wav1`: A matrix of delayed waveforms for the first receiver.
- `wav2`: A matrix of delayed waveforms for the second receiver.

# Description
This function computes the travel times from each source location to two receivers, generates delayed waveforms based on these travel times, and stacks the waveforms into matrices for each receiver.
"""
function generate_seismograms!(seis1, seis2, rec1, rec2, srcloc, wavelet, nt, dt, velocity)
    @views for i in eachindex(srcloc)
        delaywav!(seis1[:, i], wavelet, get_traveltime(rec1, srcloc[i], velocity), dt)
        delaywav!(seis2[:, i], wavelet, get_traveltime(rec2, srcloc[i], velocity), dt)
    end
    return seis1, seis2
end

# ╔═╡ 32be7b92-0b34-4ad0-8968-c77a13225383
function generate_seismograms(rec1, rec2, srcloc, wavelet, nt, dt, velocity)
    nsrc = length(srcloc)
    seis1 = zeros(Float32, nt, nsrc)
    seis2 = zeros(Float32, nt, nsrc)
    generate_seismograms!(seis1, seis2, rec1, rec2, srcloc, wavelet, nt, dt, velocity)
    return seis1, seis2
end

# ╔═╡ 2fd6cc16-519d-4d0f-ba42-874e58a32124
function normalized_inner_product_range(a, b, idxs)
    adotb = 0.0
    anorm2 = 0.0
    bnorm2 = 0.0
    @inbounds for i in idxs
        av = a[i]
        bv = b[i]
        adotb += av * bv
        anorm2 += abs2(av)
        bnorm2 += abs2(bv)
    end
    anorm2 > 0 && bnorm2 > 0 ? adotb / sqrt(anorm2 * bnorm2) : 0.0
end

# ╔═╡ e4f79b31-aa1f-4132-bc8e-b30d79380a3b
function source_contributions(cross, ref, tgrid)
    nsrc = size(cross, 2)
    contributions = zeros(Float64, nsrc)
    nsrc == 0 && return contributions

    causal = searchsortedfirst(tgrid, 0.0):length(tgrid)
    acausal = 1:searchsortedlast(tgrid, 0.0)

    @views for i in 1:nsrc
        c = cross[:, i]
        score_causal = normalized_inner_product_range(c, ref, causal)
        score_acausal = normalized_inner_product_range(c, ref, acausal)
        contributions[i] = max(score_causal, score_acausal, 0.0)
    end

    return contributions
end

# ╔═╡ d483d881-11c4-495f-a103-2d2dd6371ca3
md"## Appendix"

# ╔═╡ e059662c-69af-460c-aada-408a365d0ff7
default_plotly_template(:plotly_dark)

# ╔═╡ 366783f3-3a48-4e0d-a850-6c8cb7275377
function gaussian_wavelet(f0, bw_perc, dt)
    sigma_t = sqrt(log(2)) / (2 * π * f0 * (bw_perc / 100))
    N = max(51, 2 * ceil(Int, 6 * sigma_t / dt) + 1)
    t = dt .* (0:N-1) .- dt * (N ÷ 2)
    w = @. exp(-t^2 / (2 * sigma_t^2)) * cos(2 * π * f0 * t)
    w ./= maximum(abs, w)
    return Float32.(w)
end

# ╔═╡ 6a8aeed2-173a-4523-afb3-f9d77c252816
nt = 1028

# ╔═╡ b2c0af32-b19d-47b0-86f3-66e8fd7d7269
dt = 0.25

# ╔═╡ dafaabbb-4263-41ce-b17f-0634adc01b8f
delaywav(w, t) = delaywav(w, t, nt, dt)

# ╔═╡ b11959f8-3cf2-4488-bbc1-6597220ea599
gaussian_wavelet(f0, bw_perc) = gaussian_wavelet(f0, bw_perc, dt)

# ╔═╡ af973cc3-21b3-4dca-91ad-0c310e1214bc
tgrid_xcorr = dt .* (-(nt - 1):(nt - 1))

# ╔═╡ 09bd50ae-ab80-4cde-bdb7-ccc7cd7041d4
tgrid = dt .* (0:(nt - 1))

# ╔═╡ d9727185-8f69-49a5-a98f-7514562529ab
md"""
### Plots
"""

# ╔═╡ 1293b18a-15d5-4d54-8396-59899c89a172
function seis_heatmap(tgrid, r, title, ytitle, xtitle)
    m = maximum(abs, r)
    plot(heatmap(y=tgrid, z=r, colorscale="seismic", zmin=-m, zmax=m), Layout(title=title, yaxis_autorange="reversed", height=350, width=600, yaxis=attr(title=ytitle), xaxis=attr(title=xtitle)))
end

# ╔═╡ 8b213ed6-b4ec-44a9-8e3f-8ad03ee9f270
function randobs_safe(a, n)
    size(a, 2) == 0 ? zeros(eltype(a), size(a, 1), 1) : randobs(a, min(n, size(a, 2)))
end

# ╔═╡ ccf893d9-5b6e-4a79-b09d-3da957bc612c
function plot_line(tgrid, tr; title="", names=fill(" ", length(tr)))

    fig = PlutoPlotly.Plot(Layout(height=300, width=600,
        title=title, xaxis=attr(title="cross-correlation lag (s)"), yaxis=attr(title="amplitude"), font=attr(
            size=10), legend=attr(yanchor="bottom", y=-0.75, xanchor="left", x=0.0)))
    opacities = [1, 0.5]
    widths = [1, 0.5]
    for (T, width, opacity, name) in zip(tr, opacities, widths, names)
        add_trace!(fig, PlutoPlotly.scatter(x=tgrid, y=T, width=width, opacity=opacity, mode="lines", name=name, showlegend=true))
    end
    PlutoPlotly.plot(fig)
end

# ╔═╡ 34f81006-dead-4687-aff1-ed2753f27469
function xcorr(x, y; padmode=:longest)
    n = length(x) + length(y) - 1
    X = fft(vcat(x, zeros(eltype(x), length(y) - 1)))
    Y = fft(vcat(y, zeros(eltype(y), length(x) - 1)))
    return fftshift(real.(ifft(conj(X) .* Y)))
end

# ╔═╡ 4ce9cc5d-6808-4977-98bb-8e9cacf171cc
function true_virtual_response_for_traveltime(t0, wavelet, tgrid, dt)
    nlag = length(tgrid)
    lag0 = first(tgrid)
    acorr_wav = xcorr(wavelet, wavelet)     # length 2N-1; index length(wavelet) = lag 0
    wavelet_center = length(wavelet)
    ref = zeros(Float32, nlag)
    for τ in (t0, -t0)
        lag_center = (τ - lag0) / dt + 1
        for i in eachindex(ref)
            j = i - lag_center + wavelet_center
            j0 = floor(Int, j)
            a = j - j0
            if 1 <= j0 < length(acorr_wav)
                ref[i] += (1 - a) * acorr_wav[j0] + a * acorr_wav[j0 + 1]
            elseif j0 == length(acorr_wav)
                ref[i] += acorr_wav[j0]
            end
        end
    end
    ref
end

# ╔═╡ 15a556df-3ed3-4068-86ca-d3c130cb81d9
function xcorr!(out, x, y)
    out .= xcorr(x, y)
    return out
end

# ╔═╡ cd3a586a-b737-451e-921f-5463cc286f73
function cross_correlations!(cross, seis1, seis2)
    @views for i in axes(seis1, 2)
        xcorr!(cross[:, i], seis1[:, i], seis2[:, i])
    end
    return cross
end

# ╔═╡ 0d2266c7-312c-4c56-8f7b-8a98fba81fc6
function cross_correlations(seis1, seis2)
    nsrc = size(seis1, 2)
    nlag = size(seis1, 1) + size(seis2, 1) - 1
    cross = zeros(Float32, nlag, nsrc)
    cross_correlations!(cross, seis1, seis2)
    return cross
end

# ╔═╡ 7fe32a18-b484-4179-b0a8-d6c49059ed59
begin
    struct CanvasSourceInput
        default_pts::Vector{Vector{Int}}
        d::Float64
        f0::Float64
        bw_perc::Int
        vel::Float64
        subsample_threshold::Float64
        source_distance_scale::Float64
    end

    function CanvasSourceInput(vel; d=50.0, f0=0.1, bw_perc=40, subsample_threshold=0.0, source_distance_scale=1.0)
        pts = [[round(Int, 320 + 255*cos(a)), round(Int, 320 + 255*sin(a))]
               for a in range(0, 2π, length=501)[1:end-1]]
        CanvasSourceInput(pts, Float64(d), Float64(f0), Int(bw_perc), Float64(vel), Float64(subsample_threshold), Float64(source_distance_scale))
    end

    function colorscheme_stops_js(scheme, n=9)
        stops = map(range(0, 1, length=n)) do x
            c = get(scheme, x)
            rgb = round.(Int, 255 .* (c.r, c.g, c.b))
            "[$(round(x, digits=4)),[$(rgb[1]),$(rgb[2]),$(rgb[3])]]"
        end
        "[" * join(stops, ",") * "]"
    end

    Base.get(w::CanvasSourceInput) = Dict{String,Any}(
        "sources" => w.default_pts,
        "d" => w.d,
        "f0" => w.f0,
        "bw_perc" => w.bw_perc,
        "subsample_threshold" => w.subsample_threshold,
        "source_distance_scale" => w.source_distance_scale
    )

    function Base.show(io::IO, ::MIME"text/html", w::CanvasSourceInput)
        write(io, """
<div id="siwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    #siwidget .si-workspace{display:flex;gap:12px;align-items:flex-start;justify-content:center;width:100%}
    #siwidget .si-metrics{width:260px;min-height:640px;box-sizing:border-box;background:#050505;border:1px solid #374151;border-radius:6px;padding:12px;font:12px/1.35 sans-serif;color:#d1d5db}
    #siwidget .si-metrics h3{margin:0 0 10px 0;font-size:20px;color:#e5e7eb}
    #siwidget .si-metric-card{border:1px solid #1f2937;border-radius:6px;background:#0b0b0b;padding:10px;margin-bottom:10px}
    #siwidget .si-metric-title{font-weight:700;color:#f3f4f6;margin-bottom:7px}
    #siwidget .si-metric-row{display:grid;grid-template-columns:42px 1fr;gap:8px;margin:4px 0}
    #siwidget .si-metric-label{color:#9ca3af}
    #siwidget .si-metric-value{color:#e5e7eb}
    #siwidget .si-percent{color:#93c5fd}
    #siwidget .si-controls{width:min(914px,100%);margin-top:8px;display:grid;grid-template-columns:repeat(2,minmax(300px,1fr));gap:8px;font:12px sans-serif}
    #siwidget .si-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:9px 10px}
    #siwidget .si-control-title{font-weight:700;color:#e5e7eb;margin-bottom:7px;font-size:20px}
    #siwidget .si-control-row{display:grid;grid-template-columns:120px minmax(120px,1fr) 56px;gap:8px;align-items:center;margin:6px 0}
    #siwidget .si-control-row.wide{grid-template-columns:190px minmax(120px,1fr) 56px}
    #siwidget .si-control-row input[type=range]{width:100%;vertical-align:middle}
    #siwidget .si-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums}
    #siwidget .si-actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    #siwidget .si-presets{display:flex;gap:5px;align-items:center;flex-wrap:wrap}
    #siwidget button{border-radius:4px;border:1px solid #6b7280;background:#606060;color:#f3f4f6;padding:4px 9px}
    @media (max-width: 980px){
      #siwidget .si-workspace{flex-direction:column;align-items:center}
      #siwidget .si-metrics{width:640px;max-width:100%;min-height:0}
      #siwidget .si-controls{grid-template-columns:1fr;width:640px;max-width:100%}
    }
  </style>
  <canvas id="cmpcvs" width="914" height="140"
    style="cursor:default;background:#000;border:1px solid #374151;border-radius:4px;display:block;margin-bottom:4px"></canvas>
  <div class="si-workspace">
    <canvas id="srccvs" width="640" height="640"
      style="cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block"></canvas>
    <div id="metricspanel" class="si-metrics">
      <h3>Branch metrics</h3>
      <div style="color:#9ca3af">Waiting for comparison data...</div>
    </div>
  </div>
  <div class="si-controls">
    <div class="si-control-group">
      <div class="si-control-title">Wave / Geometry</div>
      <label class="si-control-row"><span>d (km)</span><input type="range" id="dist" min="10" max="100" step="5" value="$(w.d)"><span id="distv" class="si-value">$(round(w.d, digits=1))</span></label>
      <label class="si-control-row"><span>T (s)</span><input type="range" id="period" min="0" max="$(round(log10(25), digits=6))" step="0.001" value="$(round(log10(1/w.f0), digits=6))"><span id="periodv" class="si-value">$(round(Int, 1/w.f0)) s</span></label>
      <label class="si-control-row"><span>Bandwidth</span><input type="range" id="bw" min="5" max="100" step="5" value="$(w.bw_perc)"><span id="bwv" class="si-value">$(w.bw_perc)%</span></label>
      <label class="si-control-row"><span>Phase branch</span><input type="range" id="phasebranch" min="-20" max="20" step="1" value="0"><span id="phasebranchv" class="si-value">0</span></label>
    </div>
    <div class="si-control-group">
      <div class="si-control-title">Source distribution</div>
      <label class="si-control-row wide"><span>Coherent source subsampling</span><input type="range" id="subsamp" min="0" max="1" step="0.01" value="$(w.subsample_threshold)"><span id="subsampv" class="si-value">$(round(w.subsample_threshold, digits=2))</span></label>
      <label class="si-control-row wide"><span>Source distance scale</span><input type="range" id="srcscale" min="1" max="4" step="1" value="$(round(Int, w.source_distance_scale))"><span id="srcscalev" class="si-value">$(round(Int, w.source_distance_scale))x</span></label>
      <label class="si-control-row wide"><span>Spray width</span><input type="range" id="spraywidth" min="0" max="15" step="1" value="15"><span id="spraywidthv" class="si-value">15 km</span></label>
    </div>
    <div class="si-control-group">
      <div class="si-control-title">Actions</div>
      <div class="si-actions">
        <label><input type="checkbox" id="phasecorr" style="vertical-align:middle"> Apply phase correction</label>
        <label><input type="checkbox" id="dragsources" style="vertical-align:middle"> Drag sources</label>
        <button id="clrbtn">Clear</button>
        <span id="cnt" class="si-value">Sources: $(length(w.default_pts))</span>
      </div>
    </div>
    <div class="si-control-group">
      <div class="si-control-title">Presets</div>
      <div class="si-presets">
        <button class="preset" data-preset="circular">Circular</button>
        <button class="preset" data-preset="sectional">Sectional</button>
        <button class="preset" data-preset="leftarc">Left Arc</button>
        <button class="preset" data-preset="rightarc">Right Arc</button>
        <button class="preset" data-preset="special">Special</button>
        <button class="preset" data-preset="bimodal">Twin Beams</button>
        <button class="preset" data-preset="uniformdisk">Uniform Disk</button>
        <button class="preset" data-preset="sourcecloud">Source Cloud</button>
        <button class="preset" data-preset="twoclouds">Two Clouds</button>
      </div>
    </div>
  </div>
</div>
<script>
  const PIX=640, CMP_W=914, CMP_H=140, MID=320, SCALE=1.5, DEFAULT_SPREAD=15, PER_PT=15, RAD=200, RING=170, N=1000
  const DPR=Math.min(window.devicePixelRatio || 1, 2)
  const VEL=$(w.vel), DT=$(dt)
  const T_MIN=1.0, T_MAX=25.0
  const COLOR_STOPS=$(colorscheme_stops_js(ColorSchemes.viridis))
  const STATE_KEY='seismic_interferometry_canvas_v6'
  const par = currentScript.previousElementSibling
  const cvs = par.querySelector('#srccvs')
  const ctx = cvs.getContext('2d')
  const cmpcvs = par.querySelector('#cmpcvs')
  const cmpctx = cmpcvs.getContext('2d')
  const metricsPanel = par.querySelector('#metricspanel')
  const lbl = par.querySelector('#cnt')

  function setupHiDPICanvas(canvas, context, width, height) {
    canvas.width = Math.round(width * DPR)
    canvas.height = Math.round(height * DPR)
    canvas.style.width = width + 'px'
    canvas.style.height = height + 'px'
    context.setTransform(DPR, 0, 0, DPR, 0, 0)
  }
  setupHiDPICanvas(cvs, ctx, PIX, PIX)
  setupHiDPICanvas(cmpcvs, cmpctx, CMP_W, CMP_H)

  const toP = (d, flip) => MID + (flip ? -1 : 1) * d * SCALE

  function makeArc(a1, a2, n) {
    const out=[]
    for(let i=0;i<n;i++){
      const a=a1+(a2-a1)*i/Math.max(n-1,1)
      out.push([Math.round(MID+RING*SCALE*Math.cos(a)),Math.round(MID+RING*SCALE*Math.sin(a))])
    }
    return out
  }
  function makeUniformDisk(n) {
    const out=[]
    for(let i=0;i<n;i++){
      const r=RAD*SCALE*Math.sqrt(Math.random())
      const a=2*Math.PI*Math.random()
      out.push([Math.round(MID+r*Math.cos(a)),Math.round(MID+r*Math.sin(a))])
    }
    return out
  }
  function clampToSourceCircle(px, py) {
    const dx=px-MID, dy=py-MID
    const maxr=RAD*SCALE
    const r=Math.sqrt(dx*dx+dy*dy)
    if(r <= maxr) return [Math.round(px), Math.round(py)]
    const s=maxr/r
    return [Math.round(MID+dx*s), Math.round(MID+dy*s)]
  }
  function randn() {
    let u=0, v=0
    while(u === 0) u = Math.random()
    while(v === 0) v = Math.random()
    return Math.sqrt(-2*Math.log(u)) * Math.cos(2*Math.PI*v)
  }
  function makeSourceCloud(n, cx, cy, sigmaKm) {
    const out=[]
    for(let i=0;i<n;i++){
      const px=MID + (cx + sigmaKm*randn()) * SCALE
      const py=MID - (cy + sigmaKm*randn()) * SCALE
      out.push(clampToSourceCircle(px, py))
    }
    return out
  }
  function makeTwoClouds(n) {
    const n1=Math.floor(n/2)
    return makeSourceCloud(n1, 135, 52, 12).concat(makeSourceCloud(n-n1, 145, -54, 12))
  }
  const PRESETS = {
    circular:  makeArc(0, 2*Math.PI, N),
    sectional: makeArc(-0.1*Math.PI,0.1*Math.PI,N/2).concat(makeArc(0.9*Math.PI,1.1*Math.PI,N/2)),
    leftarc:   makeArc(Math.PI/3,-Math.PI/3,N),
    rightarc:  makeArc(2*Math.PI/3,4*Math.PI/3,N),
    special:   makeArc(Math.PI/3,2*Math.PI/3,N/2).concat(makeArc(4*Math.PI/3,5*Math.PI/3,N/2)),
    bimodal:   makeArc(0,Math.PI/6,N/2).concat(makeArc(Math.PI,7*Math.PI/6,N/2)),
  }

  // Restore persisted state or start fresh
  let saved = null
  try { saved = JSON.parse(sessionStorage.getItem(STATE_KEY)) } catch(e) {}
  let pts = saved ? saved.pts : PRESETS.circular.slice()
  let d = saved && Number.isFinite(saved.d) ? saved.d : $(w.d)
  let period = saved && Number.isFinite(saved.period) ? saved.period : $(round(Int, 1/w.f0))
  let bwPerc = saved && Number.isFinite(saved.bw_perc) ? saved.bw_perc : $(w.bw_perc)
  let phaseBranch = saved && Number.isFinite(saved.phase_branch) ? Math.round(saved.phase_branch) : 0
  let subsampThreshold = saved && Number.isFinite(saved.subsamp_threshold) ? saved.subsamp_threshold : $(w.subsample_threshold)
  let sourceScale = saved && Number.isFinite(saved.source_distance_scale) ? saved.source_distance_scale : $(w.source_distance_scale)
  let sprayWidth = saved && Number.isFinite(saved.spray_width) ? saved.spray_width : DEFAULT_SPREAD
  let applyPhaseCorrection = saved && typeof saved.apply_phase_correction === 'boolean' ? saved.apply_phase_correction : false
  let dragSourcesMode = saved && typeof saved.drag_sources === 'boolean' ? saved.drag_sources : false
  let showAverageTrace = true
  let showReferenceTrace = true
  let comparisonLegendHits = []
  let currentWeights = (window._seismic_weights && window._seismic_weights.length) ? window._seismic_weights : []
  let currentComparison = window._seismic_comparison || null

  let drawing = false
  let draggingReceiver = false
  let draggingSource = false
  let selectedSourceIndices = new Set()
  let boxSelecting = false
  let selectionStart = null
  let selectionRect = null
  let dragStart = null
  let dragInitialPositions = []

  const r1 = () => [-d/2, 0]
  const r2 = () => [ d/2, 0]
  const receiverPixel = r => [toP(r[0],false), toP(r[1],true)]
  const clampPeriod = T => Math.max(T_MIN, Math.min(T_MAX, T))
  const periodToSlider = T => Math.log10(clampPeriod(T))
  const sliderToPeriod = x => clampPeriod(Math.pow(10, x))
  const formatPeriod = T => T < 10 ? T.toFixed(2) : T.toFixed(1)

  function enforcePeriodConstraint() {
    period = clampPeriod(period)
  }

  function syncControls() {
    enforcePeriodConstraint()
    par.querySelector('#dist').value = d
    par.querySelector('#distv').textContent = d.toFixed(1)
    par.querySelector('#period').value = periodToSlider(period)
    par.querySelector('#periodv').textContent = formatPeriod(period) + ' s'
    par.querySelector('#bw').value = bwPerc
    par.querySelector('#bwv').textContent = bwPerc + '%'
    par.querySelector('#phasebranch').value = phaseBranch
    par.querySelector('#phasebranchv').textContent = phaseBranch.toFixed(0)
    par.querySelector('#subsamp').value = subsampThreshold
    par.querySelector('#subsampv').textContent = subsampThreshold.toFixed(2)
    par.querySelector('#srcscale').value = sourceScale
    par.querySelector('#srcscalev').textContent = sourceScale.toFixed(0) + 'x'
    par.querySelector('#spraywidth').value = sprayWidth
    par.querySelector('#spraywidthv').textContent = sprayWidth === 0 ? 'single' : sprayWidth.toFixed(0) + ' km'
    par.querySelector('#phasecorr').checked = applyPhaseCorrection
    par.querySelector('#dragsources').checked = dragSourcesMode
  }

  function drawAxes() {
    ctx.strokeStyle='#4b5563'; ctx.lineWidth=0.8
    ctx.beginPath(); ctx.moveTo(MID,8); ctx.lineTo(MID,PIX-8); ctx.stroke()
    ctx.beginPath(); ctx.moveTo(8,MID); ctx.lineTo(PIX-8,MID); ctx.stroke()
    ctx.fillStyle='#6b7280'; ctx.font='12px sans-serif'; ctx.textAlign='center'
    for(let v=-200;v<=200;v+=100){
      if(v===0) continue
      const px=toP(v,false), py=toP(v,true)
      ctx.strokeStyle='#4b5563'; ctx.lineWidth=0.5
      ctx.beginPath(); ctx.moveTo(px,MID-4); ctx.lineTo(px,MID+4); ctx.stroke()
      ctx.beginPath(); ctx.moveTo(MID-4,py); ctx.lineTo(MID+4,py); ctx.stroke()
      ctx.fillStyle='#6b7280'
      ctx.fillText(v,px,MID+20)
      ctx.textAlign='right'; ctx.fillText(v,MID-8,py+4); ctx.textAlign='center'
    }
    ctx.fillStyle='#9ca3af'; ctx.font='11px sans-serif'
    ctx.fillText('x',PIX-10,MID-6)
    ctx.textAlign='left'; ctx.fillText('y',MID+6,16); ctx.textAlign='center'
  }

  function drawCanvasTitle() {
    ctx.fillStyle='#e5e7eb'
    ctx.font='bold 20px sans-serif'
    ctx.textAlign='left'
    ctx.fillText('Interferometry experiment setup',12,25)
    ctx.textAlign='center'
  }

  function drawReceiver(r, color) {
    const px=toP(r[0],false), py=toP(r[1],true)
    ctx.fillStyle=color
    ctx.beginPath(); ctx.moveTo(px,py-15); ctx.lineTo(px-10,py+8); ctx.lineTo(px+10,py+8); ctx.closePath(); ctx.fill()
  }

  function drawInfo() {
    const lam = VEL * period
    const dol = d / lam
    ctx.fillStyle='rgba(0,0,0,0.82)'
    ctx.fillRect(6,PIX-38,400,30)
    ctx.fillStyle='#e5e7eb'; ctx.font='13px monospace'; ctx.textAlign='left'
    ctx.fillText('T='+formatPeriod(period)+' s  λ='+lam.toFixed(1)+' km  d='+d.toFixed(1)+' km  d/λ='+dol.toFixed(2)+'  src×='+sourceScale.toFixed(0),12,PIX-18)
    ctx.textAlign='center'
  }

  function colorSchemeScale(w, alpha=0.85) {
    const x = Math.max(0, Math.min(1, w))
    let lo = COLOR_STOPS[0], hi = COLOR_STOPS[COLOR_STOPS.length-1]
    for(let i=1;i<COLOR_STOPS.length;i++){
      if(x <= COLOR_STOPS[i][0]) { lo = COLOR_STOPS[i-1]; hi = COLOR_STOPS[i]; break }
    }
    const a = hi[0] === lo[0] ? 0 : (x - lo[0]) / (hi[0] - lo[0])
    const r = Math.round(lo[1][0] + a * (hi[1][0] - lo[1][0]))
    const g = Math.round(lo[1][1] + a * (hi[1][1] - lo[1][1]))
    const b = Math.round(lo[1][2] + a * (hi[1][2] - lo[1][2]))
    return 'rgba('+r+','+g+','+b+','+alpha+')'
  }

  function fmtVelocity(v) {
    return Number.isFinite(v) ? v.toFixed(2) + ' km/s' : '--'
  }
  function fmtError(v) {
    return Number.isFinite(v) ? '<span class="si-percent">' + v.toFixed(1) + '%</span>' : '--'
  }
  function fmtResidual(v, pct) {
    return Number.isFinite(v) && Number.isFinite(pct) ? v.toFixed(2) + ' km/s, <span class="si-percent">' + pct.toFixed(1) + '%</span>' : '--'
  }
  function fmtPhase(v) {
    if(!Number.isFinite(v)) return '--'
    if(Math.abs(v) < 0.005) return '0'
    const sign = v < 0 ? '-' : ''
    const denom = 1 / Math.abs(v)
    if(Math.abs(denom - 1) < 0.05) return sign + 'π'
    return sign + 'π/' + denom.toFixed(2)
  }
  function fmtModeResidual(v, pct) {
    return Number.isFinite(v) && Number.isFinite(pct) ? v.toFixed(2) + ' km/s, <span class="si-percent">' + pct.toFixed(1) + '%</span>' : '--'
  }
  function metricCard(title, U, Uerr, c, cerr, phasePi, r3U, r3Upct, r3c, r3cpct, rmode, rmodePct) {
    return '<div class="si-metric-card">' +
      '<div class="si-metric-title">' + title + '</div>' +
      '<div class="si-metric-row"><span class="si-metric-label">U</span><span class="si-metric-value">' + fmtVelocity(U) + ', err ' + fmtError(Uerr) + '</span></div>' +
      '<div class="si-metric-row"><span class="si-metric-label">c</span><span class="si-metric-value">' + fmtVelocity(c) + ', err ' + fmtError(cerr) + '</span></div>' +
      '<div class="si-metric-row"><span class="si-metric-label">Δφ</span><span class="si-metric-value">' + fmtPhase(phasePi) + '</span></div>' +
      '<div class="si-metric-row"><span class="si-metric-label">R3U</span><span class="si-metric-value">' + fmtResidual(r3U, r3Upct) + '</span></div>' +
      '<div class="si-metric-row"><span class="si-metric-label">R3c</span><span class="si-metric-value">' + fmtResidual(r3c, r3cpct) + '</span></div>' +
      '<div class="si-metric-row"><span class="si-metric-label">Rmode</span><span class="si-metric-value">' + fmtModeResidual(rmode, rmodePct) + '</span></div>' +
      '</div>'
  }
  function renderMetricsPanel() {
    if(!currentComparison || !currentComparison.velocity){
      metricsPanel.innerHTML = '<h3>Branch metrics</h3><div style="color:#9ca3af">Waiting for comparison data...</div>'
      return
    }
    const velocity = currentComparison.velocity || {}
    const phaseText = velocity.phase_corrected ? '<div style="color:#a7f3d0;margin-bottom:8px">phase corrected</div>' : ''
    metricsPanel.innerHTML = '<h3>Branch metrics</h3>' + phaseText +
      metricCard('Causal', velocity.U_causal, velocity.U_causal_error_pct, velocity.c_causal, velocity.c_causal_error_pct, velocity.phase_causal_pi, velocity.R3_U_causal, velocity.R3_U_causal_pct, velocity.R3_c_causal, velocity.R3_c_causal_pct, velocity.R_mode_causal, velocity.R_mode_causal_pct) +
      metricCard('Acausal', velocity.U_acausal, velocity.U_acausal_error_pct, velocity.c_acausal, velocity.c_acausal_error_pct, velocity.phase_acausal_pi, velocity.R3_U_acausal, velocity.R3_U_acausal_pct, velocity.R3_c_acausal, velocity.R3_c_acausal_pct, velocity.R_mode_acausal, velocity.R_mode_acausal_pct)
  }

  function drawComparisonPanel() {
    const cw=CMP_W, ch=CMP_H, padL=48, padR=14, padT=42, padB=24
    comparisonLegendHits = []
    cmpctx.fillStyle='#000'
    cmpctx.fillRect(0,0,cw,ch)
    cmpctx.fillStyle='#e5e7eb'; cmpctx.font='bold 18px sans-serif'; cmpctx.textAlign='left'
    cmpctx.fillText('Interstation response comparison',10,23)

    if(!currentComparison || !currentComparison.avg || !currentComparison.ref || currentComparison.avg.length < 2){
      cmpctx.fillStyle='#9ca3af'; cmpctx.font='11px sans-serif'
      cmpctx.fillText('Waiting for comparison data…',padL,ch/2+4)
      return
    }

    const avg=currentComparison.avg, ref=currentComparison.ref, t=currentComparison.t || []
    const t0=currentComparison.t0 || 0
    const velocity=currentComparison.velocity || {}
    const n=Math.min(avg.length, ref.length, t.length)
    if(velocity.phase_corrected){
      cmpctx.font='10px sans-serif'
      cmpctx.fillStyle='#a7f3d0'
      cmpctx.fillText('phase corrected', 10, 28)
    }

    // Fixed x-axis: ±3·t0 (at least 5 s so axis is never degenerate)
    const xlim = Math.max(3 * t0, 5)
    const ptx = tv => padL + (cw-padL-padR) * (tv + xlim) / (2*xlim)

    let maxAbs=0
    for(let i=0;i<n;i++){
      const a=showAverageTrace ? Math.abs(avg[i]) : 0
      const r=showReferenceTrace ? Math.abs(ref[i]) : 0
      if(Number.isFinite(a) && a>maxAbs) maxAbs=a
      if(Number.isFinite(r) && r>maxAbs) maxAbs=r
    }
    if(maxAbs <= 0) maxAbs = 1

    const py = v => padT+(ch-padT-padB)*(0.5 - 0.45*v/maxAbs)
    const yMid = py(0)

    // Horizontal zero-amplitude line
    cmpctx.strokeStyle='#374151'; cmpctx.lineWidth=0.8
    cmpctx.beginPath(); cmpctx.moveTo(padL,yMid); cmpctx.lineTo(cw-padR,yMid); cmpctx.stroke()

    // Vertical marker lines
    function vline(tv, color, dash) {
      const x = ptx(tv)
      if(x < padL || x > cw-padR) return
      cmpctx.save()
      cmpctx.strokeStyle=color; cmpctx.lineWidth=1
      if(dash) cmpctx.setLineDash([4,3])
      cmpctx.beginPath(); cmpctx.moveTo(x,padT); cmpctx.lineTo(x,ch-padB); cmpctx.stroke()
      cmpctx.setLineDash([])
      cmpctx.restore()
    }
    vline(0,   '#6b7280', false)    // lag = 0
    vline( t0, '#a3a3a3', true)     // +t0
    vline(-t0, '#a3a3a3', true)     // -t0

    // Traces
    function trace(values, color, width) {
      cmpctx.strokeStyle=color; cmpctx.lineWidth=width
      cmpctx.beginPath()
      let started=false
      for(let i=0;i<n;i++){
        const v=values[i]
        if(!Number.isFinite(v)) continue
        const x=ptx(t[i]), y=py(v)
        if(started) cmpctx.lineTo(x,y); else { cmpctx.moveTo(x,y); started=true }
      }
      cmpctx.stroke()
    }
    if(showReferenceTrace) trace(ref,'#38bdf8',3)
    if(showAverageTrace) trace(avg,'#f97316',1)

    // Legend
    cmpctx.font='15px sans-serif'; cmpctx.textAlign='left'
    function legendItem(key, label, x, color, active) {
      cmpctx.fillStyle=active ? color : 'rgba(156,163,175,0.45)'
      cmpctx.fillText(label,x,15)
      const w = cmpctx.measureText(label).width
      comparisonLegendHits.push({key:key, x0:x-4, y0:1, x1:x+w+4, y1:21})
    }
    legendItem('avg', 'averaged', cw-150, '#f97316', showAverageTrace)
    legendItem('ref', 'true', cw-78, 'rgba(56,189,248,0.75)', showReferenceTrace)

    // Axis labels: left edge, right edge, centre
    cmpctx.fillStyle='#9ca3af'; cmpctx.font='10px sans-serif'
    cmpctx.textAlign='left';   cmpctx.fillText((-xlim).toFixed(1)+' s', padL, ch-7)
    cmpctx.textAlign='right';  cmpctx.fillText((+xlim).toFixed(1)+' s', cw-padR, ch-7)
    cmpctx.textAlign='center'; cmpctx.fillText('zero lag', padL+(cw-padL-padR)/2, ch-7)

    // t0 tick labels
    if(t0 > 0){
      const xp=ptx(t0), xn=ptx(-t0)
      cmpctx.fillStyle='#6b7280'; cmpctx.textAlign='center'
      if(xp < cw-padR-2) cmpctx.fillText('+'+t0.toFixed(1), xp, ch-7)
      if(xn > padL+2)    cmpctx.fillText('-'+t0.toFixed(1), xn, ch-7)
    }
  }

  function comparisonLegendHit(px, py) {
    return comparisonLegendHits.find(hit => px >= hit.x0 && px <= hit.x1 && py >= hit.y0 && py <= hit.y1)
  }

  function redraw() {
    ctx.clearRect(0,0,PIX,PIX)
    ctx.strokeStyle='#374151'; ctx.lineWidth=1
    ctx.beginPath(); ctx.arc(MID,MID,RAD*SCALE,0,2*Math.PI); ctx.stroke()
    drawAxes()
    drawCanvasTitle()
    const hasWeights = currentWeights.length === pts.length
    pts.forEach(([px,py], i) => {
      const isDiscarded = hasWeights && subsampThreshold > 0 && currentWeights[i] <= subsampThreshold
      ctx.fillStyle = hasWeights ? (isDiscarded ? 'rgba(209,213,219,0.25)' : colorSchemeScale(currentWeights[i], 0.85)) : 'rgba(255,255,255,0.78)'
      ctx.beginPath(); ctx.arc(px,py,3,0,2*Math.PI); ctx.fill()
      if(selectedSourceIndices.has(i)){
        ctx.strokeStyle='#f8fafc'
        ctx.lineWidth=1.5
        ctx.beginPath(); ctx.arc(px,py,7,0,2*Math.PI); ctx.stroke()
      }
    })
    if(selectionRect) drawSelectionRect(selectionRect)
    drawReceiver(r1(),'#ef4444')
    drawReceiver(r2(),'#ef4444')
    drawInfo()
    lbl.textContent='Sources: '+pts.length
  }

  function emit() {
    enforcePeriodConstraint()
    try { sessionStorage.setItem(STATE_KEY, JSON.stringify({pts,d,period,bw_perc:bwPerc,phase_branch:phaseBranch,subsamp_threshold:subsampThreshold,source_distance_scale:sourceScale,spray_width:sprayWidth,apply_phase_correction:applyPhaseCorrection,drag_sources:dragSourcesMode})) } catch(e) {}
    par.value={sources:pts, d:d, f0:1/period, bw_perc:bwPerc, phase_branch:phaseBranch, subsample_threshold:subsampThreshold, source_distance_scale:sourceScale, apply_phase_correction:applyPhaseCorrection}
    par.dispatchEvent(new CustomEvent('input'))
  }

  function insideSourceCircle(px, py) {
    const dx=(px-MID)/SCALE, dy=(py-MID)/SCALE
    return dx*dx+dy*dy<=RAD*RAD
  }

  function nearestSourceIndex(px, py, radiusPx=8) {
    let best=-1, bestd2=radiusPx*radiusPx
    for(let i=pts.length-1;i>=0;i--){
      const dx=pts[i][0]-px, dy=pts[i][1]-py
      const d2=dx*dx+dy*dy
      if(d2 <= bestd2){
        best=i
        bestd2=d2
      }
    }
    return best
  }

  function selectedSourceHit(px, py, radiusPx=8) {
    let best=-1, bestd2=radiusPx*radiusPx
    selectedSourceIndices.forEach(i => {
      if(i < 0 || i >= pts.length) return
      const dx=pts[i][0]-px, dy=pts[i][1]-py
      const d2=dx*dx+dy*dy
      if(d2 <= bestd2){
        best=i
        bestd2=d2
      }
    })
    return best
  }

  function sourceInSelectionRect(pt, rect) {
    if(!rect) return false
    const x0=Math.min(rect.x0, rect.x1), x1=Math.max(rect.x0, rect.x1)
    const y0=Math.min(rect.y0, rect.y1), y1=Math.max(rect.y0, rect.y1)
    return pt[0] >= x0 && pt[0] <= x1 && pt[1] >= y0 && pt[1] <= y1
  }

  function drawSelectionRect(rect) {
    if(!rect) return
    const x=Math.min(rect.x0, rect.x1), y=Math.min(rect.y0, rect.y1)
    const w=Math.abs(rect.x1-rect.x0), h=Math.abs(rect.y1-rect.y0)
    ctx.fillStyle='rgba(147,197,253,0.16)'
    ctx.strokeStyle='rgba(147,197,253,0.9)'
    ctx.lineWidth=1
    ctx.fillRect(x,y,w,h)
    ctx.strokeRect(x,y,w,h)
  }

  function clampGroupDelta(indices, dx, dy) {
    let scale=1
    const maxr=RAD*SCALE
    indices.forEach(i => {
      if(i < 0 || i >= dragInitialPositions.length) return
      const p=dragInitialPositions[i]
      const tx=p[0]+dx, ty=p[1]+dy
      const vx=tx-MID, vy=ty-MID
      const a=dx*dx + dy*dy
      const b=2*((p[0]-MID)*dx + (p[1]-MID)*dy)
      const c=(p[0]-MID)*(p[0]-MID) + (p[1]-MID)*(p[1]-MID) - maxr*maxr
      if(vx*vx + vy*vy <= maxr*maxr || a <= Number.EPSILON) return
      const disc=b*b-4*a*c
      if(disc < 0) {
        scale=0
      } else {
        const t=(-b + Math.sqrt(disc))/(2*a)
        scale=Math.min(scale, Math.max(0, Math.min(1, t)))
      }
    })
    return [dx*scale, dy*scale]
  }

  function startGroupDrag(px, py) {
    draggingSource=true
    boxSelecting=false
    selectionRect=null
    dragStart=[px, py]
    dragInitialPositions=pts.map(p => [p[0], p[1]])
    cvs.style.cursor='grabbing'
  }

  function moveSelectedSources(px, py) {
    if(!dragStart) return
    let dx=px-dragStart[0], dy=py-dragStart[1]
    const clampedDelta = clampGroupDelta(selectedSourceIndices, dx, dy)
    dx = clampedDelta[0]
    dy = clampedDelta[1]
    selectedSourceIndices.forEach(i => {
      if(i < 0 || i >= pts.length || i >= dragInitialPositions.length) return
      pts[i]=[Math.round(dragInitialPositions[i][0]+dx), Math.round(dragInitialPositions[i][1]+dy)]
    })
    redraw()
  }

  function deleteSelectedSources() {
    if(!dragSourcesMode || selectedSourceIndices.size === 0) return false
    const selected = selectedSourceIndices
    pts = pts.filter((_, i) => !selected.has(i))
    selectedSourceIndices.clear()
    draggingSource=false
    boxSelecting=false
    selectionStart=null
    selectionRect=null
    dragStart=null
    dragInitialPositions=[]
    currentWeights=[]
    window._seismic_weights=[]
    redraw()
    emit()
    return true
  }

  function addPts(px,py) {
    let added=false
    if(sprayWidth === 0){
      if(insideSourceCircle(px, py)){
        pts.push([Math.round(px),Math.round(py)])
        added=true
      }
    } else {
      for(let i=0;i<PER_PT;i++){
        const dx=(Math.random()-0.5)*2*sprayWidth*SCALE
        const dy=(Math.random()-0.5)*2*sprayWidth*SCALE
        const nx=px+dx, ny=py+dy
        if(insideSourceCircle(nx, ny)){
          pts.push([Math.round(nx),Math.round(ny)])
          added=true
        }
      }
    }
    redraw()
    return added
  }

  function hitReceiver(px, py) {
    const hits = [receiverPixel(r1()), receiverPixel(r2())]
    return hits.some(([rx,ry]) => (px-rx)*(px-rx) + (py-ry)*(py-ry) <= 18*18)
  }

  function setDistanceFromPointer(px) {
    const x = (px - MID) / SCALE
    d = Math.max(10, Math.min(100, 2 * Math.abs(x)))
    d = Math.round(d / 5) * 5
    syncControls()
    redraw()
  }

  par.querySelector('#dist').addEventListener('input', e => {
    d = +e.target.value
    syncControls(); redraw(); emit()
  })
  par.querySelector('#period').addEventListener('input', e => {
    period = sliderToPeriod(+e.target.value)
    syncControls(); redraw(); emit()
  })
  par.querySelector('#bw').addEventListener('input', e => {
    bwPerc = +e.target.value
    syncControls(); emit()
  })
  par.querySelector('#phasebranch').addEventListener('input', e => {
    phaseBranch = Math.round(+e.target.value)
    syncControls(); emit()
  })
  par.querySelector('#subsamp').addEventListener('input', e => {
    subsampThreshold = +e.target.value
    syncControls(); redraw(); emit()
  })
  par.querySelector('#srcscale').addEventListener('input', e => {
    sourceScale = +e.target.value
    syncControls(); redraw(); emit()
  })
  par.querySelector('#spraywidth').addEventListener('input', e => {
    sprayWidth = +e.target.value
    syncControls()
    emit()
  })
  par.querySelector('#phasecorr').addEventListener('change', e => {
    applyPhaseCorrection = e.target.checked
    syncControls()
    emit()
  })
  par.querySelector('#dragsources').addEventListener('change', e => {
    dragSourcesMode = e.target.checked
    draggingSource=false
    boxSelecting=false
    selectionStart=null
    selectionRect=null
    dragStart=null
    dragInitialPositions=[]
    if(!dragSourcesMode) selectedSourceIndices.clear()
    syncControls(); redraw(); emit()
  })

  cvs.addEventListener('mousedown',e=>{
    if(hitReceiver(e.offsetX,e.offsetY)){
      draggingReceiver=true
      cvs.style.cursor='grabbing'
      setDistanceFromPointer(e.offsetX)
    } else if(dragSourcesMode) {
      const selectedHit = selectedSourceHit(e.offsetX, e.offsetY)
      const hit = selectedHit >= 0 ? selectedHit : nearestSourceIndex(e.offsetX, e.offsetY)
      if(hit >= 0){
        if(!selectedSourceIndices.has(hit)){
          selectedSourceIndices.clear()
          selectedSourceIndices.add(hit)
        }
        startGroupDrag(e.offsetX, e.offsetY)
      } else {
        draggingSource=false
        selectedSourceIndices.clear()
        boxSelecting=true
        selectionStart=[e.offsetX, e.offsetY]
        selectionRect={x0:e.offsetX, y0:e.offsetY, x1:e.offsetX, y1:e.offsetY}
        cvs.style.cursor='crosshair'
        redraw()
      }
    } else {
      if(sprayWidth === 0){
        if(addPts(e.offsetX,e.offsetY)) emit()
      } else {
        drawing=true
        addPts(e.offsetX,e.offsetY)
      }
    }
  })
  cvs.addEventListener('mousemove',e=>{
    if(draggingReceiver){
      setDistanceFromPointer(e.offsetX)
    } else if(draggingSource) {
      moveSelectedSources(e.offsetX, e.offsetY)
    } else if(boxSelecting) {
      selectionRect={x0:selectionStart[0], y0:selectionStart[1], x1:e.offsetX, y1:e.offsetY}
      redraw()
    } else if(drawing && sprayWidth > 0) {
      addPts(e.offsetX,e.offsetY)
    } else {
      if(hitReceiver(e.offsetX,e.offsetY)){
        cvs.style.cursor = 'ew-resize'
      } else if(dragSourcesMode && nearestSourceIndex(e.offsetX, e.offsetY) >= 0){
        cvs.style.cursor = 'grab'
      } else {
        cvs.style.cursor = dragSourcesMode ? 'default' : 'crosshair'
      }
    }
  })
  window.addEventListener('mouseup',()=>{
    if(draggingReceiver){
      draggingReceiver=false
      cvs.style.cursor=dragSourcesMode ? 'default' : 'crosshair'
      emit()
    } else if(draggingSource){
      draggingSource=false
      dragStart=null
      dragInitialPositions=[]
      cvs.style.cursor=dragSourcesMode ? 'default' : 'crosshair'
      emit()
    } else if(boxSelecting){
      const w=selectionRect ? Math.abs(selectionRect.x1-selectionRect.x0) : 0
      const h=selectionRect ? Math.abs(selectionRect.y1-selectionRect.y0) : 0
      selectedSourceIndices.clear()
      if(w >= 4 && h >= 4){
        pts.forEach((pt, i) => {
          if(sourceInSelectionRect(pt, selectionRect)) selectedSourceIndices.add(i)
        })
      }
      boxSelecting=false
      selectionStart=null
      selectionRect=null
      cvs.style.cursor=dragSourcesMode ? 'default' : 'crosshair'
      redraw()
    } else if(drawing){
      drawing=false
      emit()
    } else {
      drawing=false
    }
  })

  window.addEventListener('keydown', e => {
    if(e.key !== 'Delete' && e.key !== 'Backspace') return
    const target = e.target
    const tag = target && target.tagName ? target.tagName.toLowerCase() : ''
    if(tag === 'input' || tag === 'textarea' || tag === 'select' || (target && target.isContentEditable)) return
    if(deleteSelectedSources()) e.preventDefault()
  })

  cmpcvs.addEventListener('click', e => {
    const hit = comparisonLegendHit(e.offsetX, e.offsetY)
    if(!hit) return
    if(hit.key === 'avg') {
      showAverageTrace = !showAverageTrace
      if(!showAverageTrace && !showReferenceTrace) showAverageTrace = true
    } else if(hit.key === 'ref') {
      showReferenceTrace = !showReferenceTrace
      if(!showAverageTrace && !showReferenceTrace) showReferenceTrace = true
    }
    drawComparisonPanel()
  })
  cmpcvs.addEventListener('mousemove', e => {
    cmpcvs.style.cursor = comparisonLegendHit(e.offsetX, e.offsetY) ? 'pointer' : 'default'
  })

  par.querySelector('#clrbtn').addEventListener('click',()=>{
    pts=[]; selectedSourceIndices.clear(); currentWeights=[]; window._seismic_weights=[]; redraw(); emit()
  })
  par.querySelectorAll('.preset').forEach(btn=>btn.addEventListener('click',()=>{
    if(btn.dataset.preset === 'uniformdisk') {
      pts=makeUniformDisk(2*N)
    } else if(btn.dataset.preset === 'sourcecloud') {
      pts=makeSourceCloud(N, 135, 50, 14)
    } else if(btn.dataset.preset === 'twoclouds') {
      pts=makeTwoClouds(N)
    } else {
      pts=PRESETS[btn.dataset.preset].slice()
    }
    selectedSourceIndices.clear(); currentWeights=[]; window._seismic_weights=[]; redraw(); emit()
  }))

  // Listen for contribution weights fired by the reactive signal cell
  window.addEventListener('seismic-weights', e => {
    currentWeights = e.detail ? e.detail.split(',').map(Number) : []
    window._seismic_weights = currentWeights
    redraw()
  })
  window.addEventListener('seismic-comparison', e => {
    currentComparison = e.detail || null
    window._seismic_comparison = currentComparison
    drawComparisonPanel()
    renderMetricsPanel()
  })

  syncControls()
  emit()
  redraw()
  drawComparisonPanel()
  renderMetricsPanel()
</script>
        """)
    end

    const _canvas_ready = true
end

# ╔═╡ d59342b1-cb2a-4472-9717-b82fa12335f9
begin
    _canvas_ready
    WideCell(@bind _canvas_widget CanvasSourceInput(vel))
end

# ╔═╡ ca1b63bf-c868-4200-80f0-f8d6e40d9a6c
begin
    f0              = _canvas_widget["f0"]
    bw_perc         = _canvas_widget["bw_perc"]
    source_distance = stationary_phase_source_distance
end

# ╔═╡ 1ea5bea2-0c9d-4095-a395-24ad4f67ba05
function wavetype()
    return gaussian_wavelet(f0, bw_perc, dt)
end

# ╔═╡ 92057a3a-e880-46ce-86cc-7cc702085384
function generate_seismograms(rec1, rec2, srcloc)
    generate_seismograms(rec1, rec2, srcloc, wavetype(), nt, dt, vel)
end

# ╔═╡ c71bc2d6-cc33-4cae-8195-5fe26ec39d7f
begin
    srcloc_pixels = _canvas_widget["sources"]
    d_receiver = Float64(_canvas_widget["d"])
    source_distance_scale = Float64(get(_canvas_widget, "source_distance_scale", 1.0))
    rec1 = Float64[-d_receiver / 2, 0.0]
    rec_origin = Float64[0.0, 0.0]
    rec2 = Float64[ d_receiver / 2, 0.0]
    rad = 200
    srcloc = map(srcloc_pixels) do I
        x = (I[1] - 320) / 1.5
        y = -(I[2] - 320) / 1.5
        [x, y]
    end
    srcloc_for_noise = map(p -> source_distance_scale .* p, srcloc)
    xs = isempty(srcloc) ? Float64[] : getindex.(srcloc, 1)
    ys = isempty(srcloc) ? Float64[] : getindex.(srcloc, 2)
end

# ╔═╡ 280b6269-a022-4846-9071-61481649b008
let
    Th = range(0, 2pi, length=100)
    plot(rad2deg.(Th), [get_traveltime(rec1, source_distance .*[cos(th), sin(th)]) - get_traveltime(rec2, source_distance.*[cos(th), sin(th)]) for th in Th], Layout(title="Stationary phase", xaxis=attr(title="source backazimuth (deg)"), yaxis=attr(title="traveltime difference")))
end

# ╔═╡ 0279d32c-d884-4990-a1a6-68b43fae2f92
begin
    d_interstation = sqrt((rec2[1] - rec1[1])^2 + (rec2[2] - rec1[2])^2)
    d_ab = sqrt((rec_origin[1] - rec1[1])^2 + (rec_origin[2] - rec1[2])^2)
    d_bc = sqrt((rec2[1] - rec_origin[1])^2 + (rec2[2] - rec_origin[2])^2)
    lambda = vel / f0
    d_over_lambda = d_interstation / lambda
end

# ╔═╡ 3c5fd716-5b57-4e6a-aa8a-022a1cb476c4
md"""
| | |
|:---|:---|
| **λ** (wavelength) | $(round(lambda, digits=2)) km |
| **d** (inter-station distance) | $(round(d_interstation, digits=2)) km |
| **d/λ** | $(round(d_over_lambda, digits=2)) |
"""

# ╔═╡ a0f471bd-8139-4f42-bf0e-c019e6d42832
traveltime_between_receivers = get_traveltime(rec1, rec2)

# ╔═╡ d5305327-7086-4e2a-a4b7-7e3456864b8a
true_virtual_response = true_virtual_response_for_traveltime(traveltime_between_receivers, wavetype(), tgrid_xcorr, dt)

# ╔═╡ 97289274-85f2-46f7-b9bf-916f0a420301
seis1, seis2 = generate_seismograms(rec1, rec2, srcloc_for_noise, wavetype(), nt, dt, vel);

# ╔═╡ 617fd9b1-c6be-476e-9622-b6482602c7c7
let

    seis_heatmap(tgrid, randobs_safe(seis1, 200), "Noise at receiver 1", "time (s)", "Noise source index")
end

# ╔═╡ 64a67ad0-2a5f-49c4-8a90-8fe8fb55342f
let
    seis_heatmap(tgrid, randobs_safe(seis2, 200), "Noise at receiver 2", "time (s)", "Noise source index")
end

# ╔═╡ 2bad7a74-1b12-499f-b435-7f99c7d64ae9
cross = cross_correlations(seis1, seis2);

# ╔═╡ e69db63a-52b4-45c0-8b17-5239468cf5ad
let
    seis_heatmap(tgrid_xcorr, randobs_safe(cross, 1000), "Cross-correlated noise", "time lag (s)", "Source Index")
end

# ╔═╡ 3651d79a-e1f8-4960-a058-00ee5b609598
contributions = source_contributions(cross, true_virtual_response, tgrid_xcorr)

# ╔═╡ f3a2b1c0-5d4e-4f60-9012-345678901bcd
let _ = contributions
    csv = isempty(contributions) ? "" : join(round.(contributions, digits=4), ",")
    weights_js = isempty(contributions) ? "[]" : "[" * csv * "]"
    HTML("""<script>
      window._seismic_weights = $(weights_js);
      window.dispatchEvent(new CustomEvent('seismic-weights', {detail: '$(csv)'}));
    </script>""")
end

# ╔═╡ 77d2d0d5-05e8-44e4-910b-b8bc8cf3137f
function plot_point(x1, x2; weights=nothing)
    layout = Layout(xaxis_range=(-rad, rad), yaxis_range=(-rad, rad), width=600, height=400, title="Experiment setup", xaxis=attr(title="distance"), yaxis=attr(title="distance"))
    fig = Plot(layout)

    src_marker = if isnothing(weights) || isempty(weights)
        attr(size=5, color="red")
    else
        attr(size=5, color=weights, colorscale="Viridis", cmin=0, cmax=1,
             colorbar=attr(title="Contribution", thickness=12))
    end
    add_trace!(fig, scatter(x=x1, y=x2, mode="markers", showlegend=false, marker=src_marker))
    add_trace!(fig, scatter(x=[rec1[1]], y=[rec1[2]], mode="markers", name="Receiver 1", marker=attr(size=15, color="red", symbol="triangle-down")))
    add_trace!(fig, scatter(x=[rec2[1]], y=[rec2[2]], mode="markers", name="Receiver 2", marker=attr(size=15, color="red", symbol="triangle-down")))
    return PlutoPlotly.plot(fig)
end

# ╔═╡ b933216e-203a-47a1-b4f1-901fa5b38f70
plot_point(xs, ys; weights=length(contributions) == length(xs) ? contributions : nothing)

# ╔═╡ d4b78a1f-a6f7-44ef-ac09-fd0eb329bb2c
normalize_trace(x) = begin
    m = maximum(abs, x)
    m > 0 ? x ./ m : x
end

# ╔═╡ a0c7145b-cb4c-4648-a989-b6a5f88fd32a
begin
    _subsamp_threshold = Float64(get(_canvas_widget, "subsample_threshold", 0.0))
    gren = average_selected_cross_correlations(cross, contributions, _subsamp_threshold)
    plot_line(tgrid_xcorr, [normalize_trace(gren[:, 1]), normalize_trace(true_virtual_response)], names=["Averaged cross-correlation", "True virtual source response"], title="C---")
end

# ╔═╡ 8f5d89d2-ae71-4b43-bf24-9c7ce404071a
function analytic_signal_fft(x)
    n = length(x)
    n == 0 && return ComplexF64[]

    h = zeros(Float64, n)
    h[1] = 1.0
    if iseven(n)
        h[2:(n ÷ 2)] .= 2.0
        h[n ÷ 2 + 1] = 1.0
    else
        h[2:((n + 1) ÷ 2)] .= 2.0
    end

    ifft(fft(Float64.(x)) .* h)
end

# ╔═╡ 7c3bf8c7-7c8b-4a35-b2b4-f0557d9a64d6
function arrival_window_half_width(tgrid, distance, vref, f0)
    t0 = distance / vref
    dt_local = length(tgrid) > 1 ? abs(tgrid[2] - tgrid[1]) : 0.0
    min(0.75 * t0, max(1.5 / f0, 10 * dt_local))
end

# ╔═╡ 125d1e24-b6fe-4771-941d-dba92dc0e10e
function global_trace_peak(trace)
    isempty(trace) && return NaN
    peak = maximum(abs, trace)
    isfinite(peak) && peak > eps(Float64) ? peak : NaN
end

# ╔═╡ 40f06dca-4a81-471e-913a-28cde7b6f5a8
function branch_trace_peak(trace, tgrid, center, half_width)
    peak = 0.0
    has_sample = false
    @inbounds for i in eachindex(tgrid)
        if abs(tgrid[i] - center) <= half_width
            peak = max(peak, abs(trace[i]))
            has_sample = true
        end
    end
    has_sample && isfinite(peak) ? peak : NaN
end

# ╔═╡ 65820b94-ff46-4060-aacb-d5c997c86d41
function branch_is_valid(trace, tgrid, center, half_width, global_peak; threshold=0.1)
    isfinite(global_peak) || return false
    peak = branch_trace_peak(trace, tgrid, center, half_width)
    isfinite(peak) && peak >= threshold * global_peak
end

# ╔═╡ c24b2a79-6138-47e6-bd8f-d8023880bb6c
function branch_phase_differences(avg, ref, tgrid, distance, vref, f0)
    z_avg = analytic_signal_fft(avg)
    z_ref = analytic_signal_fft(ref)
    t0 = distance / vref
    half_width = arrival_window_half_width(tgrid, distance, vref, f0)
    avg_peak = global_trace_peak(avg)

    phase_for(idxs) = begin
        s = zero(ComplexF64)
        avg_norm2 = 0.0
        ref_norm2 = 0.0
        @inbounds for i in idxs
            avg_norm2 += abs2(avg[i])
            ref_norm2 += abs2(ref[i])
            s += conj(z_avg[i]) * z_ref[i]
        end
        avg_norm2 > eps(Float64) && ref_norm2 > eps(Float64) && abs(s) > eps(Float64) ? angle(s) : NaN
    end

    causal = findall(t -> abs(t - t0) <= half_width, tgrid)
    acausal = findall(t -> abs(t + t0) <= half_width, tgrid)
    phi_causal = branch_is_valid(avg, tgrid, t0, half_width, avg_peak) ? phase_for(causal) : NaN
    phi_acausal = branch_is_valid(avg, tgrid, -t0, half_width, avg_peak) ? phase_for(acausal) : NaN

    (; phase_causal=phi_causal,
     phase_acausal=phi_acausal,
     phase_causal_pi=isfinite(phi_causal) ? phi_causal / π : NaN,
     phase_acausal_pi=isfinite(phi_acausal) ? phi_acausal / π : NaN)
end

# ╔═╡ f34c5f90-6ac4-4ca6-b6a9-294a42421cd6
function phase_correct_trace(avg, ref, tgrid, distance, vref, f0)
    phases = branch_phase_differences(avg, ref, tgrid, distance, vref, f0)
    z_avg = analytic_signal_fft(avg)
    corrected = copy(z_avg)

    if isfinite(phases.phase_causal)
        rot = cis(phases.phase_causal)
        @inbounds for i in searchsortedfirst(tgrid, 0.0):length(tgrid)
            corrected[i] = z_avg[i] * rot
        end
    end

    if isfinite(phases.phase_acausal)
        rot = cis(phases.phase_acausal)
        @inbounds for i in 1:searchsortedlast(tgrid, 0.0)
            corrected[i] = z_avg[i] * rot
        end
    end

    (; trace=real.(corrected),
     phase_causal_pi=phases.phase_causal_pi,
     phase_acausal_pi=phases.phase_acausal_pi)
end

# ╔═╡ 57a034c4-9ff5-4edc-9c54-11759e02eed9
function branch_velocity_measurements(trace, tgrid, distance, vref, f0, phase_branch)
        invalid = (; U=NaN, U_error_pct=NaN, c=NaN, c_error_pct=NaN)
        isempty(trace) && return (; U_causal=NaN, U_causal_error_pct=NaN, U_acausal=NaN, U_acausal_error_pct=NaN,
            c_causal=NaN, c_causal_error_pct=NaN, c_acausal=NaN, c_acausal_error_pct=NaN)

        z = analytic_signal_fft(trace)
        env = abs.(z)
        t0 = distance / vref
        half_width = arrival_window_half_width(tgrid, distance, vref, f0)
        trace_peak = global_trace_peak(trace)
        omega = 2π * f0

        measure_branch(center) = begin
            peak_amp = -Inf
            ig = 0
            branch_is_valid(trace, tgrid, center, half_width, trace_peak) || return invalid
            @inbounds for i in eachindex(tgrid)
                if abs(tgrid[i] - center) <= half_width
                    if env[i] > peak_amp
                        peak_amp = env[i]
                        ig = i
                    end
                end
            end
            ig == 0 && return invalid
            (!isfinite(peak_amp) || peak_amp <= eps(Float64)) && return invalid

            t_group = tgrid[ig]
            abs(t_group) <= eps(Float64) && return invalid

            U = distance / abs(t_group)
            phase = angle(z[ig])
            tau = t_group - (phase + 2π * phase_branch) / omega
            c = abs(tau) > eps(Float64) ? distance / abs(tau) : NaN
            isfinite(c) && c > 0 || return invalid

            (; U=U,
             U_error_pct=100 * abs(U - vref) / vref,
             c=c,
             c_error_pct=100 * abs(c - vref) / vref)
        end

        causal = measure_branch(t0)
        acausal = measure_branch(-t0)

        (; U_causal=causal.U,
         U_causal_error_pct=causal.U_error_pct,
         U_acausal=acausal.U,
         U_acausal_error_pct=acausal.U_error_pct,
         c_causal=causal.c,
         c_causal_error_pct=causal.c_error_pct,
         c_acausal=acausal.c,
         c_acausal_error_pct=acausal.c_error_pct)
    end

# ╔═╡ 9fb5c14c-cb15-440e-aaf2-116fe146c70c
function three_station_residual(v_ac, v_ab, v_bc, d_ac, d_ab, d_bc, vref)
        if !(isfinite(v_ac) && isfinite(v_ab) && isfinite(v_bc))
            return (; value=NaN, pct=NaN)
        end

        denom = d_ab / v_ab + d_bc / v_bc
        if !(isfinite(denom) && denom > eps(Float64))
            return (; value=NaN, pct=NaN)
        end

        v_abc = d_ac / denom
        value = abs(v_abc - v_ac)
        (; value=value, pct=100 * value / vref)
    end

# ╔═╡ ee8752a4-a9ad-46dd-bbdd-92fcf52174b2
function three_station_residuals(ac, ab, bc, d_ac, d_ab, d_bc, vref)
        U_causal = three_station_residual(ac.U_causal, ab.U_causal, bc.U_causal, d_ac, d_ab, d_bc, vref)
        U_acausal = three_station_residual(ac.U_acausal, ab.U_acausal, bc.U_acausal, d_ac, d_ab, d_bc, vref)
        c_causal = three_station_residual(ac.c_causal, ab.c_causal, bc.c_causal, d_ac, d_ab, d_bc, vref)
        c_acausal = three_station_residual(ac.c_acausal, ab.c_acausal, bc.c_acausal, d_ac, d_ab, d_bc, vref)

        (; R3_U_causal=U_causal.value,
         R3_U_causal_pct=U_causal.pct,
         R3_U_acausal=U_acausal.value,
         R3_U_acausal_pct=U_acausal.pct,
         R3_c_causal=c_causal.value,
         R3_c_causal_pct=c_causal.pct,
         R3_c_acausal=c_acausal.value,
         R3_c_acausal_pct=c_acausal.pct)
    end

# ╔═╡ fcc4e475-112a-458d-b484-e75893af35c2
function gaussian_period_filter_trace(trace, dt, period, bw_perc)
        if isempty(trace) || !(isfinite(dt) && dt > 0 && isfinite(period) && period > eps(Float64) && isfinite(bw_perc) && bw_perc > 0)
            return Float64.(trace)
        end

        n = length(trace)
        center_freq = 1 / period
        sigma_freq = center_freq * (bw_perc / 100) / sqrt(log(2))
        spectrum = fft(Float64.(trace))
        weights = similar(spectrum, Float64)
        @inbounds for k in 1:n
            freq = (k - 1 <= n ÷ 2) ? (k - 1) / (n * dt) : (k - 1 - n) / (n * dt)
            weights[k] = exp(-((abs(freq) - center_freq)^2) / (2 * sigma_freq^2))
        end

        real.(ifft(spectrum .* weights))
    end

# ╔═╡ 09792337-3d9e-45c8-85fd-675b87596436
function phase_velocity_at_arrival(trace, tgrid, dt, distance, vref, period, bw_perc, branch_sign, phase_branch)
        if !(isfinite(period) && period > eps(Float64) && isfinite(branch_sign) && branch_sign != 0)
            return NaN
        end
        isempty(trace) && return NaN

        filtered = gaussian_period_filter_trace(trace, dt, period, bw_perc)
        z = analytic_signal_fft(filtered)
        env = abs.(z)
        center = branch_sign * distance / vref
        half_width = arrival_window_half_width(tgrid, distance, vref, 1 / period)
        peak_amp = -Inf
        ig = 0
        @inbounds for i in eachindex(tgrid)
            if abs(tgrid[i] - center) <= half_width && env[i] > peak_amp
                peak_amp = env[i]
                ig = i
            end
        end
        (ig == 0 || !isfinite(peak_amp) || peak_amp <= eps(Float64)) && return NaN

        t_group = tgrid[ig]
        abs(t_group) <= eps(Float64) && return NaN
        omega = 2π / period
        phase = angle(z[ig])
        tau = t_group - (phase + 2π * phase_branch) / omega
        c = abs(tau) > eps(Float64) ? distance / abs(tau) : NaN

        isfinite(c) && c > 0 ? c : NaN
    end

# ╔═╡ 5e98ea10-78ff-43e3-89d9-08bab29d1943
function mode_residual(U, c, trace, tgrid, dt, distance, vref, period, bw_perc, branch_sign, phase_branch)
        if !(isfinite(U) && U > 0 && isfinite(c) && c > 0 && isfinite(period) && period > 0)
            return (; value=NaN, pct=NaN)
        end

        c0 = phase_velocity_at_arrival(trace, tgrid, dt, distance, vref, period, bw_perc, branch_sign, phase_branch)
        isfinite(c0) || return (; value=NaN, pct=NaN)

        dT = 0.1
        Tlo = max(eps(Float64), period - dT)
        Thi = period + dT
        c_lo = Tlo < period ? phase_velocity_at_arrival(trace, tgrid, dt, distance, vref, Tlo, bw_perc, branch_sign, phase_branch) : NaN
        c_hi = phase_velocity_at_arrival(trace, tgrid, dt, distance, vref, Thi, bw_perc, branch_sign, phase_branch)

        dc_dT = if isfinite(c_lo) && isfinite(c_hi) && Thi > Tlo
            (c_hi - c_lo) / (Thi - Tlo)
        elseif isfinite(c_hi) && Thi > period
            (c_hi - c0) / (Thi - period)
        elseif isfinite(c_lo) && period > Tlo
            (c0 - c_lo) / (period - Tlo)
        else
            NaN
        end

        denom = 1 + (period / c) * dc_dT
        if !(isfinite(denom) && abs(denom) > eps(Float64))
            return (; value=NaN, pct=NaN)
        end

        value = abs(U - c / denom)
        return (; value=value, pct=100 * value / vref)
    end

# ╔═╡ dc58d268-fbfa-4e7b-9fd9-c7073cae7b35
function mode_residuals(measurements, trace, tgrid, dt, distance, vref, f0, bw_perc, phase_branch)
        isempty(trace) && return (; R_mode_causal=NaN, R_mode_causal_pct=NaN, R_mode_acausal=NaN, R_mode_acausal_pct=NaN)
        period = 1 / f0
        causal = mode_residual(measurements.U_causal, measurements.c_causal, trace, tgrid, dt, distance, vref, period, bw_perc, 1.0, phase_branch)
        acausal = mode_residual(measurements.U_acausal, measurements.c_acausal, trace, tgrid, dt, distance, vref, period, bw_perc, -1.0, phase_branch)

        (; R_mode_causal=causal.value,
         R_mode_causal_pct=causal.pct,
         R_mode_acausal=acausal.value,
         R_mode_acausal_pct=acausal.pct)
    end

# ╔═╡ a57e060f-73b6-4989-a742-553631cdd873
begin
    apply_phase_correction = Bool(get(_canvas_widget, "apply_phase_correction", false))
    phase_branch = Int(round(Float64(get(_canvas_widget, "phase_branch", 0))))
    _phase_corrected = phase_correct_trace(gren[:, 1], true_virtual_response, tgrid_xcorr, d_interstation, vel, f0)
    comparison_trace = apply_phase_correction ? _phase_corrected.trace : gren[:, 1]
    ac_measurements = branch_velocity_measurements(comparison_trace, tgrid_xcorr, d_interstation, vel, f0, phase_branch)

    _three_station_wavelet = wavetype()
    _true_ab = true_virtual_response_for_traveltime(get_traveltime(rec1, rec_origin), _three_station_wavelet, tgrid_xcorr, dt)
    _true_bc = true_virtual_response_for_traveltime(get_traveltime(rec_origin, rec2), _three_station_wavelet, tgrid_xcorr, dt)

    _seis_ab_1, _seis_ab_2 = generate_seismograms(rec1, rec_origin, srcloc_for_noise, _three_station_wavelet, nt, dt, vel)
    _cross_ab = cross_correlations(_seis_ab_1, _seis_ab_2)
    _gren_ab = average_selected_cross_correlations(_cross_ab, contributions, _subsamp_threshold)
    _phase_ab = phase_correct_trace(_gren_ab[:, 1], _true_ab, tgrid_xcorr, d_ab, vel, f0)
    _trace_ab = apply_phase_correction ? _phase_ab.trace : _gren_ab[:, 1]
    _measurements_ab = branch_velocity_measurements(_trace_ab, tgrid_xcorr, d_ab, vel, f0, phase_branch)

    _seis_bc_1, _seis_bc_2 = generate_seismograms(rec_origin, rec2, srcloc_for_noise, _three_station_wavelet, nt, dt, vel)
    _cross_bc = cross_correlations(_seis_bc_1, _seis_bc_2)
    _gren_bc = average_selected_cross_correlations(_cross_bc, contributions, _subsamp_threshold)
    _phase_bc = phase_correct_trace(_gren_bc[:, 1], _true_bc, tgrid_xcorr, d_bc, vel, f0)
    _trace_bc = apply_phase_correction ? _phase_bc.trace : _gren_bc[:, 1]
    _measurements_bc = branch_velocity_measurements(_trace_bc, tgrid_xcorr, d_bc, vel, f0, phase_branch)

    _r3_measurements = three_station_residuals(ac_measurements, _measurements_ab, _measurements_bc, d_interstation, d_ab, d_bc, vel)
    _mode_measurements = mode_residuals(ac_measurements, comparison_trace, tgrid_xcorr, dt, d_interstation, vel, f0, bw_perc, phase_branch)
    velocity_measurements = merge(
        ac_measurements,
        (; phase_causal_pi=_phase_corrected.phase_causal_pi,
         phase_acausal_pi=_phase_corrected.phase_acausal_pi,
         phase_corrected=apply_phase_correction),
        _r3_measurements,
        _mode_measurements
    )
end

# ╔═╡ 07f887de-ea77-4f39-9df4-6754734d8c5d
let
    jsnum(x) = isfinite(x) ? string(Float64(x)) : "null"
    t0 = traveltime_between_receivers
    wav = wavetype()
    win = max(5 * t0, 5 * dt * length(wav))
    mask = abs.(tgrid_xcorr) .<= win
    t_win   = tgrid_xcorr[mask]
    avg_win = normalize_trace(comparison_trace[mask])
    ref_win = normalize_trace(true_virtual_response[mask])
    t_js   = join(Float64.(t_win),   ",")
    avg_js = join(Float64.(avg_win), ",")
    ref_js = join(Float64.(ref_win), ",")
    HTML("""<script>
      window._seismic_comparison = {
        t: [$(t_js)],
        avg: [$(avg_js)],
        ref: [$(ref_js)],
        t0: $(Float64(t0)),
        velocity: {
          U_causal: $(jsnum(velocity_measurements.U_causal)),
          U_causal_error_pct: $(jsnum(velocity_measurements.U_causal_error_pct)),
          U_acausal: $(jsnum(velocity_measurements.U_acausal)),
          U_acausal_error_pct: $(jsnum(velocity_measurements.U_acausal_error_pct)),
          c_causal: $(jsnum(velocity_measurements.c_causal)),
          c_causal_error_pct: $(jsnum(velocity_measurements.c_causal_error_pct)),
          c_acausal: $(jsnum(velocity_measurements.c_acausal)),
          c_acausal_error_pct: $(jsnum(velocity_measurements.c_acausal_error_pct)),
          phase_causal_pi: $(jsnum(velocity_measurements.phase_causal_pi)),
          phase_acausal_pi: $(jsnum(velocity_measurements.phase_acausal_pi)),
          R3_U_causal: $(jsnum(velocity_measurements.R3_U_causal)),
          R3_U_causal_pct: $(jsnum(velocity_measurements.R3_U_causal_pct)),
          R3_U_acausal: $(jsnum(velocity_measurements.R3_U_acausal)),
          R3_U_acausal_pct: $(jsnum(velocity_measurements.R3_U_acausal_pct)),
          R3_c_causal: $(jsnum(velocity_measurements.R3_c_causal)),
          R3_c_causal_pct: $(jsnum(velocity_measurements.R3_c_causal_pct)),
          R3_c_acausal: $(jsnum(velocity_measurements.R3_c_acausal)),
          R3_c_acausal_pct: $(jsnum(velocity_measurements.R3_c_acausal_pct)),
          R_mode_causal: $(jsnum(velocity_measurements.R_mode_causal)),
          R_mode_causal_pct: $(jsnum(velocity_measurements.R_mode_causal_pct)),
          R_mode_acausal: $(jsnum(velocity_measurements.R_mode_acausal)),
          R_mode_acausal_pct: $(jsnum(velocity_measurements.R_mode_acausal_pct)),
          phase_corrected: $(velocity_measurements.phase_corrected)
        }
      };
      window.dispatchEvent(new CustomEvent('seismic-comparison', {detail: window._seismic_comparison}));
    </script>""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FourierTools = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MLUtils = "f1d291b0-491e-4a28-83b9-f70985020b54"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
ColorSchemes = "~3.31.0"
FFTW = "~1.10.0"
FourierTools = "~0.5.0"
LaTeXStrings = "~1.4.0"
MLUtils = "~0.4.11"
PlutoPlotly = "~0.6.6"
PlutoUI = "~0.7.83"
StatsBase = "~0.34.12"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "d9e9d0423b6adf95c3674a561532e89d1e5d829f"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractNFFTs]]
deps = ["LinearAlgebra", "Printf", "Requires", "ScopedValues"]
git-tree-sha1 = "27aa55535187ed2e62aff2883399a443944d0f55"
uuid = "7f219486-4aa7-41d6-80a7-e08ef20ceed7"
version = "0.9.1"
weakdeps = ["ChainRulesCore"]

    [deps.AbstractNFFTs.extensions]
    AbstractNFFTsChainRulesCoreExt = "ChainRulesCore"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.45"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "daa72978cd7a624246e894a4f4f067706d4e17e2"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

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

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "b8651b2eb5796a386b0398a20b519a6a6150f75c"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "1.1.3"

    [deps.Atomix.extensions]
    AtomixCUDAExt = "CUDA"
    AtomixMetalExt = "Metal"
    AtomixOpenCLExt = "OpenCL"
    AtomixoneAPIExt = "oneAPI"

    [deps.Atomix.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    OpenCL = "08131aa3-fb12-5dee-8b74-c09406e224a2"
    oneAPI = "8f75cd03-7ff8-4ecb-9b8f-daf728133b1b"

[[deps.BangBang]]
deps = ["Accessors", "ConstructionBase", "InitialValues", "LinearAlgebra"]
git-tree-sha1 = "cceb62468025be98d42a5dc581b163c20896b040"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.4.9"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTablesExt = "Tables"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BasicInterpolators]]
deps = ["LinearAlgebra", "Memoize", "Random"]
git-tree-sha1 = "3f7be532673fc4a22825e7884e9e0e876236b12a"
uuid = "26cce99e-4866-4b6d-ab74-862489e035e0"
version = "0.7.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChunkSplitters]]
git-tree-sha1 = "1c52c8e2673edc030191177ff1aee42d25149acb"
uuid = "ae650224-84b6-46f8-82ea-d812ca08434e"
version = "3.2.0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "REPL", "UUIDs"]
git-tree-sha1 = "cfb7a2e89e245a9d5016b70323db412b3a7438d5"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "3.0.2"

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
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "6fb53a69613a0b2b68a0d12671717d307ab8b24e"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.5"

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

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

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

[[deps.FourierTools]]
deps = ["ChainRulesCore", "FFTW", "IndexFunArrays", "LinearAlgebra", "MutableShiftedArrays", "NDTools", "NFFT", "Reexport"]
git-tree-sha1 = "5ba2dec7e09d79f5965879645e80f83a879bc99c"
uuid = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
version = "0.5.0"

    [deps.FourierTools.extensions]
    CUDASupportExt_FT = ["CUDA", "Adapt"]

    [deps.FourierTools.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

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

[[deps.IndexFunArrays]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f78703c7a4ba06299cddd8694799c91de0157ac"
uuid = "613c443e-d742-454e-bfc6-1d7f8dd76566"
version = "0.2.7"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

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

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "MacroTools", "PrecompileTools", "Requires", "StaticArrays", "UUIDs"]
git-tree-sha1 = "a5b87110fa95d711355af44832497745aa93fb52"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.42"

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"
    LinearAlgebraExt = "LinearAlgebra"
    SparseArraysExt = "SparseArrays"

    [deps.KernelAbstractions.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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

[[deps.MLCore]]
deps = ["DataAPI", "SimpleTraits", "Tables"]
git-tree-sha1 = "c4ab44fe709638fda6f2c0cbfea2c114932d6c2f"
uuid = "c2834f40-e789-41da-a90e-33b280584a8c"
version = "1.1.0"

    [deps.MLCore.extensions]
    MLCorePythonCallExt = "PythonCall"

    [deps.MLCore.weakdeps]
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"

[[deps.MLUtils]]
deps = ["ChainRulesCore", "CodeTracking", "Compat", "DataAPI", "DelimitedFiles", "Distributed", "InteractiveUtils", "MLCore", "NNlib", "Random", "ShowCases", "SimpleTraits", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "5409200d8edbc8db329620b525cf10d25ead00c4"
uuid = "f1d291b0-491e-4a28-83b9-f70985020b54"
version = "0.4.11"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

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

[[deps.MutableShiftedArrays]]
git-tree-sha1 = "f2edd1c03452462e351ceb8d65477519a064917a"
uuid = "d3d30c82-a38e-471c-a45a-3d24d2f4d22d"
version = "0.3.2"

    [deps.MutableShiftedArrays.extensions]
    CUDASupportExt = ["CUDA", "Adapt"]

    [deps.MutableShiftedArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.NDTools]]
deps = ["LinearAlgebra", "MutableShiftedArrays", "Random", "Statistics"]
git-tree-sha1 = "160a53a9d98f8394d9afbc8387c4f3a54eee0d7b"
uuid = "98581153-e998-4eef-8d0d-5ec2c052313d"
version = "0.8.2"

[[deps.NFFT]]
deps = ["AbstractNFFTs", "BasicInterpolators", "Distributed", "FFTW", "LinearAlgebra", "OhMyThreads", "PrecompileTools", "Printf", "Random", "Reexport", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "9a1ecb7604507a6c0c77723b93dba335b5029cb1"
uuid = "efe261a4-0d2b-5849-be55-fc731d526b0d"
version = "0.14.3"

    [deps.NFFT.extensions]
    NFFTGPUArraysExt = ["Adapt", "GPUArrays"]

    [deps.NFFT.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Random", "ScopedValues", "Statistics"]
git-tree-sha1 = "2acc74264095c06beb62c6f69691a3277b5ac3d0"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.9.37"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"
    NNlibCUDACUDNNExt = ["CUDA", "cuDNN"]
    NNlibCUDAExt = "CUDA"
    NNlibEnzymeCoreExt = "EnzymeCore"
    NNlibFFTWExt = "FFTW"
    NNlibForwardDiffExt = "ForwardDiff"
    NNlibMetalExt = "Metal"
    NNlibMooncakeCUDAExt = ["Mooncake", "CUDA"]
    NNlibSpecialFunctionsExt = "SpecialFunctions"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OhMyThreads]]
deps = ["BangBang", "ChunkSplitters", "ScopedValues", "StableTasks", "TaskLocalValues"]
git-tree-sha1 = "9a07c25c438110500d871fd5309649ec6791ef57"
uuid = "67456a42-1dca-4109-a031-0a68de7e3ad5"
version = "0.8.6"

    [deps.OhMyThreads.extensions]
    MarkdownExt = "Markdown"
    ProgressMeterExt = "ProgressMeter"

    [deps.OhMyThreads.weakdeps]
    Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
    ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

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

[[deps.ShowCases]]
git-tree-sha1 = "7f534ad62ab2bd48591bdeac81994ea8c445e4a5"
uuid = "605ecd9f-84a6-4c9e-81e2-4798472b76a3"
version = "0.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "6547cbdd8ce32efba0d21c5a40fa96d1a3548f9f"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableTasks]]
git-tree-sha1 = "c4f6610f85cb965bee5bfafa64cbeeda55a4e0b2"
uuid = "91464d47-22a1-43fe-8b7f-2d57ee82463f"
version = "0.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

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

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

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

[[deps.UnsafeAtomics]]
git-tree-sha1 = "0f30765c32d66d58e41f4cb5624d4fc8a82ec13b"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.3.1"

    [deps.UnsafeAtomics.extensions]
    UnsafeAtomicsLLVM = ["LLVM"]

    [deps.UnsafeAtomics.weakdeps]
    LLVM = "929cbde3-209d-540e-8aea-75f648917ca0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

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
"""

# ╔═╡ Cell order:
# ╠═f986cb42-7cbe-11ee-22de-cbc3714ed81d
# ╟─bd04a8af-7666-47b7-87cc-a212ec8adcbd
# ╟─d59342b1-cb2a-4472-9717-b82fa12335f9
# ╟─3c5fd716-5b57-4e6a-aa8a-022a1cb476c4
# ╟─4a592f7e-fa42-4724-8281-1534333816d0
# ╠═b8671b46-d45b-4ec7-99f2-b1a2cb0aee12
# ╠═b933216e-203a-47a1-b4f1-901fa5b38f70
# ╠═a0c7145b-cb4c-4648-a989-b6a5f88fd32a
# ╟─a57e060f-73b6-4989-a742-553631cdd873
# ╟─07f887de-ea77-4f39-9df4-6754734d8c5d
# ╠═617fd9b1-c6be-476e-9622-b6482602c7c7
# ╟─64a67ad0-2a5f-49c4-8a90-8fe8fb55342f
# ╟─1bc0ce25-4de0-4628-a308-8106544fc0ec
# ╟─e69db63a-52b4-45c0-8b17-5239468cf5ad
# ╠═5427df02-9278-41c9-9b3b-b5d1ab5b4c24
# ╠═280b6269-a022-4846-9071-61481649b008
# ╟─e0ae191c-b674-46c7-bf05-bfecade313a5
# ╠═9a74489f-8e6d-45ce-ab57-464ee36d3316
# ╠═f0e27c51-233b-45ef-be1e-34b624b5a4b2
# ╠═dafaabbb-4263-41ce-b17f-0634adc01b8f
# ╠═1ea5bea2-0c9d-4095-a395-24ad4f67ba05
# ╟─577d411c-88b5-495d-a49a-d2379336f533
# ╠═c71bc2d6-cc33-4cae-8195-5fe26ec39d7f
# ╟─8d45d812-fb70-45c5-a123-4f93844fef6d
# ╠═ca1b63bf-c868-4200-80f0-f8d6e40d9a6c
# ╠═3b144853-0d82-46e2-a4ee-301516f3ce10
# ╠═0279d32c-d884-4990-a1a6-68b43fae2f92
# ╠═4e20077c-3687-47a1-8ab0-048b99fb0a46
# ╠═e04a82d2-9eda-419f-a6b3-0af28431e554
# ╠═a0f471bd-8139-4f42-bf0e-c019e6d42832
# ╠═4ce9cc5d-6808-4977-98bb-8e9cacf171cc
# ╠═d5305327-7086-4e2a-a4b7-7e3456864b8a
# ╠═1ce924ce-5a35-4a84-835a-3bdf16b876dc
# ╠═32be7b92-0b34-4ad0-8968-c77a13225383
# ╠═92057a3a-e880-46ce-86cc-7cc702085384
# ╠═97289274-85f2-46f7-b9bf-916f0a420301
# ╠═2bad7a74-1b12-499f-b435-7f99c7d64ae9
# ╠═2fd6cc16-519d-4d0f-ba42-874e58a32124
# ╠═e4f79b31-aa1f-4132-bc8e-b30d79380a3b
# ╠═3651d79a-e1f8-4960-a058-00ee5b609598
# ╟─f3a2b1c0-5d4e-4f60-9012-345678901bcd
# ╟─d483d881-11c4-495f-a103-2d2dd6371ca3
# ╠═290b50eb-0faf-4a5f-86bf-36628e61ff06
# ╠═ec238280-1edc-442f-837e-d3c807ea5fc6
# ╠═7ede27ab-f524-4f20-8bcb-5be70600c892
# ╠═e059662c-69af-460c-aada-408a365d0ff7
# ╠═366783f3-3a48-4e0d-a850-6c8cb7275377
# ╠═b11959f8-3cf2-4488-bbc1-6597220ea599
# ╠═6a8aeed2-173a-4523-afb3-f9d77c252816
# ╠═b2c0af32-b19d-47b0-86f3-66e8fd7d7269
# ╠═af973cc3-21b3-4dca-91ad-0c310e1214bc
# ╠═09bd50ae-ab80-4cde-bdb7-ccc7cd7041d4
# ╟─d9727185-8f69-49a5-a98f-7514562529ab
# ╠═77d2d0d5-05e8-44e4-910b-b8bc8cf3137f
# ╠═1293b18a-15d5-4d54-8396-59899c89a172
# ╠═8b213ed6-b4ec-44a9-8e3f-8ad03ee9f270
# ╠═ccf893d9-5b6e-4a79-b09d-3da957bc612c
# ╠═34f81006-dead-4687-aff1-ed2753f27469
# ╠═15a556df-3ed3-4068-86ca-d3c130cb81d9
# ╠═cd3a586a-b737-451e-921f-5463cc286f73
# ╠═0d2266c7-312c-4c56-8f7b-8a98fba81fc6
# ╠═7fe32a18-b484-4179-b0a8-d6c49059ed59
# ╠═d4b78a1f-a6f7-44ef-ac09-fd0eb329bb2c
# ╠═8f5d89d2-ae71-4b43-bf24-9c7ce404071a
# ╠═7c3bf8c7-7c8b-4a35-b2b4-f0557d9a64d6
# ╠═125d1e24-b6fe-4771-941d-dba92dc0e10e
# ╠═40f06dca-4a81-471e-913a-28cde7b6f5a8
# ╠═65820b94-ff46-4060-aacb-d5c997c86d41
# ╠═c24b2a79-6138-47e6-bd8f-d8023880bb6c
# ╠═f34c5f90-6ac4-4ca6-b6a9-294a42421cd6
# ╠═57a034c4-9ff5-4edc-9c54-11759e02eed9
# ╠═9fb5c14c-cb15-440e-aaf2-116fe146c70c
# ╠═ee8752a4-a9ad-46dd-bbdd-92fcf52174b2
# ╠═fcc4e475-112a-458d-b484-e75893af35c2
# ╠═09792337-3d9e-45c8-85fd-675b87596436
# ╠═5e98ea10-78ff-43e3-89d9-08bab29d1943
# ╠═dc58d268-fbfa-4e7b-9fd9-c7073cae7b35
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
