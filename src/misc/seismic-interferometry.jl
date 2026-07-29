### A Pluto.jl notebook ###
# v0.1.0

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
    using PlutoUI, Statistics, PlutoPlotly, FourierTools, FFTW, StatsBase, MLUtils
    import PlutoUI: combine
end

# ╔═╡ 7ede27ab-f524-4f20-8bcb-5be70600c892
using LaTeXStrings

# ╔═╡ f986cb42-7cbe-11ee-22de-cbc3714ed81d
TableOfContents(include_definitions=true)

# ╔═╡ bd04a8af-7666-47b7-87cc-a212ec8adcbd
md"""
## Simple Seismic Interferometry Demo
We generate wavelet at many sources in a homogeneous medium. There are 2 receivers that are recording the wavelets from all the sources. The wavelet from a source recorded at one of the receivers is cross-correlated with the wavelet recorded at the other receiver. The cross-correlogram of all the sources are stacked to get the approximate Green's function between the 2 receivers.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India

"""

# ╔═╡ 1083c0c6-beb2-42ba-8c99-e8d8d77ce83a
@bind _wavelet_params PlutoUI.combine() do Child
    md"""
    ## Source Wavelet
    | Parameter | Value | Notes |
    |:---|:---|:---|
    | Central frequency f₀ (Hz) | $(Child("f0", Slider(1.0:0.5:50.0, default=20.0, show_value=true))) | peak of Gaussian spectrum |
    | Bandwidth (%) | $(Child("bw_perc", Slider(5:5:100, default=40, show_value=true))) | FTAN-style half-power width |
    """
end

# ╔═╡ 1bc0ce25-4de0-4628-a308-8106544fc0ec
md"""
## Cross-Correlation
"""

# ╔═╡ 3c5fd716-5b57-4e6a-aa8a-022a1cb476c4
md"""
| | |
|:---|:---|
| **λ** (wavelength) | $(round(lambda, digits=4)) distance units |
| **d** (inter-station distance) | $(round(d_interstation, digits=3)) distance units |
| **d/λ** | $(round(d_over_lambda, digits=2)) |
"""

# ╔═╡ e0ae191c-b674-46c7-bf05-bfecade313a5
md"""
## Generate Wavelets
"""

# ╔═╡ 0279d32c-d884-4990-a1a6-68b43fae2f92
begin
    dt_wav = 0.002
    d_interstation = sqrt((rec2[1] - rec1[1])^2 + (rec2[2] - rec1[2])^2)
    lambda = vel / (f0 * dt_wav)
    d_over_lambda = d_interstation / lambda
end

# ╔═╡ e1f1fa09-87ef-4cdb-94b1-0d2e0e483323
function delaywav(w, t)
    ts = w
    ts = cat(zeros(Float32, 200), ts, zeros(Float32, 1000), dims=1)
    td = shift(ts, tuple(t - 300))
    return td[1:1000]

end

# ╔═╡ 577d411c-88b5-495d-a49a-d2379336f533
md"""
## Generate Random Source Locations
"""

# ╔═╡ c71bc2d6-cc33-4cae-8195-5fe26ec39d7f
begin
    rad = 20
    srcloc = map(srcloc_pixels) do I
        x = (I[1] - 200) / 10.0
        y = -(I[2] - 200) / 10.0
        [x, y]
    end
    xs = isempty(srcloc) ? Float64[] : getindex.(srcloc, 1)
    ys = isempty(srcloc) ? Float64[] : getindex.(srcloc, 2)
end

# ╔═╡ 8d45d812-fb70-45c5-a123-4f93844fef6d
md"""
## Receivers
"""

# ╔═╡ 3b144853-0d82-46e2-a4ee-301516f3ce10
vel = 0.035

# ╔═╡ 3cf54d17-5931-461c-9e35-bad1fc8d9ac3
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
function get_traveltime(rec, src)
    dis = sqrt((rec[1] - src[1])^2 + (rec[2] - src[2])^2)
    tm = dis / vel
    return tm
end

# ╔═╡ d483d881-11c4-495f-a103-2d2dd6371ca3
md"## Appendix"

# ╔═╡ e059662c-69af-460c-aada-408a365d0ff7
default_plotly_template(:plotly_dark)

# ╔═╡ 41870289-3874-43ca-9837-3cead43d0200
function gaussian_wavelet(f0, bw_perc)
    dt = 0.002
    N = 350
    t = dt .* (0:N-1) .- dt * (N ÷ 2)
    sigma_t = sqrt(log(2)) / (2 * π * f0 * (bw_perc / 100))
    w = @. exp(-t^2 / (2 * sigma_t^2)) * cos(2 * π * f0 * t)
    w ./= maximum(abs, w)
    return Float32.(w)
end

# ╔═╡ 1ea5bea2-0c9d-4095-a395-24ad4f67ba05
function wavetype()
    return gaussian_wavelet(f0, bw_perc)
end

# ╔═╡ 09c19be2-867a-4c16-ace5-6e4a9389472c
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
function generate_seismograms(rec1, rec2, srcloc)
    tm1 = map(1:length(srcloc)) do i
        get_traveltime(rec1, srcloc[i])
    end
    tm2 = map(1:length(srcloc)) do i
        get_traveltime(rec2, srcloc[i])
    end
    wavel = wavetype()
    wav1 = map(1:length(tm1)) do h
        delaywav(wavel, tm1[h])
    end

    wav2 = map(1:length(tm2)) do h
        delaywav(wavel, tm2[h])
    end

    wav1 = stack(wav1, dims=2)
    wav2 = stack(wav2, dims=2)
    return wav1, wav2
end

# ╔═╡ 6a8aeed2-173a-4523-afb3-f9d77c252816
nt = 1999

# ╔═╡ af973cc3-21b3-4dca-91ad-0c310e1214bc
tgrid_xcorr = range(-100, stop=100, length=nt)

# ╔═╡ 09bd50ae-ab80-4cde-bdb7-ccc7cd7041d4
tgrid = range(0, stop=100, length=nt)

# ╔═╡ d9727185-8f69-49a5-a98f-7514562529ab
md"""
### Plots
"""

# ╔═╡ 1293b18a-15d5-4d54-8396-59899c89a172
function seis_heatmap(tgrid, r, title, ytitle, xtitle)
    m = maximum(abs, r)
    plot(heatmap(y=tgrid, z=r, colorscale="seismic", zmin=-m, zmax=m), Layout(title=title, yaxis_autorange="reversed", height=350, width=600, yaxis=attr(title=ytitle), xaxis=attr(title=xtitle)))
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

# ╔═╡ a153664e-24bc-4fea-b0f4-f12360eaa688
function xcorr(x, y; padmode=:longest)
    n = length(x) + length(y) - 1
    X = fft(vcat(x, zeros(eltype(x), length(y) - 1)))
    Y = fft(vcat(y, zeros(eltype(y), length(x) - 1)))
    return real.(ifft(conj(X) .* Y))
end

# ╔═╡ 5427df02-9278-41c9-9b3b-b5d1ab5b4c24
@bind _analysis_params PlutoUI.combine() do Child
    md"""
    ## Analysis
    | Parameter | Value | Notes |
    |:---|:---|:---|
    | Lag window (s) | $(Child("lag_window", Slider(0.5:0.5:20.0, default=2.0, show_value=true))) | window for stationary source highlighting |
    | Stationary phase distance | $(Child("source_distance", Slider(range(1, 10, length=100), default=5.0, show_value=true))) | source distance for stationary phase plot |
    """
end

# ╔═╡ c4339754-42eb-4969-9045-010ca1d23ef8
@bind _receiver_params confirm(PlutoUI.combine() do Child
    md"""
    ## Receiver Positions
    | Receiver | x | y |
    |:---|:---|:---|
    | R1 | $(Child("r1x", Slider(-20:1:20, default=-5, show_value=true))) | $(Child("r1y", Slider(-20:1:20, default=0, show_value=true))) |
    | R2 | $(Child("r2x", Slider(-20:1:20, default=5,  show_value=true))) | $(Child("r2y", Slider(-20:1:20, default=0, show_value=true))) |
    """
end)

# ╔═╡ ca1b63bf-c868-4200-80f0-f8d6e40d9a6c
begin
    f0              = _wavelet_params.f0
    bw_perc         = _wavelet_params.bw_perc
    lag_window      = _analysis_params.lag_window
    source_distance = _analysis_params.source_distance
    rec1 = [_receiver_params.r1x, _receiver_params.r1y]
    rec2 = [_receiver_params.r2x, _receiver_params.r2y]
end

# ╔═╡ 280b6269-a022-4846-9071-61481649b008
let
    Th = range(0, 2pi, length=100)
    plot(rad2deg.(Th), [get_traveltime(rec1, source_distance .*[cos(th), sin(th)]) - get_traveltime(rec2, source_distance.*[cos(th), sin(th)]) for th in Th], Layout(title="Stationary phase", xaxis=attr(title="source backazimuth (deg)"), yaxis=attr(title="traveltime difference")))
end

# ╔═╡ a0f471bd-8139-4f42-bf0e-c019e6d42832
traveltime_between_receivers = get_traveltime(rec1, rec2)

# ╔═╡ 6098be3e-fe12-4f11-9749-b8c66c6106a8
true_virtual_response = delaywav(wavetype(), traveltime_between_receivers)

# ╔═╡ 97289274-85f2-46f7-b9bf-916f0a420301
seis1, seis2 = generate_seismograms(rec1, rec2, srcloc);

# ╔═╡ 617fd9b1-c6be-476e-9622-b6482602c7c7
let

    seis_heatmap(tgrid, randobs(seis1, 200), "Noise at receiver 1", "time (s)", "Noise source index")
end

# ╔═╡ 64a67ad0-2a5f-49c4-8a90-8fe8fb55342f
let
    seis_heatmap(tgrid, randobs(seis2, 200), "Noise at receiver 2", "time (s)", "Noise source index")
end

# ╔═╡ 2bad7a74-1b12-499f-b435-7f99c7d64ae9
begin
    cross = map(1:size(seis1, 2)) do i
        xcorr(seis1[:, i], seis2[:, i], padmode=:longest)
    end
    cross = stack(cross, dims=2)
end;

# ╔═╡ 46c9b13c-088c-47e7-8fd3-b56aee4a8923
begin
    gren = mean(cross, dims=2)
    plot_line(tgrid_xcorr, [gren[:, 1], cat(reverse(true_virtual_response), true_virtual_response, dims=1)], names=["Averaged cross-correlation", "True virtual source response"], title="Comparision")
end

# ╔═╡ e69db63a-52b4-45c0-8b17-5239468cf5ad
let
    seis_heatmap(tgrid_xcorr, randobs(cross, 1000), "Cross-correlated noise", "time lag (s)", "Source Index")
end

# ╔═╡ 77d2d0d5-05e8-44e4-910b-b8bc8cf3137f
function plot_point(x1, x2; weights=nothing)
    layout = Layout(xaxis_range=(-rad, rad), yaxis_range=(-rad, rad), width=600, height=400, title="Experiment setup", xaxis=attr(title="distance"), yaxis=attr(title="distance"))
    fig = Plot(layout)

    src_marker = if isnothing(weights) || isempty(weights)
        attr(size=5, color="red")
    else
        attr(size=5, color=weights, colorscale="RdYlGn", cmin=0, cmax=1,
             colorbar=attr(title="Contribution", thickness=12))
    end
    add_trace!(fig, scatter(x=x1, y=x2, mode="markers", showlegend=false, marker=src_marker))
    add_trace!(fig, scatter(x=[rec1[1]], y=[rec1[2]], mode="markers", name="Receiver 1", marker=attr(size=15, color="blue", symbol="triangle-down")))
    add_trace!(fig, scatter(x=[rec2[1]], y=[rec2[2]], mode="markers", name="Receiver 2", marker=attr(size=15, color="green", symbol="triangle-down")))
    return PlutoPlotly.plot(fig)
end

# ╔═╡ b933216e-203a-47a1-b4f1-901fa5b38f70
plot_point(xs, ys; weights=length(contributions) == length(xs) ? contributions : nothing)

# ╔═╡ 7fe32a18-b484-4179-b0a8-d6c49059ed59
begin
    struct CanvasSourceInput
        default_pts::Vector{Vector{Int}}
    end

    function CanvasSourceInput()
        pts = [[round(Int, 200 + 170*cos(a)), round(Int, 200 + 170*sin(a))]
               for a in range(0, 2π, length=501)[1:end-1]]
        CanvasSourceInput(pts)
    end

    Base.get(w::CanvasSourceInput) = w.default_pts

    function Base.show(io::IO, ::MIME"text/html", w::CanvasSourceInput)
        write(io, """
<div style="display:inline-block">
  <canvas id="srccvs" width="400" height="400"
    style="cursor:crosshair;background:#111827;border-radius:6px;display:block"></canvas>
  <div style="margin-top:4px;display:flex;gap:6px;align-items:center;flex-wrap:wrap">
    <button id="clrbtn">Clear</button>
    <span id="cnt" style="color:#9ca3af;font-size:0.85em">Sources: $(length(w.default_pts))</span>
  </div>
  <div style="margin-top:4px;display:flex;gap:4px;flex-wrap:wrap">
    <span style="color:#9ca3af;font-size:0.8em;align-self:center">Preset:</span>
    <button class="preset" data-preset="circular">Circular</button>
    <button class="preset" data-preset="sectional">Sectional</button>
    <button class="preset" data-preset="leftarc">Left Arc</button>
    <button class="preset" data-preset="rightarc">Right Arc</button>
    <button class="preset" data-preset="special">Special</button>
    <button class="preset" data-preset="bimodal">Bimodal</button>
  </div>
</div>
<script>
  const SCALE=10, SPREAD=1.5, PER_PT=15, RAD=20, RING=17, N=500
  const par = currentScript.previousElementSibling
  const cvs = par.querySelector('#srccvs')
  const ctx = cvs.getContext('2d')
  const lbl = par.querySelector('#cnt')

  function makeArc(a1, a2, n) {
    const out = []
    for(let i=0;i<n;i++){
      const a = a1 + (a2-a1)*i/Math.max(n-1,1)
      out.push([Math.round(200+RING*SCALE*Math.cos(a)), Math.round(200+RING*SCALE*Math.sin(a))])
    }
    return out
  }

  const PRESETS = {
    circular:  makeArc(0, 2*Math.PI, N),
    sectional: makeArc(-0.1*Math.PI, 0.1*Math.PI, N/2).concat(makeArc(0.9*Math.PI, 1.1*Math.PI, N/2)),
    leftarc:   makeArc(Math.PI/3, -Math.PI/3, N),
    rightarc:  makeArc(2*Math.PI/3, 4*Math.PI/3, N),
    special:   makeArc(Math.PI/3, 2*Math.PI/3, N/2).concat(makeArc(4*Math.PI/3, 5*Math.PI/3, N/2)),
    bimodal:   makeArc(0, Math.PI/6, N/2).concat(makeArc(Math.PI, 7*Math.PI/6, N/2)),
  }

  let pts = PRESETS.circular.slice()
  let drawing = false

  function redraw() {
    ctx.clearRect(0,0,400,400)
    ctx.strokeStyle='#374151'; ctx.lineWidth=1
    ctx.beginPath(); ctx.arc(200,200,RAD*SCALE,0,2*Math.PI); ctx.stroke()
    ctx.strokeStyle='#1f2937'; ctx.lineWidth=0.5
    ctx.beginPath(); ctx.moveTo(200,0); ctx.lineTo(200,400); ctx.stroke()
    ctx.beginPath(); ctx.moveTo(0,200); ctx.lineTo(400,200); ctx.stroke()
    ctx.fillStyle='rgba(239,68,68,0.5)'
    pts.forEach(([px,py])=>{ ctx.beginPath(); ctx.arc(px,py,2,0,2*Math.PI); ctx.fill() })
    lbl.textContent = 'Sources: ' + pts.length
  }

  function emit() { par.value=pts; par.dispatchEvent(new CustomEvent('input')) }

  function addPts(px,py) {
    for(let i=0;i<PER_PT;i++){
      const dx=(Math.random()-0.5)*2*SPREAD*SCALE
      const dy=(Math.random()-0.5)*2*SPREAD*SCALE
      const nx=px+dx, ny=py+dy
      const dx2=(nx-200)/SCALE, dy2=(ny-200)/SCALE
      if(dx2*dx2+dy2*dy2<=RAD*RAD) pts.push([Math.round(nx),Math.round(ny)])
    }
    redraw()
  }

  cvs.addEventListener('mousedown', e=>{drawing=true; addPts(e.offsetX,e.offsetY)})
  cvs.addEventListener('mousemove', e=>{if(drawing) addPts(e.offsetX,e.offsetY)})
  window.addEventListener('mouseup', ()=>{ if(drawing){ drawing=false; emit() } else drawing=false })

  par.querySelector('#clrbtn').addEventListener('click',()=>{ pts=[]; redraw(); emit() })

  par.querySelectorAll('.preset').forEach(btn => btn.addEventListener('click', ()=>{
    pts = PRESETS[btn.dataset.preset].slice()
    redraw(); emit()
  }))

  par.value = pts
  par.dispatchEvent(new CustomEvent('input'))
  redraw()
</script>
        """)
    end
end

# ╔═╡ d59342b1-cb2a-4472-9717-b82fa12335f9
@bind srcloc_pixels CanvasSourceInput()

# ╔═╡ 8939c552-a1d9-4941-aa2e-a81a761574eb
@bind compute_contributions CounterButton("Highlight stationary sources")

# ╔═╡ 1c9caca3-375a-4ee6-9586-1fe99d1f73dc
contributions = let _ = compute_contributions
    if compute_contributions == 0 || isempty(srcloc)
        Float64[]
    else
        t0 = traveltime_between_receivers
        inwindow = abs.(tgrid_xcorr .- t0) .<= lag_window
        w = [sum(abs.(cross[inwindow, i])) for i in 1:size(cross, 2)]
        m = maximum(w)
        m > 0 ? w ./ m : w
    end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FourierTools = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
MLUtils = "f1d291b0-491e-4a28-83b9-f70985020b54"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
DSP = "~0.8.5"
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
project_hash = "5fc6e362100f69ce4bd9d05ed65b8db9562b1a74"

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

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

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

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d335b2929e1b6067951a1250df247cc5fab7d40e"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.5"

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

    [deps.DSP.weakdeps]
    OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

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

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

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

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "Setfield", "SparseArrays"]
git-tree-sha1 = "2d99b4c8a7845ab1342921733fa29366dae28b24"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.1"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"
    PolynomialsRecipesBaseExt = "RecipesBase"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

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

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

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
# ╟─1083c0c6-beb2-42ba-8c99-e8d8d77ce83a
# ╟─5427df02-9278-41c9-9b3b-b5d1ab5b4c24
# ╟─c4339754-42eb-4969-9045-010ca1d23ef8
# ╟─d59342b1-cb2a-4472-9717-b82fa12335f9
# ╟─8939c552-a1d9-4941-aa2e-a81a761574eb
# ╠═b933216e-203a-47a1-b4f1-901fa5b38f70
# ╠═3c5fd716-5b57-4e6a-aa8a-022a1cb476c4
# ╟─46c9b13c-088c-47e7-8fd3-b56aee4a8923
# ╟─617fd9b1-c6be-476e-9622-b6482602c7c7
# ╟─64a67ad0-2a5f-49c4-8a90-8fe8fb55342f
# ╟─1bc0ce25-4de0-4628-a308-8106544fc0ec
# ╟─e69db63a-52b4-45c0-8b17-5239468cf5ad
# ╠═280b6269-a022-4846-9071-61481649b008
# ╟─e0ae191c-b674-46c7-bf05-bfecade313a5
# ╠═e1f1fa09-87ef-4cdb-94b1-0d2e0e483323
# ╠═1ea5bea2-0c9d-4095-a395-24ad4f67ba05
# ╟─577d411c-88b5-495d-a49a-d2379336f533
# ╠═c71bc2d6-cc33-4cae-8195-5fe26ec39d7f
# ╟─8d45d812-fb70-45c5-a123-4f93844fef6d
# ╠═ca1b63bf-c868-4200-80f0-f8d6e40d9a6c
# ╠═3b144853-0d82-46e2-a4ee-301516f3ce10
# ╠═0279d32c-d884-4990-a1a6-68b43fae2f92
# ╠═3cf54d17-5931-461c-9e35-bad1fc8d9ac3
# ╠═a0f471bd-8139-4f42-bf0e-c019e6d42832
# ╠═6098be3e-fe12-4f11-9749-b8c66c6106a8
# ╠═09c19be2-867a-4c16-ace5-6e4a9389472c
# ╠═97289274-85f2-46f7-b9bf-916f0a420301
# ╠═2bad7a74-1b12-499f-b435-7f99c7d64ae9
# ╠═1c9caca3-375a-4ee6-9586-1fe99d1f73dc
# ╟─d483d881-11c4-495f-a103-2d2dd6371ca3
# ╠═290b50eb-0faf-4a5f-86bf-36628e61ff06
# ╠═7ede27ab-f524-4f20-8bcb-5be70600c892
# ╠═e059662c-69af-460c-aada-408a365d0ff7
# ╠═41870289-3874-43ca-9837-3cead43d0200
# ╠═6a8aeed2-173a-4523-afb3-f9d77c252816
# ╠═af973cc3-21b3-4dca-91ad-0c310e1214bc
# ╠═09bd50ae-ab80-4cde-bdb7-ccc7cd7041d4
# ╟─d9727185-8f69-49a5-a98f-7514562529ab
# ╠═77d2d0d5-05e8-44e4-910b-b8bc8cf3137f
# ╠═1293b18a-15d5-4d54-8396-59899c89a172
# ╠═ccf893d9-5b6e-4a79-b09d-3da957bc612c
# ╠═a153664e-24bc-4fea-b0f4-f12360eaa688
# ╠═7fe32a18-b484-4179-b0a8-d6c49059ed59
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
