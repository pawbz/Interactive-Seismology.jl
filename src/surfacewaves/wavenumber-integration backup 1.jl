### A Pluto.jl notebook ###
# v0.20.19

#> [frontmatter]
#> title = "Wavenumber Integration Synthetics"
#> tags = ["surfacewaves", "synthetics"]
#> layout = "layout.jlhtml"
#> description = "Layered wavenumber-integration synthetics with P/SV/SH partitions."

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 01
begin
    using LinearAlgebra
    using FFTW
    using PlutoUI
    using PlutoPlotly
    using LaTeXStrings
end

# ╔═╡ 02
ChooseDisplayMode()

# ╔═╡ 03
PlutoUI.TableOfContents()

# ╔═╡ 04
md"""
# Wavenumber Integration Synthetics (Layered Reflectivity)

This notebook now uses a **layered reflectivity / propagator approach** per horizontal wavenumber `k` to build full elastodynamic synthetics (P, SV, SH) and projects them to radial/vertical components at the free surface.

**Scope & caveats (pedagogical):**
- 1-D isotropic layers with multiple scattering handled per mode; P↔SV uses Zoeppritz interface solutions with propagator recursion, SH uses scalar reflectivity.
- Free surface is enforced via the layered reflectivity seen from the top layer; evanescent parts are tapered out.
- Ricker source in time; simple takeoff radiation factors per mode (vertical or radial force surrogate).
"""

# ╔═╡ 05
md"""
## Model & Source/Receiver Controls
"""

# ╔═╡ 06
begin
    @bind preset Select(["Sediment + Crust + Mantle" => 1, "Crust only" => 2, "Half-space" => 3], default=1)
    @bind layer_text PlutoUI.Textarea(rows=6, cols=60, default="", placeholder="thickness_km,vp_km_s,vs_km_s,rho_g_cc\n1.0,2.5,1.2,2.0\n4.0,6.0,3.4,2.7\nInf,8.0,4.5,3.3")
    @bind zs PlutoUI.Slider(0.1:0.05:20.0, default=2.0, show_value=true)  # km
    @bind r PlutoUI.Slider(0.5:0.1:50.0, default=10.0, show_value=true)   # km offset
    @bind f0 PlutoUI.Slider(0.2:0.05:5.0, default=1.0, show_value=true)  # Hz
    @bind tmax PlutoUI.Slider(5.0:0.5:80.0, default=40.0, show_value=true) # s
    @bind dt PlutoUI.Slider(0.0025:0.0025:0.05, default=0.01, show_value=true) # s
    @bind nk PlutoUI.Slider(64:16:384, default=192, show_value=true)
    @bind source Select(["Vertical force" => :vf, "Radial force" => :hf], default=:vf)
end

# ╔═╡ 07
md"""
### Helper functions
"""

# ╔═╡ 08
struct Layer
    thickness::Float64   # km; use Inf for half-space
    vp::Float64          # km/s
    vs::Float64          # km/s
    rho::Float64         # g/cc
end

const DEFAULT_LAYER_TEXT = Dict(
    1 => "thickness_km,vp_km_s,vs_km_s,rho_g_cc\n1.0,2.5,1.2,2.0\n4.0,6.0,3.4,2.7\nInf,8.0,4.5,3.3",
    2 => "thickness_km,vp_km_s,vs_km_s,rho_g_cc\n10.0,6.2,3.6,2.7\nInf,7.8,4.5,3.1",
    3 => "thickness_km,vp_km_s,vs_km_s,rho_g_cc\nInf,3.5,2.0,2.7"
)

layer_table(layers) = HTML("""
<table>
  <tr><th>Layer</th><th>Thickness (km)</th><th>Vp (km/s)</th><th>Vs (km/s)</th><th>ρ (g/cc)</th></tr>
  $(join(["<tr><td>$(i)</td><td>$(isfinite(l.thickness) ? l.thickness : "∞")</td><td>$(l.vp)</td><td>$(l.vs)</td><td>$(l.rho)</td></tr>" for (i, l) in enumerate(layers)], ""))
</table>
""")

function parse_layer_text(txt)
    lines = split(coalesce(txt, ""), '\n')
    layers = Layer[]
    for ln in lines
        ln = strip(ln)
        isempty(ln) && continue
        occursin("thickness", lowercase(ln)) && continue
        parts = split(ln, ",")
        length(parts) < 4 && continue
        ths = strip(parts[1]); vp = parse(Float64, strip(parts[2])); vs = parse(Float64, strip(parts[3])); rho = parse(Float64, strip(parts[4]))
        th = lowercase(ths) == "inf" ? Inf : parse(Float64, ths)
        push!(layers, Layer(th, vp, vs, rho))
    end
    if isempty(layers)
        layers = parse_layer_text(DEFAULT_LAYER_TEXT[preset])
    end
    if !isinf(layers[end].thickness)
        last = layers[end]
        push!(layers, Layer(Inf, last.vp, last.vs, last.rho))
    end
    return layers
end

function current_layers()
    txt = isempty(layer_text) ? DEFAULT_LAYER_TEXT[preset] : layer_text
    return parse_layer_text(txt)
end

ricker(t, f0) = (1 .- 2 .* (π .* f0 .* t).^2) .* exp.(-(π .* f0 .* t).^2)

# ╔═╡ 09
function vertical_slowness(v, p)
    arg = (1 / v^2) - p^2
    return arg >= 0 ? sqrt(arg) : im * sqrt(-arg)
end

impedance(v, rho, q) = rho * v^2 * q

function scalar_reflectivity_layered(layers::Vector{Layer}, p, ω)
    R = 0 + 0im
    for i in length(layers)-1:-1:1
        top = layers[i]
        bot = layers[i+1]
        q2 = vertical_slowness(bot.vs, p)
        phase = isfinite(bot.thickness) ? exp(2im * ω * q2 * bot.thickness) : 1 + 0im
        Rprop = phase * R

        q1 = vertical_slowness(top.vs, p)
        Z1 = impedance(top.vs, top.rho, q1)
        Z2 = impedance(bot.vs, bot.rho, q2)
        r = (Z2 - Z1) / (Z2 + Z1)
        t = 2Z2 / (Z2 + Z1)
        R = r + t^2 * Rprop / (1 - r * Rprop)
    end
    return R
end

function partition_polars(k, ω, v)
    p = k / ω
    sinθ = clamp(real(p * v), -1, 1)
    cosθ = sqrt(complex(1 - sinθ^2))
    return sinθ, cosθ
end

function wavefield_components(p, qP, qS, λ, μ, mode::Symbol, direction::Symbol)
    sign = direction == :down ? 1.0 : -1.0
    if mode == :P
        ux = 1im * p
        uz = 1im * sign * qP
        sxz = -2 * μ * p * qP * sign
        szz = -(λ * (p^2 + qP^2) + 2 * μ * qP^2)
    else
        ux = -1im * sign * qS
        uz = 1im * p
        sxz = μ * (qS^2 - p^2)
        szz = -2 * μ * sign * p * qS
    end
    return ux, uz, sxz, szz
end

function solve_interface(l1::Layer, l2::Layer, p, incident_mode::Symbol, incident_side::Symbol)
    α1, β1, ρ1 = l1.vp, l1.vs, l1.rho
    α2, β2, ρ2 = l2.vp, l2.vs, l2.rho
    qP1, qS1 = vertical_slowness(α1, p), vertical_slowness(β1, p)
    qP2, qS2 = vertical_slowness(α2, p), vertical_slowness(β2, p)
    λ1, μ1 = ρ1 * (α1^2 - 2β1^2), ρ1 * β1^2
    λ2, μ2 = ρ2 * (α2^2 - 2β2^2), ρ2 * β2^2

    # incident side setup
    inc_layer, ref_layer, tran_layer = incident_side == :top ? (l1, l1, l2) : (l2, l2, l1)
    qP_inc, qS_inc = incident_side == :top ? (qP1, qS1) : (qP2, qS2)
    λ_inc, μ_inc = incident_side == :top ? (λ1, μ1) : (λ2, μ2)
    qP_tr, qS_tr = incident_side == :top ? (qP2, qS2) : (qP1, qS1)
    λ_tr, μ_tr = incident_side == :top ? (λ2, μ2) : (λ1, μ1)
    dir_inc = incident_side == :top ? :down : :up
    dir_ref = incident_side == :top ? :up : :down
    dir_tr = incident_side == :top ? :down : :up

    inc = wavefield_components(p, qP_inc, qS_inc, λ_inc, μ_inc, incident_mode, dir_inc)
    refP = wavefield_components(p, qP_inc, qS_inc, λ_inc, μ_inc, :P, dir_ref)
    refS = wavefield_components(p, qP_inc, qS_inc, λ_inc, μ_inc, :SV, dir_ref)
    trP = wavefield_components(p, qP_tr, qS_tr, λ_tr, μ_tr, :P, dir_tr)
    trS = wavefield_components(p, qP_tr, qS_tr, λ_tr, μ_tr, :SV, dir_tr)

    A = zeros(ComplexF64, 4, 4)
    b = -collect(inc)
    for (row, comps) in enumerate((refP, refS, trP, trS))
        # columns: RP, RS, TP, TS
    end
    A[:, 1] = collect(refP)
    A[:, 2] = collect(refS)
    A[:, 3] = -collect(trP)
    A[:, 4] = -collect(trS)

    x = A \ b
    R = x[1:2]
    T = x[3:4]
    return R, T
end

function zoeppritz_interface(l1::Layer, l2::Layer, p)
    Rdown = zeros(ComplexF64, 2, 2)
    Tdown = zeros(ComplexF64, 2, 2)
    Rup = zeros(ComplexF64, 2, 2)
    Tup = zeros(ComplexF64, 2, 2)
    for (j, mode) in enumerate((:P, :SV))
        r1, t1 = solve_interface(l1, l2, p, mode, :top)
        r2, t2 = solve_interface(l1, l2, p, mode, :bottom)
        Rdown[:, j] = r1
        Tdown[:, j] = t1
        Rup[:, j] = r2
        Tup[:, j] = t2
    end
    return Rdown, Tdown, Rup, Tup
end

phase_matrix(layer::Layer, ω, p) = Diagonal(exp.(2im .* ω .* [vertical_slowness(layer.vp, p), vertical_slowness(layer.vs, p)] .* (isfinite(layer.thickness) ? layer.thickness : 0.0)))

# ╔═╡ 10.5
md"""
### Layer stack in use
"""

# ╔═╡ 10.6
layer_summary = layer_table(current_layers())

# ╔═╡ 10.7
md"""
Paste or edit CSV lines as `thickness_km,vp_km_s,vs_km_s,rho_g_cc`. Leave blank to use the preset.
"""

# ╔═╡ 10
function wavenumber_spectrum_layered(layers, zs, r, f0; tmax=40.0, dt=0.01, nk=192)
    nt = Int(ceil(tmax / dt))
    t = collect(0:dt:(nt-1) * dt)
    nf = nt
    df = 1 / (nt * dt)
    freqs = collect(0:df:df * (nf - 1))
    ω = 2π .* freqs

    s = ricker.(t .- 1.5 / f0, f0)
    Sω = fft(s)

    vmax = maximum(getfield.(layers, :vp))
    kmax = maximum(ω) / vmax * 1.5
    k = range(0, kmax, length = nk)

    Gω_r = zeros(ComplexF64, nf)
    Gω_z = zeros(ComplexF64, nf)
    Gω_parts = Dict(:P => zeros(ComplexF64, nf), :SV => zeros(ComplexF64, nf), :SH => zeros(ComplexF64, nf))

    for (ii, ωi) in enumerate(ω)
        if ωi == 0
            continue
        end
        for (jj, kj) in enumerate(k)
            dk = jj == 1 ? k[2] - k[1] : k[jj] - k[jj-1]
            p = kj / ωi

            # skip evanescent beyond top vp
            if real(p * vmax) > 1.2
                continue
            end

            # build P-SV reflection matrix via Zoeppritz recursion
            Rpsv = zeros(ComplexF64, 2, 2)
            for jjlayer in length(layers)-1:-1:1
                bot = layers[jjlayer+1]
                Rprop = Rpsv
                if isfinite(bot.thickness)
                    W = phase_matrix(bot, ωi, p)
                    Rprop = W * Rprop * W
                end
                Rd, Td, Ru, Tu = zoeppritz_interface(layers[jjlayer], bot, p)
                denom = I - Ru * Rprop
                Rpsv = Rd + Td * Rprop * (denom \ Tu)
            end

            qP = vertical_slowness(layers[1].vp, p)
            qS = vertical_slowness(layers[1].vs, p)
            sinP, cosP = partition_polars(kj, ωi, layers[1].vp)
            sinS, cosS = partition_polars(kj, ωi, layers[1].vs)

            w = source == :vf ? [cosP, -sinS] : [sinP, cosS]
            phase_up = exp.(im .* ωi .* [qP, qS] .* zs)
            amp_vec = phase_up .* (w .+ Rpsv * (phase_up .* w))
            ampP, ampS = amp_vec[1], amp_vec[2]

            # SH (scalar) with layered propagation
            RSH = scalar_reflectivity_layered(layers, p, ωi)
            qSH = qS
            phaseSH = exp(im * ωi * qSH * zs)
            wSH = source == :hf ? 1.0 : 0.3
            ampSH = phaseSH * (wSH + RSH * phaseSH * wSH) / (2qSH + 1e-12)

            geom = exp(im * kj * r)
            Gω_r[ii] += geom * (ampP * sinP + ampS * cosS) * dk
            Gω_z[ii] += geom * (ampP * cosP - ampS * sinS) * dk
            Gω_parts[:P][ii] += geom * ampP * dk
            Gω_parts[:SV][ii] += geom * ampS * dk
            Gω_parts[:SH][ii] += geom * ampSH * dk
        end
    end

    ur = real(ifft(Gω_r .* Sω))
    uz = real(ifft(Gω_z .* Sω))
    parts_time = Dict(key => real(ifft(val .* Sω)) for (key, val) in Gω_parts)

    return t, ur, uz, freqs, abs.(Gω_r), abs.(Gω_z), parts_time
end

# ╔═╡ 11
md"""
### Compute synthetic
"""

# ╔═╡ 12
begin
    layers = current_layers()
    t, ur, uz, freqs, Gmag_r, Gmag_z, parts = wavenumber_spectrum_layered(layers, zs, r, f0; tmax=tmax, dt=dt, nk=nk)
end

# ╔═╡ 13
md"""
## Plots
"""

# ╔═╡ 14
begin
    fig1 = Plot(title="Surface displacement components", xaxis_title="Time (s)", yaxis_title="Amplitude", width=760, height=320)
    add_trace!(fig1, PlutoPlotly.scatter(x=t, y=ur, mode="lines", name="Radial", line=attr(color="aqua")))
    add_trace!(fig1, PlutoPlotly.scatter(x=t, y=uz, mode="lines", name="Vertical", line=attr(color="orange")))
    PlutoPlotly.plot(fig1)
end

# ╔═╡ 15
begin
    fig2 = Plot(title="|G(ω)| amplitude (radial vs vertical)", xaxis_title="Frequency (Hz)", yaxis_title="|G|", width=760, height=320)
    add_trace!(fig2, PlutoPlotly.scatter(x=freqs, y=Gmag_r, mode="lines", name="Radial", line=attr(color="aqua")))
    add_trace!(fig2, PlutoPlotly.scatter(x=freqs, y=Gmag_z, mode="lines", name="Vertical", line=attr(color="orange")))
    PlutoPlotly.plot(fig2)
end

# ╔═╡ 15.5
begin
    fig3 = Plot(title="Mode partitions (time)", xaxis_title="Time (s)", yaxis_title="Amplitude", width=760, height=320)
    add_trace!(fig3, PlutoPlotly.scatter(x=t, y=parts[:P], mode="lines", name="P", line=attr(color="lime")))
    add_trace!(fig3, PlutoPlotly.scatter(x=t, y=parts[:SV], mode="lines", name="SV", line=attr(color="fuchsia")))
    add_trace!(fig3, PlutoPlotly.scatter(x=t, y=parts[:SH], mode="lines", name="SH", line=attr(color="gold")))
    PlutoPlotly.plot(fig3)
end

# ╔═╡ 16
md"""
### Notes
- P–SV coupling now uses Zoeppritz interface solutions plus propagator recursion; SH still uses scalar reflectivity with layered multiples.
- Radiation weights are simple surrogates for vertical or radial forces; extend if you need double-couple or moment-tensor sources.
- Evanescent components (k > ω/v) are skipped; increase `nk` and reduce `dt` for smoother spectra and higher max frequency.
- The interface solve is numeric per `(ω, k)`; for heavy models consider reducing `nk` or `tmax` to speed up.
"""

# ╔═╡ 99
default_plotly_template(:plotly_dark)
