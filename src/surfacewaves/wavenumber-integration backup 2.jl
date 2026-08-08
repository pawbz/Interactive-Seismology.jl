### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "Wavenumber Integration Synthetics"
#> tags = ["surfacewaves", "synthetics"]
#> layout = "layout.jlhtml"
#> description = "Layered wavenumber-integration synthetics with P/SV/SH partitions."

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

# ╔═╡ 207242a6-da21-11f0-b867-df46d1dd611e
begin
    using LinearAlgebra
    using FFTW
    using PlutoUI
    using PlutoPlotly
    using LaTeXStrings
    using HypertextLiteral: @htl
    using SpecialFunctions: besselj0, besselj1
end

# ╔═╡ 2072528c-da21-11f0-937e-d1687ff3f5bb
PlutoUI.TableOfContents()

# ╔═╡ 207252c8-da21-11f0-af4c-e19a10c55b26
md"""
# Wavenumber Integration Synthetics (Layered Reflectivity)

This notebook now uses a **layered reflectivity / propagator approach** per horizontal wavenumber `k` to build full elastodynamic synthetics (P, SV, SH) and projects them to radial/vertical components at the free surface.

**Scope & caveats (pedagogical):**
- 1-D isotropic layers with multiple scattering handled per mode; P↔SV uses Zoeppritz interface solutions with propagator recursion, SH uses scalar reflectivity.
- Free surface is enforced via the layered reflectivity seen from the top layer; evanescent parts are tapered out.
- Ricker source in time; simple takeoff radiation factors per mode (vertical or radial force surrogate).
"""

# ╔═╡ 2072543a-da21-11f0-b6ff-796f0819f12f
md"""
## Model & Source/Receiver Controls
"""

# ╔═╡ e629fb17-8918-4f58-acd7-87af8e767b94
@bind zs PlutoUI.Slider(0.1:0.05:20.0, default=2.0, show_value=true)  # km

# ╔═╡ 5a44cc46-f2c6-4621-a421-9557b6957890
@bind r PlutoUI.Slider(0.5:0.1:50.0, default=10.0, show_value=true)   # km offset

# ╔═╡ 972ee19b-1347-4541-b89d-7700e1fe27c1
@bind f0 PlutoUI.Slider(0.2:0.05:5.0, default=1.0, show_value=true)  # Hz

# ╔═╡ d37b6619-4b8c-4e84-9f18-ec0ff6bcd2e8
# Reference frequency for causal-Q dispersion
@bind f_ref_hz PlutoUI.Slider(0.05:0.05:5.0, default=1.0, show_value=true)

# ╔═╡ f0d706d4-8c5c-4924-bf8a-c676c87c1242
@bind tmax PlutoUI.Slider(5.0:0.5:80.0, default=40.0, show_value=true) # s

# ╔═╡ 675a3895-4872-415c-9f9e-083a7ad58edf
@bind dt PlutoUI.Slider(0.0025:0.0025:0.05, default=0.01, show_value=true) # s

# ╔═╡ 9a90eb4e-dc62-49a3-8031-722ef80e2b0a
@bind nk PlutoUI.Slider(64:16:384, default=192, show_value=true)

# ╔═╡ a2a69fab-137f-49dd-853e-75942124a5d0
@bind source Select(["Vertical force" => :vf, "Radial force" => :hf], default=:vf)

# ╔═╡ ea7fd230-40c8-4c89-9112-8b29a4c1668a
@bind n_layers Slider(1:1:12, default=3, show_value=true)

# ╔═╡ 207255d4-da21-11f0-a57c-1bd673807028
md"""
### Helper functions
"""

# ╔═╡ 6e2be5f4-5f30-4c6f-be89-f0c5a5448dc2
begin
	struct Layer
	    thickness::Float64   # km
	    vp::Float64          # km/s
	    vs::Float64          # km/s
	    rho::Float64         # g/cm³
	    Qp::Float64          # P-wave quality factor
	    Qs::Float64          # S-wave quality factor
	end
	Layer(th, vp, vs, rho) = Layer(th, vp, vs, rho, 1000.0, 1000.0)
end

# causal constant-Q dispersion with reference frequency
function causal_velocity(v::Real, Q::Real, ω::Real, ω_ref::Real)
    if Q <= 0 || !isfinite(Q)
        return complex(v)
    end
    ωref = max(ω_ref, 1e-5)
    α = atan(1 / Q) / π
    mag = (abs(ω) / ωref)^α
    return v * mag * (1 - 1im / (2Q))
end

# clamp evanescent exponentials to avoid overflow
function stable_phase(q::Complex, ω::Real, h::Real; clip=80.0)
    if !isfinite(h)
        return one(q)
    end
    atten = clamp(-2 * ω * imag(q) * h, -clip, clip)
    phase = exp(2im * ω * real(q) * h)
    return phase * exp(atten)
end

# Hankel / Bessel switch with asymptotic form for large kr
function hankel_factor(order::Int, kr::Real; switch_at=30.0)
    if kr < switch_at
        return order == 0 ? besselj0(kr) : besselj1(kr)
    else
        return sqrt(2 / (π * kr)) * cos(kr - order * π / 2 - π / 4)
    end
end

# ╔═╡ 8c6209b8-8949-4aba-8e4e-19c8e5e26a90
layer_table(layers) = HTML("""
<table>
  <tr><th>Layer</th><th>Thickness (km)</th><th>Vp (km/s)</th><th>Vs (km/s)</th><th>ρ (g/cc)</th></tr>
  $(join(["<tr><td>$(i)</td><td>$(isfinite(l.thickness) ? l.thickness : "∞")</td><td>$(l.vp)</td><td>$(l.vs)</td><td>$(l.rho)</td></tr>" for (i, l) in enumerate(layers)], ""))
</table>
""")


# ╔═╡ a3390ca6-cc58-44a7-a79c-ac684a3ad92b
const GUTENBERG_MODEL = [
    Layer(19.0, 6.14, 3.55, 2.74, 1000.0, 1000.0),
    Layer(19.0, 6.58, 3.80, 3.00, 1000.0, 1000.0),
    Layer(12.0, 8.20, 4.65, 3.32, 1000.0, 1000.0),
    Layer(10.0, 8.17, 4.62, 3.34, 1000.0, 1000.0),
    Layer(10.0, 8.14, 4.57, 3.35, 1000.0, 1000.0),
    Layer(10.0, 8.10, 4.51, 3.36, 1000.0, 1000.0),
    Layer(10.0, 8.07, 4.46, 3.37, 1000.0, 1000.0),
    Layer(10.0, 8.02, 4.41, 3.38, 1000.0, 1000.0),
    Layer(25.0, 7.93, 4.37, 3.39, 1000.0, 1000.0),
    Layer(25.0, 7.85, 4.35, 3.41, 1000.0, 1000.0),
    Layer(25.0, 7.89, 4.36, 3.43, 1000.0, 1000.0),
    Layer(25.0, 7.98, 4.38, 3.46, 1000.0, 1000.0),
    Layer(25.0, 8.10, 4.42, 3.48, 1000.0, 1000.0),
    Layer(25.0, 8.21, 4.46, 3.50, 1000.0, 1000.0),
    Layer(50.0, 8.38, 4.54, 3.53, 1000.0, 1000.0),
    Layer(50.0, 8.62, 4.68, 3.58, 1000.0, 1000.0),
    Layer(50.0, 8.87, 4.85, 3.62, 1000.0, 1000.0),
    Layer(50.0, 9.15, 5.04, 3.69, 1000.0, 1000.0),
    Layer(50.0, 9.45, 5.21, 3.82, 1000.0, 1000.0),
    Layer(100.0, 9.88, 5.45, 4.01, 1000.0, 1000.0),
    Layer(100.0, 10.30, 5.76, 4.21, 1000.0, 1000.0),
    Layer(100.0, 10.71, 6.03, 4.40, 1000.0, 1000.0),
    Layer(100.0, 11.10, 6.23, 4.56, 1000.0, 1000.0),
    Layer(100.0, 11.35, 6.32, 4.63, 1000.0, 1000.0)
]

# ╔═╡ a83d06de-cd27-40c7-adf4-c82484323517
function layer_table_input(n_layers::Int, defaults::Vector{Layer}; vp_vs_ratio=1.73)
    ui = PlutoUI.combine() do Child
        header = @htl("""
        <tr style=\"text-align:center;\">
          <th style=\"padding:6px;\">#</th>
          <th style=\"padding:6px;\">Thickness (km)</th>
          <th style=\"padding:6px;\">Vp (km/s)</th>
          <th style=\"padding:6px;\">Vs (km/s)</th>
          <th style=\"padding:6px;\">Density (g/cc)</th>
        </tr>
        """)

        rows = Any[]
        for i in 1:n_layers
            layer = i <= length(defaults) ? defaults[i] : defaults[end]
            tw = i < n_layers ? Child("thickness_$i", Slider(0.25:0.25:60.0, default=layer.thickness, show_value=true)) : "∞"
            vpw = Child("vp_$i", Slider(2.0:0.05:9.0, default=layer.vp, show_value=true))
            vsw = Child("vs_$i", Slider(0.5:0.05:5.5, default=layer.vs, show_value=true))
            rw = Child("rho_$i", Slider(1.5:0.05:4.0, default=layer.rho, show_value=true))

            push!(rows, @htl("""
                <tr>
                  <td style=\"text-align:center; padding:6px;\"><b>Layer $i</b></td>
                  <td style=\"padding:6px; text-align:center;\">$tw</td>
                  <td style=\"padding:6px; text-align:center;\">$vpw</td>
                  <td style=\"padding:6px; text-align:center;\">$vsw</td>
                  <td style=\"padding:6px; text-align:center;\">$rw</td>
                </tr>
            """))
        end

        tbl = @htl("""
        <table style=\"border-collapse:collapse; border:1px solid #ddd; width:100%;\">
          <thead>$header</thead>
          <tbody>
            $(rows...)
          </tbody>
        </table>
        """)

        md"""
        ### Layer Editor (slider mode)
        Number of layers: **$n_layers**

        $tbl
        _Last layer is a half-space (∞ thickness)._"""
    end

    return PlutoUI.Experimental.transformed_value(ui) do vals
        layers = [
            begin
                if i < n_layers
                    t = vals[4*(i-1)+1]
                    vp = vals[4*(i-1)+2]
                    vs = vals[4*(i-1)+3]
                    ρ = vals[4*(i-1)+4]
                else
                    t = Inf
                    vp = vals[4*(i-1)+1]
                    vs = vals[4*(i-1)+2]
                    ρ = vals[4*(i-1)+3]
                end
                Layer(t, vp, vs, ρ)
            end
            for i in 1:n_layers
        ]
        return layers
    end
end



# ╔═╡ 3d909011-6144-4ad2-b731-b149e780b595
@bind layers_ui layer_table_input(n_layers, GUTENBERG_MODEL)

# ╔═╡ 31eb2be7-fc30-4a94-940a-076ce813892e
ricker(t, f0) = (1 .- 2 .* (π .* f0 .* t).^2) .* exp.(-(π .* f0 .* t).^2)

# ╔═╡ 920b0e59-fda3-4126-ae8d-215ebfc359f1
function vertical_slowness(v, p)
    arg = (1 / v^2) - p^2
    return arg >= 0 ? sqrt(arg) : im * sqrt(-arg)
end

# ╔═╡ f093524a-89be-49bd-999d-ea9cde2fd271
impedance(v, rho, q) = rho * v^2 * q

# ╔═╡ f32c42a2-47ff-4456-9de4-72b410abc05a
function scalar_reflectivity_layered(layers::Vector{Layer}, p, ω, ω_ref)
    R = 0 + 0im
    for i in length(layers)-1:-1:1
        top = layers[i]
        bot = layers[i+1]
        vs2 = causal_velocity(bot.vs, bot.Qs, ω, ω_ref)
        q2 = vertical_slowness(vs2, p)
        phase = stable_phase(q2, ω, bot.thickness)
        Rprop = phase * R

        vs1 = causal_velocity(top.vs, top.Qs, ω, ω_ref)
        q1 = vertical_slowness(vs1, p)
        Z1 = impedance(vs1, top.rho, q1)
        Z2 = impedance(vs2, bot.rho, q2)
        r = (Z2 - Z1) / (Z2 + Z1)
        t = 2Z2 / (Z2 + Z1)
        R = r + t^2 * Rprop / (1 - r * Rprop)
    end
    return R
end

# ╔═╡ d3b53c8d-a5da-4efd-96c2-27615e476d33
function partition_polars(k, ω, v)
    p = k / ω
    sinθ = clamp(real(p * v), -1, 1)
    cosθ = sqrt(complex(1 - sinθ^2))
    return sinθ, cosθ
end

# ╔═╡ f60f54f2-0937-4c90-83a9-27560d27e829
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

# ╔═╡ 1df9fb85-8abe-423a-b181-6fbb36fe9393
function solve_interface(l1::Layer, l2::Layer, p, incident_mode::Symbol, incident_side::Symbol, ω, ω_ref)
    α1 = causal_velocity(l1.vp, l1.Qp, ω, ω_ref)
    β1 = causal_velocity(l1.vs, l1.Qs, ω, ω_ref)
    ρ1 = l1.rho
    α2 = causal_velocity(l2.vp, l2.Qp, ω, ω_ref)
    β2 = causal_velocity(l2.vs, l2.Qs, ω, ω_ref)
    ρ2 = l2.rho
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

# ╔═╡ 940d0e21-b7dd-419e-a5e0-1c088944122a
function zoeppritz_interface(l1::Layer, l2::Layer, p, ω, ω_ref)
    Rdown = zeros(ComplexF64, 2, 2)
    Tdown = zeros(ComplexF64, 2, 2)
    Rup = zeros(ComplexF64, 2, 2)
    Tup = zeros(ComplexF64, 2, 2)
    for (j, mode) in enumerate((:P, :SV))
        r1, t1 = solve_interface(l1, l2, p, mode, :top, ω, ω_ref)
        r2, t2 = solve_interface(l1, l2, p, mode, :bottom, ω, ω_ref)
        Rdown[:, j] = r1
        Tdown[:, j] = t1
        Rup[:, j] = r2
        Tup[:, j] = t2
    end
    return Rdown, Tdown, Rup, Tup
end

# ╔═╡ a4ee6146-3264-4bd8-ad4e-885854f78f9c
function phase_matrix(layer::Layer, ω, p, ω_ref)
    vp = causal_velocity(layer.vp, layer.Qp, ω, ω_ref)
    vs = causal_velocity(layer.vs, layer.Qs, ω, ω_ref)
    qP = vertical_slowness(vp, p)
    qS = vertical_slowness(vs, p)
    if isfinite(layer.thickness)
        return Diagonal([
            stable_phase(qP, ω, layer.thickness),
            stable_phase(qS, ω, layer.thickness)
        ])
    else
        return Diagonal([one(qP), one(qS)])
    end
end

# ╔═╡ 20725dfe-da21-11f0-9f21-cdbf80481f9a
md"""
### Layer stack in use
"""

# ╔═╡ 20725ea8-da21-11f0-9ad6-3f27a5989be9
md"""
Use the slider table below to configure the layered half-space (thickness in km, velocities in km/s, density in g/cc).
"""

# ╔═╡ 20725ef0-da21-11f0-bb8c-132d195acf15
function wavenumber_spectrum_layered(layers, zs, r, f0; tmax=40.0, dt=0.01, nk=192, f_ref=1.0)
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

            # skip evanescent beyond top vp (taper) to avoid blow-up
            if real(p * vmax) > 1.2
                continue
            end

            # build P-SV reflection matrix via Zoeppritz recursion
            Rpsv = zeros(ComplexF64, 2, 2)
            for jjlayer in length(layers)-1:-1:1
                bot = layers[jjlayer+1]
                Rprop = Rpsv
                if isfinite(bot.thickness)
                    W = phase_matrix(bot, ωi, p, f_ref)
                    Rprop = W * Rprop * W
                end
                Rd, Td, Ru, Tu = zoeppritz_interface(layers[jjlayer], bot, p, ωi, f_ref)
                denom = I - Ru * Rprop
                Rpsv = Rd + Td * Rprop * (denom \ Tu)
            end

            vp_top = causal_velocity(layers[1].vp, layers[1].Qp, ωi, f_ref)
            vs_top = causal_velocity(layers[1].vs, layers[1].Qs, ωi, f_ref)
            qP = vertical_slowness(vp_top, p)
            qS = vertical_slowness(vs_top, p)
            sinP, cosP = partition_polars(kj, ωi, vp_top)
            sinS, cosS = partition_polars(kj, ωi, vs_top)

            w = source == :vf ? [cosP, -sinS] : [sinP, cosS]
            phase_up = exp.(im .* ωi .* [qP, qS] .* zs)
            amp_vec = phase_up .* (w .+ Rpsv * (phase_up .* w))
            ampP, ampS = amp_vec[1], amp_vec[2]

            # SH (scalar) with layered propagation
            RSH = scalar_reflectivity_layered(layers, p, ωi, f_ref)
            qSH = qS
            phaseSH = exp(im * ωi * qSH * zs)
            wSH = source == :hf ? 1.0 : 0.3
            ampSH = phaseSH * (wSH + RSH * phaseSH * wSH) / (2qSH + 1e-12)

            kr = kj * r
            weight = kj * dk
            J0 = hankel_factor(0, kr)
            J1 = hankel_factor(1, kr)
            geom_r = J1
            geom_z = J0
            taper = exp(- (kr / 120)^2)  # asymptotic taper to limit noise

            Gω_r[ii] += taper * geom_r * (ampP * sinP + ampS * cosS) * weight
            Gω_z[ii] += taper * geom_z * (ampP * cosP - ampS * sinS) * weight
            Gω_parts[:P][ii] += taper * geom_z * ampP * weight
            Gω_parts[:SV][ii] += taper * geom_z * ampS * weight
            Gω_parts[:SH][ii] += taper * geom_z * ampSH * weight
        end
    end

    ur = real(ifft(Gω_r .* Sω))
    uz = real(ifft(Gω_z .* Sω))
    parts_time = Dict(key => real(ifft(val .* Sω)) for (key, val) in Gω_parts)

    return t, ur, uz, freqs, abs.(Gω_r), abs.(Gω_z), parts_time
end

# ╔═╡ 207262ae-da21-11f0-9083-3faa73ae3ee3
md"""
### Compute synthetic
"""

# ╔═╡ 207262e0-da21-11f0-98b0-ed37d7dc065b
begin
    ω_ref = 2π * f_ref_hz
    t, ur, uz, freqs, Gmag_r, Gmag_z, parts = wavenumber_spectrum_layered(layers_ui, zs, r, f0; tmax=tmax, dt=dt, nk=nk, f_ref=ω_ref)
end

# ╔═╡ 207263b2-da21-11f0-b168-a3abcfad3d26
begin
    fig1 = Plot(title="Surface displacement components", xaxis_title="Time (s)", yaxis_title="Amplitude", width=760, height=320)
    add_trace!(fig1, PlutoPlotly.scatter(x=t, y=ur, mode="lines", name="Radial", line=attr(color="aqua")))
    add_trace!(fig1, PlutoPlotly.scatter(x=t, y=uz, mode="lines", name="Vertical", line=attr(color="orange")))
    PlutoPlotly.plot(fig1)
end

# ╔═╡ 20726388-da21-11f0-9014-3ddaaac11fc5
md"""
## Plots
"""

# ╔═╡ 20726484-da21-11f0-adc1-f5f9668ccf2f
begin
    fig2 = Plot(title="|G(ω)| amplitude (radial vs vertical)", xaxis_title="Frequency (Hz)", yaxis_title="|G|", width=760, height=320)
    add_trace!(fig2, PlutoPlotly.scatter(x=freqs, y=Gmag_r, mode="lines", name="Radial", line=attr(color="aqua")))
    add_trace!(fig2, PlutoPlotly.scatter(x=freqs, y=Gmag_z, mode="lines", name="Vertical", line=attr(color="orange")))
    PlutoPlotly.plot(fig2)
end

# ╔═╡ 2072652e-da21-11f0-bd03-cf2dad4b91c3
begin
    fig3 = Plot(title="Mode partitions (time)", xaxis_title="Time (s)", yaxis_title="Amplitude", width=760, height=320)
    add_trace!(fig3, PlutoPlotly.scatter(x=t, y=parts[:P], mode="lines", name="P", line=attr(color="lime")))
    add_trace!(fig3, PlutoPlotly.scatter(x=t, y=parts[:SV], mode="lines", name="SV", line=attr(color="fuchsia")))
    add_trace!(fig3, PlutoPlotly.scatter(x=t, y=parts[:SH], mode="lines", name="SH", line=attr(color="gold")))
    PlutoPlotly.plot(fig3)
end

# ╔═╡ 20726626-da21-11f0-80ab-1bf86491b592
md"""
### Notes
- P–SV uses Zoeppritz plus a compound propagator with overflow-safe phases; SH uses scalar reflectivity with the same causal-Q dispersion.
- Causal Q is applied via a reference frequency slider; set `f_ref_hz` to your dispersion reference.
- Hankel/Bessel kernels (J0/J1) are used with asymptotic switching for large `kr`, including a high-kr taper.
- Evanescent components (k > ω/v) are skipped; increase `nk` and reduce `dt` for smoother spectra and higher max frequency.
- Radiation weights are simple surrogates for vertical or radial forces; extend if you need double-couple or moment-tensor sources.
- The interface solve is numeric per `(ω, k)`; for heavy models consider reducing `nk` or `tmax` to speed up.
"""

# ╔═╡ 207266f0-da21-11f0-8464-e9994f93792f
default_plotly_template(:plotly_dark)

# ╔═╡ Cell order:
# ╠═207242a6-da21-11f0-b867-df46d1dd611e
# ╠═2072528c-da21-11f0-937e-d1687ff3f5bb
# ╠═207252c8-da21-11f0-af4c-e19a10c55b26
# ╠═2072543a-da21-11f0-b6ff-796f0819f12f
# ╠═207263b2-da21-11f0-b168-a3abcfad3d26
# ╠═e629fb17-8918-4f58-acd7-87af8e767b94
# ╠═5a44cc46-f2c6-4621-a421-9557b6957890
# ╠═972ee19b-1347-4541-b89d-7700e1fe27c1
# ╠═d37b6619-4b8c-4e84-9f18-ec0ff6bcd2e8
# ╠═f0d706d4-8c5c-4924-bf8a-c676c87c1242
# ╠═675a3895-4872-415c-9f9e-083a7ad58edf
# ╠═9a90eb4e-dc62-49a3-8031-722ef80e2b0a
# ╠═a2a69fab-137f-49dd-853e-75942124a5d0
# ╠═ea7fd230-40c8-4c89-9112-8b29a4c1668a
# ╠═207255d4-da21-11f0-a57c-1bd673807028
# ╠═6e2be5f4-5f30-4c6f-be89-f0c5a5448dc2
# ╠═8c6209b8-8949-4aba-8e4e-19c8e5e26a90
# ╠═a3390ca6-cc58-44a7-a79c-ac684a3ad92b
# ╠═a83d06de-cd27-40c7-adf4-c82484323517
# ╠═3d909011-6144-4ad2-b731-b149e780b595
# ╠═31eb2be7-fc30-4a94-940a-076ce813892e
# ╠═920b0e59-fda3-4126-ae8d-215ebfc359f1
# ╠═f093524a-89be-49bd-999d-ea9cde2fd271
# ╠═f32c42a2-47ff-4456-9de4-72b410abc05a
# ╠═d3b53c8d-a5da-4efd-96c2-27615e476d33
# ╠═f60f54f2-0937-4c90-83a9-27560d27e829
# ╠═1df9fb85-8abe-423a-b181-6fbb36fe9393
# ╠═940d0e21-b7dd-419e-a5e0-1c088944122a
# ╠═a4ee6146-3264-4bd8-ad4e-885854f78f9c
# ╠═20725dfe-da21-11f0-9f21-cdbf80481f9a
# ╠═20725ea8-da21-11f0-9ad6-3f27a5989be9
# ╠═20725ef0-da21-11f0-bb8c-132d195acf15
# ╠═207262ae-da21-11f0-9083-3faa73ae3ee3
# ╠═207262e0-da21-11f0-98b0-ed37d7dc065b
# ╠═20726388-da21-11f0-9014-3ddaaac11fc5
# ╠═20726484-da21-11f0-adc1-f5f9668ccf2f
# ╠═2072652e-da21-11f0-bd03-cf2dad4b91c3
# ╠═20726626-da21-11f0-80ab-1bf86491b592
# ╠═207266f0-da21-11f0-8464-e9994f93792f
