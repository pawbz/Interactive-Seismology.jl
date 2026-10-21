### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> tags = ["introduction"]
#> title = "Wave Dispersion: Phase and Group Velocity"
#> description = "Why a wave packet and the crests inside it travel at different speeds, and what dispersion does to a pulse."
#> layout = "layout.jlhtml"

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

# ╔═╡ 59057ef4-796a-4f06-b1f6-5a74ab535b1d
using PlutoUI

# ╔═╡ 7a1fa0a3-b612-40da-8f9d-f1494877560c
TableOfContents()

# ╔═╡ bef462ee-988d-49b5-acc8-599d4c401067
# App shell: the Expt / Theory / Code / Help frame shared by every notebook, loaded
# from src/assets/app-shell/. Without that folder (a lone downloaded .jl) this cell just
# shows the help text and the notebook stays a plain notebook.
let
	help = md"""
	- **Medium** chooses the dispersion relation: non-dispersive, deep water, a cutoff medium, or a Love wave in a layer.
	- **k₀** slides the carrier along the ω(k) curve; **Bandwidth** sets how many wavenumbers the packet contains.
	- **Time** scrubs the packet; **Play** animates it; **Receiver** moves the station that records the seismogram.
	- Watch the red envelope (group velocity) and the crests (phase velocity) separate.

	**Try this:** pick the Love wave, play the time, and see the long waves arrive at the receiver before the short ones.
	"""
	dir = joinpath(@__DIR__, "..", "assets", "app-shell")
	if isfile(joinpath(dir, "app-shell.js"))
		HTML(string(
			"<style>", read(joinpath(dir, "app-shell.css"), String), "</style>",
			"<template class=\"app-help\">", repr(MIME"text/html"(), help), "</template>",
			"<script>", read(joinpath(dir, "app-shell.js"), String), "</script>"))
	else
		help
	end
end

# ╔═╡ 7e344ce5-8f4f-4a52-bc21-dea0934d900d
md"""
# Wave Dispersion: Phase and Group Velocity

A surface wave from a distant earthquake does not arrive as one clean pulse. The long
periods come first and the short periods straggle in behind them, so the record looks
like a stretched-out chirp. The cause is **dispersion**: waves of different wavelength
travel at different speeds, so a pulse that starts out compact spreads as it travels.

Dispersion also splits the notion of "the wave speed" in two. The crests inside a packet
move at the **phase velocity** `` v_p ``, while the packet as a whole, and the energy it
carries, moves at the **group velocity** `` v_g ``. In a non-dispersive medium the two
are equal; in a dispersive one they are not.

Use the experiment. Choose a medium, pick the wavenumber `` k_0 `` of the carrier and the
bandwidth of the packet, and scrub or play the time. The top-left panel shows the
dispersion relation `` \omega(k) `` with the two slopes that define the velocities; the
top-right panel shows the packet itself; the bottom panel is a seismogram at a receiver.

A companion notebook, *Fourier Series and the Fourier Transform*, explains why a packet
is a sum of frequencies; here we let each of those frequencies travel at its own speed.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ aab4ff38-2775-4cff-b45f-b189eb16fbf7
md"""
## How to read the experiment

- **Dispersion relation** (top left). The curve is `` \omega(k) ``. The **blue dashed
  line** joins the origin to the carrier point `` (k_0,\omega_0) ``: its slope is the phase
  velocity `` v_p=\omega_0/k_0 ``. The **red line** is the tangent to the curve at that
  point: its slope is the group velocity `` v_g=d\omega/dk ``. The shaded band is the range
  of wavenumbers in the packet. In the non-dispersive medium the curve is a straight line
  and the two slopes coincide. Move `` k_0 `` to slide the point along the curve and watch
  the two slopes change.
- **Packet** (top right). The blue curve is the wave, the red dashed curve is its
  envelope. The red vertical line marks `` x=v_g t `` and the blue vertical line marks
  `` x=v_p t ``. Play the time: the envelope follows the red line, while individual crests
  slip through it at the speed difference `` v_p-v_g ``.
- **Seismogram** (bottom). The wave recorded at the receiver `` x_r ``, with the arrival
  time `` x_r/v_g `` marked. In a dispersive medium the record is stretched and its
  oscillation period changes through the arrival.
"""

# ╔═╡ aad6728c-ee0f-4063-9f8a-00a1f86d94f4
md"""
## Two velocities from two waves

Add two plane waves of nearly equal wavenumber and frequency:

```math
\cos(k_1x-\omega_1t)+\cos(k_2x-\omega_2t)
= 2\cos\!\Big(\tfrac{\Delta k}{2}x-\tfrac{\Delta\omega}{2}t\Big)\,
\cos\!\Big(\bar k\,x-\bar\omega\,t\Big),
```

with `` \Delta k=k_2-k_1 ``, `` \Delta\omega=\omega_2-\omega_1 ``, and `` \bar k,\bar\omega ``
the averages. The second factor is a fast carrier moving at `` \bar\omega/\bar k ``. The
first factor is a slow envelope, the beat pattern, moving at `` \Delta\omega/\Delta k ``.
Letting the two waves approach each other gives the two velocities:

```math
v_p=\frac{\omega}{k},
\qquad
v_g=\frac{d\omega}{dk}.
```

Only the **dispersion relation** `` \omega(k) `` is needed to compute both. If
`` \omega=ck `` with constant `` c ``, then `` v_p=v_g=c `` and nothing disperses. If the
curve bends, the two slopes differ.
"""

# ╔═╡ d3234123-6b2d-47fe-bd11-b71671f2f080
md"""
## A packet is a Fourier sum, and dispersion reshapes it

A packet is a superposition of plane waves with a smooth distribution of wavenumbers
around `` k_0 ``, here a Gaussian of standard deviation `` \sigma_k ``:

```math
u(x,t)=\int A(k)\,e^{\,i\,[kx-\omega(k)t]}\,dk .
```

Expand the dispersion relation about the carrier:

```math
\omega(k)\approx\omega_0+v_g\,(k-k_0)+\tfrac12\,\omega''\,(k-k_0)^2 ,
\qquad \omega''=\frac{d^2\omega}{dk^2}.
```

- The **constant** term gives the carrier, which moves at `` v_p=\omega_0/k_0 ``.
- The **linear** term shifts the envelope rigidly at `` v_g ``, with no change of shape.
- The **quadratic** term, the *group-velocity dispersion*, is what changes the shape. For a
  Gaussian packet the envelope width at time `` t `` is

```math
\sigma_x(t)=\sigma_{x0}\sqrt{1+\left(\frac{\omega''\,t}{\sigma_{x0}^{2}}\right)^{2}},
\qquad \sigma_{x0}=\frac{1}{\sigma_k}.
```

A wide band of wavenumbers (small `` \sigma_{x0} ``) spreads fastest, and no spreading
occurs at all if `` \omega''=0 ``. This is the Fourier time-bandwidth trade-off at work.
The widget computes `` \omega'' `` and the predicted widening for the chosen medium.

**Normal and anomalous dispersion.** If `` v_g<v_p `` (phase velocity decreasing with
wavenumber, so long waves outrun short ones) the dispersion is called *normal*; if
`` v_g>v_p `` it is *anomalous*. Deep-water waves and Love waves are normal; the cutoff
medium has `` v_g<c<v_p ``, with `` v_p v_g=c^2 ``.
"""

# ╔═╡ 3305138e-f288-4376-932b-50ce5817204c
md"""
## Dispersion in seismology

Surface waves are guided by the layering of the Earth, and their phase velocity depends on
period because long waves sample deeper, faster rock. The **Love wave in a layer over a
half-space** is the simplest model: it is the medium called *Love wave in a layer* in the
widget (layer thickness `` h=1 ``, layer shear speed `` \beta_1=1 ``, half-space shear speed
`` \beta_2=1.6 ``). At small wavenumber the wave feels mostly the fast half-space
(`` v_p\to\beta_2 ``); at large wavenumber it is trapped in the slow layer
(`` v_p\to\beta_1 ``). Its group velocity is smaller than both, with a minimum, the famous
*Airy phase*, where energy piles up.

Measuring group velocity versus period from a seismogram, and inverting it for Earth
structure, is the basis of surface-wave tomography. The Love and Rayleigh dispersion
notebooks in this collection compute these curves for realistic models.
"""

# ╔═╡ e980ddf7-6806-4c76-9b58-b8476043ef58
md"""
## Appendix
How the experiment computes what it shows. The dispersion relations, the two velocities and the whole space-time wavefield are computed in Julia; the browser only draws them and steps through time.
"""

# ╔═╡ 38d9bc70-436f-4b9b-b733-97f6502c8a0f
md"""
### Dispersion relations
`dispersion_relation(medium)` returns ``\omega(k)`` for each medium: ``\omega = k`` (non-dispersive), ``\omega = \sqrt{k}`` (deep water, in units with ``g = 1``), ``\omega = \sqrt{1 + k^2}`` (a cutoff medium with ``c = 1``), and the fundamental Love mode. `love_omega` finds the Love-mode phase velocity ``c`` between ``\beta_1`` and ``\beta_2`` by bisection on the layer-over-half-space dispersion equation, then returns ``\omega = ck``. `velocities(ω, k0)` gives ``v_p = \omega_0/k_0``, ``v_g`` and ``\omega''`` by central finite differences of ``\omega(k)``.
"""

# ╔═╡ aedc58c8-bd48-4968-bc93-495846772259
begin
    """
        love_omega(k; h=1.0, β1=1.0, β2=1.6, μ_ratio=3.072)

    Angular frequency of the fundamental Love mode of a layer of thickness `h` and
    shear speed `β1` over a half-space of shear speed `β2`, at wavenumber `k`.
    `μ_ratio` is the rigidity ratio `μ2/μ1` (here density ratio 1.2 times
    `(β2/β1)^2`). The phase velocity `c` solves
    `k c h q1 = atan(μ_ratio q2 / q1)` with `q1 = √(1/β1² - 1/c²)` and
    `q2 = √(1/c² - 1/β2²)`, which has exactly one root in `(β1, β2)` for the
    fundamental mode; it is found by bisection.
    """
    function love_omega(k::Float64; h::Float64=1.0, β1::Float64=1.0, β2::Float64=1.6,
        μ_ratio::Float64=3.072)
        function mismatch(c)
            q1 = sqrt(1 / β1^2 - 1 / c^2)
            q2 = sqrt(1 / c^2 - 1 / β2^2)
            return k * c * h * q1 - atan(μ_ratio * q2 / q1)
        end
        lo = β1 * (1 + 1e-12)
        hi = β2 * (1 - 1e-12)
        for _ in 1:80
            mid = (lo + hi) / 2
            if mismatch(mid) < 0
                lo = mid
            else
                hi = mid
            end
        end
        return k * (lo + hi) / 2
    end

    """
        dispersion_relation(medium)

    Return the function `ω(k)` of the named medium, in natural units where the
    characteristic speed is 1: `"nondispersive"` (`ω = k`), `"water"` (deep-water
    gravity waves, `ω = √k`), `"cutoff"` (`ω = √(1 + k²)`, the Klein–Gordon relation
    with a low-frequency cutoff) and `"love"` (fundamental Love mode in a layer).
    """
    function dispersion_relation(medium::String)
        if medium == "nondispersive"
            k -> Float64(k)
        elseif medium == "water"
            k -> sqrt(Float64(k))
        elseif medium == "cutoff"
            k -> sqrt(1 + Float64(k)^2)
        else
            k -> love_omega(Float64(k))
        end
    end

    """
        medium_label(medium)

    One-line description of the medium and its dispersion relation, for the widget.
    """
    function medium_label(medium::String)
        if medium == "nondispersive"
            "Non-dispersive medium: ω = c k  (c = 1)"
        elseif medium == "water"
            "Deep-water gravity waves: ω = √(g k)  (g = 1)"
        elseif medium == "cutoff"
            "Cutoff medium (Klein–Gordon): ω = √(ω_c² + c² k²)  (ω_c = c = 1)"
        else
            "Fundamental Love mode: layer h = 1, β₁ = 1 over half-space β₂ = 1.6"
        end
    end

    """
        velocities(ω, k0)

    Phase velocity `ω(k0)/k0`, group velocity `dω/dk` and group-velocity dispersion
    `d²ω/dk²` at `k0`, from centred finite differences of the dispersion relation.
    """
    function velocities(ω, k0::Float64)
        δ = 1e-3
        w0 = ω(k0)
        vg = (ω(k0 + δ) - ω(k0 - δ)) / (2δ)
        gvd = (ω(k0 + 10δ) - 2w0 + ω(k0 - 10δ)) / (10δ)^2
        return (omega0=w0, vp=w0 / k0, vg=vg, gvd=gvd)
    end
end

# ╔═╡ 7abad04c-51be-4ed1-84ff-0d16cd987fca
md"""
### The packet as a Fourier sum
`dispersive_packet_field` evaluates the packet integral as a sum of plane waves with Gaussian weights ``A(k)`` around ``k_0``, each advanced with its own ``\omega(k)``. It returns the wave and its exact envelope ``|u|`` on a space-time grid. The amplitude is not renormalised, so spreading shows up as a drop in peak height. `envelope_rms_width` measures the root-mean-square width of an envelope.
"""

# ╔═╡ b67615bc-454a-41e6-9286-e7df1f19823b
begin
    """
        dispersive_packet_field(ω, k0, σk, xs, ts; ncomp=61)

    Gaussian wave packet in a medium with dispersion relation `ω(k)`: the sum over
    `ncomp` plane waves with wavenumbers within `3σk` of `k0`. Returns the real
    wavefield and the exact envelope `|A|` on the grid `xs × ts`, as matrices
    indexed `[ix, it]`. The envelope comes from the complex amplitude relative to the
    carrier, so it is not an approximation. The packet starts at `x = 0` with unit
    peak amplitude, and the amplitude is *not* renormalised afterwards, so the decay
    caused by spreading is preserved.
    """
    function dispersive_packet_field(ω, k0::Float64, σk::Float64,
        xs::Vector{Float64}, ts::Vector{Float64}; ncomp::Int=61)
        ks = collect(range(max(k0 - 3σk, 1e-3), k0 + 3σk, length=ncomp))
        ωk = [ω(k) for k in ks]
        ω0 = ω(k0)
        w = @. exp(-0.5 * ((ks - k0) / σk)^2)
        w ./= sum(w)
        packet = Matrix{Float64}(undef, length(xs), length(ts))
        envelope = similar(packet)
        for (it, t) in enumerate(ts), (ix, x) in enumerate(xs)
            re = 0.0
            im = 0.0
            for j in eachindex(ks)
                θ = (ks[j] - k0) * x - (ωk[j] - ω0) * t
                re += w[j] * cos(θ)
                im += w[j] * sin(θ)
            end
            A = complex(re, im)
            packet[ix, it] = real(A * cis(k0 * x - ω0 * t))
            envelope[ix, it] = abs(A)
        end
        return (packet=packet, envelope=envelope)
    end

    """
        envelope_rms_width(xs, envelope)

    Root-mean-square width of an envelope about its centroid, weighting by
    `envelope²`. Only ratios of these widths are meaningful.
    """
    function envelope_rms_width(xs::Vector{Float64}, envelope::AbstractVector{Float64})
        e2 = envelope .^ 2
        xbar = sum(xs .* e2) / sum(e2)
        return sqrt(sum((xs .- xbar) .^ 2 .* e2) / sum(e2))
    end
end

# ╔═╡ 8ffde99b-5fb0-4595-9821-9245ee5a685a
md"""
### From parameters to the widget
`wave_dispersion_payload(medium, k0, dk)` evaluates the dispersion curve, the velocities and the predicted widening, then precomputes the packet on a fixed grid (`WD_X0` to `WD_X1` in space, 0 to `WD_TMAX` in time) and packs everything into one JSON string. The widget receives it once per change of medium, wavenumber or bandwidth. Scrubbing or playing the time, and moving the receiver, read from this precomputed field without returning to Julia.
"""

# ╔═╡ d36f7477-249c-4c79-9794-1432053ce04b
begin
    const WD_X0 = -15.0
    const WD_X1 = 85.0
    const WD_NX = 321
    const WD_TMAX = 60.0
    const WD_NT = 241
    const WD_KMAX = 5.0

    json_number(x; digits::Int=5) = isfinite(x) ? string(round(Float64(x), digits=digits)) : "0"
    json_array(v; digits::Int=5) = "[" * join(json_number.(v; digits=digits), ",") * "]"
    json_string(s) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"") * "\""

    """
        wave_dispersion_payload(medium, k0, dk)

    Everything the browser needs, as a JSON string: the dispersion curve, the two
    velocities and the group-velocity dispersion at `k0`, and a precomputed space-time
    wavefield (`WD_NT` frames on `WD_NX` points) so that scrubbing and playing the
    time need no further round trip to Julia. `dk` is the full width at half maximum
    of the wavenumber distribution.
    """
    function wave_dispersion_payload(medium::String, k0::Float64, dk::Float64)
        ω = dispersion_relation(medium)
        σk = dk / 2.355
        xs = collect(range(WD_X0, WD_X1, length=WD_NX))
        ts = collect(range(0.0, WD_TMAX, length=WD_NT))
        field = dispersive_packet_field(ω, k0, σk, xs, ts)
        kc = collect(range(0.0, WD_KMAX, length=121))
        wc = [ω(k) for k in kc]
        v = velocities(ω, k0)
        kind = abs(v.vp - v.vg) < 1e-3 * v.vp ? "non-dispersive: v_g = v_p" :
               v.vg < v.vp ? "normal dispersion: v_g < v_p (long waves are faster)" :
               "anomalous dispersion: v_g > v_p (short waves are faster)"
        widen = sqrt(1 + (v.gvd * WD_TMAX * σk^2)^2)
        info = [
            "v_p = ω₀/k₀ = " * string(round(v.vp, digits=3)) * "     v_g = dω/dk = " * string(round(v.vg, digits=3)) * "     ω'' = d²ω/dk² = " * string(round(v.gvd, digits=3)),
            kind,
            "Predicted envelope widening by t = " * string(round(WD_TMAX, digits=0)) * ":  × " * string(round(widen, digits=2)),
        ]
        ie = 1:4:WD_NX
        return string(
            "{\"x0\":", json_number(WD_X0), ",\"x1\":", json_number(WD_X1), ",\"nx\":", string(WD_NX),
            ",\"nt\":", string(WD_NT), ",\"tmax\":", json_number(WD_TMAX), ",\"nxe\":", string(length(ie)),
            ",\"k0\":", json_number(k0), ",\"w0\":", json_number(v.omega0), ",\"dk\":", json_number(dk),
            ",\"vp\":", json_number(v.vp), ",\"vg\":", json_number(v.vg), ",\"kmax\":", json_number(WD_KMAX),
            ",\"wmax\":", json_number(maximum(wc)), ",\"label\":", json_string(medium_label(medium)),
            ",\"curve_k\":", json_array(kc), ",\"curve_w\":", json_array(wc),
            ",\"packet\":", json_array(vec(field.packet); digits=3),
            ",\"env\":", json_array(vec(field.envelope[ie, :]); digits=3),
            ",\"info\":[", join(json_string.(info), ","), "]}",
        )
    end
end

# ╔═╡ 8c12d133-39a7-4469-80d9-254bb3d6f399
md"""
### Verifying the implementation
The checks below compare the computed velocities with their analytic values: ``v_g = v_p/2`` for deep water, ``v_p v_g = 1`` for the cutoff medium, and ``v_p = v_g`` without dispersion. They check that the Love-mode phase velocity approaches ``\beta_2`` at small ``k`` and ``\beta_1`` at large ``k``. Finally, they check the envelope widening of a packet against the Gaussian formula above, using an exactly quadratic ``\omega(k)`` so that the formula holds without approximation.
"""

# ╔═╡ 9f3810e4-d5cd-4366-9069-716a32703113
let
    fd(ω, k) = (ω(k + 1e-4) - ω(k - 1e-4)) / 2e-4

    # 1. Velocities against their analytic values.
    k = 1.7
    water = dispersion_relation("water")
    cutoff = dispersion_relation("cutoff")
    plain = dispersion_relation("nondispersive")
    @assert isapprox(velocities(water, k).vg, 0.5 * velocities(water, k).vp; rtol=1e-4)
    @assert isapprox(velocities(cutoff, k).vp * velocities(cutoff, k).vg, 1.0; rtol=1e-4)
    @assert isapprox(velocities(plain, k).vp, velocities(plain, k).vg; rtol=1e-6)

    # 2. The Love phase velocity stays between the two shear speeds and approaches them
    #    at the two ends of the spectrum.
    love = dispersion_relation("love")
    c_low = love(0.05) / 0.05
    c_high = love(30.0) / 30.0
    @assert 1.0 < c_low < 1.6 && 1.0 < c_high < 1.6
    @assert abs(c_low - 1.6) < 0.05 && abs(c_high - 1.0) < 0.05

    # 3. Envelope spreading against the analytic Gaussian-packet formula, using an exactly
    #    quadratic dispersion relation so the formula holds without approximation.
    k0, σk, ω2, vg0 = 2.0, 0.3, 0.3, 0.5
    quad = k -> 2.0 + vg0 * (k - k0) + 0.5 * ω2 * (k - k0)^2
    xs = collect(range(-60.0, 200.0, length=2001))
    # Many components keep the periodic images of the packet far outside this wide
    # window; the rms width is very sensitive to anything in the far tails.
    f0 = dispersive_packet_field(quad, k0, σk, xs, [0.0]; ncomp=201)
    f40 = dispersive_packet_field(quad, k0, σk, xs, [40.0]; ncomp=201)
    ratio = envelope_rms_width(xs, f40.envelope[:, 1]) / envelope_rms_width(xs, f0.envelope[:, 1])
    predicted = sqrt(1 + (ω2 * 40.0 * σk^2)^2)
    @assert abs(ratio - predicted) / predicted < 0.02

    (water_vg_over_vp = round(velocities(water, k).vg / velocities(water, k).vp, digits=4),
        cutoff_vp_times_vg = round(velocities(cutoff, k).vp * velocities(cutoff, k).vg, digits=4),
        love_c_low_k = round(c_low, digits=3), love_c_high_k = round(c_high, digits=3),
        envelope_widening_measured = round(ratio, digits=3),
        envelope_widening_predicted = round(predicted, digits=3))
end

# ╔═╡ dc0c3d8e-9f21-416d-bf2b-96cfa728410c
md"""
### The Interactive Widget
`WaveDispersionInput` sends three values back to Julia: the medium, the carrier wavenumber ``k_0`` and the bandwidth. The cell after the widget coerces them into concrete types, with defaults when no browser is attached. The cell after that calls `wave_dispersion_payload` and pushes the result with a `wd-results` event. The time slider, the Play button and the receiver position are local to the browser.
"""

# ╔═╡ d0867100-ead8-4ba7-81e4-7e2204bdbc55
begin
    struct WaveDispersionInput
        medium::String
        center_k::Float64
        bandwidth::Float64
    end

    WaveDispersionInput(; medium="water", center_k=1.0, bandwidth=0.8) =
        WaveDispersionInput(String(medium), Float64(center_k), Float64(bandwidth))

    Base.get(w::WaveDispersionInput) = Dict{String,Any}(
        "medium" => w.medium,
        "center_k" => w.center_k,
        "bandwidth" => w.bandwidth,
    )

    function Base.show(io::IO, ::MIME"text/html", w::WaveDispersionInput)
        write(io, """
        <div id="wd-widget" class="app-widget">
          <style>
            #wd-widget{box-sizing:border-box;width:100%;max-width:1400px;margin:0 auto;color:#e5e7eb;font:14px system-ui,sans-serif}
            #wd-widget *{box-sizing:border-box}
            #wd-widget [hidden]{display:none!important}
            #wd-widget .wd-title{margin-bottom:10px;padding:10px 14px;text-align:center;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px}
            #wd-widget .wd-title-desc{font-size:17px;font-weight:700}
            #wd-widget .wd-title-hint{margin-top:3px;color:#9ca3af;font-size:13px}
            #wd-widget .wd-workspace{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}
            #wd-widget .wd-panel{min-width:0;padding:10px;background:#050505;border:1px solid #2f3744;border-radius:6px}
            #wd-widget .wd-wide{grid-column:1 / -1}
            #wd-widget .wd-panel-title{margin:0 0 5px;font-size:15px;font-weight:700;color:#f3f4f6}
            #wd-widget .wd-caption{min-height:34px;color:#9ca3af;font-size:13px;line-height:1.3}
            #wd-widget canvas{display:block;width:100%;height:310px;background:#000;border:1px solid #374151;border-radius:4px}
            #wd-widget canvas.wd-short{height:230px}
            #wd-widget .wd-legend{display:flex;gap:13px;flex-wrap:wrap;margin-top:7px;color:#d1d5db;font-size:12px}
            #wd-widget .wd-key{display:inline-flex;gap:5px;align-items:center}
            #wd-widget .wd-swatch{width:18px;height:3px;border-radius:2px;background:currentColor}
            #wd-widget .wd-dash{background:repeating-linear-gradient(90deg,currentColor 0 6px,transparent 6px 10px)}
            #wd-widget .wd-info{display:flex;flex-direction:column;gap:3px;margin-top:10px;padding:9px 12px;min-height:80px;background:#050505;border:1px solid #2f3744;border-radius:6px;color:#e5e7eb;font-size:14px}
            #wd-widget .wd-controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:8px;margin-top:10px}
            #wd-widget .wd-control-group{min-width:0;padding:10px;background:#050505;border:1px solid #2f3744;border-radius:6px}
            #wd-widget .wd-control-title{margin-bottom:7px;font-size:16px;font-weight:700;color:#f3f4f6}
            #wd-widget .wd-control-row{display:grid;grid-template-columns:minmax(74px,126px) minmax(70px,1fr) minmax(42px,62px);gap:7px;align-items:center;margin:7px 0;color:#d1d5db}
            #wd-widget input,#wd-widget select{min-width:0;width:100%;accent-color:#f59e0b}
            #wd-widget select{padding:4px 6px;color:#f3f4f6;background:#111827;border:1px solid #4b5563;border-radius:4px}
            #wd-widget .wd-value{color:#fbbf24;text-align:right;font-variant-numeric:tabular-nums}
            #wd-widget .wd-actions{display:flex;gap:8px;flex-wrap:wrap;margin-top:8px}
            #wd-widget button{border:1px solid #9ca3af;border-radius:4px;padding:6px 12px;background:#606060;color:#f3f4f6;font-size:14px;cursor:pointer}
            #wd-widget button:hover{background:#737373}
            @media(max-width:760px){#wd-widget .wd-workspace{grid-template-columns:1fr}#wd-widget canvas{height:260px}}
            /* App shell Expt view: the three panels fill the stage and their canvases its height. */
            html[data-sa-mode="expt"] #wd-widget .wd-workspace{grid-template-rows:minmax(0,1.35fr) minmax(0,1fr);align-items:stretch!important;justify-content:stretch!important;align-content:stretch}
            html[data-sa-mode="expt"] #wd-widget .wd-panel{display:flex;flex-direction:column;min-height:0}
            html[data-sa-mode="expt"] #wd-widget canvas{flex:1 1 auto;height:auto;min-height:0}
            @media(max-width:820px){html[data-sa-mode="expt"] #wd-widget .wd-workspace{grid-template-rows:none}html[data-sa-mode="expt"] #wd-widget canvas{flex:none;height:240px}}
          </style>
          <div class="wd-title app-title">
            <div class="wd-title-desc">Watch a wave packet and the crests inside it travel at different speeds.</div>
            <div class="wd-title-hint">Choose a medium and a carrier wavenumber &middot; play the time &middot; compare the packet with the seismogram</div>
          </div>
          <div class="wd-workspace app-stage">
            <section class="wd-panel">
              <div class="wd-panel-title">Dispersion relation &omega;(k)</div>
              <div class="wd-caption" id="wd-disp-cap"></div>
              <canvas id="wd-disp"></canvas>
              <div class="wd-legend"><span class="wd-key" style="color:#e5e7eb"><i class="wd-swatch"></i>&omega;(k)</span><span class="wd-key" style="color:#60a5fa"><i class="wd-swatch wd-dash"></i>secant: slope = v<sub>p</sub></span><span class="wd-key" style="color:#ef4444"><i class="wd-swatch"></i>tangent: slope = v<sub>g</sub></span></div>
            </section>
            <section class="wd-panel">
              <div class="wd-panel-title">Wave packet in space</div>
              <div class="wd-caption">The red envelope moves at v<sub>g</sub>; the crests inside it move at v<sub>p</sub>. The triangle marks the receiver.</div>
              <canvas id="wd-packet"></canvas>
              <div class="wd-legend"><span class="wd-key" style="color:#60a5fa"><i class="wd-swatch"></i>wave</span><span class="wd-key" style="color:#ef4444"><i class="wd-swatch wd-dash"></i>envelope</span><span class="wd-key" style="color:#ef4444"><i class="wd-swatch"></i>x = v<sub>g</sub>t</span><span class="wd-key" style="color:#60a5fa"><i class="wd-swatch"></i>x = v<sub>p</sub>t</span></div>
            </section>
            <section class="wd-panel wd-wide">
              <div class="wd-panel-title" id="wd-seis-title">Seismogram at the receiver</div>
              <div class="wd-caption">Wave recorded at the receiver as a function of time. The red dashed line is the group arrival time x<sub>r</sub>/v<sub>g</sub>; the yellow line is the current time.</div>
              <canvas id="wd-seis" class="wd-short"></canvas>
            </section>
          </div>
          <div class="wd-info" id="wd-info"></div>
          <div class="wd-controls app-controls">
            <section class="wd-control-group">
              <div class="wd-control-title">Medium</div>
              <label class="wd-control-row"><span>Medium</span><select id="wd-medium"><option value="nondispersive" $(w.medium == "nondispersive" ? "selected" : "")>Non-dispersive</option><option value="water" $(w.medium == "water" ? "selected" : "")>Deep-water waves</option><option value="cutoff" $(w.medium == "cutoff" ? "selected" : "")>Cutoff medium</option><option value="love" $(w.medium == "love" ? "selected" : "")>Love wave in a layer</option></select><span></span></label>
              <label class="wd-control-row"><span>Carrier k&#8320;</span><input id="wd-k0" type="range" min="0.5" max="3" step="0.1" value="$(w.center_k)"><span class="wd-value" id="wd-k0-v"></span></label>
              <label class="wd-control-row"><span>Bandwidth</span><input id="wd-dk" type="range" min="0.2" max="1" step="0.05" value="$(w.bandwidth)"><span class="wd-value" id="wd-dk-v"></span></label>
            </section>
            <section class="wd-control-group">
              <div class="wd-control-title">Time and receiver</div>
              <label class="wd-control-row"><span>Time</span><input id="wd-time" type="range" min="0" max="60" step="0.25" value="0"><span class="wd-value" id="wd-time-v"></span></label>
              <label class="wd-control-row"><span>Receiver x</span><input id="wd-rx" type="range" min="5" max="80" step="1" value="20"><span class="wd-value" id="wd-rx-v"></span></label>
              <div class="wd-actions"><button id="wd-play" type="button">Play</button><button id="wd-reset" type="button">Reset experiment</button></div>
            </section>
          </div>
        </div>
        <script>
        const par = currentScript.previousElementSibling;
        const byId = id => par.querySelector('#' + id);
        const state = {medium:"$(w.medium)", center_k:$(w.center_k), bandwidth:$(w.bandwidth)};
        const DEFAULTS = {medium:'water', center_k:1.0, bandwidth:0.8};
        const KEYS = {'wd-medium':'medium', 'wd-k0':'center_k', 'wd-dk':'bandwidth'};
        const LOCAL = ['wd-time', 'wd-rx'];
        let result = null, tcur = 0, rx = 20, playing = false, lastNow = 0;
        const M = {l:50, r:14, t:14, b:42};
        function fmt(v){ return String(parseFloat(v.toFixed(3))); }
        function niceTicks(lo, hi, n){
          const raw = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag;
          const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
          const out = [];
          for(let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9; v += step) out.push(parseFloat(v.toFixed(10)));
          return out;
        }
        function frame(canvas, xmin, xmax, ymin, ymax, xlabel, ylabel){
          const rect = canvas.getBoundingClientRect(), ratio = window.devicePixelRatio || 1;
          canvas.width = Math.max(1, Math.round(rect.width * ratio));
          canvas.height = Math.max(1, Math.round(rect.height * ratio));
          const ctx = canvas.getContext('2d');
          ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
          const W = rect.width, H = rect.height;
          ctx.fillStyle = '#000'; ctx.fillRect(0, 0, W, H);
          const px = x => M.l + (x - xmin) / (xmax - xmin) * (W - M.l - M.r);
          const py = y => H - M.b - (y - ymin) / (ymax - ymin) * (H - M.t - M.b);
          ctx.font = '12px sans-serif';
          ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
          for(const v of niceTicks(ymin, ymax, 5)){
            const y = py(v);
            ctx.strokeStyle = v === 0 ? '#4b5563' : '#1f2937'; ctx.lineWidth = 1;
            ctx.beginPath(); ctx.moveTo(M.l, y); ctx.lineTo(W - M.r, y); ctx.stroke();
            ctx.fillStyle = '#9ca3af'; ctx.fillText(fmt(v), M.l - 6, y);
          }
          ctx.textAlign = 'center'; ctx.textBaseline = 'top';
          for(const v of niceTicks(xmin, xmax, 8)){
            const x = px(v);
            ctx.strokeStyle = '#1f2937';
            ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, H - M.b); ctx.stroke();
            ctx.fillStyle = '#9ca3af'; ctx.fillText(fmt(v), x, H - M.b + 5);
          }
          ctx.fillStyle = '#d1d5db'; ctx.font = '13px sans-serif'; ctx.textBaseline = 'alphabetic';
          ctx.fillText(xlabel, M.l + (W - M.l - M.r) / 2, H - 6);
          ctx.save(); ctx.translate(13, M.t + (H - M.t - M.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(ylabel, 0, 0); ctx.restore();
          ctx.textAlign = 'left';
          return {ctx:ctx, px:px, py:py, W:W, H:H};
        }
        function clipPlot(f){ f.ctx.beginPath(); f.ctx.rect(M.l, M.t, f.W - M.l - M.r, f.H - M.t - M.b); f.ctx.clip(); }
        function line(f, x0, y0, x1, y1, color, width, dash){
          const ctx = f.ctx;
          ctx.save(); clipPlot(f);
          ctx.beginPath(); ctx.moveTo(f.px(x0), f.py(y0)); ctx.lineTo(f.px(x1), f.py(y1));
          ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash); ctx.stroke();
          ctx.restore();
        }
        function polyline(f, xs, ys, color, width, dash){
          const ctx = f.ctx;
          ctx.save(); clipPlot(f);
          ctx.beginPath();
          for(let i = 0; i < xs.length; i++){ const x = f.px(xs[i]), y = f.py(ys[i]); if(i) ctx.lineTo(x, y); else ctx.moveTo(x, y); }
          ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash); ctx.stroke();
          ctx.restore();
        }
        function text(f, s, x, y, color){ f.ctx.fillStyle = color; f.ctx.font = '12px sans-serif'; f.ctx.textAlign = 'left'; f.ctx.textBaseline = 'alphabetic'; f.ctx.fillText(s, x, y); }
        function frameIndex(){ return Math.max(0, Math.min(result.nt - 1, Math.round(tcur / result.tmax * (result.nt - 1)))); }
        function drawDispersion(){
          const R = result;
          const f = frame(byId('wd-disp'), 0, R.kmax, 0, R.wmax * 1.05, 'wavenumber k', 'angular frequency \\u03C9');
          const ctx = f.ctx;
          ctx.save(); clipPlot(f);
          ctx.fillStyle = 'rgba(245,158,11,0.18)';
          const kl = Math.max(0, R.k0 - R.dk / 2), kr = R.k0 + R.dk / 2;
          ctx.fillRect(f.px(kl), M.t, f.px(kr) - f.px(kl), f.H - M.t - M.b);
          ctx.restore();
          polyline(f, R.curve_k, R.curve_w, '#e5e7eb', 2.2, []);
          line(f, 0, 0, R.kmax, R.vp * R.kmax, '#60a5fa', 1.8, [6, 4]);
          const span = 1.6;
          line(f, R.k0 - span, R.w0 - R.vg * span, R.k0 + span, R.w0 + R.vg * span, '#ef4444', 2.0, []);
          ctx.fillStyle = '#f59e0b'; ctx.beginPath(); ctx.arc(f.px(R.k0), f.py(R.w0), 5, 0, 2 * Math.PI); ctx.fill();
          text(f, 'v_p = ' + R.vp.toFixed(2), f.px(R.kmax) - 78, f.py(R.vp * R.kmax) + 16, '#60a5fa');
          text(f, 'v_g = ' + R.vg.toFixed(2), f.px(R.k0) + 10, f.py(R.w0) + 18, '#ef4444');
        }
        function drawPacket(){
          const R = result, j = frameIndex(), nx = R.nx, nxe = R.nxe;
          const f = frame(byId('wd-packet'), R.x0, R.x1, -1.15, 1.15, 'position x', 'amplitude');
          const xs = new Array(nx), ys = new Array(nx), xe = new Array(nxe), ye = new Array(nxe), yn = new Array(nxe);
          for(let i = 0; i < nx; i++){ xs[i] = R.x0 + i * (R.x1 - R.x0) / (nx - 1); ys[i] = R.packet[j * nx + i]; }
          for(let i = 0; i < nxe; i++){ xe[i] = R.x0 + (4 * i) * (R.x1 - R.x0) / (nx - 1); ye[i] = R.env[j * nxe + i]; yn[i] = -ye[i]; }
          polyline(f, xe, ye, '#ef4444', 1.5, [6, 4]);
          polyline(f, xe, yn, '#ef4444', 1.5, [6, 4]);
          polyline(f, xs, ys, '#60a5fa', 2.1, []);
          const t = tcur;
          line(f, R.vg * t, -1.15, R.vg * t, 1.15, '#ef4444', 2, [5, 4]);
          line(f, R.vp * t, -1.15, R.vp * t, 1.15, '#60a5fa', 2, [5, 4]);
          if(R.vg * t <= R.x1) text(f, 'v_g t', f.px(R.vg * t) + 4, M.t + 12, '#ef4444');
          if(R.vp * t <= R.x1) text(f, 'v_p t', f.px(R.vp * t) + 4, M.t + 26, '#60a5fa');
          text(f, 't = ' + t.toFixed(1), M.l + 8, M.t + 14, '#e5e7eb');
          const ctx = f.ctx, tx = f.px(rx), ty = M.t + 4;
          ctx.fillStyle = '#facc15'; ctx.strokeStyle = '#000'; ctx.lineWidth = 1;
          ctx.beginPath(); ctx.moveTo(tx - 6, ty); ctx.lineTo(tx + 6, ty); ctx.lineTo(tx, ty + 11); ctx.closePath(); ctx.fill(); ctx.stroke();
        }
        function drawSeismogram(){
          const R = result, nx = R.nx;
          const ix = Math.max(0, Math.min(nx - 1, Math.round((rx - R.x0) / (R.x1 - R.x0) * (nx - 1))));
          const ts = new Array(R.nt), us = new Array(R.nt);
          for(let j = 0; j < R.nt; j++){ ts[j] = j * R.tmax / (R.nt - 1); us[j] = R.packet[j * nx + ix]; }
          const f = frame(byId('wd-seis'), 0, R.tmax, -1.15, 1.15, 'time t', 'amplitude');
          polyline(f, ts, us, '#60a5fa', 1.8, []);
          const arrival = rx / R.vg;
          if(arrival <= R.tmax){ line(f, arrival, -1.15, arrival, 1.15, '#ef4444', 2, [5, 4]); text(f, 'x_r / v_g', f.px(arrival) + 4, M.t + 12, '#ef4444'); }
          line(f, tcur, -1.15, tcur, 1.15, '#facc15', 2, []);
        }
        function draw(){
          if(!result) return;
          drawDispersion(); drawPacket(); drawSeismogram();
        }
        function updateLabels(){
          byId('wd-k0-v').textContent = state.center_k.toFixed(1);
          byId('wd-dk-v').textContent = state.bandwidth.toFixed(2);
          byId('wd-time-v').textContent = tcur.toFixed(1);
          byId('wd-rx-v').textContent = rx.toFixed(0);
          byId('wd-seis-title').textContent = 'Seismogram at the receiver, x = ' + rx.toFixed(0);
        }
        function showInfo(){
          const box = byId('wd-info'); box.textContent = '';
          if(!result) return;
          byId('wd-disp-cap').textContent = result.label;
          for(const s of result.info){ const d = document.createElement('div'); d.textContent = s; box.appendChild(d); }
        }
        function publish(){ par.value = Object.assign({}, state); par.dispatchEvent(new CustomEvent('input')); }
        function onControl(event){
          const id = event.target.id;
          if(event.target === par) return;
          if(LOCAL.indexOf(id) >= 0){
            event.stopImmediatePropagation();
            if(id === 'wd-time') tcur = Number(event.target.value);
            else rx = Number(event.target.value);
            updateLabels(); draw();
            return;
          }
          if(!(id in KEYS)) return;
          event.stopImmediatePropagation();
          if(event.type === 'change' && event.target.type === 'range') return;
          state[KEYS[id]] = event.target.type === 'range' ? Number(event.target.value) : event.target.value;
          updateLabels(); publish();
        }
        par.addEventListener('input', onControl, true);
        par.addEventListener('change', onControl, true);
        function tick(now){
          if(!playing) return;
          if(result){
            tcur += (now - lastNow) / 1000 * result.tmax / 10;
            if(tcur >= result.tmax) tcur = 0;
            byId('wd-time').value = tcur; updateLabels(); draw();
          }
          lastNow = now; requestAnimationFrame(tick);
        }
        byId('wd-play').addEventListener('click', () => {
          playing = !playing; byId('wd-play').textContent = playing ? 'Pause' : 'Play';
          if(playing){ lastNow = performance.now(); requestAnimationFrame(tick); }
        });
        byId('wd-reset').addEventListener('click', () => {
          Object.assign(state, DEFAULTS); tcur = 0; rx = 20;
          byId('wd-medium').value = state.medium; byId('wd-k0').value = state.center_k; byId('wd-dk').value = state.bandwidth;
          byId('wd-time').value = 0; byId('wd-rx').value = 20;
          updateLabels(); draw(); publish();
        });
        function onResults(event){ result = JSON.parse(event.detail); showInfo(); draw(); }
        window.addEventListener('wd-results', onResults);
        // Redraw whenever the widget changes size (app-shell tab switches included). Clean up
        // only once the widget has left the page: Pluto can resolve `invalidation` for a bond
        // output while keeping its DOM.
        const ro = new ResizeObserver(() => {
          if(!par.isConnected){ playing = false; window.removeEventListener('wd-results', onResults); ro.disconnect(); return; }
          draw();
        });
        ro.observe(par);
        updateLabels();
        if(window.__wdLast){ result = JSON.parse(window.__wdLast); showInfo(); draw(); }
        </script>
        """)
    end

    const _wd_ready = true
end

# ╔═╡ d376eb0c-798c-49bd-a9ce-ea75f7614231
begin
    # `WaveDispersionInput` is defined in the Appendix. The bare reference creates the
    # dependency Pluto needs on a cold notebook load.
    _wd_ready
    PlutoUI.WideCell(@bind _wd WaveDispersionInput(); max_width=1400)
end

# ╔═╡ dda8eb81-7f06-46ed-9aef-75617aa7d080
begin
    # Convert the untyped bond value once, at the widget boundary. Every downstream
    # calculation receives concrete scalar values.
    wd_safe = _wd isa AbstractDict ? _wd : Dict{String,Any}()
    wd_medium = String(get(wd_safe, "medium", "water"))
    wd_medium in ("nondispersive", "water", "cutoff", "love") || (wd_medium = "water")
    wd_center_k = clamp(Float64(get(wd_safe, "center_k", 1.0)), 0.5, 3.0)
    wd_bandwidth = clamp(Float64(get(wd_safe, "bandwidth", 0.8)), 0.2, 1.0)
    nothing
end

# ╔═╡ 5d3cc1bb-abf9-48e9-891c-2265d5b3ad21
let
    payload = wave_dispersion_payload(wd_medium, wd_center_k, wd_bandwidth)
    HTML("""<script>
      window.__wdLast = $(repr(payload));
      window.dispatchEvent(new CustomEvent('wd-results', {detail: window.__wdLast}));
    </script>""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
FFTW = "~1.10.0"
HypertextLiteral = "~0.9.5"
PlutoPlotly = "~0.6.5"
PlutoUI = "~0.7.72"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.1"
manifest_format = "2.1"
project_hash = "12cd0b4046b366044235ef63a21b66829bc04db5"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
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
registries = "General"
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
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
registries = "General"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
registries = "General"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
registries = "General"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
registries = "General"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
registries = "General"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.1+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
registries = "General"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
registries = "General"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "28278bb0053da0fd73537be94afd1682cc5a0a83"
registries = "General"
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
registries = "General"
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
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "c3b2323466378a2ba15bea4b2f73b081e022f473"
registries = "General"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
registries = "General"
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
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
registries = "General"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
registries = "General"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.2+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╟─7a1fa0a3-b612-40da-8f9d-f1494877560c
# ╟─bef462ee-988d-49b5-acc8-599d4c401067
# ╟─7e344ce5-8f4f-4a52-bc21-dea0934d900d
# ╟─d376eb0c-798c-49bd-a9ce-ea75f7614231
# ╟─dda8eb81-7f06-46ed-9aef-75617aa7d080
# ╟─5d3cc1bb-abf9-48e9-891c-2265d5b3ad21
# ╟─aab4ff38-2775-4cff-b45f-b189eb16fbf7
# ╟─aad6728c-ee0f-4063-9f8a-00a1f86d94f4
# ╟─d3234123-6b2d-47fe-bd11-b71671f2f080
# ╟─3305138e-f288-4376-932b-50ce5817204c
# ╟─e980ddf7-6806-4c76-9b58-b8476043ef58
# ╟─38d9bc70-436f-4b9b-b733-97f6502c8a0f
# ╠═aedc58c8-bd48-4968-bc93-495846772259
# ╟─7abad04c-51be-4ed1-84ff-0d16cd987fca
# ╠═b67615bc-454a-41e6-9286-e7df1f19823b
# ╟─8ffde99b-5fb0-4595-9821-9245ee5a685a
# ╠═d36f7477-249c-4c79-9794-1432053ce04b
# ╟─8c12d133-39a7-4469-80d9-254bb3d6f399
# ╠═9f3810e4-d5cd-4366-9069-716a32703113
# ╟─dc0c3d8e-9f21-416d-bf2b-96cfa728410c
# ╠═d0867100-ead8-4ba7-81e4-7e2204bdbc55
# ╠═59057ef4-796a-4f06-b1f6-5a74ab535b1d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
