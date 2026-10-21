### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> tags = ["introduction"]
#> title = "Fourier Series and the Fourier Transform"
#> description = "Build periodic waveforms from harmonics, then see how the duration of a pulse and its spectrum trade off."
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

# ╔═╡ 49067b3a-8eec-4ad0-a118-8ad350e6cb5a
using PlutoUI, FFTW

# ╔═╡ 6b768613-5b89-404a-968c-e29b81ec589f
TableOfContents()

# ╔═╡ 7be64d4b-3f1b-422c-97a5-50ec45393ce3
# App shell: the Expt / Theory / Code / Help frame shared by every notebook, loaded
# from src/assets/app-shell/. Without that folder (a lone downloaded .jl) this cell just
# shows the help text and the notebook stays a plain notebook.
let
	help = md"""
	- **Mode** switches between a Fourier *series* (a periodic signal built from harmonics) and a Fourier *transform* (a pulse and its spectrum).
	- **Series:** choose the waveform and drag **Terms** to add harmonics; tick *show individual terms* to see the first few.
	- **Pulse:** choose the shape, then drag **Duration** (and **Centre f** for the Gabor pulse) and watch the spectrum widen or narrow.
	- The readout under the plots gives the overshoot (series) or the time and frequency widths and their product (pulse).

	**Try this:** halve the duration of the Gaussian pulse and check that its spectral width doubles while Δt·Δf stays at 0.44.
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

# ╔═╡ b51929b4-3c08-4dd0-98f6-6f90f004f2be
md"""
# Fourier Series and the Fourier Transform

Every waveform can be described in two complementary ways: by how its value changes
with time (or position), and by *which frequencies it contains and how strongly*. The
**Fourier series** does this for a periodic signal, and the **Fourier transform** does
it for a signal of finite duration: a pulse, a wavelet, a seismic arrival. This notebook
lets you move back and forth between the two descriptions and see that they are the
same signal.

Start with the experiment. In **series** mode, add harmonics one at a time and watch a
square, sawtooth or triangle wave emerge, together with a bar chart of the harmonics you
have included. In **pulse** mode, change the shape and duration of a pulse and watch its
spectrum respond: squeeze the pulse in time and its spectrum spreads out.

A companion notebook, *Wave Dispersion*, uses this idea: a wave packet is a sum of
frequencies, and dispersion makes each one travel at its own speed.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ fd1f8079-b58e-46e2-b076-84fe24857701
md"""
## How to read the experiment

**Series mode.** The dashed curve is the target periodic function; the orange curve is
the sum of the first `N` sine harmonics. Switch on *show individual harmonics* to see the
first few terms that make up the sum. The bar chart on the right is the signal's
*spectrum*: one bar per harmonic, orange when it is included in the sum. Compare the three
waveforms: the bars of the square and sawtooth waves shrink slowly, the bars of the
triangle wave shrink fast, and the reconstruction follows the same pattern.

**Pulse mode.** The left panel is the pulse `y(t)`; the right panel is its amplitude
spectrum `|F(f)|`, computed with an FFT. Halve the duration and the spectrum doubles in
width. The readout under the panels reports the pulse duration `Δt`, the spectral width
`Δf`, and their product.
"""

# ╔═╡ 8720cda1-13b7-4130-a1fc-29eadc87fab8
md"""
## Fourier series: a periodic signal as a sum of harmonics

A signal that repeats every `` 2\pi `` can be written as a sum of sines and cosines whose
frequencies are whole multiples of the repetition frequency:

```math
f(x) = \frac{a_0}{2} + \sum_{n=1}^{\infty}\left[a_n\cos(nx) + b_n\sin(nx)\right],
\qquad
b_n = \frac{1}{\pi}\int_{-\pi}^{\pi} f(x)\,\sin(nx)\,dx .
```

All three target waveforms are odd about `` x=0 ``, so only the sine terms survive:

| waveform | `` b_n `` |
|:--|:--|
| square | `` 4/(\pi n) `` for odd `` n ``, zero for even `` n `` |
| sawtooth | `` 2(-1)^{n+1}/(\pi n) `` |
| triangle | `` 8(-1)^{(n-1)/2}/(\pi^2 n^2) `` for odd `` n ``, zero for even `` n `` |

**Smoothness sets the decay.** The square and sawtooth waves have jumps, and their
coefficients fall off as `` 1/n ``. The triangle wave has only a corner, and its
coefficients fall off as `` 1/n^2 ``. The smoother the signal, the faster its spectrum
dies away, so fewer harmonics are needed.

**The Gibbs overshoot does not go away.** Near a jump, the partial sum overshoots the
target. Adding more terms squeezes the overshoot toward the jump but does not shrink it:
the peak stays about `` 8.95\% `` of the jump height. The widget reports the overshoot of
the current partial sum. It is a property of the series, not a plotting artifact.
"""

# ╔═╡ 2459a91e-5d54-4af6-b73d-12ff8cb7c08a
md"""
## Fourier transform: from lines to a continuum

A pulse is not periodic, so it cannot be built from harmonics of one repetition frequency.
It needs *every* frequency, with a continuous amplitude distribution:

```math
F(f) = \int_{-\infty}^{\infty} y(t)\,e^{-2\pi i f t}\,dt,
\qquad
y(t) = \int_{-\infty}^{\infty} F(f)\,e^{2\pi i f t}\,df .
```

You can picture it as the Fourier series of a signal whose period `` T `` grows without
limit: the spacing `` 1/T `` between harmonics shrinks to zero and the bars merge into a
curve. Four pulses are available in the widget:

- **Gaussian** `` y=e^{-t^2/2s^2} ``. Its transform is also a Gaussian,
  `` |F| = s\sqrt{2\pi}\,e^{-2\pi^2 s^2 f^2} ``: a Gaussian is its own transform up to scale.
- **Ricker wavelet** `` y=(1-2\pi^2 f_0^2 t^2)\,e^{-\pi^2 f_0^2 t^2} ``, the standard
  zero-phase source wavelet in seismic modelling. It has no zero-frequency content and a
  spectrum peaked at `` f_0 ``.
- **Boxcar** (a signal switched on and off). Its sharp edges give a `` \mathrm{sinc} ``
  spectrum with slowly decaying side lobes, the same Gibbs-type behaviour as the series.
- **Gabor** (a Gaussian-windowed cosine), the building block of a wave packet. Its spectrum
  is a Gaussian centred on the carrier frequency `` f_0 ``.

**Time and frequency widths trade off.** A short pulse needs a wide band of frequencies,
and a narrow band gives a long pulse: scaling a pulse's duration by a factor scales its
bandwidth by the inverse factor, so their product stays fixed. For a Gaussian, with both
widths measured as full widths at half maximum, that product is exactly

```math
\Delta t\,\Delta f = \frac{2\ln 2}{\pi}\approx 0.44 .
```

Other shapes give other values, which the widget reports. The Ricker wavelet's is
smaller because its main lobe is narrow while part of its energy sits in side lobes, so a
half-maximum width does not capture all of it; the strict bound uses root-mean-square
widths. Either way the lesson for seismology is the same: a source with a broad spectrum
resolves thin layers, while a narrow-band source cannot separate two close arrivals.
"""

# ╔═╡ a9b7872a-6a4e-422d-b498-241cf667ca50
md"""
## Appendix
How the experiment computes what it shows. Every number in the widget comes from the Julia functions below: the series coefficients from their closed forms, the pulse spectra from an FFT. The browser only draws.
"""

# ╔═╡ 68e6eeb1-1573-4332-92e8-12fc390fd25c
md"""
### Signals and their spectra
`fourier_coefficient(waveform, n)` returns the sine coefficient ``b_n`` from the table above, and `series_target(waveform, x)` the exact periodic function it approximates. `fourier_series_panel(waveform, nterms)` sums the first `nterms` harmonics on a grid of ``x``, keeps the first few individual terms for display, lists every ``|b_n|`` for the bar spectrum, and measures the overshoot near the jump.

`pulse_signal(shape, duration, f0, t)` evaluates the four pulses. `pulse_spectrum` samples the chosen pulse finely over a long window, takes a real FFT (`rfft`, so the spectrum is computed rather than drawn from a formula), and normalises it to its peak. `half_max_width` measures a full width at half maximum, used for both ``\Delta t`` (on the pulse energy) and ``\Delta f`` (on the power spectrum).
"""

# ╔═╡ 793fa999-d11c-407b-ab91-6ae434f90d3c
begin
    """
        fourier_coefficient(waveform, n)

    Sine-series coefficient `b_n` of one period of the target waveform on
    `(-π, π)`. The three targets are odd functions, so their cosine coefficients
    vanish and these textbook sine coefficients describe them completely.
    """
    function fourier_coefficient(waveform::String, n::Int)
        if waveform == "square"
            isodd(n) ? 4 / (π * n) : 0.0
        elseif waveform == "sawtooth"
            2 * (-1)^(n + 1) / (π * n)
        else
            isodd(n) ? 8 * (-1)^((n - 1) ÷ 2) / (π^2 * n^2) : 0.0
        end
    end

    """
        series_target(waveform, x)

    Value of the target periodic function at `x` (one period on `(-π, π)`).
    """
    function series_target(waveform::String, x::Float64)
        if waveform == "square"
            sign(sin(x))
        elseif waveform == "sawtooth"
            x / π
        else
            (2 / π) * asin(sin(x))
        end
    end

    """
        fourier_series_panel(waveform, nterms; show_harmonics=false, points=801, nmax=40)

    Partial Fourier sum with `nterms` harmonics, the target waveform, up to six of the
    individual terms (when `show_harmonics`), the first `nmax` coefficients for the
    spectrum bars, and two diagnostics: the peak overshoot as a fraction of the jump
    height (zero for the continuous triangle wave) and the RMS misfit.
    """
    function fourier_series_panel(waveform::String, nterms::Int;
        show_harmonics::Bool=false, points::Int=801, nmax::Int=40)
        x = collect(range(-π, π, length=points))
        target = [series_target(waveform, xi) for xi in x]
        approx = zeros(Float64, points)
        harmonics = Vector{Vector{Float64}}()
        for n in 1:nterms
            b = fourier_coefficient(waveform, n)
            term = @. b * sin(n * x)
            approx .+= term
            if show_harmonics && b != 0.0 && length(harmonics) < 6
                push!(harmonics, term)
            end
        end
        spectrum = [fourier_coefficient(waveform, n) for n in 1:nmax]
        jump = waveform == "triangle" ? 0.0 : 2.0
        overshoot = jump == 0.0 ? 0.0 : max((maximum(approx) - 1.0) / jump, 0.0)
        rms = sqrt(sum(abs2, target .- approx) / points)
        ymax = max(1.15, 1.08 * maximum(abs, approx))
        return (x=x, target=target, approx=approx, harmonics=harmonics, spectrum=spectrum,
            overshoot=overshoot, rms=rms, ymax=ymax)
    end
end

# ╔═╡ 485fae00-47f0-4090-9407-71566c06e543
begin
    """
        pulse_signal(shape, duration, f0, t)

    Sample a unit-peak pulse on the time axis `t`. `shape` is `"gaussian"` (standard
    deviation `duration`), `"ricker"` (peak frequency `f0`), `"boxcar"` (half-width
    `duration`) or `"gabor"` (Gaussian window of standard deviation `duration` times a
    cosine at `f0`). Returns the signal and its amplitude envelope.
    """
    function pulse_signal(shape::String, duration::Float64, f0::Float64, t::Vector{Float64})
        if shape == "gaussian"
            y = @. exp(-0.5 * (t / duration)^2)
            (signal=y, envelope=copy(y))
        elseif shape == "ricker"
            a = (π * f0)^2
            y = @. (1 - 2 * a * t^2) * exp(-a * t^2)
            (signal=y, envelope=abs.(y))
        elseif shape == "boxcar"
            y = @. ifelse(abs(t) <= duration, 1.0, 0.0)
            (signal=y, envelope=copy(y))
        else
            env = @. exp(-0.5 * (t / duration)^2)
            y = @. env * cos(2π * f0 * t)
            (signal=y, envelope=env)
        end
    end

    """
        half_max_width(x, y)

    Interpolated left and right edges of the contiguous region around the global peak
    of `y` where `y` is at least half its maximum.
    """
    function half_max_width(x::AbstractVector{Float64}, y::AbstractVector{Float64})
        i0 = argmax(y)
        half = y[i0] / 2
        lo = i0
        while lo > firstindex(y) && y[lo-1] >= half
            lo -= 1
        end
        hi = i0
        while hi < lastindex(y) && y[hi+1] >= half
            hi += 1
        end
        xl = lo > firstindex(y) ? x[lo-1] + (half - y[lo-1]) / (y[lo] - y[lo-1]) * (x[lo] - x[lo-1]) : x[lo]
        xh = hi < lastindex(y) ? x[hi] + (half - y[hi]) / (y[hi+1] - y[hi]) * (x[hi+1] - x[hi]) : x[hi]
        return xl, xh
    end

    """
        pulse_spectrum(shape, duration, f0)

    Amplitude spectrum `|F(f)|` of a pulse, computed with an FFT on a long window
    (80 s, 8192 samples) so the frequency spacing is fine enough to resolve a 2 s
    pulse. Returns time and frequency axes, the signal, the spectrum, and the two
    widths `Δt`, `Δf` measured as full widths at half maximum of the energy and the
    power spectrum.
    """
    function pulse_spectrum(shape::String, duration::Float64, f0::Float64)
        N = 8192
        dt = 80.0 / N
        t = collect((-N ÷ 2):(N ÷ 2 - 1)) .* dt
        sig = pulse_signal(shape, duration, f0, t)
        f = collect(rfftfreq(N, 1 / dt))
        amp = abs.(rfft(sig.signal)) .* dt
        tl, th = half_max_width(t, sig.envelope .^ 2)
        pl, ph = half_max_width(f, amp .^ 2)
        Δt = th - tl
        Δf = argmax(amp) == 1 ? 2 * ph : ph - pl
        return (t=t, f=f, signal=sig.signal, envelope=sig.envelope, amplitude=amp, Δt=Δt, Δf=Δf)
    end
end

# ╔═╡ 5be197bd-bc1b-4d8b-b00f-9ca828a51d0c
begin
    json_number(x) = isfinite(x) ? string(round(Float64(x), digits=5)) : "0"
    json_array(v) = "[" * join(json_number.(v), ",") * "]"
    json_string(s) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"") * "\""

    """
        fourier_pair_payload(mode, waveform, nterms, show_harmonics, shape, duration, f0)

    Everything the browser needs to draw both panels, as a JSON string. The widget
    only draws these arrays: the series, the FFT and the width measurements are all
    evaluated in Julia.
    """
    function fourier_pair_payload(mode::String, waveform::String, nterms::Int,
        show_harmonics::Bool, shape::String, duration::Float64, f0::Float64)
        if mode == "pulse"
            p = pulse_spectrum(shape, duration, f0)
            it = findall(x -> abs(x) <= 4.0, p.t)
            iff = findall(x -> x <= 6.0, p.f)
            curves = String[]
            shape == "gabor" && push!(curves,
                "{\"role\":\"envelope\",\"y\":" * json_array(p.envelope[it]) * "}")
            push!(curves, "{\"role\":\"pulse\",\"y\":" * json_array(p.signal[it]) * "}")
            spec = p.amplitude[iff] ./ maximum(p.amplitude)
            product = p.Δt * p.Δf
            info = [
                "Duration Δt = " * string(round(p.Δt, digits=2)) * " s (full width at half maximum of the energy)",
                "Bandwidth Δf = " * string(round(p.Δf, digits=2)) * " Hz (full width at half maximum of the power spectrum)",
                "Time–bandwidth product Δt·Δf = " * string(round(product, digits=2)) * "  (a Gaussian gives 0.44)",
            ]
            return string(
                "{\"mode\":\"pulse\",\"left\":{\"xlabel\":\"time t (s)\",\"ylabel\":\"y(t)\",\"ymin\":-1.15,\"ymax\":1.15,\"x\":",
                json_array(p.t[it]), ",\"curves\":[", join(curves, ","), "]},",
                "\"right\":{\"kind\":\"line\",\"xlabel\":\"frequency f (Hz)\",\"ylabel\":\"|F(f)| (normalised)\",\"ymin\":0,\"ymax\":1.1,\"x\":",
                json_array(p.f[iff]), ",\"y\":", json_array(spec), "},",
                "\"info\":[", join(json_string.(info), ","), "]}",
            )
        else
            s = fourier_series_panel(waveform, nterms; show_harmonics=show_harmonics)
            ticks = [(-π, "-π"), (-π / 2, "-π/2"), (0.0, "0"), (π / 2, "π/2"), (π, "π")]
            xticks = "[" * join(["[" * json_number(v) * "," * json_string(l) * "]" for (v, l) in ticks], ",") * "]"
            curves = String["{\"role\":\"target\",\"y\":" * json_array(s.target) * "}"]
            for h in s.harmonics
                push!(curves, "{\"role\":\"harmonic\",\"y\":" * json_array(h) * "}")
            end
            push!(curves, "{\"role\":\"sum\",\"y\":" * json_array(s.approx) * "}")
            bars = abs.(s.spectrum)
            info = if waveform == "triangle"
                ["Continuous waveform: no jump, so no Gibbs overshoot.",
                    "RMS misfit with " * string(nterms) * " terms: " * string(round(s.rms, digits=3))]
            else
                ["Peak overshoot " * string(round(100 * 2 * s.overshoot, digits=1)) * " % of the target = " *
                 string(round(100 * s.overshoot, digits=1)) * " % of the jump (Gibbs limit 8.95 %)",
                    "RMS misfit with " * string(nterms) * " terms: " * string(round(s.rms, digits=3))]
            end
            return string(
                "{\"mode\":\"series\",\"left\":{\"xlabel\":\"x (radians)\",\"ylabel\":\"f(x)\",\"ymin\":",
                json_number(-s.ymax), ",\"ymax\":", json_number(s.ymax), ",\"xticks\":", xticks,
                ",\"x\":", json_array(s.x), ",\"curves\":[", join(curves, ","), "]},",
                "\"right\":{\"kind\":\"bars\",\"xlabel\":\"harmonic number n\",\"ylabel\":\"|b_n|\",\"ymin\":0,\"ymax\":",
                json_number(1.1 * maximum(bars)), ",\"included\":", string(nterms),
                ",\"x\":", json_array(1:length(bars)), ",\"y\":", json_array(bars), "},",
                "\"info\":[", join(json_string.(info), ","), "]}",
            )
        end
    end
end

# ╔═╡ c756a536-9d22-468b-835b-dc0e41ee640c
let
    # 1. Gibbs: with many terms the square-wave overshoot approaches Si(π)/π - 1/2.
    gibbs_limit = 0.08949
    gibbs_measured = fourier_series_panel("square", 39).overshoot
    @assert abs(gibbs_measured - gibbs_limit) < 0.004

    # 2. A Gaussian pulse has an analytic transform: s*sqrt(2π)*exp(-2π²s²f²).
    s = 0.5
    g = pulse_spectrum("gaussian", s, 2.0)
    analytic = @. s * sqrt(2π) * exp(-2 * π^2 * s^2 * g.f^2)
    keep = analytic .> 1e-3 * maximum(analytic)
    fft_error = maximum(abs.(g.amplitude[keep] .- analytic[keep]) ./ analytic[keep])
    @assert fft_error < 1e-6

    # 3. Gaussian time-bandwidth product (FWHM of energy x FWHM of power) = 2 ln2 / π.
    product = g.Δt * g.Δf
    @assert abs(product - 2 * log(2) / π) / (2 * log(2) / π) < 5e-3

    (square_overshoot_39_terms = round(gibbs_measured, digits=4),
        gibbs_limit = gibbs_limit,
        gaussian_fft_max_relative_error = fft_error,
        gaussian_time_bandwidth_product = round(product, digits=4),
        theory = round(2 * log(2) / π, digits=4))
end

# ╔═╡ d3196bba-fd83-4e16-add0-380a7b285004
md"""
### The Interactive Widget
`FourierPairInput` sends the mode, waveform, number of terms, harmonics toggle, pulse shape, duration and centre frequency back to Julia. The cell after the widget coerces these into concrete types (`fp_mode`, `fp_harmonics`, …, with defaults when no browser is attached), and the cell after that calls `fourier_pair_payload` and pushes the result. Changing any control therefore re-runs only those two cells.
"""

# ╔═╡ 70c4f4fc-cddb-4309-a418-d482b47e4eed
begin
    struct FourierPairInput
        mode::String
        waveform::String
        harmonics::Int
        show_harmonics::Bool
        shape::String
        duration::Float64
        center_frequency::Float64
    end

    FourierPairInput(; mode="series", waveform="square", harmonics=7, show_harmonics=true,
        shape="gaussian", duration=0.5, center_frequency=2.0) = FourierPairInput(
        String(mode), String(waveform), Int(harmonics), Bool(show_harmonics),
        String(shape), Float64(duration), Float64(center_frequency),
    )

    Base.get(w::FourierPairInput) = Dict{String,Any}(
        "mode" => w.mode,
        "waveform" => w.waveform,
        "harmonics" => w.harmonics,
        "show_harmonics" => w.show_harmonics,
        "shape" => w.shape,
        "duration" => w.duration,
        "center_frequency" => w.center_frequency,
    )

    function Base.show(io::IO, ::MIME"text/html", w::FourierPairInput)
        write(io, """
        <div id="fp-widget" class="app-widget">
          <style>
            #fp-widget{box-sizing:border-box;width:100%;max-width:1400px;margin:0 auto;color:#e5e7eb;font:14px system-ui,sans-serif}
            #fp-widget *{box-sizing:border-box}
            #fp-widget [hidden]{display:none!important}
            #fp-widget .fp-title{margin-bottom:10px;padding:10px 14px;text-align:center;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px}
            #fp-widget .fp-title-desc{font-size:17px;font-weight:700}
            #fp-widget .fp-title-hint{margin-top:3px;color:#9ca3af;font-size:13px}
            #fp-widget .fp-workspace{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}
            #fp-widget .fp-panel{min-width:0;padding:10px;background:#050505;border:1px solid #2f3744;border-radius:6px}
            #fp-widget .fp-panel-title{margin:0 0 5px;font-size:15px;font-weight:700;color:#f3f4f6}
            #fp-widget .fp-caption{min-height:34px;color:#9ca3af;font-size:13px;line-height:1.3}
            #fp-widget canvas{display:block;width:100%;height:320px;background:#000;border:1px solid #374151;border-radius:4px}
            #fp-widget .fp-legend{display:flex;gap:13px;flex-wrap:wrap;margin-top:7px;color:#d1d5db;font-size:12px}
            #fp-widget .fp-key{display:inline-flex;gap:5px;align-items:center}
            #fp-widget .fp-swatch{width:18px;height:3px;border-radius:2px;background:currentColor}
            #fp-widget .fp-dash{background:repeating-linear-gradient(90deg,currentColor 0 6px,transparent 6px 10px)}
            #fp-widget .fp-info{display:flex;flex-direction:column;gap:3px;margin-top:10px;padding:9px 12px;min-height:60px;background:#050505;border:1px solid #2f3744;border-radius:6px;color:#e5e7eb;font-size:14px}
            #fp-widget .fp-controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:8px;margin-top:10px}
            #fp-widget .fp-control-group{min-width:0;padding:10px;background:#050505;border:1px solid #2f3744;border-radius:6px}
            #fp-widget .fp-control-title{margin-bottom:7px;font-size:16px;font-weight:700;color:#f3f4f6}
            #fp-widget .fp-control-row{display:grid;grid-template-columns:minmax(74px,126px) minmax(70px,1fr) minmax(42px,62px);gap:7px;align-items:center;margin:7px 0;color:#d1d5db}
            #fp-widget .fp-control-row.fp-off{opacity:.4}
            #fp-widget input,#fp-widget select{min-width:0;width:100%;accent-color:#f59e0b}
            #fp-widget input[type=checkbox]{width:auto}
            #fp-widget select{padding:4px 6px;color:#f3f4f6;background:#111827;border:1px solid #4b5563;border-radius:4px}
            #fp-widget .fp-value{color:#fbbf24;text-align:right;font-variant-numeric:tabular-nums}
            #fp-widget button{border:1px solid #9ca3af;border-radius:4px;padding:6px 12px;background:#606060;color:#f3f4f6;font-size:14px;cursor:pointer}
            #fp-widget button:hover{background:#737373}
            @media(max-width:760px){#fp-widget .fp-workspace{grid-template-columns:1fr}#fp-widget canvas{height:260px}}
            /* App shell Expt view: the two panels fill the stage and their canvases its height. */
            html[data-sa-mode="expt"] #fp-widget .fp-workspace{grid-template-rows:minmax(0,1fr);align-items:stretch!important;justify-content:stretch!important;align-content:stretch}
            html[data-sa-mode="expt"] #fp-widget .fp-panel{display:flex;flex-direction:column;min-height:0}
            html[data-sa-mode="expt"] #fp-widget canvas{flex:1 1 auto;height:auto;min-height:0}
            @media(max-width:820px){html[data-sa-mode="expt"] #fp-widget .fp-workspace{grid-template-rows:none}html[data-sa-mode="expt"] #fp-widget canvas{flex:none;height:260px}}
          </style>
          <div class="fp-title app-title">
            <div class="fp-title-desc">Build a signal from frequencies, and see how its duration and its spectrum trade off.</div>
            <div class="fp-title-hint">Choose series or pulse mode &middot; change the waveform, the number of terms, or the pulse shape and duration</div>
          </div>
          <div class="fp-workspace app-stage">
            <section class="fp-panel">
              <div class="fp-panel-title" id="fp-left-title"></div>
              <div class="fp-caption" id="fp-left-cap"></div>
              <canvas id="fp-left"></canvas>
              <div class="fp-legend" id="fp-left-legend"></div>
            </section>
            <section class="fp-panel">
              <div class="fp-panel-title" id="fp-right-title"></div>
              <div class="fp-caption" id="fp-right-cap"></div>
              <canvas id="fp-right"></canvas>
              <div class="fp-legend" id="fp-right-legend"></div>
            </section>
          </div>
          <div class="fp-info" id="fp-info"></div>
          <div class="fp-controls app-controls">
            <section class="fp-control-group">
              <div class="fp-control-title">Experiment</div>
              <label class="fp-control-row"><span>Mode</span><select id="fp-mode"><option value="series" $(w.mode == "series" ? "selected" : "")>Series: periodic signal</option><option value="pulse" $(w.mode == "pulse" ? "selected" : "")>Transform: pulse</option></select><span></span></label>
              <div><button id="fp-reset" type="button">Reset experiment</button></div>
            </section>
            <section class="fp-control-group" id="fp-group-series">
              <div class="fp-control-title">Fourier series</div>
              <label class="fp-control-row"><span>Waveform</span><select id="fp-waveform"><option value="square" $(w.waveform == "square" ? "selected" : "")>Square</option><option value="sawtooth" $(w.waveform == "sawtooth" ? "selected" : "")>Sawtooth</option><option value="triangle" $(w.waveform == "triangle" ? "selected" : "")>Triangle</option></select><span></span></label>
              <label class="fp-control-row"><span>Terms</span><input id="fp-harmonics" type="range" min="1" max="40" step="1" value="$(w.harmonics)"><span class="fp-value" id="fp-harmonics-v"></span></label>
              <label class="fp-control-row"><span>Harmonics</span><span><input id="fp-show-h" type="checkbox" $(w.show_harmonics ? "checked" : "")> show individual terms</span><span></span></label>
            </section>
            <section class="fp-control-group" id="fp-group-pulse">
              <div class="fp-control-title">Fourier transform</div>
              <label class="fp-control-row"><span>Pulse shape</span><select id="fp-shape"><option value="gaussian" $(w.shape == "gaussian" ? "selected" : "")>Gaussian</option><option value="ricker" $(w.shape == "ricker" ? "selected" : "")>Ricker wavelet</option><option value="boxcar" $(w.shape == "boxcar" ? "selected" : "")>Boxcar</option><option value="gabor" $(w.shape == "gabor" ? "selected" : "")>Gabor (windowed cosine)</option></select><span></span></label>
              <label class="fp-control-row" id="fp-row-duration"><span>Duration</span><input id="fp-duration" type="range" min="0.1" max="2" step="0.05" value="$(w.duration)"><span class="fp-value" id="fp-duration-v"></span></label>
              <label class="fp-control-row" id="fp-row-fc"><span>Centre f</span><input id="fp-fc" type="range" min="0.5" max="5" step="0.1" value="$(w.center_frequency)"><span class="fp-value" id="fp-fc-v"></span></label>
            </section>
          </div>
        </div>
        <script>
        const par = currentScript.previousElementSibling;
        const byId = id => par.querySelector('#' + id);
        const state = {mode:"$(w.mode)", waveform:"$(w.waveform)", harmonics:$(w.harmonics), show_harmonics:$(w.show_harmonics), shape:"$(w.shape)", duration:$(w.duration), center_frequency:$(w.center_frequency)};
        const DEFAULTS = {mode:'series', waveform:'square', harmonics:7, show_harmonics:true, shape:'gaussian', duration:0.5, center_frequency:2.0};
        const KEYS = {'fp-mode':'mode', 'fp-waveform':'waveform', 'fp-harmonics':'harmonics', 'fp-show-h':'show_harmonics', 'fp-shape':'shape', 'fp-duration':'duration', 'fp-fc':'center_frequency'};
        const ROLE = {target:{color:'#9ca3af', width:1.6, dash:[6,4]}, sum:{color:'#f59e0b', width:2.4, dash:[]}, harmonic:{color:'#38bdf8', width:1.2, dash:[]}, pulse:{color:'#60a5fa', width:2.4, dash:[]}, envelope:{color:'#ef4444', width:1.6, dash:[6,4]}};
        const HCOL = ['#38bdf8', '#34d399', '#a78bfa', '#f472b6', '#fbbf24', '#fb923c'];
        const TEXT = {
          series:{lt:'Signal and its Fourier reconstruction', lc:'Dashed: target periodic function. Orange: sum of the first N sine harmonics.', rt:'Spectrum: harmonic amplitudes', rc:'One bar per harmonic. Orange bars are included in the sum, grey bars are not.',
                  ll:[['#9ca3af','target',true],['#f59e0b','partial sum',false],['#38bdf8','individual terms',false]], rl:[['#f59e0b','included',false],['#4b5563','not included',false]]},
          pulse:{lt:'Pulse in time', lc:'The signal y(t). For the Gabor pulse the dashed red curve is its Gaussian envelope.', rt:'Amplitude spectrum', rc:'Computed with an FFT in Julia and normalised to its peak.',
                 ll:[['#60a5fa','pulse y(t)',false],['#ef4444','envelope (Gabor only)',true]], rl:[['#34d399','|F(f)|',false]]}
        };
        let result = null;
        const M = {l:50, r:14, t:14, b:42};
        function fmt(v){ return String(parseFloat(v.toFixed(3))); }
        function niceTicks(lo, hi, n){
          const raw = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag;
          const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
          const out = [];
          for(let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9; v += step) out.push(parseFloat(v.toFixed(10)));
          return out;
        }
        function frame(canvas, xmin, xmax, ymin, ymax, xticks, xlabel, ylabel){
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
          for(const t of xticks){
            const x = px(t[0]);
            ctx.strokeStyle = '#1f2937';
            ctx.beginPath(); ctx.moveTo(x, M.t); ctx.lineTo(x, H - M.b); ctx.stroke();
            ctx.fillStyle = '#9ca3af'; ctx.fillText(t[1], x, H - M.b + 5);
          }
          ctx.fillStyle = '#d1d5db'; ctx.font = '13px sans-serif'; ctx.textBaseline = 'alphabetic';
          ctx.fillText(xlabel, M.l + (W - M.l - M.r) / 2, H - 6);
          ctx.save(); ctx.translate(13, M.t + (H - M.t - M.b) / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(ylabel, 0, 0); ctx.restore();
          ctx.textAlign = 'left';
          return {ctx:ctx, px:px, py:py, W:W, H:H};
        }
        function clipPlot(f){ f.ctx.beginPath(); f.ctx.rect(M.l, M.t, f.W - M.l - M.r, f.H - M.t - M.b); f.ctx.clip(); }
        function polyline(f, xs, ys, style){
          const ctx = f.ctx;
          ctx.save(); clipPlot(f);
          ctx.beginPath();
          for(let i = 0; i < xs.length; i++){ const x = f.px(xs[i]), y = f.py(ys[i]); if(i) ctx.lineTo(x, y); else ctx.moveTo(x, y); }
          ctx.strokeStyle = style.color; ctx.lineWidth = style.width; ctx.setLineDash(style.dash); ctx.stroke();
          ctx.restore();
        }
        function autoTicks(lo, hi){ return niceTicks(lo, hi, 8).map(v => [v, fmt(v)]); }
        function draw(){
          if(!result) return;
          const L = result.left, R = result.right;
          const lx = L.x, rx = R.x;
          let f = frame(byId('fp-left'), lx[0], lx[lx.length - 1], L.ymin, L.ymax, L.xticks || autoTicks(lx[0], lx[lx.length - 1]), L.xlabel, L.ylabel);
          let hi = 0;
          for(const c of L.curves){
            let style = ROLE[c.role];
            if(c.role === 'harmonic'){ style = {color:HCOL[hi % HCOL.length], width:style.width, dash:style.dash}; hi++; }
            polyline(f, lx, c.y, style);
          }
          if(R.kind === 'bars'){
            f = frame(byId('fp-right'), 0.5, rx.length + 0.5, R.ymin, R.ymax, autoTicks(1, rx.length), R.xlabel, R.ylabel);
            const bw = Math.max(2, (f.px(2) - f.px(1)) * 0.7);
            for(let i = 0; i < rx.length; i++){
              f.ctx.fillStyle = (i + 1) <= R.included ? '#f59e0b' : '#4b5563';
              const x = f.px(rx[i]), y = f.py(R.y[i]);
              f.ctx.fillRect(x - bw / 2, y, bw, f.py(0) - y);
            }
          } else {
            f = frame(byId('fp-right'), rx[0], rx[rx.length - 1], R.ymin, R.ymax, autoTicks(rx[0], rx[rx.length - 1]), R.xlabel, R.ylabel);
            const ctx = f.ctx;
            ctx.save(); clipPlot(f);
            ctx.beginPath(); ctx.moveTo(f.px(rx[0]), f.py(0));
            for(let i = 0; i < rx.length; i++) ctx.lineTo(f.px(rx[i]), f.py(R.y[i]));
            ctx.lineTo(f.px(rx[rx.length - 1]), f.py(0)); ctx.closePath();
            ctx.fillStyle = 'rgba(52,211,153,0.18)'; ctx.fill();
            ctx.restore();
            polyline(f, rx, R.y, {color:'#34d399', width:2.2, dash:[]});
          }
        }
        function legendHtml(items){
          return items.map(it => '<span class="fp-key" style="color:' + it[0] + '"><i class="fp-swatch' + (it[2] ? ' fp-dash' : '') + '"></i>' + it[1] + '</span>').join('');
        }
        function applyMode(){
          const T = TEXT[state.mode];
          byId('fp-left-title').textContent = T.lt; byId('fp-left-cap').textContent = T.lc;
          byId('fp-right-title').textContent = T.rt; byId('fp-right-cap').textContent = T.rc;
          byId('fp-left-legend').innerHTML = legendHtml(T.ll); byId('fp-right-legend').innerHTML = legendHtml(T.rl);
          byId('fp-group-series').hidden = state.mode !== 'series';
          byId('fp-group-pulse').hidden = state.mode !== 'pulse';
          const used = {gaussian:[true,false], ricker:[false,true], boxcar:[true,false], gabor:[true,true]}[state.shape];
          byId('fp-duration').disabled = !used[0]; byId('fp-row-duration').classList.toggle('fp-off', !used[0]);
          byId('fp-fc').disabled = !used[1]; byId('fp-row-fc').classList.toggle('fp-off', !used[1]);
        }
        function updateLabels(){
          byId('fp-harmonics-v').textContent = state.harmonics;
          byId('fp-duration-v').textContent = state.duration.toFixed(2) + ' s';
          byId('fp-fc-v').textContent = state.center_frequency.toFixed(1) + ' Hz';
        }
        function showInfo(){
          const box = byId('fp-info'); box.textContent = '';
          for(const s of (result ? result.info : [])){ const d = document.createElement('div'); d.textContent = s; box.appendChild(d); }
        }
        function publish(){ par.value = Object.assign({}, state); par.dispatchEvent(new CustomEvent('input')); }
        function readControl(el){
          const key = KEYS[el.id];
          if(!key) return false;
          if(key === 'show_harmonics') state[key] = el.checked;
          else if(el.type === 'range') state[key] = Number(el.value);
          else state[key] = el.value;
          return true;
        }
        function onControl(event){
          if(event.target === par || !(event.target.id in KEYS)) return;
          event.stopImmediatePropagation();
          if(event.type === 'change' && event.target.type === 'range') return;
          readControl(event.target); applyMode(); updateLabels(); publish();
        }
        par.addEventListener('input', onControl, true);
        par.addEventListener('change', onControl, true);
        byId('fp-reset').addEventListener('click', () => {
          Object.assign(state, DEFAULTS);
          for(const id in KEYS){ const el = byId(id), key = KEYS[id]; if(key === 'show_harmonics') el.checked = state[key]; else el.value = state[key]; }
          applyMode(); updateLabels(); publish();
        });
        function onResults(event){ result = JSON.parse(event.detail); showInfo(); draw(); }
        window.addEventListener('fp-results', onResults);
        // Redraw whenever the widget changes size (app-shell tab switches included). Clean up
        // only once the widget has left the page: Pluto can resolve `invalidation` for a bond
        // output while keeping its DOM.
        const ro = new ResizeObserver(() => {
          if(!par.isConnected){ window.removeEventListener('fp-results', onResults); ro.disconnect(); return; }
          draw();
        });
        ro.observe(par);
        applyMode(); updateLabels();
        if(window.__fpLast){ result = JSON.parse(window.__fpLast); showInfo(); draw(); }
        </script>
        """)
    end

    const _fp_ready = true
end

# ╔═╡ 19bbf319-4410-4247-b0bc-67dda062b9d0
begin
    # `FourierPairInput` is defined in the Appendix. The bare reference creates the
    # dependency Pluto needs on a cold notebook load.
    _fp_ready
    PlutoUI.WideCell(@bind _fp FourierPairInput(); max_width=1400)
end

# ╔═╡ bde1db20-c268-4a80-8a09-4c8e95f0c9d2
begin
    # Convert the untyped bond value once, at the widget boundary. Every downstream
    # calculation receives concrete scalar values.
    fp_safe = _fp isa AbstractDict ? _fp : Dict{String,Any}()
    fp_mode = String(get(fp_safe, "mode", "series")) == "pulse" ? "pulse" : "series"
    fp_waveform = String(get(fp_safe, "waveform", "square"))
    fp_waveform in ("square", "sawtooth", "triangle") || (fp_waveform = "square")
    fp_harmonics = clamp(round(Int, Float64(get(fp_safe, "harmonics", 7))), 1, 40)
    fp_show_harmonics = get(fp_safe, "show_harmonics", true) === true
    fp_shape = String(get(fp_safe, "shape", "gaussian"))
    fp_shape in ("gaussian", "ricker", "boxcar", "gabor") || (fp_shape = "gaussian")
    fp_duration = clamp(Float64(get(fp_safe, "duration", 0.5)), 0.1, 2.0)
    fp_center_frequency = clamp(Float64(get(fp_safe, "center_frequency", 2.0)), 0.5, 5.0)
    nothing
end

# ╔═╡ 7d3a73c0-d343-401a-a50e-4fd2b0d56a1a
let
    payload = fourier_pair_payload(fp_mode, fp_waveform, fp_harmonics, fp_show_harmonics,
        fp_shape, fp_duration, fp_center_frequency)
    HTML("""<script>
      window.__fpLast = $(repr(payload));
      window.dispatchEvent(new CustomEvent('fp-results', {detail: window.__fpLast}));
    </script>""")
end

# ╔═╡ 29f0e0bc-ba1f-4268-aeaf-805db1652d3e
md"""
### From parameters to the widget
`fourier_pair_payload` runs the series or the pulse branch for the current settings and packs the curves, the bars and the readouts into one JSON string. That string is pushed to the browser with a `fp-results` event; the widget never recomputes any of it.
"""

# ╔═╡ 27302cd1-e37c-4bbc-9064-d229c59f3e7a
md"""
### Verifying the implementation
The checks below assert three independent facts. With 39 terms the square-wave overshoot is within 0.004 of the Gibbs limit ``\mathrm{Si}(\pi)/\pi - 1/2 \approx 0.0895``. The FFT spectrum of a Gaussian matches its analytic transform ``s\sqrt{2\pi}\,e^{-2\pi^2 s^2 f^2}`` to better than one part in a million. And the measured Gaussian time-bandwidth product agrees with ``2\ln 2/\pi`` to 0.5%.
"""

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
# ╟─6b768613-5b89-404a-968c-e29b81ec589f
# ╟─7be64d4b-3f1b-422c-97a5-50ec45393ce3
# ╟─b51929b4-3c08-4dd0-98f6-6f90f004f2be
# ╟─19bbf319-4410-4247-b0bc-67dda062b9d0
# ╟─bde1db20-c268-4a80-8a09-4c8e95f0c9d2
# ╟─7d3a73c0-d343-401a-a50e-4fd2b0d56a1a
# ╟─fd1f8079-b58e-46e2-b076-84fe24857701
# ╟─8720cda1-13b7-4130-a1fc-29eadc87fab8
# ╟─2459a91e-5d54-4af6-b73d-12ff8cb7c08a
# ╟─a9b7872a-6a4e-422d-b498-241cf667ca50
# ╟─68e6eeb1-1573-4332-92e8-12fc390fd25c
# ╠═793fa999-d11c-407b-ab91-6ae434f90d3c
# ╠═485fae00-47f0-4090-9407-71566c06e543
# ╟─29f0e0bc-ba1f-4268-aeaf-805db1652d3e
# ╠═5be197bd-bc1b-4d8b-b00f-9ca828a51d0c
# ╟─27302cd1-e37c-4bbc-9064-d229c59f3e7a
# ╠═c756a536-9d22-468b-835b-dc0e41ee640c
# ╟─d3196bba-fd83-4e16-add0-380a7b285004
# ╠═70c4f4fc-cddb-4309-a418-d482b47e4eed
# ╠═49067b3a-8eec-4ad0-a118-8ad350e6cb5a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
