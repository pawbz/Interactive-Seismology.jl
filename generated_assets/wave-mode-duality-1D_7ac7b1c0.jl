### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Normal Modes 1D"
#> layout = "layout.jlhtml"
#> tags = ["normalmodes"]
#> description = "Paint a layered 1-D string, watch its normal modes one at a time or summed together, and see how superposition turns standing waves into a traveling pulse."

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

# ╔═╡ ec847573-27e7-4845-9375-e92539d7aff6
begin
    using PlutoUI
    using LinearAlgebra
    using FFTW
end

# ╔═╡ 46d9831d-2b46-4030-bb15-cdc732d72e68
TableOfContents()

# ╔═╡ 319925c0-1512-459b-8be2-1d2650678f35
md"""
# Wave-Mode Duality

A vibrating string -- or a 1-D Earth -- can be described two completely equivalent ways.
**Propagate** it: solve the wave equation forward in time from a source, and watch a
disturbance travel, reflect, and interfere. Or **decompose** it: find the string's own
natural shapes of vibration (its *normal modes*), each oscillating forever at its own fixed
frequency, and write the motion as a sum of those shapes. Neither view is more "real" than
the other -- they're the same physics, seen through two different bases. This notebook lets
you paint a layered 1-D string and watch both pictures side by side: look at **one normal
mode in isolation** (a pure standing wave), or watch a **superposition** of many modes
sharpen, as you add more of them, into a localized pulse that bounces between the string's
fixed ends.

The transverse displacement of the string can be expanded as
```math
u(x,t) = \sum_{n} c_n\,U_n(x)\cos(\omega_n t),
```
with normal modes ``U_n(x)`` and their eigenfrequencies ``\omega_n``. The source-free
elastic wave equation for a 1-D medium of density ``\rho(x)`` and shear modulus ``\mu(x)``,
written as a coupled first-order system,
```math
\partial_t\sigma - \mu\,\partial_x\dot u = 0, \qquad \rho\,\partial_t\dot u - \partial_x\sigma = 0,
```
combines into the familiar second-order wave equation
```math
\rho\,\partial_t^2 u = \partial_x(\mu\,\partial_x u).
```
Separation of variables (below) turns this into a Sturm–Liouville eigenvalue problem for
``U_n(x)`` and ``\omega_n``. Two different ways of representing displacements --
normal-mode summation, and propagating waves -- fall out of the same equation.

Seismologists use exactly this duality for real Earth structure: one can extrapolate depth
with a shooting method to match boundary conditions, or solve the same problem as an
eigenvalue problem for the Earth's own free oscillations (its normal modes, directly
observable after large earthquakes -- see below).

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)
Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 75617746-5345-4074-bbdf-7bd3e8505aae
md"""
## From the Wave Equation to Normal Modes

Look for a solution that oscillates in time at a single frequency ``\omega``,
``u(x,t) = U(x)\cos(\omega t)``. Substituting into ``\rho\,\partial_t^2 u =
\partial_x(\mu\,\partial_x u)`` and dividing out ``\cos(\omega t)`` gives an ordinary
differential equation in ``x`` alone -- a **Sturm–Liouville eigenvalue problem**:
```math
\frac{d}{dx}\!\left(\mu(x)\frac{dU}{dx}\right) = -\rho(x)\,\omega^2\,U(x).
```
Taking density uniform (``\rho=1``, so ``\mu(x)=c(x)^2``) and the string clamped
(``U=0``) at both ends turns this into a genuine eigenvalue problem: only special values
``\omega_n`` admit a nontrivial solution ``U_n(x)``, exactly like the resonant frequencies
of a guitar string.

!!! note "Why the ends being fixed matters"
	A free string (no clamping) can move as a whole with **any** frequency -- there's
	nothing to select discrete values. It's the boundary condition, not the wave equation
	alone, that turns a continuum of possible frequencies into the discrete ladder
	``\omega_1<\omega_2<\dots`` you'll see in the mode-spectrum panel below.
"""

# ╔═╡ 9a5e731a-4e70-4938-b0a9-b7fd796f4e59
md"""
### Homogeneous String: A Closed-Form Answer

If ``c(x)=c`` is constant over the whole string (length ``L``, both ends fixed), the
Sturm–Liouville equation collapses to ``U''=-(\omega/c)^2U``, whose general solution is
``U(x)=A\sin(\omega x/c)+B\cos(\omega x/c)``. The boundary conditions pin down
everything: fixing the left end forces ``B=0`` (up to where you put the origin), and
fixing the right end forces ``\sin(\omega L/c)=0``, i.e.
```math
\omega_n = \frac{n\pi c}{L}, \qquad U_n(x) = \sin\!\left(\frac{n\pi x}{L}\right), \qquad n=1,2,3,\dots
```
This is the textbook vibrating-string spectrum, and it's the case the widget's own
numerical solver is checked against below (machine precision).

### Layered String: Match Segments, Don't Force Smoothness

A real (seismological) medium isn't one homogeneous rod -- it's layers with different
wave speeds. Inside each layer the solution is still the same ``A\sin(\omega x/c_i)+
B\cos(\omega x/c_i)`` (constant ``c_i`` there), but at an interface the physics only
demands **displacement `` U `` and stress `` \mu U' ``** stay continuous -- *not* that
``U'`` itself is continuous, since ``\mu`` jumps. Chaining a `` 2\times2 `` "transfer
matrix" through each layer, from the fixed left end to the fixed right end, turns the
right-end boundary condition into a single scalar equation in ``\omega``; its roots are
the layered string's exact eigenfrequencies -- see [`find_eigenfrequencies`](@ref) in the
Appendix. This is exact for any number of layers, with no discretization at all: no grid
to refine, no kink in ``U'`` to resolve, because the piecewise-sinusoidal solution is
matched analytically across every interface.
"""

# ╔═╡ 1d0767df-945c-45cf-9f59-9ba0e3accefe
md"""
## The Two Pictures, Side by Side

**Click** a single stem in the Mode Spectrum panel to isolate that mode:
``U_n(x)\cos(\omega_n t)`` on its own -- a pure standing wave, frozen shape, oscillating
amplitude, `` n-1 `` fixed nodes never moving. This is the "resonance" picture -- the
string's own natural way of vibrating, with nothing to do with where you'd actually
strike it.

**Drag** across a range of stems instead, and the widget sums every mode in that range,
each weighted by how strongly a source at your chosen position excites it
(``c_n=U_n(x_{\text{src}})`` for a point/delta source, since the normalized modes satisfy
`` \int U_mU_n\,dx=\delta_{mn} `` -- see "Verifying the Solver"). Drag a narrow range
starting at mode 1 and the result looks like a smooth, spread-out wobble; widen it toward
mode 40 and the sum sharpens into a **localized pulse** that reflects back and forth
between the fixed ends -- the same physics as directly time-stepping the wave equation from
that source, just built entirely out of standing waves. Grow the selection and watch this
happen: that convergence *is* the duality.
"""

# ╔═╡ c17c4841-10ae-4e5c-84eb-762f7ce72af9
md"""
### Real Earth Normal Modes
"""

# ╔═╡ 235eff94-9dca-40e6-a37f-e906d04eb094
PlutoUI.LocalResource("../assets/images/normal-modes1.jpg")

# ╔═╡ 04f765cb-826e-4939-800e-9a4952e895bb
md"""[[image source]](https://www.gwrinstruments.com/low-frequency-seismology-monitoring.html)
Spectra of normal modes below `1 mHz` were measured using a Superconducting Gravimeter (SG).
Normal modes from 3 major earthquakes are compared: Blue - Japan (11 Mar. 2011, Mw=9.0),
Red - Chile (27 Feb. 2010, Mw=8.8), and Black - Sumatra (26 Dec. 2004, Mw=9.3). Data sets
are 60 hours long starting 4 hours after the onset of each event and were recorded by SG
C021 at Membach, Belgium. The same idea as the receiver-record panel in the widget above,
at planetary scale: a giant earthquake excites the whole Earth's normal-mode spectrum at
once, and a long-running gravimeter record's own Fourier spectrum shows those
eigenfrequencies directly as spectral peaks.
"""

# ╔═╡ 44fb45cd-448b-4927-856c-0bf7799bd40e
md"## Appendix"

# ╔═╡ 8de036a2-2116-4d76-975a-b87f37100a36
md"""
### Exact Layered-String Eigensolver

The real workhorse. Piecewise-constant wavespeed means the *exact* solution inside layer
`` i `` is `` A_i\sin(\omega x/c_i)+B_i\cos(\omega x/c_i) `` -- no discretization at
all. [`chain_transfer_matrix`](@ref) propagates `` (A_i,B_i) `` across each interface by
matching displacement and stress; [`find_eigenfrequencies`](@ref) turns the right-end
boundary condition into a scalar root-finding problem in `` \omega ``;
[`eigenfunction`](@ref) evaluates the resulting piecewise-sinusoidal mode shape at any `x`.
"""

# ╔═╡ e4a8ceca-2cd9-47b6-bca7-146b70c6ede4
"""
	chain_transfer_matrix(ω, velocities, widths)

Propagate the layer-local coefficients `` (A_i,B_i) `` of ``U_i(s)=A_i\\sin(k_is)+B_i
\\cos(k_is) `` (`` s `` measured from that layer's own left edge, `` k_i=\\omega/c_i ``)
across every interface, starting from `` A_1=1,\\,B_1=0 `` (displacement fixed to zero at
the string's left end, amplitude normalized arbitrarily -- rescaled later). Matching
`` U `` and stress `` \\mu U'=c^2U' `` at each interface gives the two update equations
below. Returns the full per-layer `` (A_i,B_i) `` lists (not just the final layer's), so
[`eigenfunction`](@ref) can evaluate the mode shape anywhere, not only at the right end.
"""
function chain_transfer_matrix(ω, velocities, widths)
    n = length(velocities)
    As, Bs = zeros(n), zeros(n)
    As[1], Bs[1] = 1.0, 0.0
    for i in 1:n-1
        k = ω / velocities[i]
        s, cph = sincos(k * widths[i])
        As[i+1] = (velocities[i] / velocities[i+1]) * (As[i] * cph - Bs[i] * s)
        Bs[i+1] = As[i] * s + Bs[i] * cph
    end
    return As, Bs
end

# ╔═╡ 2731c578-e51d-4fb0-8ca7-3216d49f649a
"""
	find_eigenfrequencies(velocities, boundaries; nmodes, ωmax, dω)

The `nmodes` lowest positive roots of `` U(x{=}1)=0 `` (the right-end fixed boundary
condition), found by scanning `` \\omega \\in (0,\\omega_{\\max}] `` in steps of `` d\\omega ``
for sign changes and bisecting each one -- robust and simple, and cheap enough at this
problem size that a smarter root-finder isn't worth the complexity. `boundaries` (length
`length(velocities)-1`) are the layer-interface positions in `(-1,1)`, increasing left to
right.
"""
function find_eigenfrequencies(velocities, boundaries; nmodes=40, ωmax=80.0, dω=0.002)
    xs = vcat(-1.0, boundaries, 1.0)
    widths = diff(xs)
    function residual(ω)
        As, Bs = chain_transfer_matrix(ω, velocities, widths)
        k = ω / velocities[end]
        s, cph = sincos(k * widths[end])
        return As[end] * s + Bs[end] * cph
    end
    ws = 1.0e-6:dω:ωmax
    roots = Float64[]
    prev = residual(ws[1])
    for i in 2:length(ws)
        cur = residual(ws[i])
        if sign(cur) != sign(prev) && !isnan(cur) && !isnan(prev)
            a, b, fa = ws[i-1], ws[i], prev
            for _ in 1:80
                m = (a + b) / 2
                fm = residual(m)
                sign(fm) == sign(fa) ? ((a, fa) = (m, fm)) : (b = m)
            end
            push!(roots, (a + b) / 2)
            length(roots) >= nmodes && break
        end
        prev = cur
    end
    return roots
end

# ╔═╡ 7056f891-2b79-4da5-95e9-b133fa86043f
"""
	eigenfunction(xsample, ω, velocities, boundaries)

Evaluate the (un-normalized, `` A_1{=}1 ``) mode shape at eigenfrequency `ω` on the sample
points `xsample`, by locating which layer each point falls in and evaluating that layer's
own `` A_i\\sin(k_is)+B_i\\cos(k_is) `` there.
"""
function eigenfunction(xsample, ω, velocities, boundaries)
    xs = vcat(-1.0, boundaries, 1.0)
    widths = diff(xs)
    As, Bs = chain_transfer_matrix(ω, velocities, widths)
    u = similar(xsample)
    for (j, x) in enumerate(xsample)
        i = clamp(searchsortedlast(xs, x), 1, length(velocities))
        s = x - xs[i]
        k = ω / velocities[i]
        u[j] = As[i] * sin(k * s) + Bs[i] * cos(k * s)
    end
    return u
end

# ╔═╡ 29ca2faf-d4e0-4e9d-ba07-525aba0b65c5
md"""
### Normalized Modes, Source Projection, and Summation

Sturm–Liouville theory guarantees the *continuous* problem's eigenfunctions are orthogonal;
normalizing each to unit energy, `` \int U_n^2\,dx=1 ``, makes `` \{U_n\} `` a genuine
orthonormal basis (checked numerically below). A delta-function source at `` x_s `` then
excites mode `` n `` with amplitude exactly `` c_n=\int\delta(x-x_s)U_n(x)\,dx=U_n(x_s) ``
-- no matrix inverse required, unlike projecting onto a *non-orthogonal* numerical basis
(which is what naively transposing a non-symmetric matrix's eigenvectors would silently
get wrong).
"""

# ╔═╡ c728bc06-abc9-4a60-97a4-598b721e78bc
begin
    """
    	normalized_eigenfunctions(xsample, ωs, velocities, boundaries)

    [`eigenfunction`](@ref) for every `ω` in `ωs`, each rescaled to unit `` L^2 `` norm
    (`` \\int U_n^2\\,dx=1 ``, trapezoidal quadrature on `xsample`) so they form an
    orthonormal set. Returns an `(length(xsample), length(ωs))` matrix, one mode per column.
    """
    function normalized_eigenfunctions(xsample, ωs, velocities, boundaries)
        dx = xsample[2] - xsample[1]
        U = reduce(hcat, [eigenfunction(xsample, ω, velocities, boundaries) for ω in ωs])
        norms = sqrt.(sum(abs2, U; dims=1) .* dx)
        return U ./ norms
    end

    """
    	mode_series(U, ωs, coeffs, tgrid)

    The displacement field `` u(x,t)=\\sum_n c_n\\,U_n(x)\\cos(\\omega_n t) `` at every time
    in `tgrid`, for normalized modes `U` (columns), eigenfrequencies `ωs`, and per-mode
    coefficients `coeffs` (pass a one-hot vector to isolate a single mode; pass
    `U[isrc, :]` -- the source-position row -- for the source-weighted superposition).
    Returns an `(length(tgrid), size(U,1))` array, matching `plot_up_frame`'s layout.
    """
    function mode_series(U, ωs, coeffs, tgrid)
        nt, nx = length(tgrid), size(U, 1)
        up = Matrix{Float64}(undef, nt, nx)
        weighted = U .* coeffs'   # (nx, nmodes), each mode's spatial shape pre-scaled by c_n
        for (it, t) in enumerate(tgrid)
            up[it, :] .= weighted * cos.(ωs .* t)
        end
        return up
    end
end

# ╔═╡ b6a4c9d1-70e1-4e2a-9c1d-3a5b9e0d2f11
md"""
### Medium Setup
"""

# ╔═╡ 94568f52-b993-4f62-8170-776e83945438
"""
	nearest_index(rng, value)

The index into the range/vector `rng` closest to `value` -- used to snap a continuous
source/receiver position to the sample grid.
"""
nearest_index(rng, value) = argmin(abs.(rng .- value))

# ╔═╡ 8f6fb806-7877-48d5-989b-9c9944e35ce5
md"""
### Verifying the Solver
"""

# ╔═╡ d7db01d4-f552-4c65-b122-ba392943f3b3
let
    c0 = 2.0
    L = 2.0
    ω_exact = find_eigenfrequencies([c0], Float64[]; nmodes=10)
    ω_analytic = [n * pi * c0 / L for n in 1:10]
    residual = maximum(abs.(ω_exact .- ω_analytic) ./ ω_analytic)
    @assert residual < 1.0e-9
    md"""
    !!! correct "Homogeneous string: exact solver matches the closed form"
    	Fixed-fixed string, `` c=2 ``\ km/s, length `` L=2 ``\ km:
    	`` \omega_n=n\pi c/L ``. Maximum relative error over the first 10 modes,
    	transfer-matrix solver vs. the closed form: **$(round(residual, sigdigits=3))** --
    	machine precision.
    """
end

# ╔═╡ ca1d876b-0b59-47c5-944c-4ff45e40ff6d
let
    xsample = range(-1, 1; length=2001)
    velocities, boundaries = [1.8, 2.6, 3.1], [-0.3, 0.4]
    ωs = find_eigenfrequencies(velocities, boundaries; nmodes=8)
    U = normalized_eigenfunctions(collect(xsample), ωs, velocities, boundaries)
    dx = step(xsample)
    G = (U' * U) .* dx
    offdiag_max = maximum(abs.(G .- I(size(G, 1))))
    @assert offdiag_max < 1.0e-4
    md"""
    !!! correct "Layered string (3 layers): normalized modes are orthonormal"
    	`` \int U_mU_n\,dx `` for the first 8 modes of a 3-layer string (velocities
    	1.8/2.6/3.1 km/s), compared to the identity matrix: maximum deviation
    	**$(round(offdiag_max, sigdigits=3))** -- confirms the source-projection formula
    	`` c_n=U_n(x_s) `` above is valid, and that these modes really are the physically
    	correct orthogonal basis Sturm–Liouville theory promises.
    """
end

# ╔═╡ dbf086b6-8b44-4f7e-ad99-7a31612c299a
md"""
### Pushing Results to the Widget

`ModeExplorerPush` does no physics -- it takes the already-computed displacement field,
mode spectrum, and receiver record, and hands them to the *already-rendered*
[`ModeExplorerInput`](@ref) widget via a `CustomEvent`, the same pattern this repo's other
widgets use.
"""

# ╔═╡ b5b41af3-2699-43d9-8f60-2cac6df2d040
_nm_flatten(v) = join(v, ",")

# ╔═╡ 8eb5c988-5471-44df-b0a2-38196a3fe828
begin
    struct ModeExplorerInput
        selLo::Int               # 1-based, first mode index in the selected range
        selHi::Int               # 1-based, last mode index; selLo==selHi means a single mode
        velocities::Vector{Float64}
        boundaries::Vector{Float64}
        srcX::Float64
        recX::Float64
    end
    ModeExplorerInput(; selLo=1, selHi=1, velocities=[2.0],
        boundaries=Float64[], srcX=-0.4, recX=0.55) =
        ModeExplorerInput(selLo, selHi, velocities, boundaries, srcX, recX)

    Base.get(w::ModeExplorerInput) = Dict{String,Any}(
        "selLo" => w.selLo, "selHi" => w.selHi,
        "velocities" => w.velocities, "boundaries" => w.boundaries,
        "srcX" => w.srcX, "recX" => w.recX)

    function Base.show(io::IO, ::MIME"text/html", w::ModeExplorerInput)
        write(io, """
        <div id="nmwidget">
        <style>
        pluto-cell:has(#nmwidget) { width: min(80vw, 1500px) !important;
          margin-left: calc((100% - min(80vw, 1500px)) / 2) !important; }
        #nmwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #nmwidget .nm-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #nmwidget .nm-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #nmwidget .nm-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #nmwidget .nm-primary{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start;margin-bottom:14px}
        #nmwidget .nm-controls-row{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:stretch;margin-bottom:14px}
        #nmwidget .nm-secondary{display:flex;justify-content:center}
        #nmwidget .nm-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #nmwidget .nm-panel-title{font-size:14px;font-weight:700;color:#e5e7eb;text-align:center;margin-bottom:6px}
        #nmwidget .nm-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:4px}
        #nmwidget canvas{display:block;cursor:crosshair}
        #nmwidget .nm-control-group{flex:1 1 220px;min-width:200px;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #nmwidget .nm-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #nmwidget .nm-control-row{display:grid;grid-template-columns:minmax(60px,90px) minmax(70px,1fr) minmax(40px,58px);gap:6px;align-items:center;margin:5px 0}
        #nmwidget .nm-control-row label{font-size:13px;color:#9ca3af}
        #nmwidget .nm-control-row input[type=range]{width:100%;min-width:0}
        #nmwidget .nm-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #nmwidget .nm-actions{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
        #nmwidget select{background:#0b0b0b;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:3px;width:100%}
        #nmwidget .nm-vinput{position:absolute;width:52px;font-size:11px;text-align:center;pointer-events:auto;
          background:#0b0b0bcc;color:#e5e7eb;border-width:1.5px;border-style:solid;border-radius:4px;padding:1px 2px}
        #nmwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:13px;cursor:pointer}
        #nmwidget button.active{background:#2563eb;border-color:#93c5fd}
        </style>

        <div class="nm-title">
          <div class="nm-title-desc">Watch one normal mode oscillate, or sum many of them into a traveling pulse.</div>
          <div class="nm-title-hint">click a mode below to isolate it, or drag across several to sum them &middot; drag the star (source) or triangle (receiver) &middot; Play to animate</div>
        </div>

        <div class="nm-primary">
          <div>
            <div class="nm-panel-title">1-D Earth Displacement</div>
            <div class="nm-panel" style="position:relative">
              <canvas id="nm-string"></canvas>
              <div id="nm-vlabels" style="position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none"></div>
            </div>
            <div class="nm-caption" id="nm-string-caption"></div>
          </div>
          <div>
            <div class="nm-panel-title">Mode Spectrum</div>
            <div class="nm-panel"><canvas id="nm-spectrum"></canvas></div>
            <div class="nm-caption" id="nm-spectrum-caption"></div>
          </div>
        </div>

        <div class="nm-controls-row">
          <div class="nm-control-group">
            <div class="nm-control-title">Playback</div>
            <div class="nm-caption" style="margin:0 0 10px 0;text-align:left">click a mode in the Mode Spectrum panel to view it alone, or drag across a range to sum them.</div>
            <div class="nm-actions"><button id="nm-play" type="button">Play</button><button id="nm-reset" type="button">Reset</button></div>
          </div>
          <div class="nm-control-group">
            <div class="nm-control-title">Medium</div>
            <div class="nm-control-row"><label>layers</label><select id="nm-nlayers"><option value="1">1</option><option value="2">2</option><option value="3">3</option><option value="4">4</option><option value="5">5</option></select><span></span></div>
            <div class="nm-caption" style="margin:8px 0 0 0;text-align:left">edit a layer's speed directly on the plot above &middot; drag a dashed boundary line to move it</div>
          </div>
        </div>

        <div class="nm-secondary">
          <div>
            <div class="nm-panel-title">Receiver Record</div>
            <div class="nm-panel"><canvas id="nm-record"></canvas></div>
            <div class="nm-caption">time series (left) vs. amplitude spectrum with eigenfrequencies marked (right)</div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        const XMIN=-1, XMAX=1, MAXMODES=40;
        const LAYER_COLORS = ['#3b82f6','#f97316','#22c55e','#a855f7','#eab308'];

        let selLo = $(w.selLo), selHi = $(w.selHi);
        let velocities = [$(join(w.velocities, ","))];
        let boundaries = [$(join(w.boundaries, ","))];
        let srcX = $(w.srcX), recX = $(w.recX);

        const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1500);
        const GAP = 16, PANEL_OVERHEAD = 18;
        const SEC_W = Math.max(280, Math.floor((availW - GAP - 2*PANEL_OVERHEAD) / 2));
        const SEC_H = Math.round(SEC_W * 0.42);
        const DPR = window.devicePixelRatio || 1;

        function hidpi(canvas, context, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          context.setTransform(DPR,0,0,DPR,0,0);
        }

        const cv = par.querySelector('#nm-string');
        const ctx = cv.getContext('2d');
        hidpi(cv, ctx, SEC_W, SEC_H);

        const spCv = par.querySelector('#nm-spectrum');
        const spCtx = spCv.getContext('2d');
        hidpi(spCv, spCtx, SEC_W, SEC_H);

        const RECORD_W = availW - PANEL_OVERHEAD, RECORD_H = Math.round((availW - PANEL_OVERHEAD) * 0.28);
        const rCv = par.querySelector('#nm-record');
        const rCtx = rCv.getContext('2d');
        hidpi(rCv, rCtx, RECORD_W, RECORD_H);

        function toScreen(x){ return (x-XMIN)/(XMAX-XMIN)*SEC_W; }
        function toWorld(px){ return XMIN + px/SEC_W*(XMAX-XMIN); }

        function drawXAxis(c, W, H, xmin, xmax, step, unit){
          const bandH = 16;
          c.fillStyle = 'rgba(0,0,0,0.55)'; c.fillRect(0, H-bandH, W, bandH);
          c.strokeStyle = '#4b5563'; c.beginPath(); c.moveTo(0,H-bandH+0.5); c.lineTo(W,H-bandH+0.5); c.stroke();
          c.fillStyle = '#9ca3af'; c.font = '9px sans-serif'; c.textBaseline = 'top';
          for(let x=xmin; x<=xmax+1e-9; x+=step){
            const px = Math.min(W, (x-xmin)/(xmax-xmin)*W);
            c.strokeStyle = '#6b7280'; c.beginPath(); c.moveTo(px,H-bandH); c.lineTo(px,H-bandH+4); c.stroke();
            const label = (Math.round(x*100)/100)+(unit||'');
            const tw = c.measureText(label).width;
            const tx = Math.max(2, Math.min(W - tw - 2, px - tw/2));
            c.textAlign = 'left';
            c.fillText(label, tx, H-bandH+5);
          }
        }
        function drawYAxis(c, H, ymin, ymax, step, unit, digits){
          const bandW = 34;
          c.fillStyle = 'rgba(0,0,0,0.55)'; c.fillRect(0,0,bandW,H);
          c.strokeStyle = '#4b5563'; c.beginPath(); c.moveTo(bandW+0.5,0); c.lineTo(bandW+0.5,H); c.stroke();
          c.fillStyle = '#9ca3af'; c.font = '9px sans-serif'; c.textAlign = 'left';
          for(let y=ymin; y<=ymax+1e-9; y+=step){
            const py = H - Math.min(H, (y-ymin)/(ymax-ymin)*H);
            c.strokeStyle = '#6b7280'; c.beginPath(); c.moveTo(bandW-4,py); c.lineTo(bandW,py); c.stroke();
            const label = y.toFixed(digits||0)+(unit||'');
            const tw = c.measureText(label).width;
            const tx = Math.max(2, bandW-4-tw);
            const ty = Math.max(1, Math.min(H-1, py));
            c.textBaseline = py<=1 ? 'top' : (py>=H-1 ? 'bottom' : 'middle');
            c.fillText(label, tx, ty);
          }
        }

        function drawStarMarker(c, cx, cy, r, fill, stroke){
          const spikes = 5, rOuter = r, rInner = r * 0.45;
          c.beginPath();
          for(let i=0; i<spikes*2; i++){
            const rad = i % 2 === 0 ? rOuter : rInner;
            const ang = -Math.PI/2 + i*Math.PI/spikes;
            const x = cx + rad*Math.cos(ang), y = cy + rad*Math.sin(ang);
            i===0 ? c.moveTo(x,y) : c.lineTo(x,y);
          }
          c.closePath(); c.fillStyle = fill; c.fill(); c.strokeStyle = stroke; c.lineWidth = 1; c.stroke();
        }
        function drawTriangleDownMarker(c, cx, cy, r, fill, stroke){
          c.beginPath();
          for(let i=0; i<3; i++){
            const ang = Math.PI/2 + i*2*Math.PI/3;
            const x = cx + r*Math.cos(ang), y = cy + r*Math.sin(ang);
            i===0 ? c.moveTo(x,y) : c.lineTo(x,y);
          }
          c.closePath(); c.fillStyle = fill; c.fill(); c.strokeStyle = stroke; c.lineWidth = 1.5; c.stroke();
        }

        let pushed = null; // {xgrid,up,nt,nx,omegas,record,spectrum,freqaxis,upmax,recordmax,specmax}
        let playing = false, rafId = null, tPhase = 0, lastTs = null;
        let dragMode = null; // 'src' | 'rec'
        let specAnchor = null; // mode index where a spectrum click/drag started

        function drawString(){
          ctx.clearRect(0,0,SEC_W,SEC_H);
          const midY = SEC_H*0.5;
          // one colored band per layer, so a painted velocity structure reads visually
          // on the string itself -- each band's own speed is editable in place via the
          // number field positioned over it (see positionVelocityInputs), not a slider
          const edges = [XMIN, ...boundaries, XMAX];
          for(let i=0; i<velocities.length; i++){
            const x0 = toScreen(edges[i]), x1 = toScreen(edges[i+1]);
            ctx.fillStyle = LAYER_COLORS[i % LAYER_COLORS.length];
            ctx.globalAlpha = 0.16;
            ctx.fillRect(x0, 0, x1-x0, SEC_H);
            ctx.globalAlpha = 1;
          }
          // layer-boundary lines -- draggable (see mousedown/mousemove below)
          ctx.strokeStyle = '#9ca3af'; ctx.setLineDash([4,3]); ctx.lineWidth = 1.5;
          for(const b of boundaries){ const px = toScreen(b); ctx.beginPath(); ctx.moveTo(px,0); ctx.lineTo(px,SEC_H); ctx.stroke(); }
          ctx.setLineDash([]); ctx.lineWidth = 1;
          positionVelocityInputs(edges);
          if(pushed){
            // interpolate between the two nearest precomputed time samples rather than
            // rounding to the nearest one -- at the slowed-down playback speed a fast
            // mode needs (see currentSimSpeed), tPhase can sit on the same rounded frame
            // for many real animation frames in a row, which reads as stutter/lag even
            // though nothing is actually slow; a fractional blend moves continuously.
            const frameF = Math.min(pushed.nt-1, Math.max(0, tPhase/pushed.tmax*(pushed.nt-1)));
            const f0 = Math.floor(frameF), f1 = Math.min(pushed.nt-1, f0+1), frac = frameF-f0;
            const base0 = f0*pushed.nx, base1 = f1*pushed.nx;
            const mx = pushed.upmax;
            const amp = SEC_H*0.38/mx;
            ctx.beginPath();
            for(let ix=0; ix<pushed.nx; ix++){
              const x = toScreen(pushed.xgrid[ix]);
              const y = midY - (pushed.up[base0+ix]*(1-frac) + pushed.up[base1+ix]*frac)*amp;
              ix===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
            }
            ctx.strokeStyle = '#38bdf8'; ctx.lineWidth = 2; ctx.stroke();
          }
          ctx.strokeStyle = '#4b5563'; ctx.beginPath(); ctx.moveTo(0,midY); ctx.lineTo(SEC_W,midY); ctx.stroke();
          drawStarMarker(ctx, toScreen(srcX), midY-SEC_H*0.42, 10, '#facc15', '#0a0f18');
          drawTriangleDownMarker(ctx, toScreen(recX), midY-SEC_H*0.42, 8, '#f5f3ef', '#0a0f18');
          drawXAxis(ctx, SEC_W, SEC_H, XMIN, XMAX, 0.5, ' km');
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1; ctx.strokeRect(0.5,0.5,SEC_W-1,SEC_H-1);
          const capEl = par.querySelector('#nm-string-caption');
          if(selLo===selHi) capEl.textContent = 'mode ' + selLo + (pushed ? ' \\u00b7 \\u03c9 = ' + pushed.omegas[selLo-1].toFixed(3) + ' rad/s' : '');
          else capEl.textContent = 'modes ' + selLo + '\\u2013' + selHi + ' summed, weighted by source position';
        }

        function drawSpectrum(){
          spCtx.clearRect(0,0,SEC_W,SEC_H);
          if(!pushed){
            spCtx.fillStyle = '#6b7280'; spCtx.font = '12px sans-serif'; spCtx.fillText('computing...', 10, 18);
            return;
          }
          const n = pushed.omegas.length;
          const wmax = Math.max(...pushed.omegas);
          const bandL = 34, bandB = 16;
          const plotW = SEC_W - bandL - 6, plotH = SEC_H - bandB - 6;
          const stemW = plotW/n;
          spCtx.fillStyle = '#facc15'; spCtx.globalAlpha = 0.14;
          spCtx.fillRect(bandL+(selLo-1)*stemW, 0, (selHi-selLo+1)*stemW, SEC_H-bandB);
          spCtx.globalAlpha = 1;
          for(let i=0;i<n;i++){
            const active = (i+1>=selLo && i+1<=selHi);
            const px = bandL + (i+0.5)*stemW;
            const py = SEC_H - bandB - (pushed.omegas[i]/wmax)*plotH;
            spCtx.strokeStyle = active ? '#facc15' : '#4b5563';
            spCtx.lineWidth = active ? 2.5 : 1.5;
            spCtx.beginPath(); spCtx.moveTo(px, SEC_H-bandB); spCtx.lineTo(px, py); spCtx.stroke();
            spCtx.fillStyle = active ? '#facc15' : '#9ca3af';
            spCtx.beginPath(); spCtx.arc(px, py, active?3.5:2, 0, 2*Math.PI); spCtx.fill();
          }
          drawYAxis(spCtx, SEC_H-bandB, 0, wmax, wmax/4, '', 1);
          spCtx.strokeStyle = '#9ca3af'; spCtx.font = '9px sans-serif'; spCtx.textAlign='center'; spCtx.textBaseline='top';
          spCtx.fillStyle = '#9ca3af';
          spCtx.fillText('mode index n', bandL+plotW/2, SEC_H-bandB+5);
          spCtx.strokeStyle = '#374151'; spCtx.lineWidth = 1; spCtx.strokeRect(0.5,0.5,SEC_W-1,SEC_H-1);
          par.querySelector('#nm-spectrum-caption').textContent = '\\u03c9 (rad/s) vs. mode index \\u00b7 click or drag to select';
        }

        function pxToMode(px){
          if(!pushed) return null;
          const n = pushed.omegas.length;
          const bandL = 34, bandB = 16;
          const plotW = SEC_W - bandL - 6;
          const i = Math.round((px-bandL)/plotW*n - 0.5);
          return Math.max(1, Math.min(n, i+1));
        }

        function drawRecord(){
          rCtx.clearRect(0,0,RECORD_W,RECORD_H);
          rCtx.strokeStyle = '#374151'; rCtx.lineWidth = 1; rCtx.strokeRect(0.5,0.5,RECORD_W-1,RECORD_H-1);
          if(!pushed){
            rCtx.fillStyle = '#6b7280'; rCtx.font = '12px sans-serif'; rCtx.fillText('computing...', 10, 18);
            return;
          }
          const halfW = Math.floor(RECORD_W/2);
          const bandL = 34, bandB = 16;
          const plotH = RECORD_H - bandB;
          const plotW = halfW - bandL - 4;
          const nt = pushed.record.length;
          // time series (left half)
          {
            const mx = pushed.recordmax, midY = plotH*0.5, amp = plotH*0.42/mx;
            rCtx.beginPath();
            for(let it=0; it<nt; it++){
              const x = bandL + it/(nt-1)*plotW, y = midY - pushed.record[it]*amp;
              it===0 ? rCtx.moveTo(x,y) : rCtx.lineTo(x,y);
            }
            rCtx.strokeStyle = '#38bdf8'; rCtx.lineWidth = 1.5; rCtx.stroke();
            rCtx.strokeStyle = '#4b5563'; rCtx.beginPath(); rCtx.moveTo(bandL,midY); rCtx.lineTo(bandL+plotW,midY); rCtx.stroke();
            drawYAxis(rCtx, plotH, -mx, mx, mx/2, '', 2);
            const tmax = pushed.tmax, tStep = tmax>3 ? 1 : (tmax>1 ? 0.5 : 0.2);
            rCtx.save(); rCtx.translate(bandL,0);
            drawXAxis(rCtx, plotW, RECORD_H, 0, tmax, tStep, ' s');
            rCtx.restore();
          }
          // amplitude spectrum (right half) with eigenfrequency markers
          {
            const nf = pushed.spectrum.length;
            const fmax = pushed.freqaxis[nf-1];
            const mx = pushed.specmax;
            rCtx.strokeStyle = '#ef4444';
            for(const om of pushed.omegas){
              const f = om/(2*Math.PI);
              if(f>fmax) continue;
              const x = halfW + bandL + f/fmax*plotW;
              rCtx.globalAlpha = 0.35; rCtx.beginPath(); rCtx.moveTo(x,0); rCtx.lineTo(x,plotH); rCtx.stroke(); rCtx.globalAlpha=1;
            }
            rCtx.beginPath();
            for(let i=0;i<nf;i++){
              const x = halfW + bandL + i/(nf-1)*plotW, y = plotH - (pushed.spectrum[i]/mx)*plotH*0.92;
              i===0 ? rCtx.moveTo(x,y) : rCtx.lineTo(x,y);
            }
            rCtx.strokeStyle = '#38bdf8'; rCtx.lineWidth = 1.5; rCtx.stroke();
            rCtx.save(); rCtx.translate(halfW,0);
            drawYAxis(rCtx, plotH, 0, mx, mx/4, '', 1);
            rCtx.restore();
            rCtx.save(); rCtx.translate(halfW+bandL,0);
            drawXAxis(rCtx, plotW, RECORD_H, 0, fmax, fmax>10?5:1, ' Hz');
            rCtx.restore();
          }
          rCtx.strokeStyle = '#4b5563'; rCtx.beginPath(); rCtx.moveTo(halfW+0.5,0); rCtx.lineTo(halfW+0.5,RECORD_H); rCtx.stroke();
          rCtx.fillStyle = '#e5e7eb'; rCtx.font = '11px sans-serif'; rCtx.textAlign='left';
          rCtx.fillText('u(t) at receiver', bandL+2, 12);
          rCtx.fillText('|U(f)|, eigenfrequencies in red', halfW+bandL+2, 12);
        }

        function draw(){ drawString(); drawSpectrum(); drawRecord(); }

        // one <input type=number> per layer, absolutely positioned over its own colored
        // band on the String Displacement canvas -- rebuilt only when the layer count
        // changes (ensureVelocityInputs), repositioned every draw (positionVelocityInputs,
        // cheap style-only updates so an in-progress edit never loses focus).
        const vlabelsEl = par.querySelector('#nm-vlabels');
        function ensureVelocityInputs(){
          while(vlabelsEl.children.length > velocities.length){ vlabelsEl.removeChild(vlabelsEl.lastChild); }
          while(vlabelsEl.children.length < velocities.length){
            const idx = vlabelsEl.children.length;
            const inp = document.createElement('input');
            inp.type = 'number'; inp.min = '1'; inp.max = '4'; inp.step = '0.1';
            inp.className = 'nm-vinput'; inp.title = 'wavespeed, km/s';
            inp.dataset.i = idx;
            inp.style.borderColor = LAYER_COLORS[idx % LAYER_COLORS.length];
            vlabelsEl.appendChild(inp);
          }
          for(let i=0;i<velocities.length;i++){ vlabelsEl.children[i].value = velocities[i].toFixed(1); }
        }
        function positionVelocityInputs(edges){
          for(let i=0;i<velocities.length;i++){
            const inp = vlabelsEl.children[i];
            if(!inp) continue;
            const x0 = toScreen(edges[i]), x1 = toScreen(edges[i+1]);
            inp.style.left = Math.round((x0+x1)/2 - 26) + 'px';
            inp.style.top = '4px';
          }
        }
        vlabelsEl.addEventListener('input', e=>{
          if(!e.target.classList.contains('nm-vinput')) return;
          const i = +e.target.dataset.i, v = +e.target.value;
          if(!isFinite(v)) return;
          velocities[i] = Math.max(1, Math.min(4, v));
          draw(); throttledCommit();
        });
        vlabelsEl.addEventListener('change', e=>{
          if(!e.target.classList.contains('nm-vinput')) return;
          commit();
        });

        function syncControls(){
          par.querySelector('#nm-nlayers').value = velocities.length;
          ensureVelocityInputs();
        }

        let commitInFlight = false;
        function commit(){
          commitInFlight = true;
          par.value = { selLo, selHi, velocities, boundaries, srcX, recX };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledCommit(){ if(!commitInFlight) commit(); }

        par.addEventListener('nm-update', e=>{
          pushed = e.detail;
          commitInFlight = false;
          draw();
        });

        function hitMarker(px, x0){
          return Math.abs(px - toScreen(x0)) < 12;
        }
        cv.addEventListener('mousedown', e=>{
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left;
          if(hitMarker(px, srcX)){ dragMode='src'; return; }
          if(hitMarker(px, recX)){ dragMode='rec'; return; }
          for(let i=0;i<boundaries.length;i++){
            if(hitMarker(px, boundaries[i])){ dragMode='bnd'+i; return; }
          }
        });
        spCv.addEventListener('mousedown', e=>{
          const rect = spCv.getBoundingClientRect();
          const m = pxToMode(e.clientX-rect.left);
          if(m===null) return;
          specAnchor = m; selLo = m; selHi = m;
          draw(); throttledCommit();
        });
        window.addEventListener('mousemove', e=>{
          if(dragMode){
            const rect = cv.getBoundingClientRect();
            const px = Math.max(0, Math.min(SEC_W, e.clientX-rect.left));
            const x = Math.max(XMIN, Math.min(XMAX, toWorld(px)));
            if(dragMode==='src'){ srcX = x; }
            else if(dragMode==='rec'){ recX = x; }
            else {
              const i = +dragMode.slice(3);
              const lo = (i===0 ? XMIN : boundaries[i-1]) + 0.05;
              const hi = (i===boundaries.length-1 ? XMAX : boundaries[i+1]) - 0.05;
              boundaries[i] = Math.max(lo, Math.min(hi, x));
            }
            draw(); throttledCommit();
          }
          if(specAnchor!==null){
            const rect = spCv.getBoundingClientRect();
            const m = pxToMode(Math.max(0, Math.min(SEC_W, e.clientX-rect.left)));
            if(m!==null){ selLo = Math.min(specAnchor, m); selHi = Math.max(specAnchor, m); draw(); throttledCommit(); }
          }
        });
        window.addEventListener('mouseup', ()=>{
          if(dragMode) commit();
          dragMode = null;
          if(specAnchor!==null) commit();
          specAnchor = null;
        });

        par.querySelector('#nm-nlayers').addEventListener('change', e=>{
          const n = +e.target.value;
          while(velocities.length < n){ velocities.push(2.0 + 0.4*velocities.length); }
          velocities.length = n;
          const defaultsB = [-0.6,-0.2,0.2,0.6];
          while(boundaries.length < n-1){ boundaries.push(defaultsB[boundaries.length]); }
          boundaries.length = Math.max(0, n-1);
          syncControls(); draw(); commit();
        });

        par.addEventListener('input', e=>{
          if(e.target===par) return;
          const id = e.target.id, v = +e.target.value;
          const cMatch = id.match(/^nm-c(\\d)\$/);
          if(cMatch){ e.stopImmediatePropagation(); const i = +cMatch[1]-1; velocities[i] = v; par.querySelector('#'+id+'-v').textContent = v.toFixed(1); throttledCommit(); return; }
          const bMatch = id.match(/^nm-b(\\d)\$/);
          if(bMatch){ e.stopImmediatePropagation(); const i = +bMatch[1]-1; boundaries[i] = v; par.querySelector('#'+id+'-v').textContent = v.toFixed(2); throttledCommit(); return; }
        }, true);

        par.addEventListener('change', e=>{
          if(e.target===par) return;
          const id = e.target.id;
          if(/^nm-c\\d\$/.test(id)||/^nm-b\\d\$/.test(id)) { e.stopImmediatePropagation(); commit(); }
        }, true);

        // Real time (SIM_SPEED=1) is a comfortable pace for the low modes, but at 40
        // modes omega_40 is ~25x omega_1 -- played back at real speed it just flickers.
        // Slow the clock down for whatever's currently selected so its fastest included
        // mode never completes more than one cycle roughly every TARGET_PERIOD seconds;
        // never speed a slow mode up past real time, only ever slow a fast one down.
        const TARGET_OMEGA = 2*Math.PI/1.6;
        function currentSimSpeed(){
          if(!pushed) return 1.0;
          const wMax = pushed.omegas[selHi-1];
          return Math.min(1.0, TARGET_OMEGA/wMax);
        }
        const playBtn = par.querySelector('#nm-play');
        function stepAnim(ts){
          if(lastTs===null) lastTs = ts;
          const dt = Math.min(0.1, (ts-lastTs)/1000);
          lastTs = ts;
          if(pushed){
            tPhase += dt*currentSimSpeed();
            if(tPhase > pushed.tmax) tPhase = 0;
          }
          draw();
          rafId = requestAnimationFrame(stepAnim);
        }
        playBtn.addEventListener('click', ()=>{
          playing = !playing;
          playBtn.textContent = playing ? 'Pause' : 'Play';
          if(playing){ lastTs=null; rafId = requestAnimationFrame(stepAnim); }
          else if(rafId){ cancelAnimationFrame(rafId); rafId=null; }
        });
        par.querySelector('#nm-reset').addEventListener('click', ()=>{
          selLo = $(w.selLo); selHi = $(w.selHi);
          velocities = [$(join(w.velocities, ","))]; boundaries = [$(join(w.boundaries, ","))];
          srcX = $(w.srcX); recX = $(w.recX); tPhase = 0;
          syncControls(); draw(); commit();
        });

        syncControls(); draw();
        }
        </script>
        """)
    end

    const _nm_ready = true
end

# ╔═╡ 27ac2670-a5b2-431d-9626-c95241a00b55
begin
    _nm_ready
    WideCell(@bind nm ModeExplorerInput(); max_width=1500)
end

# ╔═╡ 4b695e82-0860-4193-ae03-8d317b7280c7
# The bond starts as `nothing` until the widget's first real interaction in a live browser
# reports back -- fall back to the same defaults the widget itself opens with.
nm_safe = nm isa AbstractDict ? nm : Dict{String,Any}(
    "selLo" => 1, "selHi" => 1,
    "velocities" => [2.0], "boundaries" => Float64[], "srcX" => -0.4, "recX" => 0.55)

# ╔═╡ d5670bb8-aab6-421b-83fa-bed20faf7c33
let
    nl = length(nm_safe["velocities"])
    lo, hi = round(Int, nm_safe["selLo"]), round(Int, nm_safe["selHi"])
    readout = "Viewing " * (lo == hi ?
        "mode **$lo** in isolation" :
        "modes **$lo**–**$hi** summed") *
        " of a **$nl**-layer string (velocities " * join(round.(nm_safe["velocities"], digits=1), ", ") * " km/s)."
    Markdown.parse(readout)
end

# ╔═╡ 1fa011bc-69bc-4e4a-aa00-b89672bbd376
# a fixed, generous sample resolution for plotting/animating the (exact, closed-form) mode
# shapes -- the transfer-matrix solver itself has no notion of a discretization grid at
# all, so this is purely a plotting choice, not a numerics one.
const _nm_xgrid = collect(range(-1.0, 1.0; length=400))

# ╔═╡ 089009a9-41e7-406e-9eea-5faa31126c6e
_nm_velocities, _nm_boundaries = Float64.(nm_safe["velocities"]), Float64.(nm_safe["boundaries"])

# ╔═╡ 0c723ed7-79f6-4334-b575-baa15b6dc6b5
# the expensive step (a root-search over omega, then normalizing each mode) -- gated behind
# the medium only, so dragging the source/receiver or scrubbing the animation never
# re-triggers it.
_nm_omegas = find_eigenfrequencies(_nm_velocities, _nm_boundaries; nmodes=40)

# ╔═╡ 10ada4dd-c11c-4605-8b68-0e550279234f
_nm_U = normalized_eigenfunctions(_nm_xgrid, _nm_omegas, _nm_velocities, _nm_boundaries)

# ╔═╡ 476860e2-3ba3-4132-8848-09cd2e6ae52a
const _nm_tgrid = range(0.0, step=0.02, length=250)

# ╔═╡ ab6ad92c-2d8c-4cee-8643-539b71f9b462
_nm_selLo = clamp(round(Int, nm_safe["selLo"]), 1, length(_nm_omegas))

# ╔═╡ ead6d118-a3a8-4ca6-8cf5-7af2997d24f5
_nm_selHi = clamp(round(Int, nm_safe["selHi"]), _nm_selLo, length(_nm_omegas))

# ╔═╡ 6fe6ac9c-e7fa-4c75-a10b-68128fb786be
_nm_isrc = nearest_index(_nm_xgrid, nm_safe["srcX"])

# ╔═╡ 9ae40035-5296-4cf5-a4fb-be2619ddc992
_nm_irec = nearest_index(_nm_xgrid, nm_safe["recX"])

# ╔═╡ c40fa7ad-622f-402f-8444-2d26f79ac1f9
# a single selected mode plays back at unit amplitude, ignoring where the source sits (a
# pure resonance has nothing to do with excitation) -- a selected range of modes is instead
# weighted by how strongly a delta source at x_src actually excites each one.
_nm_up = if _nm_selLo == _nm_selHi
    mode_series(_nm_U[:, _nm_selLo:_nm_selLo], [_nm_omegas[_nm_selLo]], [1.0], _nm_tgrid)
else
    mode_series(_nm_U[:, _nm_selLo:_nm_selHi], _nm_omegas[_nm_selLo:_nm_selHi], _nm_U[_nm_isrc, _nm_selLo:_nm_selHi], _nm_tgrid)
end

# ╔═╡ e2d0f6a1-7a4c-4b9d-8f2e-3c6d5a1b8e97
_nm_record = _nm_up[:, _nm_irec]

# ╔═╡ f3e1a7b2-8b5d-4c0e-9f3f-4d7e6b2c9fa8
_nm_freqaxis, _nm_spectrum = rfftfreq(length(_nm_tgrid), inv(step(_nm_tgrid))), abs.(rfft(_nm_record))

# ╔═╡ a4f2b8c3-9c6e-4d1f-a04f-5e8f7c3daeb9
begin
    struct ModeExplorerPush
        xgrid::Any
        up::Any
        omegas::Any
        record::Any
        spectrum::Any
        freqaxis::Any
        nt::Int
        nx::Int
        tmax::Float64
        upmax::Float64
        recordmax::Float64
        specmax::Float64
    end
    function Base.show(io::IO, ::MIME"text/html", p::ModeExplorerPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('nmwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('nm-update', { detail: {
            xgrid: [$(_nm_flatten(p.xgrid))],
            up: [$(_nm_flatten(p.up))],
            omegas: [$(_nm_flatten(p.omegas))],
            record: [$(_nm_flatten(p.record))],
            spectrum: [$(_nm_flatten(p.spectrum))],
            freqaxis: [$(_nm_flatten(p.freqaxis))],
            nt: $(p.nt), nx: $(p.nx), tmax: $(p.tmax),
            upmax: $(p.upmax), recordmax: $(p.recordmax), specmax: $(p.specmax),
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ b8c4d2e6-1a7f-4e9c-b3d8-6f2a9c1e7db4
ModeExplorerPush(_nm_xgrid, vec(permutedims(_nm_up)), _nm_omegas, _nm_record, _nm_spectrum, _nm_freqaxis,
    length(_nm_tgrid), length(_nm_xgrid), _nm_tgrid[end],
    max(maximum(abs, _nm_up), 1.0e-12), max(maximum(abs, _nm_record), 1.0e-12), max(maximum(_nm_spectrum), 1.0e-12))

# ╔═╡ c9d5e3f7-2b8a-4f1d-a4e9-7d3b0f2c8ea5
md"""
## Additional Reading
- Weideman, J. A. C., and L. N. Trefethen. “The Eigenvalues of Second-Order Spectral Differentiation Matrices.” SIAM Journal on Numerical Analysis, vol. 25, no. 6, 1988, pp. 1279–98. [pdf](https://www.jstor.org/stable/2157485?seq=3#metadata_info_tab_contents)
- [Prof. John K Hunter; Notes from UC Davis](https://www.math.ucdavis.edu/~hunter/m280_09/ch4.pdf)
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
FFTW = "~1.10.0"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "a084269afdb10b773d44cec96b75aa461ae7f805"

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

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

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

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

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

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
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

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "5253f44481f18cd938d4559d5e44fa82198408a6"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

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
# ╠═46d9831d-2b46-4030-bb15-cdc732d72e68
# ╟─319925c0-1512-459b-8be2-1d2650678f35
# ╟─27ac2670-a5b2-431d-9626-c95241a00b55
# ╟─75617746-5345-4074-bbdf-7bd3e8505aae
# ╟─9a5e731a-4e70-4938-b0a9-b7fd796f4e59
# ╟─1d0767df-945c-45cf-9f59-9ba0e3accefe
# ╟─c17c4841-10ae-4e5c-84eb-762f7ce72af9
# ╟─235eff94-9dca-40e6-a37f-e906d04eb094
# ╟─04f765cb-826e-4939-800e-9a4952e895bb
# ╟─44fb45cd-448b-4927-856c-0bf7799bd40e
# ╠═ec847573-27e7-4845-9375-e92539d7aff6
# ╟─8de036a2-2116-4d76-975a-b87f37100a36
# ╠═e4a8ceca-2cd9-47b6-bca7-146b70c6ede4
# ╠═2731c578-e51d-4fb0-8ca7-3216d49f649a
# ╠═7056f891-2b79-4da5-95e9-b133fa86043f
# ╟─29ca2faf-d4e0-4e9d-ba07-525aba0b65c5
# ╠═c728bc06-abc9-4a60-97a4-598b721e78bc
# ╟─b6a4c9d1-70e1-4e2a-9c1d-3a5b9e0d2f11
# ╠═94568f52-b993-4f62-8170-776e83945438
# ╟─8f6fb806-7877-48d5-989b-9c9944e35ce5
# ╟─d7db01d4-f552-4c65-b122-ba392943f3b3
# ╟─ca1d876b-0b59-47c5-944c-4ff45e40ff6d
# ╟─dbf086b6-8b44-4f7e-ad99-7a31612c299a
# ╠═b5b41af3-2699-43d9-8f60-2cac6df2d040
# ╠═8eb5c988-5471-44df-b0a2-38196a3fe828
# ╠═4b695e82-0860-4193-ae03-8d317b7280c7
# ╟─d5670bb8-aab6-421b-83fa-bed20faf7c33
# ╠═1fa011bc-69bc-4e4a-aa00-b89672bbd376
# ╠═089009a9-41e7-406e-9eea-5faa31126c6e
# ╠═0c723ed7-79f6-4334-b575-baa15b6dc6b5
# ╠═10ada4dd-c11c-4605-8b68-0e550279234f
# ╠═476860e2-3ba3-4132-8848-09cd2e6ae52a
# ╠═ab6ad92c-2d8c-4cee-8643-539b71f9b462
# ╠═ead6d118-a3a8-4ca6-8cf5-7af2997d24f5
# ╠═6fe6ac9c-e7fa-4c75-a10b-68128fb786be
# ╠═9ae40035-5296-4cf5-a4fb-be2619ddc992
# ╠═c40fa7ad-622f-402f-8444-2d26f79ac1f9
# ╠═e2d0f6a1-7a4c-4b9d-8f2e-3c6d5a1b8e97
# ╠═f3e1a7b2-8b5d-4c0e-9f3f-4d7e6b2c9fa8
# ╠═a4f2b8c3-9c6e-4d1f-a04f-5e8f7c3daeb9
# ╠═b8c4d2e6-1a7f-4e9c-b3d8-6f2a9c1e7db4
# ╟─c9d5e3f7-2b8a-4f1d-a4e9-7d3b0f2c8ea5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
