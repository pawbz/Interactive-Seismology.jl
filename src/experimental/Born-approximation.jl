### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Born Approximation"
#> tags = ["imaging"]
#> layout = "layout.jlhtml"
#> description = "Paint subsurface scatterers, watch the linearized (Born) scattered wavefield propagate, and see whether migration recovers what you drew."

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

# ╔═╡ 218dc0cf-7d17-4e9b-8393-ceca9c22acaf
begin
    using FFTW
    using Bessels
    using Tullio
    using LinearAlgebra
    using PlutoUI
end

# ╔═╡ 44d58ae7-f626-4bf6-a8ac-9dae63f2d9d6
TableOfContents()

# ╔═╡ 2f429d1f-ecd1-4b97-b46d-b4714144e24f
md"""
# Born Approximation

A seismic wave crossing a heterogeneous Earth does two things at once: it propagates through the
smooth background structure, and it **scatters** off whatever local heterogeneities it encounters
-- a fracture, an ore body, a buried fault. The Born approximation is the standard way to make
that second effect tractable: treat the scattered field as a small, *linear* correction on top of
the background field, valid as long as the heterogeneity is weak enough that multiply-scattered
energy is negligible.

This notebook derives that linearization, then hands you a widget where you **paint scatterers
directly onto a subsurface model** and watch three things update live: the scattered wavefield
propagating away from what you painted, the data it produces at a receiver array, and -- the real
payoff -- the **migration image**, the cheapest possible linear inverse, which tries to recover
your painted scatterers from the recorded data alone.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)
Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 5a0cea57-1627-4c2a-841e-d3a7db5bdedc
md"""
## From the Wave Equation to a Scattering Integral

Work in the frequency domain with the scalar (acoustic) Helmholtz equation, written so the medium
parameter appears *linearly* -- the squared slowness `` m(\mathbf x) = 1/c(\mathbf x)^2 ``, not
the slowness or the velocity itself:

```math
L(u, m) \;\equiv\; \omega^2 m\,u + \nabla^2 u.
```

The **background** state has medium `` m_0(\mathbf x) `` and wavefield `` U(\mathbf x,\mathbf x_s) ``
excited by a point source `` F(\omega)\delta(\mathbf x-\mathbf x_s) ``:

```math
L(U, m_0) = F(\omega)\,\delta(\mathbf x - \mathbf x_s).
```

The **perturbed** state -- what you get once you paint scatterers -- adds a small `` \delta m(\mathbf x) ``
(what the widget's brush actually paints) and the wavefield responds with a correction
`` \delta U ``:

```math
L(U+\delta U,\ m_0+\delta m) = F(\omega)\,\delta(\mathbf x - \mathbf x_s).
```

Subtract the background equation from the perturbed one and expand:

```math
\underbrace{\omega^2 m_0\,\delta U + \nabla^2\delta U}_{L(\delta U,\, m_0)} \;+\; \omega^2\delta m\,(U+\delta U) \;=\; 0.
```

!!! danger "The Born approximation, in one line"
	The term `` \omega^2\delta m\,\delta U `` is second order -- a product of two already-small
	quantities. **Dropping it** is the entire Born approximation: it linearizes the scattering
	problem by keeping `` \delta U `` only to first order in `` \delta m ``, i.e. it accounts for a
	wave scattering *once* off the heterogeneity, never twice. What survives is a wave equation for
	`` \delta U `` sourced by the *background* field `` U `` scattering off `` \delta m ``:

```math
L(\delta U,\ m_0) \;=\; -\,\omega^2\delta m\, U. \qquad(\text{Born})
```
"""

# ╔═╡ 1ed696fe-3a60-43a7-8dc3-474b8f68c2b2
md"""
## The Scattering Integral and Its Adjoint

(Born) is a wave equation for `` \delta U `` with a source term `` -\omega^2\delta m\,U ``
spread over the whole scattering volume, not a point -- Green's-function machinery turns it into
an explicit integral. Let `` G(\mathbf x,\mathbf x_r) `` solve the *same* background operator for
an impulsive point source at the receiver, `` L(G,m_0)=\delta(\mathbf x-\mathbf x_r) ``. Combining
it with (Born) via Green's second identity, and using that a volume scatterer's outgoing waves
satisfy the Sommerfeld radiation condition (so the usual surface-integral term vanishes), gives the
**Born scattering integral equation**:

```math
\delta U(\mathbf x_r,\mathbf x_s) \;=\; -\,\omega^2\!\int_V G(\mathbf x,\mathbf x_r)\,\delta m(\mathbf x)\,U(\mathbf x,\mathbf x_s)\;dV(\mathbf x).
```

Since the background field is itself just the source convolved with the background Green's
function, `` U(\mathbf x,\mathbf x_s)=G(\mathbf x,\mathbf x_s)F(\omega) ``, this becomes

```math
\delta U(\mathbf x_r,\mathbf x_s) \;=\; -\,\omega^2 F(\omega)\!\int_V G(\mathbf x,\mathbf x_r)\,G(\mathbf x,\mathbf x_s)\,\delta m(\mathbf x)\;dV(\mathbf x). \qquad(\star)
```

Discretizing `` (\star) `` over the grid you paint on turns the integral into a matrix-vector
product `` \delta d = G\,\delta m `` -- exactly what [`get_forward_operator`](@ref) builds and
[`get_scattered_wavefield`](@ref) applies, one entry of `` G `` per (frequency, receiver, grid
cell) combination.

!!! tip "Migration: running the physics backwards"
	`` (\star) `` is *linear* in `` \delta m ``, so the cheapest possible estimate of `` \delta m ``
	from observed data `` \delta d `` is to apply the **adjoint** of `` G `` instead of solving a
	linear system: `` \delta m_{\text{image}} = G^{\!\top}\delta d ``. This is exactly
	[`get_migration_image`](@ref) -- literally the forward-modeling operator run backwards. It is
	not a true inverse (`` G^{\!\top}G\neq I `` in general), so the recovered image is smeared and
	amplitude-distorted compared to what you painted -- but it's the foundation every more
	elaborate seismic imaging algorithm builds on, and the widget's migration panel lets you see
	directly how much aperture (how many live sources and receivers) it takes to make that smearing
	tolerable.
"""

# ╔═╡ 689dd440-942d-449d-9b75-825ad3b052af
md"""
## Appendix
"""

# ╔═╡ e4de01e7-6d0a-4211-b83c-c292984c76f3
md"""
### Acquisition & Grid
"""

# ╔═╡ 7a4b1527-e09d-4393-bba1-69f7f7048e6e
begin
    const XMAX = 800.0  # m -- domain width
    const ZMAX = 400.0  # m -- domain depth
    const NX = 61        # imaging/painting grid -- also used for the migration image
    const NZ = 31
    const xgrid = range(0.0, XMAX; length=NX)
    const zgrid = range(0.0, ZMAX; length=NZ)
    const NXM = 31       # coarser grid, animated wavefield-movie panel only -- keeps the pushed
    const NZM = 16        # per-frame data small enough to stay snappy
    const xgridM = range(0.0, XMAX; length=NXM)
    const zgridM = range(0.0, ZMAX; length=NZM)
    const STANDOFF = 15.0 # m -- fixed elevation of sources/receivers above z=0, avoids the G0
    # singularity a source or receiver sitting exactly on top of a scatterer would cause
    const NSRC = 3
    const NREC = 31
    const recX = collect(range(0.0, XMAX; length=NREC))
    const TGRID = range(0.0, 1.0; length=201) # s
end

# ╔═╡ 65a857c3-d6ab-4f9a-8839-1f829c4621c4
md"""
### 2-D Green's Function
"""

# ╔═╡ 4ba10533-065c-44b0-8295-2363fb99a8ea
begin
    """
    	rad(sx, sz, rx, rz)

    Euclidean distance (m) between a point `(sx, sz)` and a point `(rx, rz)`.
    """
    rad(sx, sz, rx, rz) = sqrt(abs2(sx - rx) + abs2(sz - rz))

    """
    	G0(rx, rz, sx, sz, k, rho)

    The 2-D acoustic Green's function for a homogeneous medium: the frequency-domain response at
    `(rx, rz)` to a unit impulsive point source at `(sx, sz)`, wavenumber `k = omega/c` and density
    `rho`. Built from the zero-order Hankel function of the second kind -- an outgoing cylindrical
    wave under the `` e^{-i\\omega t} `` convention this notebook uses throughout.
    """
    function G0(rx, rz, sx, sz, k, rho)
        # H0(2) is genuinely singular at r=0 (Y0(x) -> -Inf as x -> 0), not just numerically
        # touchy -- this floor regularizes that self-term. It matters because the wavefield-movie
        # panel reuses the model's own grid points as its "receivers" (see resample_to_grid /
        # the acqM cell), which guarantees an exact-coincidence r=0 case on the grid diagonal; the
        # floor is well below the smallest real (non-coincident) grid spacing used anywhere in
        # this notebook, so it never affects a genuine source-receiver-scatterer distance.
        r = max(rad(sx, sz, rx, rz), 2.0)
        return -0.25 * rho * im * hankelh2(0, k * r)
    end
end

# ╔═╡ dff8cb60-4cbd-47bc-a803-9efbb71ee47b
md"""
### Source Spectrum
"""

# ╔═╡ db3a4cb2-a525-46b2-aaf5-9a7d9078119f
"""
	source_spectrum(tgrid, fpeak)

A Gaussian-shaped amplitude spectrum peaked at `fpeak` Hz, sampled on the real-FFT frequency grid
implied by `tgrid`. Returns `(freqgrid, Fsource)`; the DC bin is zeroed since a DC source carries
no travel-time information here.
"""
function source_spectrum(tgrid, fpeak)
    freqgrid = collect(rfftfreq(length(tgrid), inv(step(tgrid))))
    Fsource = exp.(-abs2.(freqgrid .- fpeak) .* 1.0e-3)
    Fsource[1] = 0.0
    return freqgrid, Fsource
end

# ╔═╡ e18b63b0-ae7b-4d53-81f4-f72eb00a09d2
"""
	acoustic_medium(vp0, fpeak, xg, zg; rho0=1000.0)

Bundle a homogeneous-background medium and its source spectrum into the named tuple every
forward/migration function below expects: grid coordinates `xgrid`/`zgrid`, the real-FFT frequency
grid, the matching wavenumber grid `kgrid = 2*pi*freqgrid/vp0`, the source amplitude spectrum
`Fsource`, density `rho0`, and `tgrid` itself (kept only for the final `irfft` length).
"""
function acoustic_medium(vp0, fpeak, xg, zg; rho0=1000.0)
    freqgrid, Fsource = source_spectrum(TGRID, fpeak)
    kgrid = 2 * pi * freqgrid ./ vp0
    return (; xgrid=xg, zgrid=zg, tgrid=TGRID, freqgrid, kgrid, Fsource, rho0)
end

# ╔═╡ 856ee675-dff3-4d61-9b13-efd83f85cfc3
md"""
### Forward Modeling
"""

# ╔═╡ fc2711fe-68d2-4432-9ed8-e511c7fa71fe
"""
	get_forward_operator(pa, acq, sloc_x, sloc_z)

Build the Born single-scattering sensitivity operator `G` for one source at `(sloc_x, sloc_z)`,
discretizing `` (\\star) `` above over `pa`'s spatial grid: `G[iω, ir, iz, ix]` is the response a
unit `δm` at grid cell `(iz, ix)` produces at receiver `ir`, frequency `iω`, i.e.
`` G_0(\\text{scatterer}\\!\\to\\!\\text{receiver})\\cdot G_0(\\text{source}\\!\\to\\!\\text{scatterer})
\\cdot\\omega^2 F(\\omega) ``. Flattened to `(nr*nω, nz*nx)` so scattering/migration below are plain
matrix-vector products, `` \\delta d = G\\,\\delta m `` / `` \\delta m_{\\text{image}} = G^{\\!\\top}\\delta d ``.
The DC frequency row is zeroed (the source spectrum itself carries no DC energy, and `k=0` would
be singular in `G0` anyway).
"""
function get_forward_operator(pa, acq, sloc_x, sloc_z)
    (; xgrid, zgrid, freqgrid, kgrid, Fsource, rho0) = pa
    # xgrid/zgrid are already concrete, directly-indexable ranges -- no need to `collect`
    # them into a Vector first, that would just be an extra allocation for no benefit.
    nx, nz = length(xgrid), length(zgrid)
    nr = length(acq.rlocs_x)
    nω = length(freqgrid)
    @tullio F[iω] := freqgrid[iω] * freqgrid[iω] * Fsource[iω] * 4.0 * pi * pi
    @tullio G[iω, ir, iz, ix] := G0(xgrid[ix], zgrid[iz], sloc_x, sloc_z, kgrid[iω], rho0) *
                                  G0(acq.rlocs_x[ir], acq.rlocs_z[ir], xgrid[ix], zgrid[iz], kgrid[iω], rho0) * F[iω]
    G[1, :, :, :] .= complex(0.0)
    return reshape(G, nr * nω, nx * nz)
end

# ╔═╡ 04207c33-ce75-4b16-a22d-df46cc949ef9
"""
	get_reference_wavefield(pa, acq, tgrid)

The reference (unperturbed) wavefield at every receiver in `acq`, for every source in `acq`, in
the time domain: builds `` D(\\omega) = G_0(\\text{source}\\!\\to\\!\\text{receiver})\\,F(\\omega) ``
then inverse-FFTs it. Returns an `(nt, nr, ns)` array.
"""
function get_reference_wavefield(pa, acq, tgrid)
    (; freqgrid, kgrid, Fsource, rho0) = pa
    @tullio D[iω, ir, is] := G0(acq.rlocs_x[ir], acq.rlocs_z[ir], acq.slocs_x[is], acq.slocs_z[is], kgrid[iω], rho0) * Fsource[iω]
    D[1, :, :] .= complex(0.0)
    return irfft(D, length(tgrid), 1)
end

# ╔═╡ a219febc-278d-4093-aba4-2c7dc736f6ab
"""
	get_scattered_wavefield(dm, G, acq, pa)

The Born-approximate scattered field at every receiver, in the time domain, for a squared-slowness
perturbation field `dm` (`(nz,nx)`, Julia's own natural column-major layout -- see
[`resample_to_grid`](@ref) for how the movie panel gets its own copy) and its precomputed
sensitivity operator `G` (from [`get_forward_operator`](@ref)): one matrix-vector product
`` \\delta d = G\\,\\delta m `` followed by an inverse FFT. Returns an `(nt, nr)` array.
"""
function get_scattered_wavefield(dm, G, acq, pa)
    nr = length(acq.rlocs_x)
    nω = length(pa.freqgrid)
    d = G * vec(dm)
    d = reshape(d, nω, nr)
    return irfft(d, length(pa.tgrid), 1)
end

# ╔═╡ c58bc7c7-e686-423f-94b2-5bb2ae531643
begin
    """
    	nearest_index(rng, value)

    The index into the uniform range `rng` closest to `value`, computed directly by arithmetic --
    `O(1)` and allocation-free, unlike `argmin(abs.(value .- rng))`, which scans the whole range
    and allocates a temporary array to do it. Valid because every grid this notebook resamples
    between is a uniform range (`range(...; length=...)`). Clamped to stay within `rng`'s valid
    indices.
    """
    nearest_index(rng, value) = clamp(round(Int, (value - first(rng)) / step(rng)) + 1, 1, length(rng))

    """
    	resample_to_grid(field_fine, zg_fine, xg_fine, zg_coarse, xg_coarse)

    Nearest-neighbour resample a `(length(zg_fine), length(xg_fine))` field onto a different
    rectangular grid -- used to carry the painted perturbation (defined on the fine imaging grid)
    onto the coarser grid the animated wavefield-movie panel uses. Runs once per paint stroke, so
    the allocation-free [`nearest_index`](@ref) lookup (rather than a linear `argmin` search) keeps
    painting responsive.
    """
    function resample_to_grid(field_fine, zg_fine, xg_fine, zg_coarse, xg_coarse)
        [field_fine[nearest_index(zg_fine, z), nearest_index(xg_fine, x)] for z in zg_coarse, x in xg_coarse]
    end
end

# ╔═╡ 90fcf549-8982-4bde-895d-eb9df3b32a41
md"""
### Migration (Imaging)
"""

# ╔═╡ cd5c240c-db82-4307-8d5b-702a0b49a946
"""
	get_migration_image(δd, G, acq, pa)

The migrated (adjoint-recovered) image from time-domain scattered data `δd` (`nt × nr`): forward-
FFT the data, then apply the adjoint of the same sensitivity operator `G` used to model it,
`` \\delta m_{\\text{image}} = G^{\\!\\top}\\delta d `` -- literally running the forward operator
backwards, per the "Migration" tip above. Returns the real part, reshaped to `(nz, nx)`.
"""
function get_migration_image(δd, G, acq, pa)
    (; xgrid, zgrid) = pa
    nx, nz = length(xgrid), length(zgrid)
    δdf = rfft(δd, 1)
    δm = G' * vec(δdf)
    return real.(reshape(δm, nz, nx))
end

# ╔═╡ 80f38939-748f-41cc-83ec-762dbfaa6267
md"""
### Verifying the Migration Operator
"""

# ╔═╡ c43e7a70-5895-4130-a8a2-988e5c9e6d9f
md"""
!!! correct "Checking that migration really is the adjoint of modeling"
	`get_migration_image` only makes physical sense as "running the physics backwards" if `G'` is
	genuinely the adjoint of `G` under the standard complex inner product
	`` \\langle u,v\\rangle=\\overline{u}^{\\!\\top}v ``, i.e. `` \\langle Gs,d\\rangle = \\langle s,G^{\\!\\top}d\\rangle ``
	for *any* `s`, `d` -- not just for the specific fields this notebook happens to compute. Checked
	below on a small, self-contained problem with random `s`/`d` (Julia's `'` on a complex matrix
	*is* the conjugate transpose, so this is exactly testing what `get_migration_image` uses).
"""

# ╔═╡ 51d129d0-48f2-444a-964b-dbbec58a0b38
let
    xg, zg = range(0.0, 100.0; length=4), range(0.0, 100.0; length=4)
    tg = range(0.0, 0.5; step=0.01)
    fg, Fs = source_spectrum(tg, 25.0)
    kg = 2 * pi * fg ./ 1500.0
    pa_test = (; xgrid=xg, zgrid=zg, tgrid=tg, freqgrid=fg, kgrid=kg, Fsource=Fs, rho0=1000.0)
    acq_test = (; rlocs_x=[10.0, 50.0, 90.0], rlocs_z=fill(-5.0, 3), slocs_x=[50.0], slocs_z=[-5.0])
    G = get_forward_operator(pa_test, acq_test, 50.0, -5.0)
    s = randn(ComplexF64, size(G, 2))
    d = randn(ComplexF64, size(G, 1))
    lhs = dot(G * s, d)   # ⟨Gs,d⟩
    rhs = dot(s, G' * d)  # ⟨s,Gᵀd⟩
    residual = abs(lhs - rhs) / max(abs(lhs), 1e-12)
    @assert residual < 1e-10
    "adjoint check: |⟨Gs,d⟩ - ⟨s,Gᵀd⟩| / |⟨Gs,d⟩| = $(round(residual, sigdigits=3)) ✓"
end

# ╔═╡ 08d196f8-86cb-4c0f-8050-811003182671
md"""
### Computing the Fields
"""

# ╔═╡ 38cd5bf1-203f-49ce-810b-7ed60e88c705
md"""
### The Interactive Widget
"""

# ╔═╡ 59b5207f-7af0-4508-bb22-465ac206a10d
md"""
### Pushing Results to the Widget
"""

# ╔═╡ 527bf8f5-1b91-48a1-b235-c5774c7b1c70
md"""
`BsPush` does no physics -- it takes the already-computed gather/image/movie arrays and hands them
to the *already-rendered* [`BornScatteringInput`](@ref) widget by dispatching a browser
`CustomEvent`, the same pattern this repo's other widgets use (see e.g. `PfrPush` in
`lame-theorem.jl`).
"""

# ╔═╡ 271bca15-eab3-4263-b2c3-deaedbfe12a8
begin
    _bs_flatten(v) = join(v, ",")

    """
    	flatten_grid_rowmajor(M)

    Flatten a `(NZ,NX)`-style Julia array (natural column-major, z fastest) into the row-major
    `iz*NX+ix` order the widget's JS expects, comma-joined and ready to interpolate into a pushed
    `CustomEvent`.
    """
    flatten_grid_rowmajor(M) = join(vec(permutedims(M)), ",")

    """
    	flatten_movie_rowmajor(M, NZg, NXg)

    Flatten an `(nt, NZg*NXg)` time series (Julia-natural columns, from
    [`get_reference_wavefield`](@ref)/[`get_scattered_wavefield`](@ref)) into one long row-major
    string, frame by frame, via [`flatten_grid_rowmajor`](@ref) -- what the widget's animated
    wavefield panel indexes into as `frame*NZg*NXg + iz*NXg+ix`.
    """
    function flatten_movie_rowmajor(M, NZg, NXg)
        nt = size(M, 1)
        join([flatten_grid_rowmajor(reshape(view(M, it, :), NZg, NXg)) for it in 1:nt], ",")
    end
end

# ╔═╡ fb798d80-3f07-4294-a8af-7fc7c34a61f7
begin
    struct BornScatteringInput
        vp0::Float64            # background velocity, m/s
        fpeak::Float64          # source peak frequency, Hz
        pert::Vector{Float64}   # flat NZ*NX painted perturbation, % of background m0=1/vp0^2, JS row-major (iz*NX+ix)
        srcX::Vector{Float64}   # NSRC source x positions, m
        viewSrc::Int             # 1..NSRC -- which source's gather/movie panels are shown
        recMuted::Vector{Float64} # NREC flags, 0=active 1=muted
    end
    BornScatteringInput(; vp0=2500.0, fpeak=40.0, pert=zeros(NZ * NX), srcX=[150.0, 400.0, 650.0],
        viewSrc=2, recMuted=zeros(NREC)) =
        BornScatteringInput(vp0, fpeak, pert, srcX, viewSrc, recMuted)

    Base.get(w::BornScatteringInput) = Dict{String,Any}(
        "vp0" => w.vp0, "fpeak" => w.fpeak, "pert" => w.pert,
        "srcX" => w.srcX, "viewSrc" => w.viewSrc, "recMuted" => w.recMuted)

    function Base.show(io::IO, ::MIME"text/html", w::BornScatteringInput)
        write(io, """
        <div id="bswidget">
        <style>
        pluto-cell:has(#bswidget) { width: min(80vw, 1500px) !important;
          margin-left: calc((100% - min(80vw, 1500px)) / 2) !important; }
        #bswidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #bswidget .bs-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #bswidget .bs-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #bswidget .bs-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #bswidget .bs-workspace{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start;margin-bottom:14px}
        #bswidget .bs-row2{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start}
        #bswidget .bs-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #bswidget .bs-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:4px}
        #bswidget canvas{display:block;cursor:crosshair}
        #bswidget .bs-controls{flex:0 0 540px;width:540px;display:grid;
          grid-template-columns:repeat(2, minmax(0,1fr));gap:8px;font:14px sans-serif;align-content:start}
        #bswidget .bs-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #bswidget .bs-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #bswidget .bs-control-row{display:grid;grid-template-columns:70px minmax(0,1fr) 58px;gap:6px;align-items:center;margin:5px 0}
        #bswidget .bs-control-row label{font-size:13px;color:#9ca3af}
        #bswidget .bs-control-row input[type=range]{width:100%;min-width:0}
        #bswidget .bs-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #bswidget .bs-actions{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
        #bswidget .bs-readout{font-size:13px;color:#d1d5db;line-height:1.6}
        #bswidget .bs-readout b{color:#e5e7eb}
        #bswidget select{background:#0b0b0b;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:3px;width:100%}
        #bswidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:13px;cursor:pointer}
        #bswidget button.active{background:#2563eb;border-color:#93c5fd}
        #bswidget .bs-cbar-label{font-size:11px;color:#9ca3af;margin:8px 0 2px;text-align:center}
        #bswidget .bs-cbar{cursor:default;display:block;margin:0 auto}
        </style>

        <div class="bs-title">
          <div class="bs-title-desc">Paint subsurface scatterers and see the scattered wavefield, the data it produces, and whether migration recovers it.</div>
          <div class="bs-title-hint">drag on the model to paint &middot; drag a star to move a source (click to view it) &middot; click a triangle to mute a receiver &middot; press Play to animate</div>
        </div>

        <div class="bs-workspace">
          <div>
            <div class="bs-panel"><canvas id="bs-model"></canvas></div>
            <div class="bs-caption" id="bs-model-caption"></div>
            <div class="bs-cbar-label">scattering strength, % of background</div>
            <canvas id="bs-cbar" class="bs-cbar"></canvas>
          </div>
          <div class="bs-controls">
            <div class="bs-control-group">
              <div class="bs-control-title">Paint</div>
              <div class="bs-actions">
                <button id="bs-paint-neg" type="button">Slower (red)</button>
                <button id="bs-paint-pos" type="button">Faster (blue)</button>
              </div>
              <div class="bs-control-row"><label>brush</label><input type="range" id="bs-brush" min="20" max="150" step="10" value="60"><span class="bs-value" id="bs-brush-v"></span></div>
              <div class="bs-actions"><button id="bs-clear" type="button">Clear paint</button></div>
            </div>
            <div class="bs-control-group">
              <div class="bs-control-title">Background</div>
              <div class="bs-control-row"><label>v&#8320; (m/s)</label><input type="range" id="bs-vp0" min="1500" max="4000" step="50" value="$(w.vp0)"><span class="bs-value" id="bs-vp0-v"></span></div>
              <div class="bs-control-row"><label>f peak (Hz)</label><input type="range" id="bs-fpeak" min="15" max="80" step="1" value="$(w.fpeak)"><span class="bs-value" id="bs-fpeak-v"></span></div>
            </div>
            <div class="bs-control-group">
              <div class="bs-control-title">Wavefield movie</div>
              <div class="bs-actions"><button id="bs-play" type="button">Play</button><button id="bs-reset" type="button">Reset</button></div>
              <div class="bs-actions" style="margin-top:6px">
                <button id="bs-view-ref" type="button">Reference</button>
                <button id="bs-view-tot" type="button">Total</button>
                <button id="bs-view-sca" type="button">Scattered</button>
              </div>
            </div>
            <div class="bs-control-group">
              <div class="bs-control-title">Readouts</div>
              <div class="bs-readout" id="bs-readout"></div>
            </div>
          </div>
        </div>

        <div class="bs-row2">
          <div>
            <div class="bs-panel"><canvas id="bs-gather"></canvas></div>
            <div class="bs-caption">shot gather, source #<span id="bs-gather-src"></span> &middot; reference (left) vs. scattered (right)</div>
          </div>
          <div>
            <div class="bs-panel"><canvas id="bs-image"></canvas></div>
            <div class="bs-caption">migration image (all $(NSRC) sources stacked)</div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        const NX=$(NX), NZ=$(NZ), NXM=$(NXM), NZM=$(NZM), XMAX=$(XMAX), ZMAX=$(ZMAX);
        const NSRC=$(NSRC), NREC=$(NREC), NT=$(length(TGRID)), TMAX=$(TGRID[end]);
        const recX = [$(_bs_flatten(recX))];

        let vp0=$(w.vp0), fpeak=$(w.fpeak);
        let pert = $(w.pert == zeros(NZ*NX) ? "new Array(NX*NZ).fill(0)" : "[" * join(w.pert, ",") * "]");
        let srcX = [$(_bs_flatten(w.srcX))];
        let viewSrc = $(w.viewSrc);
        let recMuted = [$(_bs_flatten(w.recMuted))];

        const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1500);
        const CONTROLS_W = 540+16; // matches .bs-controls' 2-column grid width -- keep in sync
        const SEC_W = Math.max(360, Math.min(availW - CONTROLS_W, 620));
        const SEC_H = Math.round(SEC_W * ZMAX/XMAX);
        const DPR = window.devicePixelRatio || 1;

        function hidpi(canvas, context, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          context.setTransform(DPR,0,0,DPR,0,0);
        }

        const cv = par.querySelector('#bs-model');
        const ctx = cv.getContext('2d');
        hidpi(cv, ctx, SEC_W, SEC_H);

        const cbarCv = par.querySelector('#bs-cbar');
        const cbarCtx = cbarCv.getContext('2d');
        const CBAR_W = SEC_W, CBAR_H = 34;
        hidpi(cbarCv, cbarCtx, CBAR_W, CBAR_H);

        const GATHER_W = SEC_W, GATHER_H = Math.round(SEC_W * 0.62);
        const gCv = par.querySelector('#bs-gather');
        const gCtx = gCv.getContext('2d');
        hidpi(gCv, gCtx, GATHER_W, GATHER_H);

        const IMG_W = availW - SEC_W - 16, IMG_H = Math.round(IMG_W * ZMAX/XMAX);
        const iCv = par.querySelector('#bs-image');
        const iCtx = iCv.getContext('2d');
        hidpi(iCv, iCtx, Math.max(260,IMG_W), Math.max(130,Math.round(Math.max(260,IMG_W)*ZMAX/XMAX)));

        function toScreen(x, z){ return [x/XMAX*SEC_W, z/ZMAX*SEC_H]; }
        function toWorld(px, pz){ return [px/SEC_W*XMAX, pz/SEC_H*ZMAX]; }
        function xOf(ix,n){ return ix/(n-1)*XMAX; }
        function zOf(iz,n){ return iz/(n-1)*ZMAX; }

        // diverging colormap, reused for the paint field / gather / migration image alike --
        // blue=faster/positive, red=slower/negative (this repo's usual convention).
        function velColor(v, mx){
          const t = Math.max(-1, Math.min(1, v/mx));
          if(t >= 0) return [Math.round(255*(1-t)), Math.round(255*(1-t)), 255];
          const s = -t;
          return [255, Math.round(255*(1-s)), Math.round(255*(1-s))];
        }

        const FIELD_PCT_MAX = 10; // fixed colorbar scale: painted perturbation saturates at +-10%

        function drawFieldColorbar(){
          const w = CBAR_W, h = CBAR_H;
          cbarCtx.clearRect(0,0,w,h);
          const barH = 14, barY = 2;
          for(let i=0;i<w;i++){
            const t = i/(w-1);
            const v = -FIELD_PCT_MAX + t*2*FIELD_PCT_MAX;
            const [r,g,b] = velColor(v, FIELD_PCT_MAX);
            cbarCtx.fillStyle = 'rgb('+r+','+g+','+b+')';
            cbarCtx.fillRect(i, barY, 1, barH);
          }
          cbarCtx.strokeStyle = '#4b5563'; cbarCtx.lineWidth = 1;
          cbarCtx.strokeRect(0.5, barY+0.5, w-1, barH-1);
          cbarCtx.fillStyle = '#9ca3af'; cbarCtx.font = '10px sans-serif';
          cbarCtx.textAlign = 'left'; cbarCtx.fillText('-'+FIELD_PCT_MAX+'%', 0, barY+barH+12);
          cbarCtx.textAlign = 'right'; cbarCtx.fillText('+'+FIELD_PCT_MAX+'%', w, barY+barH+12);
          cbarCtx.textAlign = 'center'; cbarCtx.fillText('0', w/2, barY+barH+12);
        }

        let pushed = null; // {dref,dscat,image,drefMovie,dscatMovie,dmax,imgmax,wavemax} from Julia
        let playing = false, rafId = null, tPhase = 0, lastTs = null;
        let paintMode = 'neg'; // 'neg' (red, +delta-m, slower) | 'pos' (blue, -delta-m, faster)
        let waveView = 'ref'; // 'ref' | 'tot' | 'sca'
        let dragMode = null, dragIdx = -1, painting = false;

        // ---- painting on the (coarse) imaging grid ----
        const PAINT_STEP_PCT = 0.6; // % of the +-10 cap added/removed per brush touch
        function paintAt(xkm, zkm, brushKm){
          const sign = paintMode==='pos' ? 1 : -1;
          for(let iz=0; iz<NZ; iz++){
            const dz = zOf(iz,NZ)-zkm;
            for(let ix=0; ix<NX; ix++){
              const dx = xOf(ix,NX)-xkm;
              const d2 = dx*dx+dz*dz;
              const r2 = brushKm*brushKm;
              if(d2 < r2*4){
                const falloff = Math.exp(-d2/(2*r2/4));
                const amp = falloff * PAINT_STEP_PCT * sign;
                const idx = iz*NX+ix;
                pert[idx] = Math.max(-FIELD_PCT_MAX, Math.min(FIELD_PCT_MAX, pert[idx] + amp));
              }
            }
          }
        }

        function drawField(){
          const img = ctx.createImageData(Math.round(SEC_W*DPR), Math.round(SEC_H*DPR));
          const wpx = img.width, hpx = img.height;
          for(let py=0; py<hpx; py++){
            const zkm = py/hpx*ZMAX;
            const iz = Math.min(NZ-1, Math.round(zkm/ZMAX*(NZ-1)));
            for(let px2=0; px2<wpx; px2++){
              const xkm = px2/wpx*XMAX;
              const ix = Math.min(NX-1, Math.round(xkm/XMAX*(NX-1)));
              const [r,g,b] = velColor(pert[iz*NX+ix], FIELD_PCT_MAX);
              const idx = (py*wpx+px2)*4;
              img.data[idx]=r; img.data[idx+1]=g; img.data[idx+2]=b; img.data[idx+3]=255;
            }
          }
          ctx.putImageData(img, 0, 0);
        }

        function drawWaveOverlay(){
          if(!pushed) return;
          const nFrame = Math.min(NT-1, Math.max(0, Math.round(tPhase/TMAX*(NT-1))));
          const base = nFrame*NZM*NXM;
          const ref = pushed.drefMovie, sca = pushed.dscatMovie;
          const mx = pushed.wavemax;
          const cw = SEC_W/NXM, ch = SEC_H/NZM;
          ctx.save(); ctx.globalAlpha = 0.62;
          for(let iz=0; iz<NZM; iz++){
            for(let ix=0; ix<NXM; ix++){
              const k = base + iz*NXM+ix;
              let v = ref[k];
              if(waveView==='sca') v = sca[k];
              else if(waveView==='tot') v = ref[k] + sca[k];
              if(Math.abs(v) < mx*0.02) continue;
              const [r,g,b] = velColor(v, mx);
              ctx.fillStyle = 'rgb('+r+','+g+','+b+')';
              ctx.fillRect(ix*cw, iz*ch, cw+0.5, ch+0.5);
            }
          }
          ctx.restore();
        }

        function drawTriangleDownMarker(cx, cy, r, fill, stroke){
          ctx.beginPath();
          for(let i=0; i<3; i++){
            const ang = Math.PI/2 + i*2*Math.PI/3;
            const x = cx + r*Math.cos(ang), y = cy + r*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1.5; ctx.stroke();
        }

        function drawStarMarker(cx, cy, r, fill, stroke){
          const spikes = 5, rOuter = r, rInner = r * 0.45;
          ctx.beginPath();
          for(let i=0; i<spikes*2; i++){
            const rad = i % 2 === 0 ? rOuter : rInner;
            const ang = -Math.PI/2 + i*Math.PI/spikes;
            const x = cx + rad*Math.cos(ang), y = cy + rad*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1; ctx.stroke();
        }

        function drawGeom(cvx){
          const zStand = -30; // px above z=0 to draw the surface markers, purely cosmetic
          for(let i=0;i<NREC;i++){
            const [px] = toScreen(recX[i], 0);
            const muted = recMuted[i]===1;
            drawTriangleDownMarker(px, 6, 6, muted ? '#4b5563' : '#f5f3ef', '#0a0f18');
          }
          for(let i=0;i<NSRC;i++){
            const [px] = toScreen(srcX[i], 0);
            const active = (i+1)===viewSrc;
            drawStarMarker(px, 18, active?11:8, active?'#facc15':'#f5f3ef', '#0a0f18');
          }
        }

        function draw(){
          drawField();
          drawWaveOverlay();
          drawGeom();
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1; ctx.strokeRect(0.5,0.5,SEC_W-1,SEC_H-1);
          par.querySelector('#bs-model-caption').textContent = 'x: 0–' + XMAX.toFixed(0) + ' m  ·  z (depth): 0–' + ZMAX.toFixed(0) + ' m';
          drawFieldColorbar();
          drawGather();
          drawImage();
          updateReadout();
        }

        function drawGather(){
          gCtx.clearRect(0,0,GATHER_W,GATHER_H);
          gCtx.strokeStyle = '#374151'; gCtx.lineWidth = 1; gCtx.strokeRect(0.5,0.5,GATHER_W-1,GATHER_H-1);
          par.querySelector('#bs-gather-src').textContent = viewSrc;
          if(!pushed){
            gCtx.fillStyle = '#6b7280'; gCtx.font = '12px sans-serif'; gCtx.fillText('computing...', 10, 18);
            return;
          }
          const halfW = Math.floor(GATHER_W/2);
          const mx = pushed.dmax;
          function panel(data, x0, w){
            const img = gCtx.createImageData(Math.round(w*DPR), Math.round(GATHER_H*DPR));
            const wpx = img.width, hpx = img.height;
            for(let py=0; py<hpx; py++){
              const it = Math.min(NT-1, Math.round(py/hpx*(NT-1)));
              for(let px2=0; px2<wpx; px2++){
                const ir = Math.min(NREC-1, Math.round(px2/wpx*(NREC-1)));
                const v = data[it + ir*NT];
                const [r,g,b] = velColor(v, mx);
                const idx = (py*wpx+px2)*4;
                img.data[idx]=r; img.data[idx+1]=g; img.data[idx+2]=b; img.data[idx+3]=255;
              }
            }
            gCtx.putImageData(img, Math.round(x0*DPR), 0);
          }
          panel(pushed.dref, 0, halfW);
          panel(pushed.dscat, halfW, GATHER_W-halfW);
          gCtx.strokeStyle = '#4b5563'; gCtx.beginPath(); gCtx.moveTo(halfW+0.5,0); gCtx.lineTo(halfW+0.5,GATHER_H); gCtx.stroke();
          gCtx.fillStyle = '#e5e7eb'; gCtx.font = '11px sans-serif'; gCtx.textAlign='left';
          gCtx.fillText('reference', 4, 12);
          gCtx.fillText('scattered', halfW+4, 12);
        }

        function drawImage(){
          const W = iCv.width/DPR, H = iCv.height/DPR;
          iCtx.clearRect(0,0,W,H);
          iCtx.strokeStyle = '#374151'; iCtx.lineWidth = 1; iCtx.strokeRect(0.5,0.5,W-1,H-1);
          if(!pushed){
            iCtx.fillStyle = '#6b7280'; iCtx.font = '12px sans-serif'; iCtx.fillText('computing...', 10, 18);
            return;
          }
          const mx = pushed.imgmax;
          const img = iCtx.createImageData(Math.round(W*DPR), Math.round(H*DPR));
          const wpx = img.width, hpx = img.height;
          for(let py=0; py<hpx; py++){
            const iz = Math.min(NZ-1, Math.round(py/hpx*(NZ-1)));
            for(let px2=0; px2<wpx; px2++){
              const ix = Math.min(NX-1, Math.round(px2/wpx*(NX-1)));
              const v = pushed.image[iz*NX+ix];
              const [r,g,b] = velColor(v, mx);
              const idx = (py*wpx+px2)*4;
              img.data[idx]=r; img.data[idx+1]=g; img.data[idx+2]=b; img.data[idx+3]=255;
            }
          }
          iCtx.putImageData(img, 0, 0);
          function toImgScreen(x,z){ return [x/XMAX*W, z/ZMAX*H]; }
          for(let i=0;i<NREC;i++){
            const [px] = toImgScreen(recX[i], 0);
            const muted = recMuted[i]===1;
            iCtx.beginPath();
            const r=4;
            for(let k=0;k<3;k++){ const ang=Math.PI/2+k*2*Math.PI/3; const x=px+r*Math.cos(ang), y=6+r*Math.sin(ang); k===0?iCtx.moveTo(x,y):iCtx.lineTo(x,y); }
            iCtx.closePath(); iCtx.fillStyle = muted ? '#4b5563' : '#f5f3ef'; iCtx.fill(); iCtx.strokeStyle='#0a0f18'; iCtx.lineWidth=1; iCtx.stroke();
          }
          for(let i=0;i<NSRC;i++){
            const [px] = toImgScreen(srcX[i], 0);
            iCtx.beginPath(); iCtx.arc(px,14,3,0,2*Math.PI); iCtx.fillStyle='#facc15'; iCtx.fill();
          }
        }

        function updateReadout(){
          const nmuted = recMuted.reduce((a,b)=>a+(b===1?1:0),0);
          par.querySelector('#bs-readout').innerHTML =
            'v&#8320; <b>' + vp0.toFixed(0) + '</b> m/s &middot; f<sub>peak</sub> <b>' + fpeak.toFixed(0) + '</b> Hz<br>' +
            '<b>' + (NREC-nmuted) + '/' + NREC + '</b> receivers active<br>' +
            'viewing source <b>#' + viewSrc + '</b><br>' +
            (pushed ? 'movie t <b>' + tPhase.toFixed(2) + '</b> s' : 'computing...');
        }

        function syncControls(){
          par.querySelector('#bs-vp0').value = vp0; par.querySelector('#bs-vp0-v').textContent = vp0.toFixed(0);
          par.querySelector('#bs-fpeak').value = fpeak; par.querySelector('#bs-fpeak-v').textContent = fpeak.toFixed(0);
          par.querySelector('#bs-paint-neg').classList.toggle('active', paintMode==='neg');
          par.querySelector('#bs-paint-pos').classList.toggle('active', paintMode==='pos');
          par.querySelector('#bs-view-ref').classList.toggle('active', waveView==='ref');
          par.querySelector('#bs-view-tot').classList.toggle('active', waveView==='tot');
          par.querySelector('#bs-view-sca').classList.toggle('active', waveView==='sca');
        }

        let commitInFlight = false;
        function commit(){
          commitInFlight = true;
          par.value = { vp0, fpeak, pert, srcX, viewSrc, recMuted };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledCommit(){ if(!commitInFlight) commit(); }

        par.addEventListener('bs-update', e=>{
          pushed = e.detail;
          commitInFlight = false;
          draw();
        });

        cv.addEventListener('mousedown', e=>{
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, pz = e.clientY-rect.top;
          let best=null, bestD=14, bestIdx=-1;
          for(let i=0;i<NSRC;i++){
            const [spx,spz] = toScreen(srcX[i], 0);
            const d = Math.hypot(px-spx, pz-18);
            if(d<bestD){bestD=d; best='src'; bestIdx=i;}
          }
          if(best){ dragMode = best; dragIdx = bestIdx; viewSrc = bestIdx+1; draw(); return; }
          for(let i=0;i<NREC;i++){
            const [rpx] = toScreen(recX[i], 0);
            const d = Math.hypot(px-rpx, pz-6);
            if(d<10){ recMuted[i] = recMuted[i]===1 ? 0 : 1; draw(); commit(); return; }
          }
          painting = true;
          const [xkm,zkm] = toWorld(px,pz);
          const brushKm = +par.querySelector('#bs-brush').value;
          paintAt(xkm, zkm, brushKm);
          draw();
        });
        window.addEventListener('mousemove', e=>{
          const rect = cv.getBoundingClientRect();
          const px = e.clientX-rect.left, pz = e.clientY-rect.top;
          if(dragMode==='src'){
            const [xkm] = toWorld(Math.max(0,Math.min(SEC_W,px)), 0);
            srcX[dragIdx] = Math.max(0, Math.min(XMAX, xkm));
            draw(); throttledCommit();
            return;
          }
          if(!painting) return;
          const [xkm,zkm] = toWorld(px,pz);
          const brushKm = +par.querySelector('#bs-brush').value;
          paintAt(xkm, zkm, brushKm);
          draw();
        });
        window.addEventListener('mouseup', ()=>{
          if(painting || dragMode) commit();
          painting = false; dragMode = null; dragIdx = -1;
        });

        par.querySelector('#bs-paint-neg').addEventListener('click', ()=>{ paintMode='neg'; syncControls(); });
        par.querySelector('#bs-paint-pos').addEventListener('click', ()=>{ paintMode='pos'; syncControls(); });
        par.querySelector('#bs-view-ref').addEventListener('click', ()=>{ waveView='ref'; syncControls(); draw(); });
        par.querySelector('#bs-view-tot').addEventListener('click', ()=>{ waveView='tot'; syncControls(); draw(); });
        par.querySelector('#bs-view-sca').addEventListener('click', ()=>{ waveView='sca'; syncControls(); draw(); });
        par.querySelector('#bs-clear').addEventListener('click', ()=>{
          pert = new Array(NX*NZ).fill(0); draw(); commit();
        });

        par.addEventListener('input', e=>{
          if(e.target===par) return;
          e.stopImmediatePropagation();
          const id = e.target.id, v = e.target.value;
          if(id==='bs-vp0'){ vp0=+v; par.querySelector('#bs-vp0-v').textContent=vp0.toFixed(0); }
          else if(id==='bs-fpeak'){ fpeak=+v; par.querySelector('#bs-fpeak-v').textContent=fpeak.toFixed(0); }
          else if(id==='bs-brush'){ par.querySelector('#bs-brush-v').textContent=v+' m'; return; }
          else return;
          updateReadout(); throttledCommit();
        }, true);

        par.addEventListener('change', e=>{
          if(e.target===par) return;
          e.stopImmediatePropagation();
          const id = e.target.id;
          if(id==='bs-vp0'||id==='bs-fpeak'){ commit(); return; }
        }, true);

        const SIM_SPEED = 0.35;
        const playBtn = par.querySelector('#bs-play');
        function stepAnim(ts){
          if(lastTs===null) lastTs = ts;
          const dt = Math.min(0.1, (ts-lastTs)/1000);
          lastTs = ts;
          tPhase += dt*SIM_SPEED;
          if(tPhase > TMAX) tPhase = 0;
          draw();
          rafId = requestAnimationFrame(stepAnim);
        }
        playBtn.addEventListener('click', ()=>{
          playing = !playing;
          playBtn.textContent = playing ? 'Pause' : 'Play';
          if(playing){ lastTs=null; rafId = requestAnimationFrame(stepAnim); }
          else if(rafId){ cancelAnimationFrame(rafId); rafId=null; }
        });

        par.querySelector('#bs-reset').addEventListener('click', ()=>{
          vp0=$(w.vp0); fpeak=$(w.fpeak);
          pert = new Array(NX*NZ).fill(0);
          srcX = [$(_bs_flatten(w.srcX))]; viewSrc = $(w.viewSrc);
          recMuted = [$(_bs_flatten(w.recMuted))];
          tPhase = 0;
          syncControls(); draw(); commit();
        });

        par.querySelector('#bs-brush-v').textContent = par.querySelector('#bs-brush').value+' m';
        syncControls(); draw();
        }
        </script>
        """)
    end

    const _bs_ready = true
end

# ╔═╡ 24fab919-51b8-4ca3-ac6c-cdf2c4dbb407
begin
    _bs_ready
    WideCell(@bind bs BornScatteringInput(); max_width=1500)
end

# ╔═╡ c70710d4-3918-4397-aaba-e3ce9c87f96a
# The bond starts as `nothing` until the widget's first real interaction in a live browser reports
# back -- fall back to the same defaults the widget itself opens with.
bs_safe = bs isa AbstractDict ? bs : Dict{String,Any}(
    "vp0" => 2500.0, "fpeak" => 40.0, "pert" => zeros(NZ * NX),
    "srcX" => [150.0, 400.0, 650.0], "viewSrc" => 2, "recMuted" => zeros(NREC))

# ╔═╡ b313a158-3384-490f-8cdf-6786cd3facd4
let
    nmuted = count(==(1.0), Float64.(bs_safe["recMuted"]))
    md"""
    Background velocity **$(round(Int, bs_safe["vp0"]))** m/s, source peak frequency
    **$(round(Int, bs_safe["fpeak"]))** Hz &middot; **$(NSRC)** sources, **$(NREC - nmuted)/$(NREC)**
    receivers active &middot; viewing source **#$(bs_safe["viewSrc"])**'s gather and wavefield movie
    below (the migration image always stacks *all* $(NSRC) sources).
    """
end

# ╔═╡ 27ef0a2c-50b9-4d13-b50e-68786de4db91
# background medium + source spectrum -- depends only on vp0/fpeak, not on what's painted, so this
# (and every G-matrix cell below) only recomputes when the user actually moves those sliders.
_bs_pa = acoustic_medium(bs_safe["vp0"], bs_safe["fpeak"], xgrid, zgrid)

# ╔═╡ a3b1c2d4-5e6f-4708-9a1b-2c3d4e5f6071
# a JS bond's whole-number values can deserialize as either Int64 or Float64 depending on how they
# were last written (see this repo's pluto_bond_value_type memory note) -- viewSrc is used as an
# array index below, and Julia refuses to index with a Float64 even when it's exactly integer-
# valued, so round+clamp it once here rather than trusting the raw bond value at each use site.
_bs_viewSrc = clamp(round(Int, bs_safe["viewSrc"]), 1, NSRC)

# ╔═╡ c006bf6d-1ca6-4349-a7f0-e5e81c1fa67b
_bs_paM = acoustic_medium(bs_safe["vp0"], bs_safe["fpeak"], xgridM, zgridM)

# ╔═╡ 649a946a-8459-4bce-9ecd-aa6f771b8e6c
_bs_acq = (; rlocs_x=recX, rlocs_z=fill(-STANDOFF, NREC),
    slocs_x=Float64.(bs_safe["srcX"]), slocs_z=fill(-STANDOFF, NSRC))

# ╔═╡ ca1a7fed-2c65-4360-9e19-dfaf798bea1e
# reference field never depends on the painted scatterers -- computed once for all sources at once
_bs_dref = get_reference_wavefield(_bs_pa, _bs_acq, TGRID)

# ╔═╡ 082041a0-4a7a-41e4-aa28-9db8c421d463
begin
    _bs_gridZ = [z for z in zgridM, x in xgridM]  # (NZM,NXM), Julia-natural (z fastest when vec'd)
    _bs_gridX = [x for z in zgridM, x in xgridM]
    _bs_acqM = (; rlocs_x=vec(_bs_gridX), rlocs_z=vec(_bs_gridZ),
        slocs_x=[Float64(bs_safe["srcX"][_bs_viewSrc])], slocs_z=[-STANDOFF])
end

# ╔═╡ ba46a17b-af4d-423e-8398-9c344d0df8bf
_bs_Gmovie = get_forward_operator(_bs_paM, _bs_acqM, _bs_acqM.slocs_x[1], -STANDOFF)

# ╔═╡ 604ef000-2a93-44bf-8373-e36a4294deb0
_bs_drefM = get_reference_wavefield(_bs_paM, _bs_acqM, TGRID)[:, :, 1]

# ╔═╡ dca30e9e-6501-44df-a086-ee2780fb9f89
# one sensitivity matrix per source (real acquisition, full imaging grid) -- the expensive step,
# gated behind background/geometry only, exactly as the "Computing the Fields" architecture intends
_bs_Greal = [get_forward_operator(_bs_pa, _bs_acq, bs_safe["srcX"][i], -STANDOFF) for i in 1:NSRC]

# ╔═╡ 022d3632-8aab-4b44-9867-a8a639bcdfa6
# painted perturbation, JS row-major (iz*NX+ix) -> Julia-natural (NZ,NX), then to physical squared-
# slowness units -- this is the one cell that changes on every brush stroke, and everything above it
# (the G matrices) is untouched by it, so repainting only re-triggers the cheap cells below.
_bs_dm = permutedims(reshape(Float64.(bs_safe["pert"]), NX, NZ)) ./ 100 ./ bs_safe["vp0"]^2

# ╔═╡ 6679f50a-7af5-455b-abf5-2f166248950b
_bs_dmM = resample_to_grid(_bs_dm, zgrid, xgrid, zgridM, xgridM)

# ╔═╡ 8d3124af-39b2-4400-9fc7-246ee7180305
_bs_dscatM = get_scattered_wavefield(_bs_dmM, _bs_Gmovie, _bs_acqM, _bs_paM)

# ╔═╡ 2d04488e-29bd-409f-8a19-5e84efc6b9b6
_bs_muted = findall(==(1.0), Float64.(bs_safe["recMuted"]))

# ╔═╡ a144eb61-7b1b-4ec5-9037-ea13bf6cd140
_bs_dscat = map(1:NSRC) do i
    d = get_scattered_wavefield(_bs_dm, _bs_Greal[i], _bs_acq, _bs_pa)
    d[:, _bs_muted] .= 0.0
    d
end

# ╔═╡ 34fb2cf5-8b2a-43c9-9bdf-c52d0583a909
# migration always stacks ALL sources, muted receivers included as zeros exactly like a real survey
# with dead channels -- this is the panel that answers "did the image recover what I painted?"
_bs_image = sum(get_migration_image(_bs_dscat[i], _bs_Greal[i], _bs_acq, _bs_pa) for i in 1:NSRC)

# ╔═╡ 5f7b9237-8a45-4df7-a92e-eae3ade48ecd
begin
    struct BsPush
        dref::Any
        dscat::Any
        image::Any
        drefMovie::Any
        dscatMovie::Any
        dmax::Float64
        imgmax::Float64
        wavemax::Float64
    end
    function Base.show(io::IO, ::MIME"text/html", p::BsPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('bswidget');
        if(w){
          w.dispatchEvent(new CustomEvent('bs-update', { detail: {
            dref: [$(p.dref)],
            dscat: [$(p.dscat)],
            image: [$(p.image)],
            drefMovie: [$(p.drefMovie)],
            dscatMovie: [$(p.dscatMovie)],
            dmax: $(p.dmax),
            imgmax: $(p.imgmax),
            wavemax: $(p.wavemax),
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 7bfd9cbb-c4d1-412a-8c13-af5db5bcca06
begin
    _bs_dref_view = _bs_dref[:, :, _bs_viewSrc]     # (nt, NREC)
    _bs_dscat_view = _bs_dscat[_bs_viewSrc]          # (nt, NREC)
    _bs_dmax = max(maximum(abs, _bs_dref_view), maximum(abs, _bs_dscat_view), 1.0e-12)
    _bs_imgmax = max(maximum(abs, _bs_image), 1.0e-12)
    _bs_wavemax = max(maximum(abs, _bs_drefM), maximum(abs, _bs_dscatM), 1.0e-12)
end

# ╔═╡ 60c5d139-4df7-45ea-a8a9-c11e13cce883
BsPush(
    _bs_flatten(vec(_bs_dref_view)), _bs_flatten(vec(_bs_dscat_view)),
    flatten_grid_rowmajor(_bs_image),
    flatten_movie_rowmajor(_bs_drefM, NZM, NXM), flatten_movie_rowmajor(_bs_dscatM, NZM, NXM),
    _bs_dmax, _bs_imgmax, _bs_wavemax)

# ╔═╡ 978a06e3-f582-4070-ae25-5deb12ca4c28
md"""
## References
- [Stanford Notes](http://sepwww.stanford.edu/public/docs/sep131/rgunther1/paper_html/node3.html)
- [SEG Wiki](https://wiki.seg.org/wiki/Born-approximate_modeling_formula)
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Tullio = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"

[compat]
Bessels = "~0.2.8"
FFTW = "~1.10.0"
PlutoUI = "~0.7.83"
Tullio = "~0.3.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "c5bd21b734aed2f6a54a15f0f931f6b4a6f87ebe"

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

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

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

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

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

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "dbd2e8cd2c1c27f0b584f6661b4309609c5a685e"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

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

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "c3ac026e735264e9bdc6a9bcbd1b1e781b36e3bc"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.8.3"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

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

[[deps.Tullio]]
deps = ["DiffRules", "LinearAlgebra", "Requires"]
git-tree-sha1 = "de0febfe1243e89f352abd4ca0e9de6c8e6190c5"
uuid = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"
version = "0.3.9"

    [deps.Tullio.extensions]
    TullioCUDAExt = "CUDA"
    TullioChainRulesCoreExt = "ChainRulesCore"
    TullioFillArraysExt = "FillArrays"
    TullioTrackerExt = "Tracker"

    [deps.Tullio.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

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
# ╠═44d58ae7-f626-4bf6-a8ac-9dae63f2d9d6
# ╟─2f429d1f-ecd1-4b97-b46d-b4714144e24f
# ╟─24fab919-51b8-4ca3-ac6c-cdf2c4dbb407
# ╟─b313a158-3384-490f-8cdf-6786cd3facd4
# ╟─5a0cea57-1627-4c2a-841e-d3a7db5bdedc
# ╟─1ed696fe-3a60-43a7-8dc3-474b8f68c2b2
# ╟─689dd440-942d-449d-9b75-825ad3b052af
# ╟─e4de01e7-6d0a-4211-b83c-c292984c76f3
# ╠═7a4b1527-e09d-4393-bba1-69f7f7048e6e
# ╟─65a857c3-d6ab-4f9a-8839-1f829c4621c4
# ╠═4ba10533-065c-44b0-8295-2363fb99a8ea
# ╟─dff8cb60-4cbd-47bc-a803-9efbb71ee47b
# ╠═db3a4cb2-a525-46b2-aaf5-9a7d9078119f
# ╠═e18b63b0-ae7b-4d53-81f4-f72eb00a09d2
# ╟─856ee675-dff3-4d61-9b13-efd83f85cfc3
# ╠═fc2711fe-68d2-4432-9ed8-e511c7fa71fe
# ╠═04207c33-ce75-4b16-a22d-df46cc949ef9
# ╠═a219febc-278d-4093-aba4-2c7dc736f6ab
# ╠═c58bc7c7-e686-423f-94b2-5bb2ae531643
# ╟─90fcf549-8982-4bde-895d-eb9df3b32a41
# ╠═cd5c240c-db82-4307-8d5b-702a0b49a946
# ╟─80f38939-748f-41cc-83ec-762dbfaa6267
# ╟─c43e7a70-5895-4130-a8a2-988e5c9e6d9f
# ╠═51d129d0-48f2-444a-964b-dbbec58a0b38
# ╟─08d196f8-86cb-4c0f-8050-811003182671
# ╠═c70710d4-3918-4397-aaba-e3ce9c87f96a
# ╠═27ef0a2c-50b9-4d13-b50e-68786de4db91
# ╠═a3b1c2d4-5e6f-4708-9a1b-2c3d4e5f6071
# ╠═c006bf6d-1ca6-4349-a7f0-e5e81c1fa67b
# ╠═649a946a-8459-4bce-9ecd-aa6f771b8e6c
# ╠═082041a0-4a7a-41e4-aa28-9db8c421d463
# ╠═dca30e9e-6501-44df-a086-ee2780fb9f89
# ╠═ba46a17b-af4d-423e-8398-9c344d0df8bf
# ╠═ca1a7fed-2c65-4360-9e19-dfaf798bea1e
# ╠═604ef000-2a93-44bf-8373-e36a4294deb0
# ╠═022d3632-8aab-4b44-9867-a8a639bcdfa6
# ╠═6679f50a-7af5-455b-abf5-2f166248950b
# ╠═2d04488e-29bd-409f-8a19-5e84efc6b9b6
# ╠═a144eb61-7b1b-4ec5-9037-ea13bf6cd140
# ╠═8d3124af-39b2-4400-9fc7-246ee7180305
# ╠═34fb2cf5-8b2a-43c9-9bdf-c52d0583a909
# ╟─38cd5bf1-203f-49ce-810b-7ed60e88c705
# ╠═fb798d80-3f07-4294-a8af-7fc7c34a61f7
# ╟─59b5207f-7af0-4508-bb22-465ac206a10d
# ╟─527bf8f5-1b91-48a1-b235-c5774c7b1c70
# ╠═271bca15-eab3-4263-b2c3-deaedbfe12a8
# ╠═5f7b9237-8a45-4df7-a92e-eae3ade48ecd
# ╠═7bfd9cbb-c4d1-412a-8c13-af5db5bcca06
# ╠═60c5d139-4df7-45ea-a8a9-c11e13cce883
# ╟─978a06e3-f582-4070-ae25-5deb12ca4c28
# ╠═218dc0cf-7d17-4e9b-8393-ceca9c22acaf
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
