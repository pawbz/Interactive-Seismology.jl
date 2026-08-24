### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Multipole Scattering and the ka-kL Diagram"
#> tags = ["scattering"]
#> layout = "layout.jlhtml"
#> description = "Derive multipole scattering off a cylinder from the Helmholtz equation, then use it to build a classification diagram for single vs. multiple scattering."

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

# ╔═╡ eb35764e-ab0a-4b36-a09c-3e7b744685b1
begin
    using Bessels
    using PlutoUI
end

# ╔═╡ fcf7c3ea-a7af-4eda-b77a-cd797302aa60
TableOfContents()

# ╔═╡ 67f3ae02-dc84-46ef-b527-6a7ceab9ed82
md"""
# Multipole Scattering and the ka–kL Diagram

Seismologists classify how a wave "sees" a heterogeneous medium -- a scatterer, a patch of
crust, the whole mantle -- using two numbers: **ka**, the heterogeneity's size relative to
the wavelength, and **kL**, the propagation distance relative to the wavelength. A plot
built from these two numbers (our own version of it below, not a reproduction of any
published figure) sorts wave propagation into three regimes: an *equivalent homogeneous
body* (too short a path to notice the heterogeneity at all), *ray theory / single
scattering* (individual scatterers act independently), and *multiple scattering* (the wave
bounces between many heterogeneities before arriving).

This notebook builds that diagram from scratch, starting from the wave equation itself. The
path is: solve the Helmholtz equation for a plane wave hitting a single cylindrical
scatterer (a **multipole expansion**), use that solution's cross-section to define a
**mean free path** for a medium full of such scatterers, and turn the mean free path into
the diagram's boundary curves.

Play with the two widgets below first -- the derivation that follows explains *why* they
work the way they do.
"""

# ╔═╡ 5b0963f4-dfa6-4739-9614-dfe60d386e71
md"""
## What Does `ka` Actually Mean?

The circle you just dragged sets `ka = k·a`: the cylinder's radius `a` measured in units of
"radians per wavelength." At small `ka` the scattered pattern is broad and weak -- the
**Rayleigh** regime. At large `ka` the pattern narrows into a forward-peaked beam and the
scattering strength grows and then saturates -- the **geometric-optics** regime, where the
cylinder starts behaving like a solid obstacle blocking a beam of light. Everything on the
ka–kL diagram's vertical axis comes from this one number. The rest of this section
derives where the scattered pattern and its cross-section σ(ka) actually come from,
starting from the wave equation.
"""

# ╔═╡ a54d1f2a-46da-4b76-af72-d5bfcc9a0f58
md"""
### From the Wave Equation to Cylindrical Waves

A time-harmonic (angular frequency `ω`) scalar wave field
``p(\mathbf{x},t) = \mathrm{Re}[p(\mathbf{x})e^{-i\omega t}]`` satisfies the **Helmholtz
equation**

```math
\nabla^2 p + k^2 p = 0, \qquad k = \omega/v,
```

the frequency-domain form of the ordinary scalar wave equation, `v` the background wave
speed. Working in cylindrical coordinates ``(r,\theta)`` and separating variables
``p=R(r)\Theta(\theta)`` turns this into two ordinary differential equations: ``\Theta``
satisfies a simple harmonic equation with solutions ``e^{in\theta}`` for integer `n`
(single-valuedness in ``\theta`` forces `n` to be an integer), and `R` satisfies **Bessel's
equation** of order `n`,

```math
r^2 R'' + rR' + (k^2r^2 - n^2)R = 0,
```

whose two independent solutions are the Bessel function `J_n(kr)` (finite at `r=0`) and the
Hankel function `H_n^{(2)}(kr)` (outgoing as ``r\to\infty``, for the ``e^{-i\omega t}``
convention used throughout this notebook -- the same convention this repo's
[`Born-approximation.jl`](Born-approximation.jl) notebook already uses for its own
point-source Green's function). Every solution of the Helmholtz equation is built from the
cylindrical harmonics ``J_n(kr)e^{in\theta}`` and ``H_n^{(2)}(kr)e^{in\theta}`` -- this is the
alphabet the rest of the notebook writes in.
"""

# ╔═╡ 7c40b8e6-f103-4a6e-954b-1a9a0489c6ff
md"""
### The Multipole Expansion of a Plane Wave

An incident plane wave ``p_{\mathrm{inc}} = e^{ikx} = e^{ikr\cos\theta}`` (traveling in the
`+x` direction) is *itself* an exact sum of these cylindrical harmonics -- the
**Jacobi–Anger expansion**:

```math
e^{ikr\cos\theta} = \sum_{n=-\infty}^{\infty} i^n J_n(kr)\, e^{in\theta}.
```

This is the multipole expansion the notebook is named for: a single plane wave, written as
an infinite sum of `J_n` "partial waves," each with its own angular pattern ``e^{in\theta}``.
It is not an approximation -- the sum converges to the plane wave exactly -- but for a
fixed `kr` only the terms up to roughly ``n\sim kr`` contribute meaningfully (`J_n(kr)` decays
rapidly once `n` exceeds its argument), which is exactly why a *finite* truncated sum is
enough to compute with. The Appendix checks this expansion numerically against a direct
plane-wave evaluation before using it for anything.
"""

# ╔═╡ 2b29f7cb-15ba-4f89-86ea-bd5d94fb6417
md"""
### Scattering Off a Cylinder

Place a cylinder of radius `a` at the origin, with interior density ``\rho_2`` and wave
speed `v_2` (exterior ``\rho_1``, `v_1`; interior wavenumber ``k_2=\omega/v_2``). Since the
incident wave is already a sum of `J_n`-harmonics, look for the *total* field, order by
order in `n`, as

```math
p = \sum_n i^n\Big[J_n(k_1r) + b_n H_n^{(2)}(k_1r)\Big]e^{in\theta} \ \ (r>a), \qquad
p = \sum_n i^n a_n J_n(k_2r)\,e^{in\theta} \ \ (r<a),
```

i.e. the exterior field is the incident wave plus an *outgoing* scattered correction
(`H_n^{(2)}`, satisfying the radiation condition at infinity), and the interior field is
regular at the origin. The unknowns `a_n` (interior amplitude) and `b_n` (scattering
coefficient) are fixed by matching pressure and normal velocity across `r=a`, order by
order:

```math
J_n(k_1a) + b_n H_n^{(2)}(k_1a) = a_n J_n(k_2a), \qquad
\frac{k_1}{\rho_1}\Big[J_n'(k_1a) + b_n H_n^{(2)\prime}(k_1a)\Big]
= \frac{k_2}{\rho_2}\,a_n J_n'(k_2a).
```

This is a 2×2 linear system per order `n` -- [`cylinder_scattering_coefficients`](@ref)
in the Appendix solves it directly (numerically, per order) rather than by hand-deriving a
closed-form ratio. Sending the density ratio ``\rho_2/\rho_1\to\infty`` forces the interior
normal velocity to vanish (a **rigid** cylinder, Neumann condition); sending it to `0`
forces the interior pressure to vanish (a **pressure-release** cylinder, Dirichlet
condition) -- both limits are checked numerically against their textbook closed forms in
the Appendix. The `ρ₂/ρ₁` and `v₂/v₁` sliders on the widget above dial continuously between
these limits, including through the impedance-matched point ``\rho_1v_1=\rho_2v_2`` where the
contrast -- and with it, all scattering -- vanishes.
"""

# ╔═╡ fe0647ae-9326-493d-9556-e3c8258052d0
md"""
### From Coefficients to a Cross Section

Far from the cylinder, the `H_n^{(2)}(k_1r)` asymptotic form factors a common
``e^{-ik_1r}/\sqrt{r}`` out of every term in the scattered sum, leaving a purely angular
**far-field amplitude** ``f(\theta) = \sum_n b_n e^{in\theta}`` (see
[`far_field_amplitude`](@ref)'s docstring for exactly how the `i^n` phase cancels). The
total **scattering cross-section** (in 2-D, a length -- how much incident power per unit
wavefront the cylinder removes) is defined here by *directly integrating* the far-field
power, ``\sigma \propto \int_0^{2\pi}|f(\theta)|^2\,d\theta``, with the proportionality
constant fixed by matching the textbook closed form ``\sigma=(4/k)\sum_n|b_n|^2`` via
Parseval's identity -- see [`scattering_cross_section`](@ref) and its self-check below, so
no closed-form prefactor is trusted from memory alone.

σ(ka) (right-hand panel above) has two clean asymptotic regimes, both checked numerically
in the Appendix: at small `ka`, ``\sigma\propto (ka)^3`` (**Rayleigh scattering** --
small heterogeneities scatter weakly, and much more weakly still at low frequency); at
large `ka`, ``\sigma`` saturates toward `4a` (**geometric optics** -- twice the cylinder's
geometric shadow width, the classic *extinction paradox*: even a perfectly sharp obstacle
scatters more than its silhouette suggests, because diffraction around the edges adds as
much again as the direct shadow). These two regimes are exactly what makes the ka–kL
diagram's boundary curves bend, in the next section.
"""

# ╔═╡ a4d829be-5596-477b-8329-e64afc9ec003
md"""
## Mean Free Path and the D=1 Boundary

Now scatter a *field* of these cylinders: a dilute 2-D random medium with `n_d` scatterers
per unit area, each of radius `a` and size parameter `ka`. A wave loses, on average, a
fraction ``n_d\cdot\sigma(ka)`` of its coherent amplitude per unit distance traveled, which
defines the **mean free path**

```math
\ell = \frac{1}{n_d\,\sigma(ka)}, \qquad \sigma(ka) = a\cdot\text{scattering\_cross\_section}(ka)
```

(the extra factor of `a` rescales the unit-radius result above to a cylinder of actual
radius `a` -- the boundary-matching problem depends only on `ka`, so this scaling is exact,
not an approximation; see [`mean_free_path`](@ref)). The **optical depth** ``D=L/\ell`` counts
how many mean free paths a wave has crossed after traveling a distance `L`: ``D\ll 1`` means
the wave has essentially only single-scattered (or not scattered at all), while
``D\gtrsim 1`` means it has very likely bounced off more than one heterogeneity -- **multiple
scattering**. Setting `D=1` and solving for `L` as a function of `ka` gives the diagram's
second boundary curve directly, in the diagram's own coordinates
([`d1_boundary_kL`](@ref)):

```math
kL^*(ka) = \frac{ka}{n_d\,a^2\,\text{scattering\_cross\_section}(ka)}.
```

The widget below draws this curve live, alongside the trivial `ka=kL` line (where the
heterogeneity size `a` equals the propagation distance `L` itself -- too short a path to
speak of "the medium" statistically at all). Drag the point around to see which regime a
given frequency/distance combination falls into.
"""

# ╔═╡ 757357c5-7636-4fd2-85ac-8f8cf78dfe71
md"""
### Why the Boundary Bends

Notice the `D=1` curve is not a straight line: it comes down steeply from the Rayleigh
regime (small `ka`, where ``\sigma`` is tiny, so a wave must travel *very* far -- many
wavelengths' worth of heterogeneity size -- before it scatters even once), reaches a
minimum somewhere around ``ka\sim O(1)``, then rises again through the geometric-optics
regime (large `ka`, where ``\sigma`` has saturated to a constant, so ``kL^*\propto ka`` grows
linearly). That bend is not an artifact of drawing two separate line segments -- it falls
directly out of ``\sigma(ka)`` having two different power laws, exactly the two regimes
explored in the first widget. A classification diagram built this way from real, validated
multipole physics reproduces the qualitative shape of the textbook ka–kL diagrams
(Aki & Richards, 1980; Sato, Fehler & Maeda) without needing to copy their specific numbers.
"""

# ╔═╡ 7e1f4e2e-7ed9-4921-b750-988da8bf2ed7
md"""
## Appendix
"""

# ╔═╡ 9bc68957-8895-4913-b6a9-c42de31d952a
md"""
## Single-Cylinder Multipole Solution
"""

# ╔═╡ e9baabb0-3e63-402d-9830-4bd71b2ec6d1
begin
    # Derivative of Jₙ, H⁽²⁾ₙ via the standard three-term recurrence -- Bessels.jl doesn't
    # expose a derivative method directly.
    dbesselj(n, x) = 0.5 * (besselj(n - 1, x) - besselj(n + 1, x))
    dhankelh2(n, x) = 0.5 * (hankelh2(n - 1, x) - hankelh2(n + 1, x))
end

# ╔═╡ bd980060-a4ce-4035-bad6-e971d727b603
"""
	cylinder_scattering_coefficients(ka; density_ratio=8.0, velocity_ratio=0.6, nmax=nothing)

Multipole (partial-wave) coefficients for a plane wave scattering off a penetrable (fluid)
cylinder of nondimensional size `ka = k₁a`, with interior/exterior density ratio
`` \\rho_2/\\rho_1 = `` `density_ratio` and velocity ratio `` v_2/v_1 = `` `velocity_ratio`.
The interior wavenumber is `k₂ = ka/velocity_ratio`.

Matches pressure and normal-velocity continuity at the boundary order by order in `n`,
solving the resulting 2×2 linear system numerically (via Cramer's rule) rather than a
hand-derived closed form. `density_ratio → ∞` recovers the rigid (Neumann) limit
`` b_n = -J_n'(ka)/H_n'^{(2)}(ka) ``; `density_ratio → 0` recovers the pressure-release
(Dirichlet) limit `` b_n = -J_n(ka)/H_n^{(2)}(ka) `` -- both checked numerically below.

Returns `(a_n, b_n)`, complex vectors indexed `1:nmax+1` for orders `n=0:nmax` (only
non-negative orders are solved; every sum over the full `-N:N` range uses the symmetry
`b₋ₙ = bₙ`, `a₋ₙ = aₙ`).
"""
function cylinder_scattering_coefficients(ka; density_ratio=8.0, velocity_ratio=0.6, nmax=nothing)
    N = nmax === nothing ? ceil(Int, ka) + 20 : nmax
    k1, k2 = ka, ka / velocity_ratio
    rho1, rho2 = 1.0, density_ratio
    an = zeros(ComplexF64, N + 1)
    bn = zeros(ComplexF64, N + 1)
    for n in 0:N
        m11 = -besselj(n, k2)
        m12 = hankelh2(n, k1)
        m21 = -(k2 / rho2) * dbesselj(n, k2)
        m22 = (k1 / rho1) * dhankelh2(n, k1)
        r1 = -besselj(n, k1)
        r2 = -(k1 / rho1) * dbesselj(n, k1)
        det = m11 * m22 - m12 * m21
        an[n+1] = (r1 * m22 - m12 * r2) / det
        bn[n+1] = (m11 * r2 - r1 * m21) / det
    end
    return an, bn
end

# ╔═╡ c22db7a4-2b03-42c5-bfc8-2b3f5d1f3ad5
"""
	far_field_amplitude(theta, b_n)

Far-field scattering amplitude `` f(\\theta) = \\sum_n b_n e^{in\\theta} `` (sum over all
integers `-N:N`), built from the non-negative-order coefficients `b_n` returned by
[`cylinder_scattering_coefficients`](@ref) using `b₋ₙ = bₙ`. This is what remains of the
scattered field `` \\sum_n i^n b_n H_n^{(2)}(k_1r)e^{in\\theta} `` once the common
`H_n^{(2)}` large-argument asymptotic
`` H_n^{(2)}(k_1r) \\sim \\sqrt{2/(\\pi k_1r)}\\,e^{-i(k_1r-n\\pi/2-\\pi/4)} `` is factored
out: the `iⁿ` prefactor and the asymptotic's own `` e^{in\\pi/2} `` phase exactly cancel
(`` i^ne^{-in\\pi/2} = i^n(-i)^n = 1 ``), leaving this clean real-argument sum. The
scattered pressure at large `r` is then
`` p_{sc}\\approx f(\\theta)\\,e^{-i(k_1r-\\pi/4)}/\\sqrt{r} `` (up to a real constant).
"""
function far_field_amplitude(theta, b_n)
    s = b_n[1]
    for n in 1:length(b_n)-1
        s += b_n[n+1] * (cis(n * theta) + cis(-n * theta))
    end
    return s
end

# ╔═╡ a168c575-f6b0-47e9-97a6-7502b991bc8d
"""
	scattering_cross_section(ka; density_ratio=8.0, velocity_ratio=0.6)

Total scattering cross-section ("scattering width") of a *unit-radius* cylinder, defined by
directly integrating the far-field power,
`` \\sigma = \\frac{2}{\\pi ka}\\int_0^{2\\pi}|f(\\theta)|^2\\,d\\theta ``. The `2/(πka)`
prefactor is *derived*, not guessed: by Parseval's identity
`` \\int|f|^2d\\theta = 2\\pi\\sum_n|b_n|^2 ``, this expression is exactly the standard
closed form `` \\sigma = (4/k)\\sum_n|b_n|^2 `` (both are checked to agree numerically
below).

For a cylinder of actual radius `a`, the physical cross-section is
`a·scattering_cross_section(ka;...)` -- the boundary-matching problem in
[`cylinder_scattering_coefficients`](@ref) depends only on the nondimensional `ka`, so this
rescaling is exact (see [`mean_free_path`](@ref)).
"""
function scattering_cross_section(ka; density_ratio=8.0, velocity_ratio=0.6)
    _, b_n = cylinder_scattering_coefficients(ka; density_ratio, velocity_ratio)
    nth = 720
    flux = 0.0
    for i in 0:nth-1
        theta = 2pi * i / nth
        flux += abs2(far_field_amplitude(theta, b_n))
    end
    flux *= (2pi / nth)
    return (2 / (pi * ka)) * flux
end

# ╔═╡ da37e582-d104-4ecf-be26-339961f735fc
md"""
### Verifying the Multipole Solution
"""

# ╔═╡ 3a30e542-560c-47f3-8575-b0896e738cb9
let
    # 1. Jacobi-Anger expansion: exp(i k r cos θ) = Σ iⁿ Jₙ(kr) e^{inθ}
    r, th = 2.3, 0.7
    series = sum(cis(n * pi / 2) * besselj(n, r) * cis(n * th) for n in -30:30)
    exact = cis(r * cos(th))
    println("Jacobi-Anger check: |series - exact| = ", abs(series - exact))
    @assert abs(series - exact) < 1e-10

    # 2. rigid-cylinder limit (density_ratio -> Inf) vs the textbook closed form
    ka_test = 5.0
    _, bn_rigid = cylinder_scattering_coefficients(ka_test; density_ratio=1e8, velocity_ratio=1.0)
    closed = [-dbesselj(n, ka_test) / dhankelh2(n, ka_test) for n in 0:5]
    maxd = maximum(abs.(bn_rigid[1:6] .- closed))
    println("Rigid-limit check (ka=$ka_test): max|bₙ - closed form| = ", maxd)
    @assert maxd < 1e-6

    # 3. cross section vs the textbook closed form (4/k)Σ|bₙ|², cross-checking the derived
    #    Parseval prefactor in scattering_cross_section
    an, bn = cylinder_scattering_coefficients(ka_test; density_ratio=8.0, velocity_ratio=0.6)
    closed_form_sigma = (4 / ka_test) * (abs2(bn[1]) + 2 * sum(abs2, bn[2:end]))
    our_sigma = scattering_cross_section(ka_test; density_ratio=8.0, velocity_ratio=0.6)
    println("Cross-section check: our σ=$our_sigma, closed-form σ=$closed_form_sigma")
    @assert abs(our_sigma - closed_form_sigma) < 1e-9

    # 4. Rayleigh-regime power-law slope (expect ~3, a strongly-scattering near-rigid contrast
    #    so the asymptote is clean)
    kas = [0.003, 0.01, 0.03]
    sig = [scattering_cross_section(k; density_ratio=1e6, velocity_ratio=1.0) for k in kas]
    slope = log(sig[3] / sig[1]) / log(kas[3] / kas[1])
    println("Rayleigh-regime slope (expect ≈3): ", slope)
    @assert 2.9 < slope < 3.1

    # 5. geometric-optics plateau (extinction paradox: σ → 4a, a=1 here)
    plateau = scattering_cross_section(300.0; density_ratio=1e6, velocity_ratio=1.0) / 4
    println("Geometric-limit σ/4a (expect → 1): ", plateau)
    @assert plateau > 0.95

    "All multipole self-checks passed"
end

# ╔═╡ 41f157c2-ec27-4c07-b1f9-4f5b1d49c790
md"""
## Mean Free Path and the D=1 Boundary
"""

# ╔═╡ 95084277-fe1a-4b8c-bd70-5c7dac8dd277
"""
	mean_free_path(ka, a, n_d; density_ratio=8.0, velocity_ratio=0.6)

Scattering mean free path `` \\ell = 1/(n_d\\sigma) `` for a dilute 2-D random medium of
number density `n_d` (scatterers per unit area) of cylinders with radius `a` and size
parameter `ka`. Uses `` \\sigma = a\\cdot `` [`scattering_cross_section`](@ref)`(ka;...)`,
the physical (dimensional) cross-section obtained by scaling the unit-radius result by `a`
-- the problem is scale-invariant in `ka` alone, so this scaling is exact, not an
approximation.
"""
mean_free_path(ka, a, n_d; density_ratio=8.0, velocity_ratio=0.6) =
    1 / (n_d * a * scattering_cross_section(ka; density_ratio, velocity_ratio))

# ╔═╡ b32c9097-3026-4863-83b9-882d0b5bd236
"""
	d1_boundary_kL(ka, a, n_d; density_ratio=8.0, velocity_ratio=0.6)

The nondimensional propagation distance `kL` at which the optical depth `D=L/ℓ` first
reaches 1, as a function of `ka`, for a medium of scatterer radius `a` and number density
`n_d`. Derived by solving `` D = n_d\\,a\\,\\sigma(ka)\\,L = 1 `` for `L` and multiplying by
`` k=(ka)/a ``:

``` math
kL^*(ka) = \\frac{ka}{n_d\\,a^2\\,\\sigma(ka)}.
```

This is the "Ray theory `D=1`" boundary curve of the ka–kL classification diagram
(see [`mean_free_path`](@ref)); it is *not* a straight line in log-log `(kL,ka)` space,
because [`scattering_cross_section`](@ref)`(ka)` itself has two different power-law
regimes (Rayleigh at small `ka`, geometric-optics plateau at large `ka`) -- the curve bends
between the two.
"""
d1_boundary_kL(ka, a, n_d; density_ratio=8.0, velocity_ratio=0.6) =
    ka / (n_d * a^2 * scattering_cross_section(ka; density_ratio, velocity_ratio))

# ╔═╡ 20fd089f-d7e8-413f-a6da-98f091787c28
md"""
## Widget A: Single-Cylinder Multipole Explorer
"""

# ╔═╡ da44d725-33ec-4c52-bf84-1f85ea291e11
begin
    struct MultipoleScatteringInput
        ka::Float64
        density_ratio::Float64
        velocity_ratio::Float64
    end
    MultipoleScatteringInput(; ka=3.0, density_ratio=8.0, velocity_ratio=0.6) =
        MultipoleScatteringInput(ka, density_ratio, velocity_ratio)

    Base.get(w::MultipoleScatteringInput) = Dict{String,Any}(
        "ka" => w.ka, "density_ratio" => w.density_ratio, "velocity_ratio" => w.velocity_ratio)

    function Base.show(io::IO, ::MIME"text/html", w::MultipoleScatteringInput)
        write(io, """
        <div id="mswidget">
        <style>
        pluto-cell:has(#mswidget) { width: min(80vw, 1150px) !important;
          margin-left: calc((100% - min(80vw, 1150px)) / 2) !important; }
        #mswidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #mswidget .ms-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #mswidget .ms-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #mswidget .ms-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #mswidget .ms-primary{display:flex;gap:14px;flex-wrap:wrap;justify-content:center;align-items:flex-start;margin-bottom:14px}
        #mswidget .ms-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #mswidget .ms-panel-title{font-size:14px;font-weight:700;color:#e5e7eb;text-align:center;margin-bottom:6px}
        #mswidget .ms-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:4px;max-width:300px}
        #mswidget canvas{display:block}
        #mswidget #ms-polar{cursor:grab}
        #mswidget .ms-controls-row{display:flex;justify-content:center}
        #mswidget .ms-control-group{flex:1 1 700px;max-width:1000px;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 14px}
        #mswidget .ms-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #mswidget .ms-control-row{display:grid;grid-template-columns:140px minmax(0,1fr) 90px;gap:8px;align-items:center;margin:6px 0}
        #mswidget .ms-control-row label{font-size:13px;color:#9ca3af}
        #mswidget .ms-control-row input[type=range]{width:100%;min-width:0}
        #mswidget .ms-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        </style>

        <div class="ms-title">
          <div class="ms-title-desc">Drag the scatterer at the center of the pattern to resize it (ka) and see the multipole expansion change character.</div>
          <div class="ms-title-hint">small ka &rarr; weak, nearly isotropic (Rayleigh) scattering &middot; large ka &rarr; strong, forward-peaked (near-geometric) scattering</div>
        </div>

        <div class="ms-primary">
          <div>
            <div class="ms-panel-title">Scattering Pattern</div>
            <div class="ms-panel"><canvas id="ms-polar"></canvas></div>
            <div class="ms-caption">|f(&theta;)|&sup2; &middot; drag the circle &middot; incident wave arrives from the left</div>
          </div>
          <div>
            <div class="ms-panel-title">Partial-Wave Coefficients</div>
            <div class="ms-panel"><canvas id="ms-bars"></canvas></div>
            <div class="ms-caption" id="ms-bars-caption"></div>
          </div>
          <div>
            <div class="ms-panel-title">Cross Section &sigma;(ka)</div>
            <div class="ms-panel"><canvas id="ms-sigma"></canvas></div>
            <div class="ms-caption">Rayleigh (&sigma;&prop;ka&sup3;) at small ka &middot; geometric plateau at large ka</div>
          </div>
        </div>

        <div class="ms-controls-row">
          <div class="ms-control-group">
            <div class="ms-control-title">Contrast</div>
            <div class="ms-control-row"><label>&rho;&#8322;/&rho;&#8321; (density)</label><input type="range" id="ms-dr" min="1" max="30" step="0.1" value="$(w.density_ratio)"><span class="ms-value" id="ms-dr-v"></span></div>
            <div class="ms-control-row"><label>v&#8322;/v&#8321; (velocity)</label><input type="range" id="ms-vr" min="0.15" max="3" step="0.01" value="$(w.velocity_ratio)"><span class="ms-value" id="ms-vr-v"></span></div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        let state = { ka: $(w.ka), density_ratio: $(w.density_ratio), velocity_ratio: $(w.velocity_ratio) };
        let pushed = null; // {bn_mag, pattern, ka_grid, sigma_vals, sigma_current} from Julia
        let commitInFlight = false;

        const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1150);
        const GAP = 14, PANEL_OVERHEAD = 18;
        const PW = Math.max(220, Math.floor((availW - 2*GAP - 3*PANEL_OVERHEAD) / 3));
        const PH = PW;
        const DPR = window.devicePixelRatio || 1;

        function hidpi(canvas, ctx, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          ctx.setTransform(DPR,0,0,DPR,0,0);
        }

        const polarCv = par.querySelector('#ms-polar'), polarCtx = polarCv.getContext('2d');
        hidpi(polarCv, polarCtx, PW, PH);
        const barsCv = par.querySelector('#ms-bars'), barsCtx = barsCv.getContext('2d');
        hidpi(barsCv, barsCtx, PW, PH);
        const sigmaCv = par.querySelector('#ms-sigma'), sigmaCtx = sigmaCv.getContext('2d');
        hidpi(sigmaCv, sigmaCtx, PW, PH);

        // the polar plot's center doubles as a draggable "scatterer" -- its radius encodes ka
        // on a log scale, dragged/scrubbed directly instead of a separate ka slider.
        const POLAR_CX = PW/2, POLAR_CY = PH/2, POLAR_R = Math.min(PW,PH)/2 - 26;
        const KA_MIN = 0.1, KA_MAX = 31.6;
        const CIRCLE_R_MIN = 10, CIRCLE_R_MAX = POLAR_R - 6;
        function kaToCircleR(ka){
          const t = (Math.log10(ka)-Math.log10(KA_MIN)) / (Math.log10(KA_MAX)-Math.log10(KA_MIN));
          return CIRCLE_R_MIN + (CIRCLE_R_MAX-CIRCLE_R_MIN) * Math.min(1,Math.max(0,t));
        }
        function circleRToKa(r){
          const t = Math.min(1, Math.max(0, (r-CIRCLE_R_MIN)/(CIRCLE_R_MAX-CIRCLE_R_MIN)));
          return Math.pow(10, Math.log10(KA_MIN) + t*(Math.log10(KA_MAX)-Math.log10(KA_MIN)));
        }
        let draggingKa = false;

        function emit(){
          commitInFlight = true;
          par.value = { ka: state.ka, density_ratio: state.density_ratio, velocity_ratio: state.velocity_ratio };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledEmit(){ if(!commitInFlight) emit(); }

        function drawPolar(){
          const ctx = polarCtx, W = PW, H = PH;
          ctx.clearRect(0,0,W,H);
          const cx = POLAR_CX, cy = POLAR_CY, R = POLAR_R;
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          for(const frac of [0.25,0.5,0.75,1.0]){
            ctx.beginPath(); ctx.arc(cx, cy, R*frac, 0, 2*Math.PI); ctx.stroke();
          }
          ctx.fillStyle = '#9ca3af'; ctx.font = '10px sans-serif'; ctx.textAlign='left';
          ctx.fillText('max', cx + R + 3, cy);
          ctx.strokeStyle = '#9ca3af'; ctx.fillStyle = '#9ca3af'; ctx.lineWidth = 1.5;
          const ax0 = cx - R - 20, ax1 = cx - R - 4;
          ctx.beginPath(); ctx.moveTo(ax0, cy); ctx.lineTo(ax1, cy); ctx.stroke();
          ctx.beginPath(); ctx.moveTo(ax1, cy); ctx.lineTo(ax1-6, cy-4); ctx.lineTo(ax1-6, cy+4); ctx.closePath(); ctx.fill();
          ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
          ctx.fillText('incident', ax0 - 10, cy - 8);
          if(pushed){
            const pattern = pushed.pattern;
            const n = pattern.length;
            const mx = Math.max(...pattern, 1e-12);
            ctx.beginPath();
            for(let i=0;i<=n;i++){
              const theta = (i%n) * 2*Math.PI/n;
              const r = R * Math.sqrt(pattern[i%n]/mx);
              const x = cx + r*Math.cos(theta), y = cy - r*Math.sin(theta);
              i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
            }
            ctx.closePath();
            ctx.strokeStyle = '#3b82f6'; ctx.lineWidth = 2; ctx.stroke();
            ctx.fillStyle = 'rgba(59,130,246,0.25)'; ctx.fill();
          }
          // draggable "scatterer" -- its radius (log scale) IS ka, dragged/scrubbed directly
          // instead of a separate slider; always drawn on top so it stays graspable.
          const scR = kaToCircleR(state.ka);
          ctx.beginPath(); ctx.arc(cx, cy, scR, 0, 2*Math.PI);
          ctx.fillStyle = 'rgba(15,20,30,0.9)'; ctx.fill();
          ctx.strokeStyle = draggingKa ? '#f3f4f6' : '#9ca3af'; ctx.lineWidth = 2; ctx.stroke();
          ctx.fillStyle = '#f3f4f6'; ctx.font = 'bold 13px sans-serif';
          ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
          ctx.fillText('ka=' + state.ka.toFixed(state.ka<1?2:1), cx, cy);
          ctx.textBaseline = 'alphabetic';
        }

        function drawBars(){
          const ctx = barsCtx, W = PW, H = PH;
          ctx.clearRect(0,0,W,H);
          const cap = par.querySelector('#ms-bars-caption');
          if(!pushed){ return; }
          const bn = pushed.bn_mag;
          const mx = Math.max(...bn, 1e-12);
          let last = 0;
          for(let i=0;i<bn.length;i++) if(bn[i] > 1e-3*mx) last = i;
          const nshow = Math.min(bn.length, Math.max(last+4, 6));
          if(cap) cap.textContent = nshow + ' of ' + bn.length + ' computed orders shown (n=0..' + (nshow-1) + ')';
          const padL=28, padB=20, padT=10, padR=8;
          const plotW = W-padL-padR, plotH = H-padT-padB;
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          ctx.beginPath(); ctx.moveTo(padL, padT); ctx.lineTo(padL, H-padB); ctx.lineTo(W-padR, H-padB); ctx.stroke();
          const bw = plotW/nshow;
          for(let i=0;i<nshow;i++){
            const h = (bn[i]/mx) * (plotH-4);
            ctx.fillStyle = '#3b82f6';
            ctx.fillRect(padL + i*bw + bw*0.15, H-padB-h, bw*0.7, h);
          }
          ctx.fillStyle = '#9ca3af'; ctx.font = '10px sans-serif'; ctx.textAlign='center';
          for(let i=0;i<nshow;i+=Math.max(1,Math.round(nshow/8))){
            ctx.fillText(String(i), padL + i*bw + bw/2, H-padB+12);
          }
          ctx.save(); ctx.translate(10, H/2); ctx.rotate(-Math.PI/2);
          ctx.fillText('|b\\u2099|', 0, 0); ctx.restore();
        }

        function drawSigma(){
          const ctx = sigmaCtx, W = PW, H = PH;
          ctx.clearRect(0,0,W,H);
          if(!pushed){ return; }
          const kaG = pushed.ka_grid, sig = pushed.sigma_vals;
          const padL=34, padB=22, padT=10, padR=10;
          const plotW = W-padL-padR, plotH = H-padT-padB;
          const kaMin = Math.log10(kaG[0]), kaMax = Math.log10(kaG[kaG.length-1]);
          let sMin = Math.min(...sig), sMax = Math.max(...sig);
          sMin = Math.max(sMin, 1e-8); sMin = Math.log10(sMin); sMax = Math.log10(Math.max(sMax,1e-7));
          if(sMax - sMin < 1) { sMax = sMin+1; }
          function X(ka){ return padL + (Math.log10(ka)-kaMin)/(kaMax-kaMin) * plotW; }
          function Y(s){ return padT + (1 - (Math.log10(Math.max(s,1e-12))-sMin)/(sMax-sMin)) * plotH; }
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          ctx.beginPath(); ctx.moveTo(padL, padT); ctx.lineTo(padL, H-padB); ctx.lineTo(W-padR, H-padB); ctx.stroke();
          ctx.beginPath();
          for(let i=0;i<kaG.length;i++){
            const x = X(kaG[i]), y = Y(sig[i]);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.strokeStyle = '#38bdf8'; ctx.lineWidth = 2; ctx.stroke();
          const mx = X(state.ka), my = Y(pushed.sigma_current);
          ctx.beginPath(); ctx.arc(mx, my, 4, 0, 2*Math.PI);
          ctx.fillStyle = '#ef4444'; ctx.fill();
          ctx.fillStyle = '#9ca3af'; ctx.font = '10px sans-serif'; ctx.textAlign='center';
          ctx.fillText('ka', W/2, H-4);
          ctx.save(); ctx.translate(10, H/2); ctx.rotate(-Math.PI/2);
          ctx.fillText('\\u03c3', 0, 0); ctx.restore();
        }

        function draw(){ drawPolar(); drawBars(); drawSigma(); }

        function syncControls(){
          par.querySelector('#ms-dr-v').textContent = state.density_ratio.toFixed(1);
          par.querySelector('#ms-vr-v').textContent = state.velocity_ratio.toFixed(2);
        }

        function onControl(event){
          const id = event.target.id;
          if(id === 'ms-dr') state.density_ratio = Number(event.target.value);
          else if(id === 'ms-vr') state.velocity_ratio = Number(event.target.value);
          else return;
          syncControls();
          throttledEmit();
        }
        par.querySelectorAll('input[type=range]').forEach(el => el.addEventListener('input', onControl));

        function polarPointer(ev){
          const r = polarCv.getBoundingClientRect();
          return [ev.clientX - r.left, ev.clientY - r.top];
        }
        function updateKaFromPointer(px, py){
          state.ka = circleRToKa(Math.hypot(px-POLAR_CX, py-POLAR_CY));
          drawPolar();
          drawSigma();
          throttledEmit();
        }
        polarCv.addEventListener('mousedown', ev => {
          draggingKa = true;
          polarCv.style.cursor = 'grabbing';
          updateKaFromPointer(...polarPointer(ev));
        });
        window.addEventListener('mousemove', ev => {
          if(!draggingKa) return;
          updateKaFromPointer(...polarPointer(ev));
        });
        window.addEventListener('mouseup', () => {
          if(!draggingKa) return;
          draggingKa = false;
          polarCv.style.cursor = 'grab';
          drawPolar();
        });

        par.addEventListener('ms-results', event => {
          pushed = event.detail || null;
          commitInFlight = false;
          draw();
        });

        syncControls();
        draw();
        }
        </script>
        """)
    end

    const _ms_ready = true
end

# ╔═╡ 501ec109-c79b-49f0-b3bc-8fee55e8653f
begin
    _ms_ready
    WideCell(@bind ms MultipoleScatteringInput(); max_width=1150)
end

# ╔═╡ 1d8bf334-c69e-4f90-af0a-28cd3af723bd
begin
    struct MsPush
        bn_mag::String
        pattern::String
        ka_grid::String
        sigma_vals::String
        sigma_current::Float64
    end
    function Base.show(io::IO, ::MIME"text/html", p::MsPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('mswidget');
        if(w){
          w.dispatchEvent(new CustomEvent('ms-results', { detail: {
            bn_mag: [$(p.bn_mag)],
            pattern: [$(p.pattern)],
            ka_grid: [$(p.ka_grid)],
            sigma_vals: [$(p.sigma_vals)],
            sigma_current: $(p.sigma_current),
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 4d3f3baa-4d61-4701-923c-dd5fd940cf50
begin
    local ms_ka = ms isa AbstractDict ? ms["ka"] : 3.0
    local ms_density_ratio = ms isa AbstractDict ? ms["density_ratio"] : 8.0
    local ms_velocity_ratio = ms isa AbstractDict ? ms["velocity_ratio"] : 0.6
    local an, bn = cylinder_scattering_coefficients(ms_ka; density_ratio=ms_density_ratio, velocity_ratio=ms_velocity_ratio)
    local nth = 180
    local pattern = [abs2(far_field_amplitude(2pi * i / nth, bn)) for i in 0:nth-1]
    local ka_grid = 10.0 .^ range(log10(0.1), log10(31.6); length=60)
    local sigma_vals = [scattering_cross_section(k; density_ratio=ms_density_ratio, velocity_ratio=ms_velocity_ratio) for k in ka_grid]
    local sigma_current = scattering_cross_section(ms_ka; density_ratio=ms_density_ratio, velocity_ratio=ms_velocity_ratio)
    MsPush(join(abs.(bn), ","), join(pattern, ","), join(ka_grid, ","), join(sigma_vals, ","), sigma_current)
end

# ╔═╡ 1b9b4706-65af-46b8-b249-68c6417e2b95
md"""
## Widget B: Scattering Regime Map
"""

# ╔═╡ 6750958f-4d00-4671-b70b-a0c7b59c3058
begin
    struct ScatteringRegimeMapInput
        a::Float64      # scatterer radius, km
        n_d::Float64    # number density, scatterers per km^2
        v::Float64      # background velocity, km/s -- display-only (converts ka,kL to f,L), not used by the Julia physics
        ka::Float64     # initial marker position
        kL::Float64
    end
    ScatteringRegimeMapInput(; a=1.0, n_d=1.0e-3, v=4.0, ka=3.0, kL=50.0) =
        ScatteringRegimeMapInput(a, n_d, v, ka, kL)

    Base.get(w::ScatteringRegimeMapInput) = Dict{String,Any}(
        "a" => w.a, "n_d" => w.n_d, "v" => w.v, "ka" => w.ka, "kL" => w.kL)

    function Base.show(io::IO, ::MIME"text/html", w::ScatteringRegimeMapInput)
        write(io, """
        <div id="srmwidget">
        <style>
        pluto-cell:has(#srmwidget) { width: min(80vw, 1050px) !important;
          margin-left: calc((100% - min(80vw, 1050px)) / 2) !important; }
        #srmwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #srmwidget .srm-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #srmwidget .srm-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #srmwidget .srm-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #srmwidget .srm-row{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start}
        #srmwidget .srm-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #srmwidget .srm-panel-title{font-size:14px;font-weight:700;color:#e5e7eb;text-align:center;margin-bottom:6px}
        #srmwidget canvas{display:block;cursor:grab}
        #srmwidget .srm-side{display:flex;flex-direction:column;gap:12px;min-width:230px}
        #srmwidget .srm-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #srmwidget .srm-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #srmwidget .srm-control-row{display:grid;grid-template-columns:64px minmax(0,1fr) 64px;gap:6px;align-items:center;margin:6px 0}
        #srmwidget .srm-control-row label{font-size:13px;color:#9ca3af}
        #srmwidget .srm-control-row input[type=range]{width:100%;min-width:0}
        #srmwidget .srm-value{font-size:12px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #srmwidget .srm-readout-row{display:grid;grid-template-columns:1fr auto;gap:6px;margin:4px 0;font-size:13px}
        #srmwidget .srm-readout-row label{color:#9ca3af}
        #srmwidget .srm-readout-row span{color:#e5e7eb;text-align:right}
        #srmwidget .srm-regime{font-size:15px;font-weight:700;text-align:center;margin-top:6px;padding:8px;border-radius:4px}
        #srmwidget .srm-legend{display:flex;gap:10px;flex-wrap:wrap;justify-content:center;font-size:11px;color:#9ca3af;margin-top:8px}
        #srmwidget .srm-swatch{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:middle}
        </style>

        <div class="srm-title">
          <div class="srm-title-desc">Drag the point to explore the ka&ndash;kL classification diagram.</div>
          <div class="srm-title-hint">horizontal = kL (propagation distance) &middot; vertical = ka (heterogeneity size), both relative to wavelength</div>
        </div>

        <div class="srm-row">
          <div>
            <div class="srm-panel-title">Regime Map</div>
            <div class="srm-panel"><canvas id="srm-map"></canvas></div>
            <div class="srm-legend">
              <span><span class="srm-swatch" style="background:#22c55e"></span>equivalent homogeneous body</span>
              <span><span class="srm-swatch" style="background:#3b82f6"></span>ray theory / single scattering</span>
              <span><span class="srm-swatch" style="background:#ef4444"></span>multiple scattering</span>
            </div>
          </div>
          <div class="srm-side">
            <div class="srm-control-group">
              <div class="srm-control-title">Medium</div>
              <div class="srm-control-row"><label>a (km)</label><input type="range" id="srm-a" min="-1" max="1" step="0.01" value="$(log10(w.a))"><span class="srm-value" id="srm-a-v"></span></div>
              <div class="srm-control-row"><label>n&#8320; (/km&sup2;)</label><input type="range" id="srm-nd" min="-5" max="-2" step="0.02" value="$(log10(w.n_d))"><span class="srm-value" id="srm-nd-v"></span></div>
              <div class="srm-control-row"><label>v (km/s)</label><input type="range" id="srm-v" min="1" max="8" step="0.1" value="$(w.v)"><span class="srm-value" id="srm-v-v"></span></div>
            </div>
            <div class="srm-control-group">
              <div class="srm-control-title">Readouts</div>
              <div class="srm-readout-row"><label>frequency</label><span id="srm-f"></span></div>
              <div class="srm-readout-row"><label>distance L</label><span id="srm-L"></span></div>
              <div class="srm-readout-row"><label>mean free path &#8467;</label><span id="srm-ell"></span></div>
              <div class="srm-readout-row"><label>optical depth D</label><span id="srm-D"></span></div>
              <div class="srm-regime" id="srm-regime"></div>
            </div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        let state = { a: $(w.a), n_d: $(w.n_d), v: $(w.v), ka: $(w.ka), kL: $(w.kL) };
        let pushed = null; // {ka_grid, kL_d1} from Julia
        let commitInFlight = false;
        let dragging = false;

        const MAPW = 560, MAPH = 560;
        const KA_LOG_MIN = -1, KA_LOG_MAX = 4;   // ka: 0.1 .. 1e4
        const KL_LOG_MIN = 0, KL_LOG_MAX = 5;    // kL: 1 .. 1e5
        const DPR = window.devicePixelRatio || 1;

        function hidpi(canvas, ctx, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          ctx.setTransform(DPR,0,0,DPR,0,0);
        }
        const mapCv = par.querySelector('#srm-map'), mapCtx = mapCv.getContext('2d');
        hidpi(mapCv, mapCtx, MAPW, MAPH);

        const padL = 46, padB = 30, padT = 10, padR = 10;
        const plotW = MAPW - padL - padR, plotH = MAPH - padT - padB;
        function X(kL){ return padL + (Math.log10(kL)-KL_LOG_MIN)/(KL_LOG_MAX-KL_LOG_MIN) * plotW; }
        function Y(ka){ return padT + (1 - (Math.log10(ka)-KA_LOG_MIN)/(KA_LOG_MAX-KA_LOG_MIN)) * plotH; }
        function invX(px){ return Math.pow(10, KL_LOG_MIN + (px-padL)/plotW*(KL_LOG_MAX-KL_LOG_MIN)); }
        function invY(py){ return Math.pow(10, KA_LOG_MIN + (1-(py-padT)/plotH)*(KA_LOG_MAX-KA_LOG_MIN)); }

        function interpKLd1(ka){
          if(!pushed) return NaN;
          const g = pushed.ka_grid, d = pushed.kL_d1;
          if(ka <= g[0]) return d[0];
          if(ka >= g[g.length-1]) return d[d.length-1];
          let lo = 0, hi = g.length-1;
          while(hi-lo>1){ const mid=(lo+hi)>>1; (g[mid] <= ka) ? lo=mid : hi=mid; }
          const t = (Math.log10(ka)-Math.log10(g[lo])) / (Math.log10(g[hi])-Math.log10(g[lo]));
          return Math.pow(10, Math.log10(d[lo]) + t*(Math.log10(d[hi])-Math.log10(d[lo])));
        }

        function emit(){
          commitInFlight = true;
          par.value = { a: state.a, n_d: state.n_d, v: state.v, ka: state.ka, kL: state.kL };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledEmit(){ if(!commitInFlight) emit(); }

        function fillRegions(){
          if(!pushed) return;
          const g = pushed.ka_grid, d = pushed.kL_d1;
          const COL_HOMOG = 'rgba(34,197,94,0.16)', COL_SINGLE = 'rgba(59,130,246,0.14)', COL_MULTI = 'rgba(239,68,68,0.16)';
          for(let i=0;i<g.length-1;i++){
            const kaLo = g[i], kaHi = g[i+1];
            const yTop = Y(kaHi), yBot = Y(kaLo);
            const b1 = kaLo, b2 = d[i];
            const homogX1 = X(Math.max(1, Math.min(b1, Math.pow(10,KL_LOG_MAX))));
            mapCtx.fillStyle = COL_HOMOG;
            mapCtx.fillRect(padL, yTop, Math.max(0,homogX1-padL), yBot-yTop);
            const hi = Math.max(b1,b2);
            const singleX0 = X(Math.max(1,Math.min(b1,Math.pow(10,KL_LOG_MAX))));
            const singleX1 = X(Math.max(1,Math.min(hi,Math.pow(10,KL_LOG_MAX))));
            mapCtx.fillStyle = COL_SINGLE;
            mapCtx.fillRect(singleX0, yTop, Math.max(0,singleX1-singleX0), yBot-yTop);
            const multiX0 = singleX1;
            mapCtx.fillStyle = COL_MULTI;
            mapCtx.fillRect(multiX0, yTop, Math.max(0,(padL+plotW)-multiX0), yBot-yTop);
          }
        }

        function drawAxes(){
          const ctx = mapCtx;
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          ctx.strokeRect(padL, padT, plotW, plotH);
          ctx.font = '10px sans-serif';
          ctx.textAlign='center';
          for(let e=KL_LOG_MIN; e<=KL_LOG_MAX; e++){
            const x = X(Math.pow(10,e));
            ctx.strokeStyle = '#1f2937'; ctx.beginPath(); ctx.moveTo(x,padT); ctx.lineTo(x,padT+plotH); ctx.stroke();
            ctx.fillStyle = '#9ca3af';
            ctx.fillText('1e'+e, x, padT+plotH+14);
          }
          ctx.textAlign='right';
          for(let e=KA_LOG_MIN; e<=KA_LOG_MAX; e++){
            const y = Y(Math.pow(10,e));
            ctx.strokeStyle = '#1f2937'; ctx.beginPath(); ctx.moveTo(padL,y); ctx.lineTo(padL+plotW,y); ctx.stroke();
            ctx.fillStyle = '#9ca3af';
            ctx.fillText('1e'+e, padL-4, y+3);
          }
          ctx.textAlign='center'; ctx.fillStyle = '#d1d5db'; ctx.font = '12px sans-serif';
          ctx.fillText('kL', padL+plotW/2, MAPH-2);
          ctx.save(); ctx.translate(10, padT+plotH/2); ctx.rotate(-Math.PI/2);
          ctx.fillText('ka', 0, 0); ctx.restore();
        }

        function drawBoundaries(){
          const ctx = mapCtx;
          ctx.strokeStyle = '#e5e7eb'; ctx.setLineDash([5,4]); ctx.lineWidth = 1.5;
          ctx.beginPath();
          ctx.moveTo(X(Math.pow(10,KA_LOG_MIN)), Y(Math.pow(10,KA_LOG_MIN)));
          ctx.lineTo(X(Math.pow(10,KA_LOG_MAX)), Y(Math.pow(10,KA_LOG_MAX)));
          ctx.stroke();
          ctx.setLineDash([]);
          if(pushed){
            const g = pushed.ka_grid, d = pushed.kL_d1;
            ctx.strokeStyle = '#facc15'; ctx.lineWidth = 2;
            ctx.beginPath();
            for(let i=0;i<g.length;i++){
              const x = X(Math.min(Math.max(d[i],1),Math.pow(10,KL_LOG_MAX)));
              const y = Y(g[i]);
              i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
            }
            ctx.stroke();
          }
        }

        function drawMarker(){
          const ctx = mapCtx;
          const x = X(Math.min(Math.max(state.kL,1),1e5)), y = Y(Math.min(Math.max(state.ka,0.1),1e4));
          ctx.beginPath(); ctx.arc(x,y,7,0,2*Math.PI);
          ctx.fillStyle = '#f3f4f6'; ctx.fill();
          ctx.strokeStyle = '#111827'; ctx.lineWidth = 2; ctx.stroke();
        }

        function fmt(x, unit){
          if(!isFinite(x)) return '&mdash;';
          if(x >= 1e4 || x < 1e-2) return x.toExponential(2) + (unit?(' '+unit):'');
          return x.toPrecision(3) + (unit?(' '+unit):'');
        }

        function updateReadouts(){
          const ka = state.ka, kL = state.kL, a = state.a, v = state.v;
          const k = ka/a;
          const f = k*v/(2*Math.PI);
          const L = kL*a/ka;
          const kLd1 = interpKLd1(ka);
          const D = kL/kLd1;
          const ell = L/D;
          par.querySelector('#srm-f').innerHTML = fmt(f,'Hz');
          par.querySelector('#srm-L').innerHTML = fmt(L,'km');
          par.querySelector('#srm-ell').innerHTML = fmt(ell,'km');
          par.querySelector('#srm-D').innerHTML = fmt(D,'');
          const reg = par.querySelector('#srm-regime');
          if(kL < ka){ reg.textContent = 'Equivalent homogeneous body'; reg.style.background='#052e16'; reg.style.color='#4ade80'; }
          else if(D < 1){ reg.textContent = 'Ray theory / single scattering'; reg.style.background='#0c1e3d'; reg.style.color='#60a5fa'; }
          else { reg.textContent = 'Multiple scattering'; reg.style.background='#3d0c0c'; reg.style.color='#f87171'; }
        }

        function draw(){
          mapCtx.clearRect(0,0,MAPW,MAPH);
          fillRegions();
          drawAxes();
          drawBoundaries();
          drawMarker();
          updateReadouts();
        }

        function syncControls(){
          par.querySelector('#srm-a-v').textContent = state.a.toFixed(2)+' km';
          par.querySelector('#srm-nd-v').textContent = state.n_d.toExponential(1);
          par.querySelector('#srm-v-v').textContent = state.v.toFixed(1)+' km/s';
        }

        function onControl(event){
          const id = event.target.id;
          if(id === 'srm-a') state.a = Math.pow(10, Number(event.target.value));
          else if(id === 'srm-nd') state.n_d = Math.pow(10, Number(event.target.value));
          else if(id === 'srm-v'){ state.v = Number(event.target.value); syncControls(); draw(); return; }
          else return;
          syncControls();
          throttledEmit();
        }
        par.querySelectorAll('input[type=range]').forEach(el => el.addEventListener('input', onControl));

        function hitMarker(px,py){
          const x = X(state.kL), y = Y(state.ka);
          return Math.hypot(px-x,py-y) < 14;
        }
        function pxFromEvent(ev){
          const r = mapCv.getBoundingClientRect();
          return [ (ev.clientX-r.left), (ev.clientY-r.top) ];
        }
        mapCv.addEventListener('mousedown', ev => {
          const [px,py] = pxFromEvent(ev);
          dragging = true;
          mapCv.style.cursor = 'grabbing';
          moveMarkerTo(px,py);
        });
        window.addEventListener('mousemove', ev => {
          if(!dragging) return;
          const [px,py] = pxFromEvent(ev);
          moveMarkerTo(px,py);
        });
        window.addEventListener('mouseup', () => { dragging = false; mapCv.style.cursor = 'grab'; });
        function moveMarkerTo(px,py){
          const cx = Math.min(Math.max(px, padL), padL+plotW);
          const cy = Math.min(Math.max(py, padT), padT+plotH);
          state.kL = invX(cx);
          state.ka = invY(cy);
          draw();
        }

        par.addEventListener('srm-results', event => {
          pushed = event.detail || null;
          commitInFlight = false;
          draw();
        });

        syncControls();
        draw();
        }
        </script>
        """)
    end

    const _srm_ready = true
end

# ╔═╡ e4c3bf1e-beb5-46d5-b538-faa044c02f72
begin
    _srm_ready
    WideCell(@bind srm ScatteringRegimeMapInput(); max_width=1050)
end

# ╔═╡ c80bd048-a116-470d-82e9-58f5a2e7759a
begin
    struct SrmPush
        ka_grid::String
        kL_d1::String
    end
    function Base.show(io::IO, ::MIME"text/html", p::SrmPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('srmwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('srm-results', { detail: {
            ka_grid: [$(p.ka_grid)],
            kL_d1: [$(p.kL_d1)],
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ fb54e27f-578d-472c-bc64-0c1946d719eb
begin
    local srm_a = srm isa AbstractDict ? srm["a"] : 1.0
    local srm_n_d = srm isa AbstractDict ? srm["n_d"] : 1.0e-3
    local ka_grid = 10.0 .^ range(log10(0.1), log10(1.0e4); length=150)
    local kL_d1 = [d1_boundary_kL(k, srm_a, srm_n_d) for k in ka_grid]
    SrmPush(join(ka_grid, ","), join(kL_d1, ","))
end

# ╔═╡ 6875a305-087a-49bb-9da2-114538b98dc8
md"""
## References
- Aki, K., & Richards, P. G. (1980 / 2002). *Quantitative Seismology*.
- Sato, H., Fehler, M. C., & Maeda, T. (2012). *Seismic Wave Propagation and Scattering in
  the Heterogeneous Earth*, 2nd ed. Springer. (Source of the ka–kL classification
  diagram this notebook builds its own version of.)
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Bessels = "~0.2.8"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "da0f5d56ff09b918d2e3ac3de0d1f5c24047003e"

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

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

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

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

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
"""

# ╔═╡ Cell order:
# ╟─eb35764e-ab0a-4b36-a09c-3e7b744685b1
# ╠═fcf7c3ea-a7af-4eda-b77a-cd797302aa60
# ╟─67f3ae02-dc84-46ef-b527-6a7ceab9ed82
# ╟─501ec109-c79b-49f0-b3bc-8fee55e8653f
# ╟─5b0963f4-dfa6-4739-9614-dfe60d386e71
# ╟─a54d1f2a-46da-4b76-af72-d5bfcc9a0f58
# ╟─7c40b8e6-f103-4a6e-954b-1a9a0489c6ff
# ╟─2b29f7cb-15ba-4f89-86ea-bd5d94fb6417
# ╟─fe0647ae-9326-493d-9556-e3c8258052d0
# ╟─4d3f3baa-4d61-4701-923c-dd5fd940cf50
# ╟─a4d829be-5596-477b-8329-e64afc9ec003
# ╟─e4c3bf1e-beb5-46d5-b538-faa044c02f72
# ╟─fb54e27f-578d-472c-bc64-0c1946d719eb
# ╟─757357c5-7636-4fd2-85ac-8f8cf78dfe71
# ╟─7e1f4e2e-7ed9-4921-b750-988da8bf2ed7
# ╟─9bc68957-8895-4913-b6a9-c42de31d952a
# ╠═e9baabb0-3e63-402d-9830-4bd71b2ec6d1
# ╠═bd980060-a4ce-4035-bad6-e971d727b603
# ╠═c22db7a4-2b03-42c5-bfc8-2b3f5d1f3ad5
# ╠═a168c575-f6b0-47e9-97a6-7502b991bc8d
# ╟─da37e582-d104-4ecf-be26-339961f735fc
# ╠═3a30e542-560c-47f3-8575-b0896e738cb9
# ╟─41f157c2-ec27-4c07-b1f9-4f5b1d49c790
# ╠═95084277-fe1a-4b8c-bd70-5c7dac8dd277
# ╠═b32c9097-3026-4863-83b9-882d0b5bd236
# ╟─20fd089f-d7e8-413f-a6da-98f091787c28
# ╠═da44d725-33ec-4c52-bf84-1f85ea291e11
# ╠═1d8bf334-c69e-4f90-af0a-28cd3af723bd
# ╟─1b9b4706-65af-46b8-b249-68c6417e2b95
# ╠═6750958f-4d00-4671-b70b-a0c7b59c3058
# ╠═c80bd048-a116-470d-82e9-58f5a2e7759a
# ╟─6875a305-087a-49bb-9da2-114538b98dc8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
