### A Pluto.jl notebook ###
# v0.2.6

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

# ╔═╡ 81b3647e-d741-11f0-b633-478f0d655f0b
using PlutoUI, Interpolations, DelimitedFiles

# ╔═╡ 98f30e60-f1e5-4102-8539-2c6600c3de90
TableOfContents(include_definitions=true)

# ╔═╡ 169e2164-b62f-4b76-bf12-bec1a6aba243
md"""
# What's Inside the Earth?

If you assume the Earth is a uniform ball of average rock, you would guess the pressure and
temperature at the center are large but unremarkable — maybe a few times what you would
feel at the bottom of the ocean, a few hundred degrees hotter than a kitchen oven. Both
guesses are wrong by more than two orders of magnitude: the real center sits near 364 GPa
and roughly 5000-6000 K, hot enough to melt iron and under more pressure than exists
anywhere else in the solar system's rocky bodies.

We cannot drill there — the deepest hole ever bored (the Kola Superdeep) reached about
12 km, barely a scratch on a 6371 km radius. Instead, this notebook builds up pressure and
temperature the way seismologists actually do: from a density profile constrained by how
seismic waves travel through the planet, combined with the physics of hydrostatic
equilibrium. **The interactive panel below lets you explore how temperature — the part we
understand least — changes the picture; the panel after it shows that the Earth's interior
is not even spherically symmetric, using a real seismic tomography model.**
"""

# ╔═╡ 98788e51-80d7-4a50-adf8-f0731d99f2b0
md"""
## Where Density Comes From: Hydrostatic Equilibrium

A parcel of mantle at radius ``r`` is squeezed by the weight of everything above it. In
equilibrium, that squeeze exactly balances the pressure gradient trying to push the parcel
outward:

```math
\frac{dP}{dr} = -\rho(r)\, g(r)
```

Two more pieces close the system. The gravitational acceleration depends on how much mass
is enclosed inside radius ``r``:

```math
g(r) = \frac{G\, M(r)}{r^2}, \qquad M(r) = 4\pi \int_0^r \rho(r')\, r'^2\, dr'
```

Given a density profile ``\rho(r)`` — from PREM, next — these three equations are all it
takes to build up gravity and pressure from the surface to the center. See
`enclosed_mass_at_r`, `g_profile`, and `pressure_profile` in the Appendix for exactly how.
"""

# ╔═╡ 3f6270d2-a8f3-404c-87eb-d544a36cbe29
md"""
## The PREM Reference Model

PREM — the **P**reliminary **R**eference **E**arth **M**odel (Dziewonski & Anderson, 1981)
— is the standard 1-D (depth-only) description of the Earth's density and seismic velocity,
built by simultaneously fitting a huge catalog of earthquake travel times, normal-mode
frequencies, and astronomical constraints (the Earth's total mass and moment of inertia).

The density profile used here is PREM's actual published values at every one of its ~160
radial sample points — parsed directly from the same reference model file this repo's
normal-mode notebook (`specnm.jl`) uses for gravito-elastic mode calculations — linearly
interpolated only *between* those real points, not reduced to a handful of anchors. Every
discontinuity — the Moho, the 410 and 660 km transition-zone boundaries, the core-mantle
boundary (CMB), and the inner-core boundary (ICB) — is a jump where mineral phase or
composition changes abruptly, exactly as PREM specifies it.

!!! note "Crust simplified to a continental column"
	PREM itself models an *oceanic* crustal column (with a water layer on top); this
	notebook substitutes a simpler continental crust for the top ~25 km, since "how thick
	is the ocean" is not the point here.
"""

# ╔═╡ 5a36b51d-da51-4ae2-95d8-227650643530
md"""
## Geotherms, Adiabats, and the Solidus

Density and pressure are almost boring: they increase steadily and predictably with depth
because rock is nearly incompressible and the overburden only grows. Temperature is the
interesting one, because heat moves through the Earth by two different mechanisms:

- **Conduction**, near the surface, where the lithosphere is too rigid and cold to convect.
  Temperature rises roughly *linearly* with depth here — the `T surface` and `gradient`
  sliders above control this.
- **Convection**, through most of the mantle, where hot rock physically rises and cold rock
  sinks. A convecting fluid's temperature profile is close to an **adiabat**: nearly
  constant potential temperature, rising only slowly with depth as pressure compresses it.
  The `potential T` and `adiabatic gradient` sliders control this piece.

Where the geotherm crosses the **solidus** — the temperature at which rock starts to melt —
is where partial melting begins: the mechanism behind mid-ocean-ridge basalts and the
low-viscosity asthenosphere. Watch for that crossing in the zoomed temperature panel above.

!!! danger "The solidus line is real physics, but only shallow"
	The solidus curve comes from a real experimental fit (Hirschmann, 2000), not a toy
	formula — but that fit is only calibrated up to about 10 GPa (roughly 300 km). It
	deliberately stops being drawn past that depth rather than extrapolating into
	meaningless territory — see `peridotite_solidus` in the Appendix.
"""

# ╔═╡ 070fe588-0e0e-4c02-8888-b724670a7f0a
md"""
## Beyond One Dimension

Everything above assumes the Earth is spherically symmetric — density, pressure, and
temperature depend only on depth, never on location. That is a good approximation for
density and pressure, but a poor one for temperature: subducting slabs are colder than
their surroundings, rising plumes are hotter, and both leave a fingerprint in how fast
seismic waves travel through them.

The panel below shows a real global seismic tomography model instead of the 1-D profile
above — but it can also show density, pressure, or seismic velocity from PREM, switched
with the **parameter** menu. Three independent tomography models are on offer there too:
**GyPSuM** (Simmons et al. 2010), **S40RTS** (Ritsema et al. 2011), and **GLAD-M25** (Lei et
al. 2020, a much higher-resolution full-waveform/adjoint model) — each built from different
data and methods, so switching between them is a real look at how much (or how little)
independent seismic inversions agree. **Drag out a line on the globe** — press to drop pin
A, drag to stretch pin B out to where you want it — or fine-tune either pin afterward by
dragging it on its own (right-drag rotates the view): the depth section on the right redraws
continuously along the great circle between them, from the surface down to whatever depth
that parameter actually has data for. For the tomography models that stops just above the
CMB — S-waves don't resolve the liquid outer core — but for PREM's profiles, which don't
vary with location, the section keeps going all the way to the literal centre of the Earth,
drawn as concentric shells instead of a line. Coastlines are drawn on the globe throughout,
purely as a geographic reference so you can tell which anomaly sits under which continent.
The **preset** menu jumps straight to a few well-known features without needing to know
coordinates — mantle plumes (Hawaii, Iceland, Yellowstone) and subducting slabs (the
Farallon/Pacific remnant, Japan, Tonga).
"""

# ╔═╡ 98799827-8c07-4706-9354-faf9ee1a3e15
md"""
## What dVs% Actually Means

The color you are seeing is **dVs**, the percent difference between the shear-wave (S-wave)
speed at a point and the global average at that depth. Positive dVs (blue here) means
seismic waves travel *faster* than average; negative (red) means *slower*.

Fast and slow are usually read as cold and hot — colder rock is stiffer and denser, so
waves move through it faster — and that is the dominant effect. Subducting slabs (fast/
blue) and mantle plumes (slow/red) are the textbook examples; look for both by rotating
the globe.

!!! warning "Fast is not always cold"
	Composition and mineral phase also change seismic speed, independent of temperature —
	a chemically distinct, iron-enriched pile near the CMB can be slow without being
	especially hot. Tomography alone cannot separate these effects; that requires
	combining it with mineral physics and, ideally, other data (density, geoid, mineral
	phase boundaries).
"""

# ╔═╡ 22a88956-3e85-4d61-a3fa-aec2da926fd9
md"""
## Appendix

Everything below is implementation: the density/gravity/pressure pipeline, the temperature
models, the tomography data loading, and finally the two interactive widgets' own code.
"""

# ╔═╡ e2493968-3488-4a46-b0c1-222010c6dcb9
md"### The PREM Density Profile"

# ╔═╡ 67a3d0fe-0a8c-4163-9d79-5e9a511ddab0
begin
    # PREM (Dziewonski & Anderson, 1981, Phys. Earth Planet. Inter. 25(4):297-356) density
    # and velocity profile, parsed directly from the same reference model file specnm.jl
    # uses for normal-mode calculations: src/assets/data/specnm_models/prem_ani, an AXISEM card-deck
    # listing (radius, rho, vpv, vsv, vph, vsh, eta, qka, qmu) at every one of PREM's real
    # radial sample points -- not an independently-typed subset. We use the vertically-
    # polarized vp/vs (vpv/vsv) as a single representative velocity at each depth, ignoring
    # PREM's mild radial anisotropy (vpv != vph in a few layers) -- fine for this
    # notebook's purpose of showing the overall velocity profile, not modeling anisotropy.
    depth_km = Float64[]
    rho = Float64[]
    vp = Float64[]
    vs = Float64[]
    open(joinpath(@__DIR__, "..", "assets", "data", "specnm_models", "prem_ani")) do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#") || !occursin(r"^[0-9.]", s)) && continue
            cols = split(s)
            r_m = parse(Float64, cols[1])
            rho_val = parse(Float64, cols[2])
            vp_val = parse(Float64, cols[3])
            vs_val = parse(Float64, cols[4])
            d_km = (6371.0e3 - r_m) / 1e3
            # Discontinuities repeat the same radius on consecutive rows (a density jump);
            # nudge the depth by a tiny amount so it stays strictly increasing, which
            # LinearInterpolation requires.
            if !isempty(depth_km) && d_km <= depth_km[end]
                d_km = depth_km[end] + 1e-4
            end
            push!(depth_km, d_km)
            push!(rho, rho_val)
            push!(vp, vp_val)
            push!(vs, vs_val)
        end
    end

    depth = depth_km .* 1e3
    R = 6371e3  # Earth's radius (m)

    interp_rho = LinearInterpolation(depth, rho; extrapolation_bc=Line())
    interp_vp = LinearInterpolation(depth, vp; extrapolation_bc=Line())
    interp_vs = LinearInterpolation(depth, vs; extrapolation_bc=Line())

    (depth, rho, vp, vs, interp_rho, interp_vp, interp_vs, R)
end

# ╔═╡ 0defc265-fe5a-42c6-a876-5e5e16cde7fb
md"### Layer 1: Mass, Gravity, and Pressure"

# ╔═╡ f2e9afae-368b-4224-bfc9-7aa681e0d4e7
begin
    Gconst = 6.67430e-11 # gravitational constant, m^3 kg^-1 s^-2

    """
    	enclosed_mass_at_r(depth, rho_interp, R)

    Mass enclosed within radius `r`, ``M(r) = 4\\pi \\int_0^r \\rho(r')\\, r'^2\\, dr'``,
    evaluated on a fine radial grid by cumulative summation. `rho_interp` gives density as
    a function of *depth* below the surface, so it is sampled at `R - r`.

    Returns `(rgrid, mass_cum)`, both sampled on the same fine radial grid.
    """
    function enclosed_mass_at_r(depth, rho_interp, R)
        rgrid = range(0.0, stop=R, length=3000)
        dens = rho_interp.(R .- rgrid)
        integrand = 4π .* rgrid .^ 2 .* dens
        mass_cum = cumsum(integrand) .* (rgrid[2] - rgrid[1])
        return rgrid, mass_cum
    end

    """
    	g_profile(rgrid, mass_cum)

    Gravitational acceleration ``g(r) = G\\,M(r)/r^2`` at every radius in `rgrid`, given the
    enclosed mass `mass_cum` from [`enclosed_mass_at_r`](@ref). `g(0) = 0` by convention —
    there is no enclosed mass at the center to pull on anything.
    """
    function g_profile(rgrid, mass_cum)
        g = similar(rgrid)
        for (i, r) in enumerate(rgrid)
            g[i] = r == 0 ? 0.0 : Gconst * mass_cum[i] / r^2
        end
        return g
    end

    """
    	pressure_profile(rgrid, mass_cum, dens)

    Hydrostatic pressure ``P(r) = \\int_r^{R} \\rho(r')\\, g(r')\\, dr'``, integrated inward
    from zero pressure at the surface using a reverse cumulative trapezoid rule.

    !!! note "Sign convention"
    	The differential form is ``dP/dr = -\\rho(r)\\,g(r)`` — pressure *decreases* outward.
    	Integrating from the surface inward, as this function does, keeps every partial sum
    	positive, which is why the loop below runs from the last index down to the first
    	rather than the more obvious-looking forward direction.
    """
    function pressure_profile(rgrid, mass_cum, dens)
        dr = rgrid[2] - rgrid[1]
        g = g_profile(rgrid, mass_cum)
        integrand = dens .* g
        P = zeros(length(rgrid))
        for i in length(rgrid):-1:1
            P[i] = i == length(rgrid) ? 0.0 : P[i+1] + 0.5 * (integrand[i] + integrand[i+1]) * dr
        end
        return P
    end

    rgrid, mass_cum = enclosed_mass_at_r(depth, interp_rho, R)
    dens_grid = interp_rho.(R .- rgrid)
    ggrid = g_profile(rgrid, mass_cum)
    Pgrid = pressure_profile(rgrid, mass_cum, dens_grid)

    # Present results referenced to depth from the surface, not radius from the center.
    depth_grid = R .- rgrid

    # depth_grid runs center -> surface (increasing r, decreasing depth); LinearInterpolation
    # needs its knots strictly increasing, so reverse both arrays together.
    interp_pressure = LinearInterpolation(reverse(depth_grid), reverse(Pgrid); extrapolation_bc=Line())

    (depth_grid, dens_grid, ggrid, Pgrid, interp_pressure)
end

# ╔═╡ bd4aff4f-a4af-4d01-967f-02d7acf2a098
md"#### Validating hydrostatic equilibrium"

# ╔═╡ 259dab7b-4d62-4a5a-aaee-de0de4fd86fb
let
    # ### Validating hydrostatic equilibrium
    # Three independent sanity checks on the mass/gravity/pressure pipeline above, each
    # against a number that doesn't depend on this notebook's code at all.
    g_surface = ggrid[end]
    total_mass = mass_cum[end]
    earth_mass_true = 5.972e24  # kg, independent reference value
    # rgrid runs center (index 1) -> surface (index end), so pressure should DECREASE
    # monotonically along it, reaching zero at the surface.
    monotonic_out = all(diff(Pgrid) .<= 0)

    md"""
    - surface gravity: **$(round(g_surface, digits=2)) m/s²** (true value ≈ 9.8 m/s²)
    - total enclosed mass: **$(round(total_mass, sigdigits=4)) kg** (true value ≈ $(earth_mass_true) kg,
      $(round(100*(total_mass-earth_mass_true)/earth_mass_true, digits=1))% off)
    - pressure decreases monotonically from center to surface: **$(monotonic_out)**
    """
end

# ╔═╡ 7e1b745c-5174-4e36-90a0-51db489d4768
md"### Layer 2: Geotherms, Adiabats, and the Solidus"

# ╔═╡ 306d2112-cab7-4fdc-b441-e9228e243cdc
begin
    """
    	linear_geotherm(depth_grid; T_surface=300.0, gradient=0.5)

    A conductive geotherm: temperature increases linearly with depth from `T_surface` (K)
    at `gradient` K/km. A reasonable model for the lithosphere, where heat moves by
    conduction — see [`adiabatic_profile`](@ref) for the convecting interior, where it
    doesn't.
    """
    function linear_geotherm(depth_grid; T_surface=300.0, gradient=0.5)
        return T_surface .+ gradient .* (depth_grid ./ 1e3)
    end

    """
    	adiabatic_profile(depth_grid; T_potential=1600.0, adi_grad=0.3)

    A simple adiabatic temperature profile for the convecting mantle: `T_potential` (K) is
    the temperature the mantle would have if brought adiabatically to the surface, and
    `adi_grad` (K/km) is the small, compression-driven increase with depth along the
    adiabat.
    """
    function adiabatic_profile(depth_grid; T_potential=1600.0, adi_grad=0.3)
        return T_potential .+ adi_grad .* (depth_grid ./ 1e3)
    end

    # Hirschmann (2000), "Mantle solidus: Experimental constraints and the effects of
    # peridotite composition," Geochem. Geophys. Geosyst. 1(10), 1042, doi:10.1029/2000GC000070.
    # The fit is calibrated against melting experiments only up to this pressure.
    const P_MAX_GPA = 10.0

    """
    	peridotite_solidus(P_pa)

    Anhydrous peridotite solidus temperature (K), from Hirschmann (2000): a quadratic fit
    in pressure, ``T(\\mathrm{°C}) = 1120.661 + 132.899\\,P - 5.104\\,P^2`` with `P` in GPa,
    calibrated against melting experiments up to about 10 GPa (`P_MAX_GPA`).

    !!! warning "Only valid below `P_MAX_GPA`"
    	Evaluated far outside its calibration range — e.g. at core-mantle-boundary pressure
    	— this fit produces physically meaningless numbers; it is not a model of the whole
    	mantle's melting behavior, only of shallow peridotite. Above `P_MAX_GPA`, this
    	function returns `missing` rather than silently extrapolating, so the solidus curve
    	simply stops being drawn once it leaves the range the experiments actually cover.
    """
    function peridotite_solidus(P_pa)
        P_gpa = P_pa ./ 1e9
        return map(P_gpa) do p
            p > P_MAX_GPA ? missing : 1120.661 + 132.899 * p - 5.104 * p^2 + 273.15
        end
    end

    (linear_geotherm, adiabatic_profile, peridotite_solidus)
end

# ╔═╡ 0fad939d-cb1e-43bf-80d1-bdd8957515f9
md"### Layer 3: The Tomography Grid"

# ╔═╡ 1faacc81-f508-4832-9a31-148aa3e3c7cc
begin
    """
        load_tomo_csv(path)

    Reads a `lon,lat,depth_km,dvs_pct` CSV (see each `data/*_dvs_coarse.csv` file for that
    model's full provenance) into `(lons, lats, depths, grid)`, where `grid` is one
    `(nlat, nlon)` matrix per depth layer -- the shape every tomography model in this
    notebook shares once loaded, regardless of each model's own native resolution.
    """
    function load_tomo_csv(path)
        raw, _header = readdlm(path, ','; comments=true, comment_char='#', header=true)
        lons = sort(unique(raw[:, 1]))
        lats = sort(unique(raw[:, 2]))
        depths = sort(unique(raw[:, 3]))
        nlo, nla = length(lons), length(lats)
        lon_idx = Dict(v => i for (i, v) in enumerate(lons))
        lat_idx = Dict(v => i for (i, v) in enumerate(lats))
        dep_idx = Dict(v => i for (i, v) in enumerate(depths))
        grid = [zeros(nla, nlo) for _ in depths]
        for row in eachrow(raw)
            lo, la, de, v = row
            grid[dep_idx[de]][lat_idx[la], lon_idx[lo]] = v
        end
        (lons, lats, depths, grid)
    end

    # GyPSuM (Simmons, Forte, Boschi & Grand, 2010, doi:10.1029/2010JB007631), downloaded
    # from IRIS EMC. All 100 native depth layers kept; lat/lon coarsened to 3x3 deg (still
    # finer than the globe's own render mesh) to keep the embedded payload a few MB instead
    # of tens of MB. See ../assets/data/gypsum_dvs_coarse.csv for full provenance.
    tomo_lons, tomo_lats, tomo_depths, tomo_grid =
        load_tomo_csv(joinpath(@__DIR__, "..", "assets", "data", "gypsum_dvs_coarse.csv"))

    # S40RTS (Ritsema, Deuss, van Heijst & Woodhouse, 2011,
    # doi:10.1111/j.1365-246X.2010.04884.x). See ../assets/data/s40rts_dvs_coarse.csv for provenance.
    s40rts_lons, s40rts_lats, s40rts_depths, s40rts_grid =
        load_tomo_csv(joinpath(@__DIR__, "..", "assets", "data", "s40rts_dvs_coarse.csv"))

    # GLAD-M25 (Lei et al., 2020, doi:10.1093/gji/ggaa253), a full-waveform/adjoint-tomography
    # model -- higher resolution than GyPSuM/S40RTS, since it comes from a much finer
    # spectral-element inversion rather than a truncated spherical-harmonic expansion. See
    # ../assets/data/gladm25_dvs_coarse.csv for full provenance, including how dVs% was derived
    # (GLAD-M25 ships absolute Vsv, not a PREM-relative perturbation like GyPSuM/S40RTS) and
    # why its usable depth range starts at 80 km rather than the surface (a real data-fill
    # bug in the shallowest layers, plus crustal structure sharp enough to blow out this
    # widget's mantle-focused color scale).
    gladm25_lons, gladm25_lats, gladm25_depths, gladm25_grid =
        load_tomo_csv(joinpath(@__DIR__, "..", "assets", "data", "gladm25_dvs_coarse.csv"))

    (tomo_grid, s40rts_grid, gladm25_grid)
end

# ╔═╡ b8d8c763-ccf3-4c29-af91-3959a0b903eb
let
    function checks(lons, lats, depths, grid)
        shape_ok = all(size(M) == (length(lats), length(lons)) for M in grid)
        range_ok = maximum(m -> maximum(abs, m), grid) < 10.0
        depths_sorted = issorted(depths)
        shape_ok && range_ok && depths_sorted
    end

    gypsum_ok = checks(tomo_lons, tomo_lats, tomo_depths, tomo_grid)
    s40rts_ok = checks(s40rts_lons, s40rts_lats, s40rts_depths, s40rts_grid)
    gladm25_ok = checks(gladm25_lons, gladm25_lats, gladm25_depths, gladm25_grid)

    md"""
    #### Verifying the loaded grids
    Each check: every depth layer shaped `(nlat, nlon)`, `|dVs| < 10%` everywhere, depths
    sorted increasing.
    - GyPSuM: **$(gypsum_ok)**
    - S40RTS: **$(s40rts_ok)**
    - GLAD-M25: **$(gladm25_ok)**
    """
end

# ╔═╡ 4bef3d37-3bfd-4c7f-83c6-31b4318f9e5a
md"### Layer 4: PREM Profiles on the Shared Depth Axis"

# ╔═╡ 7cb6d994-2002-44f2-bbc7-e73f146d1704
begin
    # The globe widget lets you pick between GyPSuM (laterally-varying, one value per
    # lat/lon/depth) and PREM's density/velocity/pressure (radially symmetric -- one value
    # per depth, the same everywhere on a given shell). To drive both through the same
    # depth slider and the same cut-face rendering code, sample the PREM/pressure profiles
    # onto GyPSuM's own depth axis (`tomo_depths`, km) rather than keeping them on their
    # native finer grid.
    tomo_depths_m = tomo_depths .* 1e3
    prem_rho_at_tomo_depths = interp_rho.(tomo_depths_m)                    # kg/m^3
    prem_vp_at_tomo_depths = interp_vp.(tomo_depths_m) ./ 1e3               # km/s
    prem_vs_at_tomo_depths = interp_vs.(tomo_depths_m) ./ 1e3               # km/s
    prem_pressure_at_tomo_depths = interp_pressure.(tomo_depths_m) ./ 1e9   # GPa

    (prem_rho_at_tomo_depths, prem_vp_at_tomo_depths, prem_vs_at_tomo_depths, prem_pressure_at_tomo_depths)
end

# ╔═╡ 5c8a2f14-6e91-4b3a-9d27-8f41c0a9e6b2
md"### Layer 4b: PREM Profiles at Full Depth Resolution (Surface to Centre)"

# ╔═╡ d3f7b920-1a44-4e6c-8b19-2c5d9a7f0e33
begin
    # The cross-section fan (in the globe widget below) lets a PREM parameter's fill reach
    # the literal centre of the Earth -- unlike GyPSuM, which only has data to the CMB --
    # so PREM needs its own depth axis running the full 0-6371 km, independent of
    # `tomo_depths` (Layer 4 above, which stays confined to GyPSuM's shallower range since
    # it also drives the globe's own depth slider). 300 points matches the resolution
    # already used for `EarthProfileInput`'s own depth axis (`epi_depth_km`) elsewhere in
    # this notebook.
    fan_depths_km = collect(range(0, R / 1e3, length=300))
    fan_depths_m = fan_depths_km .* 1e3
    prem_fan_rho = interp_rho.(fan_depths_m)                    # kg/m^3
    prem_fan_vp = interp_vp.(fan_depths_m) ./ 1e3                # km/s
    prem_fan_vs = interp_vs.(fan_depths_m) ./ 1e3                # km/s
    prem_fan_pressure = interp_pressure.(fan_depths_m) ./ 1e9    # GPa

    (fan_depths_km, prem_fan_rho, prem_fan_vp, prem_fan_vs, prem_fan_pressure)
end

# ╔═╡ 3ac914f3-56e4-4025-81d1-9461a0bfc835
md"### Layer 5: Coastlines"

# ╔═╡ 9bebe536-3dc8-4e45-b27d-160dc0b4f608
begin
    # Natural Earth 1:110m coastline (public domain), vendored offline from GeoMakie.jl's
    # bundled assets -- see ../assets/data/coastlines_110m.csv for provenance. Drawn on the globe
    # widget's outer shell purely as a geographic reference layer (it carries no data of
    # its own), so a viewer can tell "that fast anomaly is under East Asia" at a glance
    # instead of having to mentally project bare lat/lon onto a rotating sphere.
    coast_raw, _coast_header = readdlm(joinpath(@__DIR__, "..", "assets", "data", "coastlines_110m.csv"), ',';
        comments=true, comment_char='#', header=true)
    coast_line_id = Int.(coast_raw[:, 1])
    coast_lon = Float64.(coast_raw[:, 2])
    coast_lat = Float64.(coast_raw[:, 3])

    (coast_line_id, coast_lon, coast_lat)
end

# ╔═╡ 8438fcce-2d66-4693-b747-641a4c2264b7
md"### The Interactive Widgets"

# ╔═╡ a41bd2f5-bce3-4a30-8f31-5b86c19e59fb
begin
    # Downsample the fine (3000-point) mass/gravity/pressure grid to something light enough
    # to bake directly into the profile widget below, and reorder to run surface -> center
    # (increasing depth) to match how the canvas draws depth top-to-bottom.
    epi_idx = round.(Int, range(1, length(depth_grid), length=300))
    epi_depth_km = reverse(depth_grid[epi_idx]) ./ 1e3
    epi_dens = reverse(dens_grid[epi_idx])
    epi_press_gpa = reverse(Pgrid[epi_idx]) ./ 1e9

    (epi_depth_km, epi_dens, epi_press_gpa)
end

# ╔═╡ a916b96c-ad21-4cd1-ba7e-6434ce2e5396
begin
    struct EarthProfileInput
        T_surface::Float64
        geotherm_gradient::Float64
        potential_T::Float64
        adi_grad::Float64
        show_solidus::Bool
    end

    function EarthProfileInput(; T_surface=300.0, geotherm_gradient=0.5, potential_T=1600.0,
        adi_grad=0.3, show_solidus=true)
        EarthProfileInput(Float64(T_surface), Float64(geotherm_gradient), Float64(potential_T),
            Float64(adi_grad), Bool(show_solidus))
    end

    Base.get(w::EarthProfileInput) = Dict{String,Any}(
        "T_surface" => w.T_surface,
        "geotherm_gradient" => w.geotherm_gradient,
        "potential_T" => w.potential_T,
        "adi_grad" => w.adi_grad,
        "show_solidus" => w.show_solidus,
    )

    function Base.show(io::IO, ::MIME"text/html", w::EarthProfileInput)
        write(io, """
<div id="epiwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#epiwidget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #epiwidget .epi-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #epiwidget .epi-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #epiwidget .epi-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #epiwidget .epi-workspace{display:flex;gap:12px;align-items:flex-start;justify-content:center;width:100%}
    #epiwidget canvas{background:#000;border:1px solid #374151;border-radius:6px;display:block}
    #epiwidget .epi-controls{width:min(var(--totalw,960px),100%);margin-top:8px;display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:8px;font:14px sans-serif}
    #epiwidget .epi-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
    #epiwidget .epi-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:20px}
    #epiwidget .epi-control-row{display:grid;grid-template-columns:minmax(90px,150px) minmax(70px,1fr) minmax(40px,64px);gap:6px;align-items:center;margin:7px 0}
    #epiwidget .epi-control-row input[type=range]{width:100%;min-width:0;vertical-align:middle}
    #epiwidget .epi-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums}
    #epiwidget .epi-actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    #epiwidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px}
    #epiwidget label{color:#d1d5db}
    @media (max-width: 980px){
      #epiwidget .epi-workspace{flex-direction:column;align-items:center}
      #epiwidget .epi-controls{grid-template-columns:1fr;width:660px;max-width:100%}
    }
  </style>
  <div class="epi-title">
    <div class="epi-title-desc">Density and pressure grow steadily toward the center, but temperature depends on how heat moves &mdash; conduction near the surface, convection below.</div>
    <div class="epi-title-hint">drag the geotherm/adiabat sliders &middot; watch where the geotherm crosses the solidus in the zoomed temperature panel</div>
  </div>
  <div class="epi-workspace">
    <canvas id="denscvs" width="220" height="620"></canvas>
    <canvas id="presscvs" width="220" height="620"></canvas>
    <canvas id="tempcvs" width="260" height="620"></canvas>
  </div>
  <div class="epi-controls">
    <div class="epi-control-group">
      <div class="epi-control-title">Geotherm</div>
      <label class="epi-control-row"><span>T surface</span><input type="range" id="tsurf" min="200" max="600" step="10" value="$(w.T_surface)"><span id="tsurfv" class="epi-value">$(round(Int, w.T_surface)) K</span></label>
      <label class="epi-control-row"><span>gradient</span><input type="range" id="ggrad" min="0.1" max="1.5" step="0.1" value="$(w.geotherm_gradient)"><span id="ggradv" class="epi-value">$(round(w.geotherm_gradient, digits=1)) K/km</span></label>
    </div>
    <div class="epi-control-group">
      <div class="epi-control-title">Adiabat</div>
      <label class="epi-control-row"><span>potential T</span><input type="range" id="ptemp" min="1200" max="1900" step="50" value="$(w.potential_T)"><span id="ptempv" class="epi-value">$(round(Int, w.potential_T)) K</span></label>
      <label class="epi-control-row"><span>adi. gradient</span><input type="range" id="agrad" min="0.1" max="0.6" step="0.05" value="$(w.adi_grad)"><span id="agradv" class="epi-value">$(round(w.adi_grad, digits=2)) K/km</span></label>
    </div>
    <div class="epi-control-group">
      <div class="epi-control-title">Display</div>
      <div class="epi-actions">
        <label><input type="checkbox" id="showsol" $(w.show_solidus ? "checked" : "") style="vertical-align:middle"> show <b>peridotite solidus</b></label>
      </div>
      <div style="margin-top:8px;color:#9ca3af;font-size:13px">
        the solidus line stops where its calibration (&lt;10 GPa) ends &mdash; see the Appendix for why
      </div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8)
  const heightBudget = Math.max(420, window.innerHeight - 420)
  const H = Math.round(Math.min(heightBudget, 700))
  const totalW = Math.max(700, availW)
  const W1 = Math.round(totalW*0.30), W2 = Math.round(totalW*0.30), W3 = Math.round(totalW*0.34)
  par.style.setProperty('--totalw', Math.round(totalW) + 'px')

  const DPR = Math.min(window.devicePixelRatio || 1, 2)
  const DEPTH_KM = $(epi_depth_km)
  const DENS = $(epi_dens)
  const PRESS = $(epi_press_gpa)
  const MAX_DEPTH = DEPTH_KM[DEPTH_KM.length-1]
  const TEMP_MAX_DEPTH = 410

  function hidpi(cv, w, h){
    cv.width = Math.round(w*DPR); cv.height = Math.round(h*DPR)
    cv.style.width = w+'px'; cv.style.height = h+'px'
    const cx = cv.getContext('2d'); cx.setTransform(DPR,0,0,DPR,0,0); return cx
  }
  const dcvs = par.querySelector('#denscvs'), dctx = hidpi(dcvs, W1, H)
  const pcvs = par.querySelector('#presscvs'), pctx = hidpi(pcvs, W2, H)
  const tcvs = par.querySelector('#tempcvs'), tctx = hidpi(tcvs, W3, H)

  let T_surface = $(w.T_surface), geothermGradient = $(w.geotherm_gradient)
  let potentialT = $(w.potential_T), adiGrad = $(w.adi_grad)
  let showSolidus = $(w.show_solidus)

  const PAD_L = 46, PAD_R = 12, PAD_T = 30, PAD_B = 10

  function drawPanel(ctx, w, h, xmax, maxDepth, title, series, xlabel){
    ctx.fillStyle = '#000'; ctx.fillRect(0,0,w,h)
    const x0 = PAD_L, x1 = w-PAD_R, y0 = PAD_T, y1 = h-PAD_B
    const X = v => x0 + (v/xmax)*(x1-x0)
    const Y = d => y0 + (d/maxDepth)*(y1-y0)
    ctx.strokeStyle = '#1f2937'; ctx.fillStyle = '#9ca3af'; ctx.font='11px sans-serif'; ctx.textAlign='right'
    const nTicks = 6
    for(let i=0;i<=nTicks;i++){
      const d = maxDepth*i/nTicks
      const y = Y(d)
      ctx.beginPath(); ctx.moveTo(x0,y); ctx.lineTo(x1,y); ctx.stroke()
      ctx.fillText(Math.round(d)+'', x0-6, y+4)
    }
    ctx.strokeStyle = '#4b5563'; ctx.beginPath(); ctx.moveTo(x0,y0); ctx.lineTo(x0,y1); ctx.stroke()
    ctx.textAlign='center'; ctx.fillStyle='#e5e7eb'; ctx.font='13px sans-serif'
    ctx.fillText(title, (x0+x1)/2, 16)
    ctx.save(); ctx.translate(12, (y0+y1)/2); ctx.rotate(-Math.PI/2)
    ctx.fillStyle='#9ca3af'; ctx.font='11px sans-serif'; ctx.textAlign='center'
    ctx.fillText('depth (km)', 0, 0)
    ctx.restore()
    ctx.fillStyle='#9ca3af'; ctx.textAlign='center'; ctx.font='11px sans-serif'
    ctx.fillText(xlabel, (x0+x1)/2, h-2)
    for(const s of series){
      ctx.strokeStyle = s.color; ctx.lineWidth = 2
      if(s.dash) ctx.setLineDash(s.dash); else ctx.setLineDash([])
      ctx.beginPath()
      let started = false
      for(let i=0;i<s.depth.length;i++){
        const dv = s.depth[i], xv = s.x[i]
        if(xv === null || xv === undefined || Number.isNaN(xv) || dv > maxDepth) continue
        const px = X(xv), py = Y(dv)
        if(!started){ ctx.moveTo(px,py); started=true } else ctx.lineTo(px,py)
      }
      ctx.stroke()
    }
    ctx.setLineDash([])
  }

  function redraw(){
    drawPanel(dctx, W1, H, Math.max(...DENS)*1.05, MAX_DEPTH, 'Density', [{depth:DEPTH_KM, x:DENS, color:'#38bdf8'}], 'kg/m³')
    drawPanel(pctx, W2, H, Math.max(...PRESS)*1.05, MAX_DEPTH, 'Pressure', [{depth:DEPTH_KM, x:PRESS, color:'#f59e0b'}], 'GPa')

    const n = 120
    const tdepth = Array.from({length:n}, (_,i)=> TEMP_MAX_DEPTH*i/(n-1))
    const geo = tdepth.map(d => T_surface + geothermGradient*d)
    const adi = tdepth.map(d => potentialT + adiGrad*d)
    // Interpolate pressure(depth) from the fixed PRESS/DEPTH_KM arrays -- mirrors
    // peridotite_solidus in the Appendix; keep the two in sync if either changes.
    function pressAt(d){
      for(let i=1;i<DEPTH_KM.length;i++){
        if(DEPTH_KM[i] >= d){
          const t = (d-DEPTH_KM[i-1])/(DEPTH_KM[i]-DEPTH_KM[i-1])
          return PRESS[i-1] + t*(PRESS[i]-PRESS[i-1])
        }
      }
      return PRESS[PRESS.length-1]
    }
    const P_MAX_GPA = 10.0
    const sol = tdepth.map(d => {
      const p = pressAt(d)
      return p > P_MAX_GPA ? null : 1120.661 + 132.899*p - 5.104*p*p + 273.15
    })
    const tempSeries = [
      {depth: tdepth, x: geo, color: '#38bdf8'},
      {depth: tdepth, x: adi, color: '#f59e0b'},
    ]
    if(showSolidus) tempSeries.push({depth: tdepth, x: sol, color: '#22c55e', dash:[5,3]})
    const allT = geo.concat(adi).concat(showSolidus ? sol.filter(v=>v!==null) : [])
    drawPanel(tctx, W3, H, Math.max(...allT)*1.05, TEMP_MAX_DEPTH, 'Temperature (zoomed to 410 km)', tempSeries, 'K')

    tctx.font='11px sans-serif'; tctx.textAlign='left'
    tctx.fillStyle='#38bdf8'; tctx.fillText('geotherm', PAD_L, H-2-24)
    tctx.fillStyle='#f59e0b'; tctx.fillText('adiabat', PAD_L, H-2-12)
    if(showSolidus){ tctx.fillStyle='#22c55e'; tctx.fillText('solidus', PAD_L+70, H-2-24) }
  }

  redraw()

  function bindSlider(id, vid, fmt, setter){
    par.querySelector('#'+id).addEventListener('input', e=>{
      const v = parseFloat(e.target.value)
      setter(v)
      par.querySelector('#'+vid).textContent = fmt(v)
      redraw()
    })
  }
  bindSlider('tsurf','tsurfv', v=>Math.round(v)+' K', v=>T_surface=v)
  bindSlider('ggrad','ggradv', v=>v.toFixed(1)+' K/km', v=>geothermGradient=v)
  bindSlider('ptemp','ptempv', v=>Math.round(v)+' K', v=>potentialT=v)
  bindSlider('agrad','agradv', v=>v.toFixed(2)+' K/km', v=>adiGrad=v)
  par.querySelector('#showsol').addEventListener('change', e=>{ showSolidus = e.target.checked; redraw() })
</script>
""")
    end

    # Pluto's static dependency analysis doesn't reliably detect a `@bind x SomeType()`
    # expression as depending on the cell that defines `SomeType` -- on a fresh restart
    # it can run the bind cell before this one, throwing `UndefVarError`. Referencing
    # this flag as a bare statement at the top of the bind cell (see below) forces an
    # explicit dependency edge, guaranteeing this cell runs first. Same pattern as
    # `_geoid_canvas_ready` in geoid-kernel.jl.
    const _epi_ready = true
end

# ╔═╡ c8266e03-23be-43ea-9292-97359ab7aa41
begin
    _epi_ready
    @bind _epi EarthProfileInput()
end

# ╔═╡ 81868aa1-b26e-4aea-ab9a-2fc2a9e522b0
begin
    """
        TomographyGlobeInput(; depth_idx=2, param="gypsum", lonA=-155.5, latA=19.5, lonB=-155.5, latB=4.5)

    A globe paired with a great-circle depth-section "fan": drag out a line on the globe's
    surface (press to drop **A**, drag to stretch **B** out live) and the fan panel draws a
    continuous cross-section of the selected quantity along the great circle between them,
    from the surface down to either the active tomography model's data limit (~2900 km --
    S-waves don't resolve the liquid outer core) or, for PREM's radially-symmetric profiles,
    all the way to Earth's centre. Either pin can also be dragged individually to fine-tune
    it, and right-drag rotates the view. Three independent whole-mantle shear-velocity
    models are selectable from the same **parameter** menu -- GyPSuM (Simmons et al. 2010),
    S40RTS (Ritsema et al. 2011), and GLAD-M25 (Lei et al. 2020, a full-waveform/adjoint
    tomography model, higher resolution than the other two) -- each on its own native
    lon/lat/depth grid, so switching between them is a real comparison of independent
    seismic inversions, not just a re-color of the same data. Presets jump straight to
    well-known features -- mantle plumes (Hawaii, Iceland, Yellowstone) and subducting slabs
    (the Farallon/Pacific remnant, Japan, Tonga) -- without needing to know coordinates.
    The globe's outer shell still shows the selected depth's value everywhere (via the depth
    slider) -- cheap to redraw every frame since it only ever samples that one depth -- and
    coastlines are drawn as a geographic reference throughout.
    """
    struct TomographyGlobeInput
        depth_idx::Int
        param::String
        lonA::Float64
        latA::Float64
        lonB::Float64
        latB::Float64
    end

    TomographyGlobeInput(; depth_idx=2, param="gypsum", lonA=-155.5, latA=19.5, lonB=-155.5, latB=4.5) =
        TomographyGlobeInput(Int(depth_idx), param, Float64(lonA), Float64(latA), Float64(lonB), Float64(latB))

    Base.get(w::TomographyGlobeInput) = Dict{String,Any}(
        "depth_idx" => w.depth_idx,
        "param" => w.param,
        "lonA" => w.lonA, "latA" => w.latA,
        "lonB" => w.lonB, "latB" => w.latB,
    )

    # Compact [[[dvs,...],...],...] literal: layer -> lat-row -> lon-value. Built by hand
    # rather than relying on Julia's default array `show`, which uses `;`-separated rows
    # for a Matrix -- not valid JS array syntax.
    function tomo_js_literal(grid)
        "[" * join(["[" * join(["[" * join(round.(row, digits=2), ",") * "]" for row in eachrow(M)], ",") * "]" for M in grid], ",") * "]"
    end

    """
        coast_js_literal(line_id, lon, lat)

    Groups the flat `(line_id, lon, lat)` coastline arrays into one `[[lon,lat],...]`
    array per polyline and renders that as a `[[[lon,lat],...],...]` JS literal, so the
    widget's JS side can iterate "one polyline at a time" without redoing the grouping.
    """
    function coast_js_literal(line_id, lon, lat)
        segs = Dict{Int,Vector{Tuple{Float64,Float64}}}()
        for i in eachindex(line_id)
            push!(get!(() -> Tuple{Float64,Float64}[], segs, line_id[i]), (lon[i], lat[i]))
        end
        ids = sort(collect(keys(segs)))
        "[" * join(["[" * join(["[$(lo),$(la)]" for (lo, la) in segs[id]], ",") * "]" for id in ids], ",") * "]"
    end

    function Base.show(io::IO, ::MIME"text/html", w::TomographyGlobeInput)
        write(io, """
<div id="tgiwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#tgiwidget){width:min(80vw,1500px)!important;margin-left:calc((100% - min(80vw,1500px))/2)!important}
    #tgiwidget .tgi-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #tgiwidget .tgi-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #tgiwidget .tgi-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #tgiwidget .tgi-panels{display:flex;gap:14px;align-items:flex-start;justify-content:center;width:100%}
    #tgiwidget .tgi-col-globe, #tgiwidget .tgi-col-fan{flex:1 1 0;min-width:0;display:flex;flex-direction:column;align-items:center}
    #tgiwidget .tgi-panel-label{font-size:12px;color:#6b7280;margin-bottom:4px;text-transform:uppercase;letter-spacing:.04em}
    #tgiwidget .tgi-panel-info{font-size:12px;color:#9ca3af;line-height:1.4;text-align:center;min-height:34px}
    #tgiwidget canvas{cursor:crosshair;background:#000;border:1px solid #374151;border-radius:6px;display:block;max-width:100%}
    #tgiwidget .tgi-controls{display:flex;gap:10px;flex-wrap:wrap;width:100%;margin-top:12px}
    #tgiwidget .tgi-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px;flex:1 1 200px;min-width:200px}
    #tgiwidget .tgi-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #tgiwidget .tgi-control-row{display:grid;grid-template-columns:minmax(60px,90px) minmax(60px,1fr) minmax(40px,80px);gap:6px;align-items:center;margin:7px 0}
    #tgiwidget .tgi-control-row input[type=range]{width:100%;min-width:0;vertical-align:middle}
    #tgiwidget .tgi-value{color:#d1d5db;text-align:left;font-variant-numeric:tabular-nums}
    #tgiwidget label{color:#d1d5db}
    #tgiwidget select{width:100%;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:5px 6px;font-size:14px}
    #tgiwidget .tgi-coords{font-size:13px;line-height:1.9}
    #tgiwidget .tgi-cbar-label{font-size:11px;color:#9ca3af;margin:6px 0 2px}
    #tgiwidget .tgi-cbar{cursor:default;border:none;background:transparent}
    @media (max-width: 900px){
      #tgiwidget .tgi-panels{flex-direction:column;align-items:center}
      #tgiwidget .tgi-col-globe, #tgiwidget .tgi-col-fan{width:100%;max-width:640px}
    }
  </style>
  <div class="tgi-title">
    <div class="tgi-title-desc">Earth's interior isn't spherically symmetric &mdash; mantle plumes and slabs show up as lateral structure that a 1-D depth profile alone can't capture.</div>
    <div class="tgi-title-hint">drag empty globe to draw a new section &middot; drag a pin (A or B) to move it &middot; right-drag to rotate</div>
  </div>
  <div class="tgi-panels">
    <div class="tgi-col-globe">
      <div class="tgi-panel-label">Map view, at the selected depth</div>
      <div class="tgi-panel-info" id="tgGlobeInfo"></div>
      <canvas id="tgcvs" width="600" height="600"></canvas>
    </div>
    <div class="tgi-col-fan">
      <div class="tgi-panel-label">Depth cross-section, A &rarr; B</div>
      <div class="tgi-panel-info" id="tgFanInfo"></div>
      <canvas id="tgfan" width="600" height="600"></canvas>
    </div>
  </div>
  <div class="tgi-controls">
    <div class="tgi-control-group">
      <div class="tgi-control-title">Parameter</div>
      <select id="tgparam">
        <option value="gypsum" $(w.param == "gypsum" ? "selected" : "")>GyPSuM (dVs %)</option>
        <option value="s40rts" $(w.param == "s40rts" ? "selected" : "")>S40RTS (dVs %)</option>
        <option value="gladm25" $(w.param == "gladm25" ? "selected" : "")>GLAD-M25 (dVs %)</option>
        <option value="density" $(w.param == "density" ? "selected" : "")>PREM density</option>
        <option value="vp" $(w.param == "vp" ? "selected" : "")>PREM Vp</option>
        <option value="vs" $(w.param == "vs" ? "selected" : "")>PREM Vs</option>
        <option value="pressure" $(w.param == "pressure" ? "selected" : "")>Pressure</option>
      </select>
    </div>
    <div class="tgi-control-group">
      <div class="tgi-control-title">Preset</div>
      <select id="tgpreset">
        <option value="custom">Custom</option>
        <option value="hawaii" selected>Hawaii plume</option>
        <option value="iceland">Iceland plume</option>
        <option value="yellowstone">Yellowstone plume</option>
        <option value="pacific">Farallon (Pacific) slab</option>
        <option value="japan">Japan slab</option>
        <option value="tonga">Tonga slab</option>
      </select>
    </div>
    <div class="tgi-control-group">
      <div class="tgi-control-title">Depth (globe shell)</div>
      <label class="tgi-control-row"><span>depth</span><input type="range" id="tgdepth" min="0" max="$(length(tomo_depths)-1)" step="1" value="$(w.depth_idx)"><span id="tgdepthv" class="tgi-value">$(round(Int, tomo_depths[w.depth_idx+1])) km</span></label>
    </div>
    <div class="tgi-control-group">
      <div class="tgi-control-title">Section endpoints</div>
      <div class="tgi-coords"><b style="color:#facc15">A</b> <span id="tgAll"></span><br><b style="color:#f472b6">B</b> <span id="tgBll"></span></div>
    </div>
    <div class="tgi-control-group" style="flex:1 1 260px">
      <div class="tgi-control-title">Color scale</div>
      <div class="tgi-cbar-label" id="tgcbarLabel1">surface, this depth</div>
      <canvas id="tgcbar1" class="tgi-cbar" width="240" height="34"></canvas>
      <div class="tgi-cbar-label" id="tgcbarLabel2">section, all depths</div>
      <canvas id="tgcbar2" class="tgi-cbar" width="240" height="34"></canvas>
      <div id="tglegend" style="font-size:13px;line-height:1.6;margin-top:8px"></div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8)
  const heightBudget = Math.max(360, window.innerHeight - 460)
  const panelW = Math.max(300, Math.round((availW - 14)/2))
  const SEC = Math.round(Math.min(panelW, heightBudget, 600))

  const DPR = Math.min(window.devicePixelRatio || 1, 2)
  const R_OUT = Math.round(SEC*0.42)
  const CX = SEC/2, CY = SEC/2
  const R_EARTH_KM = 6371, CMB_KM = 2891, ICB_KM = 5150

  const LONS_GYPSUM = $(tomo_lons), LATS_GYPSUM = $(tomo_lats), DEPTHS_GYPSUM = $(tomo_depths)
  const TOMO_GYPSUM = $(tomo_js_literal(tomo_grid))
  const LONS_S40RTS = $(s40rts_lons), LATS_S40RTS = $(s40rts_lats), DEPTHS_S40RTS = $(s40rts_depths)
  const TOMO_S40RTS = $(tomo_js_literal(s40rts_grid))
  const LONS_GLADM25 = $(gladm25_lons), LATS_GLADM25 = $(gladm25_lats), DEPTHS_GLADM25 = $(gladm25_depths)
  const TOMO_GLADM25 = $(tomo_js_literal(gladm25_grid))

  const PREM_RHO = $(prem_rho_at_tomo_depths), PREM_VP = $(prem_vp_at_tomo_depths)
  const PREM_VS = $(prem_vs_at_tomo_depths), PREM_PRESSURE = $(prem_pressure_at_tomo_depths)
  const PREM_FAN_DEPTHS = $(fan_depths_km)
  const PREM_FAN_RHO = $(prem_fan_rho), PREM_FAN_VP = $(prem_fan_vp)
  const PREM_FAN_VS = $(prem_fan_vs), PREM_FAN_PRESSURE = $(prem_fan_pressure)
  const COAST = $(coast_js_literal(coast_line_id, coast_lon, coast_lat))

  // Max |value| across an entire lateral grid (every depth, every lat/lon) -- used once per
  // model to build the fan's fixed color scale (see colorForFan).
  function globalAbsMax(data){
    let m = 0
    for(const layer of data) for(const row of layer) for(const v of row) m = Math.max(m, Math.abs(v))
    return m || 1
  }

  // One entry per selectable parameter. The three tomography (lateral) models each carry
  // their OWN lon/lat/depth grid -- they were built independently, at different native
  // resolutions -- plus a precomputed global |dVs| max for the cross-section's fixed color
  // scale. The four PREM profiles all share one depth-only grid (no lon/lat) since PREM has
  // no lateral variation to begin with.
  const PARAMS = {
    gypsum:   {lateral:true, data:TOMO_GYPSUM,  lons:LONS_GYPSUM,  lats:LATS_GYPSUM,  depths:DEPTHS_GYPSUM,
               gmax:globalAbsMax(TOMO_GYPSUM),  label:'dVs', unit:'%',
               modelLabel:'GyPSuM (Simmons et al. 2010)', modelName:'GyPSuM'},
    s40rts:   {lateral:true, data:TOMO_S40RTS,  lons:LONS_S40RTS,  lats:LATS_S40RTS,  depths:DEPTHS_S40RTS,
               gmax:globalAbsMax(TOMO_S40RTS),  label:'dVs', unit:'%',
               modelLabel:'S40RTS (Ritsema et al. 2011)', modelName:'S40RTS'},
    gladm25:  {lateral:true, data:TOMO_GLADM25, lons:LONS_GLADM25, lats:LATS_GLADM25, depths:DEPTHS_GLADM25,
               gmax:globalAbsMax(TOMO_GLADM25), label:'dVs', unit:'%',
               modelLabel:'GLAD-M25 (Lei et al. 2020)', modelName:'GLAD-M25'},
    density:  {lateral:false, data:PREM_RHO,      fanData:PREM_FAN_RHO,      fanDepths:PREM_FAN_DEPTHS, depths:DEPTHS_GYPSUM, label:'density',  unit:'kg/m³'},
    vp:       {lateral:false, data:PREM_VP,       fanData:PREM_FAN_VP,       fanDepths:PREM_FAN_DEPTHS, depths:DEPTHS_GYPSUM, label:'Vp',       unit:'km/s'},
    vs:       {lateral:false, data:PREM_VS,       fanData:PREM_FAN_VS,       fanDepths:PREM_FAN_DEPTHS, depths:DEPTHS_GYPSUM, label:'Vs',       unit:'km/s'},
    pressure: {lateral:false, data:PREM_PRESSURE, fanData:PREM_FAN_PRESSURE, fanDepths:PREM_FAN_DEPTHS, depths:DEPTHS_GYPSUM, label:'pressure', unit:'GPa'},
  }

  // Illustrative cross-sections through well-known features -- picked for a clean,
  // recognizable section (GyPSuM is a coarsened 3° model here, not precision-migrated
  // tomography), not scientific precision.
  const PRESETS = {
    hawaii:      {label:'Hawaii plume',           lonA:-155.5, latA:19.5,  lonB:-155.5, latB:4.5},
    iceland:     {label:'Iceland plume',          lonA:-17.0,  latA:64.9,  lonB:-17.0,  latB:49.9},
    yellowstone: {label:'Yellowstone plume',      lonA:-110.5, latA:44.5,  lonB:-110.5, latB:29.5},
    pacific:     {label:'Farallon (Pacific) slab',lonA:-125.0, latA:40.0,  lonB:-75.0,  latB:38.0},
    japan:       {label:'Japan slab',             lonA:148.0,  latA:38.0,  lonB:128.0,  latB:38.0},
    tonga:       {label:'Tonga slab',             lonA:-170.0, latA:-20.0, lonB:178.0,  latB:-18.0},
  }

  const cvs = par.querySelector('#tgcvs'), ctx = cvs.getContext('2d')
  const fcvs = par.querySelector('#tgfan'), fctx = fcvs.getContext('2d')
  function hidpi(cv, cx, w, h){ cv.width=Math.round(w*DPR); cv.height=Math.round(h*DPR); cv.style.width=w+'px'; cv.style.height=h+'px'; cx.setTransform(DPR,0,0,DPR,0,0) }
  hidpi(cvs, ctx, SEC, SEC)
  hidpi(fcvs, fctx, SEC, SEC)

  // The sphere+coastlines render (20,000 quads) only actually changes with yaw, pitch,
  // param, or depthIdx -- never with A/B. Caching it to an offscreen canvas and blitting
  // that with drawImage() means dragging a pin, or drawing a brand-new section, redraws
  // only the cheap marker/arc overlay instead of the whole sphere every mousemove.
  const sphereCache = document.createElement('canvas')
  const sphereCtx = sphereCache.getContext('2d')
  hidpi(sphereCache, sphereCtx, SEC, SEC)

  const cbar1 = par.querySelector('#tgcbar1'), cbar1ctx = cbar1.getContext('2d')
  const cbar2 = par.querySelector('#tgcbar2'), cbar2ctx = cbar2.getContext('2d')
  hidpi(cbar1, cbar1ctx, 240, 34)
  hidpi(cbar2, cbar2ctx, 240, 34)

  let yaw = 0.9, pitch = 0.35, dragMode = null, lastX = 0, lastY = 0
  let depthIdx = $(w.depth_idx)
  let param = "$(w.param)"
  let lonA = $(w.lonA), latA = $(w.latA), lonB = $(w.lonB), latB = $(w.latB)
  let currentPreset = 'hawaii'

  function rot(p){
    const x = p[0]*Math.cos(yaw) - p[1]*Math.sin(yaw)
    const y = p[0]*Math.sin(yaw) + p[1]*Math.cos(yaw)
    const z = p[2]
    return [x, y*Math.cos(pitch) - z*Math.sin(pitch), y*Math.sin(pitch) + z*Math.cos(pitch)]
  }
  // Algebraic inverse of rot() (pitch undone, then yaw undone) -- turns a point already in
  // rotated/screen space back into the sphere's own unrotated coordinates, needed to figure
  // out which lat/lon a mouse click landed on.
  function unrot(P){
    const y1 = P[1]*Math.cos(pitch) + P[2]*Math.sin(pitch)
    const z1 = -P[1]*Math.sin(pitch) + P[2]*Math.cos(pitch)
    const px = P[0]*Math.cos(yaw) + y1*Math.sin(yaw)
    const py = -P[0]*Math.sin(yaw) + y1*Math.cos(yaw)
    return [px, py, z1]
  }
  function sph(theta, phi, r){ return [r*Math.sin(theta)*Math.cos(phi), r*Math.sin(theta)*Math.sin(phi), r*Math.cos(theta)] }
  // Screen X is CX - x (not CX + x): with phi = longitude directly (see sampleAt/
  // llToThetaPhi), the raw math convention x = r*sinθ*cosφ puts increasing longitude at
  // DEcreasing x -- i.e. east would render to the left of a north-up view, backwards from
  // the standard "east is right" convention for a globe viewed from outside. Negating x
  // here is the one place that needs to change to fix it.
  function proj(p){ return [CX - p[0], CY - p[2]] }
  function llToThetaPhi(lon, lat){ return [(90-lat)*Math.PI/180, lon*Math.PI/180] }

  // Great-circle helpers, all on unit vectors -- used both to draw the arc between the two
  // markers on the globe and to sample continuously along it for the fan panel.
  function llToXYZ(lon, lat){ const [theta, phi] = llToThetaPhi(lon, lat); return sph(theta, phi, 1) }
  function xyzToThetaPhi(v){
    const n = Math.hypot(v[0], v[1], v[2])
    return [Math.acos(Math.max(-1, Math.min(1, v[2]/n))), Math.atan2(v[1], v[0])]
  }
  function dot3(a, b){ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
  function greatCircleAngle(vA, vB){ return Math.acos(Math.max(-1, Math.min(1, dot3(vA, vB)))) }
  function slerp(vA, vB, ang, f){
    if(ang < 1e-9) return vA
    const s = Math.sin(ang)
    const a = Math.sin((1-f)*ang)/s, b = Math.sin(f*ang)/s
    return [a*vA[0]+b*vB[0], a*vA[1]+b*vB[1], a*vA[2]+b*vB[2]]
  }
  // Inverse of the globe's own render pipeline: given a screen click, solve for the point
  // on the R_OUT sphere facing the camera, then unrotate it back to lon/lat. Returns null
  // for a click that misses the sphere's silhouette entirely.
  function screenToLL(mx, my){
    const px = CX - mx, pz = CY - my
    const rem = R_OUT*R_OUT - px*px - pz*pz
    if(rem < 0) return null
    const py = Math.sqrt(rem)
    const v = unrot([px, py, pz])
    const n = Math.hypot(v[0], v[1], v[2])
    const theta = Math.acos(Math.max(-1, Math.min(1, v[2]/n)))
    const phi = Math.atan2(v[1], v[0])
    return [phi*180/Math.PI, 90 - theta*180/Math.PI]
  }

  // Opaque diverging ramp for GyPSuM: positive (fast/cold, dense-like) -> blue, negative
  // (slow/hot, light-like) -> red. Same convention as geoid-kernel.jl's densityColor.
  function divergingColor(v, mx){
    const t = Math.max(-1, Math.min(1, v/mx))
    const bg = [18, 22, 29]
    const hi = t >= 0 ? [59, 130, 246] : [239, 68, 68]
    const a = Math.min(1, Math.abs(t) * 1.6)
    const c = [0,1,2].map(i => Math.round(bg[i] + (hi[i]-bg[i])*a))
    return 'rgb('+c[0]+','+c[1]+','+c[2]+')'
  }
  // Sequential ramp for the PREM profiles: dark (low) -> amber (high), normalized across
  // the WHOLE profile's range (now including the fan's full surface-to-centre reach where
  // available), not per depth-layer -- these are monotonic-ish with depth, so a per-layer
  // max would make every layer look identically saturated and hide the depth trend.
  function sequentialColor(v, vmin, vmax){
    const t = vmax > vmin ? Math.max(0, Math.min(1, (v-vmin)/(vmax-vmin))) : 0
    const bg = [18, 22, 29], hi = [250, 204, 21]
    const c = [0,1,2].map(i => Math.round(bg[i] + (hi[i]-bg[i])*t))
    return 'rgb('+c[0]+','+c[1]+','+c[2]+')'
  }
  function profileRange(p){
    // Never used for lateral params (colorFor/colorForFan only consult range in the
    // sequential/PREM branch) -- and it must bail out early for them, not just skip the
    // *use* of the result: p.data for a lateral param is a 3-level nested array (depth x
    // lat x lon), so `for(const v of arr)` binds v to a whole 2D layer, and Math.min/max
    // coerce that to a number via an implicit array-to-string stringify of ~7,000+
    // elements, PER LAYER. That silent quadratic-ish cost (not the quad count, not the
    // canvas fill calls) was the actual source of GyPSuM's per-frame slowdown.
    if(p.lateral) return [0, 0]
    const arr = p.fanData || p.data
    let mn = Infinity, mx = -Infinity
    for(const v of arr){ mn = Math.min(mn, v); mx = Math.max(mx, v) }
    return [mn, mx]
  }
  function layerMax(p, di){
    let mx = 0
    for(const row of p.data[di]) for(const v of row) mx = Math.max(mx, Math.abs(v))
    return mx || 1
  }
  // Takes the layer's max PRE-COMPUTED (gmax), not (p, di) -- drawGlobe() calls this once
  // per quad (tens of thousands of times per frame while dragging), and layerMax() itself
  // is an O(lat*lon) scan; computing it fresh per-quad made every drag frame O(quads *
  // lat*lon), the actual cost behind a sluggish globe. Callers hoist layerMax(p, depthIdx)
  // once per frame instead.
  function colorFor(p, v, gmax, range){
    return p.lateral ? divergingColor(v, gmax) : sequentialColor(v, range[0], range[1])
  }
  // Fan coloring uses one FIXED diverging scale per model (p.gmax, precomputed once when
  // PARAMS is built) instead of layerMax's per-selected-layer scale -- otherwise the fan's
  // colors would rescale confusingly every time the (otherwise unrelated) depth slider
  // moved, even though the fan always shows every depth at once.
  function colorForFan(p, v, range){
    return p.lateral ? divergingColor(v, p.gmax) : sequentialColor(v, range[0], range[1])
  }

  // Each lateral (tomography) param carries its OWN lon/lat/depth grid -- GyPSuM, S40RTS,
  // and GLAD-M25 were built on three different native grids -- so sampling always reads
  // p.lons/p.lats/p.depths rather than a single shared global.
  function sampleAt(p, di, theta, phi){
    if(!p.lateral) return p.data[di]
    const latDeg = 90 - theta*180/Math.PI
    let lonDeg = (phi*180/Math.PI)
    lonDeg = ((lonDeg + 180) % 360 + 360) % 360 - 180
    const LATSp = p.lats, LONSp = p.lons
    let li = Math.round((latDeg - LATSp[0]) / (LATSp[1]-LATSp[0]))
    let lj = Math.round((lonDeg - LONSp[0]) / (LONSp[1]-LONSp[0]))
    li = Math.max(0, Math.min(LATSp.length-1, li))
    lj = Math.max(0, Math.min(LONSp.length-1, lj))
    return p.data[di][li][lj]
  }
  function sampleAtDepth(p, depthKm, theta, phi){
    const D = p.depths
    if(depthKm <= D[0]) return sampleAt(p, 0, theta, phi)
    if(depthKm >= D[D.length-1]) return sampleAt(p, D.length-1, theta, phi)
    for(let i=1;i<D.length;i++){
      if(D[i] >= depthKm){
        const t = (depthKm-D[i-1])/(D[i]-D[i-1])
        const a = sampleAt(p, i-1, theta, phi), b = sampleAt(p, i, theta, phi)
        return a + t*(b-a)
      }
    }
    return sampleAt(p, D.length-1, theta, phi)
  }
  // Fan-only sampler: GyPSuM defers to sampleAtDepth (correctly clamped at its true ~2900
  // km data limit); a PREM parameter instead interpolates its own full 0-6371 km axis, so
  // the fan can genuinely reach the centre rather than clamping at GyPSuM's shallower limit.
  function sampleForFan(p, depthKm, theta, phi){
    if(p.lateral) return sampleAtDepth(p, depthKm, theta, phi)
    const D = p.fanDepths, V = p.fanData
    if(depthKm <= D[0]) return V[0]
    if(depthKm >= D[D.length-1]) return V[D.length-1]
    for(let i=1;i<D.length;i++){
      if(D[i] >= depthKm){
        const t = (depthKm-D[i-1])/(D[i]-D[i-1])
        return V[i-1] + t*(V[i]-V[i-1])
      }
    }
    return V[V.length-1]
  }

  function drawCoastlines(tctx){
    tctx.strokeStyle = 'rgba(245,243,239,0.85)'; tctx.lineWidth = 1
    const rc_out = R_OUT * 1.004
    for(const seg of COAST){
      tctx.beginPath()
      let started = false
      for(const [lon, lat] of seg){
        const [theta, phi] = llToThetaPhi(lon, lat)
        const rp = rot(sph(theta, phi, rc_out))
        if(rp[1] < 0){ started = false; continue }   // back hemisphere -- lift the pen
        const s = proj(rp)
        if(!started){ tctx.moveTo(s[0], s[1]); started = true } else tctx.lineTo(s[0], s[1])
      }
      tctx.stroke()
    }
  }

  // The great-circle arc between A and B, drawn on the globe's own surface, plus the two
  // draggable pins themselves.
  function markerScreenPos(lon, lat){
    const [theta, phi] = llToThetaPhi(lon, lat)
    const rp = rot(sph(theta, phi, R_OUT))
    if(rp[1] < 0) return null
    return proj(rp)
  }
  function drawMarker(lon, lat, label, color){
    const [theta, phi] = llToThetaPhi(lon, lat)
    const rp = rot(sph(theta, phi, R_OUT*1.012))
    if(rp[1] < 0) return
    const s = proj(rp)
    ctx.beginPath(); ctx.arc(s[0], s[1], 6, 0, 2*Math.PI)
    ctx.fillStyle = color; ctx.fill()
    ctx.strokeStyle = '#0a0f18'; ctx.lineWidth = 1.5; ctx.stroke()
    ctx.fillStyle = '#e5e7eb'; ctx.font = 'bold 12px sans-serif'
    ctx.fillText(label, s[0]+9, s[1]-8)
  }
  function drawMarkersAndArc(){
    const vA = llToXYZ(lonA, latA), vB = llToXYZ(lonB, latB)
    const ang = greatCircleAngle(vA, vB)
    ctx.beginPath()
    let started = false
    const NA = 80
    for(let i=0;i<=NA;i++){
      const v = slerp(vA, vB, ang, i/NA)
      const rp = rot([v[0]*R_OUT*1.006, v[1]*R_OUT*1.006, v[2]*R_OUT*1.006])
      if(rp[1] < 0){ started = false; continue }
      const s = proj(rp)
      if(!started){ ctx.moveTo(s[0], s[1]); started = true } else ctx.lineTo(s[0], s[1])
    }
    ctx.strokeStyle = '#22d3ee'; ctx.lineWidth = 2; ctx.stroke()
    drawMarker(lonA, latA, 'A', '#facc15')
    drawMarker(lonB, latB, 'B', '#f472b6')
  }

  function fmtLL(lon, lat){
    const ns = lat >= 0 ? 'N' : 'S', ew = lon >= 0 ? 'E' : 'W'
    return Math.abs(lat).toFixed(1)+'°'+ns+', '+Math.abs(lon).toFixed(1)+'°'+ew
  }
  function updateReadouts(){
    par.querySelector('#tgAll').textContent = fmtLL(lonA, latA)
    par.querySelector('#tgBll').textContent = fmtLL(lonB, latB)
  }
  function setCustomPreset(){
    if(currentPreset !== 'custom'){
      currentPreset = 'custom'
      par.querySelector('#tgpreset').value = 'custom'
    }
  }
  // Rotates the globe so the midpoint of A and B faces the camera, by solving for the
  // yaw/pitch that puts that point's rotated x and z components at zero (screen center).
  function snapView(){
    const vA = llToXYZ(lonA, latA), vB = llToXYZ(lonB, latB)
    const ang = greatCircleAngle(vA, vB)
    const vm = slerp(vA, vB, ang, 0.5)
    const r = Math.hypot(vm[0], vm[1])
    yaw = Math.atan2(vm[0], vm[1])
    pitch = Math.atan2(-vm[2], r)
  }

  function updateLegend(){
    const p = PARAMS[param]
    const el = par.querySelector('#tglegend')
    if(p.lateral){
      el.innerHTML = 'blue = fast / cold (dense-like) &middot; red = slow / hot (light-like)<br>'+
        'model: '+p.modelLabel+' &middot; grid stops just above the CMB &mdash; '+
        "S-waves don't resolve the liquid outer core.<br>"+
        'drag <b style="color:#facc15">A</b> / <b style="color:#f472b6">B</b> on the globe, or pick a preset'
    } else {
      const cur = p.data[depthIdx].toFixed(p.unit === 'kg/m³' ? 0 : 2)
      el.innerHTML = 'PREM (Dziewonski &amp; Anderson 1981) &middot; '+p.label+' at this depth: <b>'+cur+' '+p.unit+'</b><br>'+
        'uniform at every location &mdash; the section (right) shows how it changes with depth, all the way to the centre'
    }
    updateColorbars()
  }

  // Renders one gradient legend strip, sampled directly from the same divergingColor/
  // sequentialColor functions the globe and fan actually paint with -- not an approximated
  // canvas gradient -- so what's shown here always matches what's on screen exactly,
  // including divergingColor's non-linear saturation curve near the extremes.
  function drawColorbar(cctx, w, h, isDiverging, mx, vmin, vmax, fmt){
    cctx.clearRect(0,0,w,h)
    const barH = 14, barY = 2
    for(let i=0;i<w;i++){
      const t = i/(w-1)
      const v = isDiverging ? (-mx + t*2*mx) : (vmin + t*(vmax-vmin))
      cctx.fillStyle = isDiverging ? divergingColor(v, mx) : sequentialColor(v, vmin, vmax)
      cctx.fillRect(i, barY, 1, barH)
    }
    cctx.strokeStyle = '#4b5563'; cctx.lineWidth = 1
    cctx.strokeRect(0.5, barY+0.5, w-1, barH-1)
    cctx.fillStyle = '#9ca3af'; cctx.font = '10px sans-serif'
    const lo = isDiverging ? -mx : vmin, hi = isDiverging ? mx : vmax
    cctx.textAlign = 'left'; cctx.fillText(fmt(lo), 0, barY+barH+12)
    cctx.textAlign = 'right'; cctx.fillText(fmt(hi), w, barY+barH+12)
    if(isDiverging){ cctx.textAlign = 'center'; cctx.fillText('0', w/2, barY+barH+12) }
    cctx.textAlign = 'left'
  }

  // Two DIFFERENT normalizations for a tomography model, by design (per the widget's whole
  // premise -- "this depth" vs "the whole section"): the map view saturates its diverging
  // scale at THIS depth layer's own |dVs| max (layerMax), so a quiet layer still shows
  // contrast; the cross-section instead uses one FIXED scale (p.gmax, the |dVs| max across
  // every depth of the active model) so colors stay comparable as you look up and down the
  // same section. PREM fields have no such distinction -- both panels already share one
  // physical range -- so only one bar is shown for them.
  function updateColorbars(){
    const p = PARAMS[param]
    const label2 = par.querySelector('#tgcbarLabel2'), c2 = par.querySelector('#tgcbar2')
    if(p.lateral){
      par.querySelector('#tgcbarLabel1').textContent = 'surface, this depth (' + Math.round(p.depths[depthIdx]) + ' km)'
      drawColorbar(cbar1ctx, 240, 34, true, layerMax(p, depthIdx), 0, 0, v => v.toFixed(1)+'%')
      label2.style.display = ''; c2.style.display = ''
      label2.textContent = 'section, all depths (A → B)'
      drawColorbar(cbar2ctx, 240, 34, true, p.gmax, 0, 0, v => v.toFixed(1)+'%')
    } else {
      par.querySelector('#tgcbarLabel1').textContent = p.label + ' (' + p.unit + '), full depth range'
      const range = profileRange(p)
      const fmt = v => p.unit === 'kg/m³' ? Math.round(v).toString() : v.toFixed(1)
      drawColorbar(cbar1ctx, 240, 34, false, 1, range[0], range[1], fmt)
      label2.style.display = 'none'; c2.style.display = 'none'
    }
  }

  // The expensive part: the outer sphere's 20,000-quad shell (one value per quad, at the
  // single currently-selected depth) plus coastlines, rendered into the offscreen cache.
  // Only depends on yaw, pitch, param, and depthIdx -- NOT on A/B -- so it only needs to
  // re-run when one of those actually changes, not on every pin-drag/draw-a-section frame.
  function drawSphereLayer(){
    sphereCtx.clearRect(0,0,SEC,SEC)
    const p = PARAMS[param]
    const range = profileRange(p)
    // Hoisted OUT of the per-quad loop below -- see the comment on colorFor().
    const gmax = p.lateral ? layerMax(p, depthIdx) : 1

    const NT = 100, NP = 200
    for(let i=0;i<NT;i++){
      const t0 = Math.PI*i/NT, t1 = Math.PI*(i+1)/NT
      for(let j=0;j<NP;j++){
        const p0 = 2*Math.PI*j/NP, p1 = 2*Math.PI*(j+1)/NP
        const pm = (p0+p1)/2
        const c = [sph(t0,p0,R_OUT), sph(t0,p1,R_OUT), sph(t1,p1,R_OUT), sph(t1,p0,R_OUT)]
        const rc = c.map(rot)
        if((rc[0][1]+rc[1][1]+rc[2][1]+rc[3][1])/4 < 0) continue
        const v = sampleAt(p, depthIdx, (t0+t1)/2, pm)
        sphereCtx.beginPath()
        const s0 = proj(rc[0]); sphereCtx.moveTo(s0[0], s0[1])
        for(let k=1;k<4;k++){ const s=proj(rc[k]); sphereCtx.lineTo(s[0], s[1]) }
        sphereCtx.closePath()
        sphereCtx.fillStyle = colorFor(p, v, gmax, range)
        sphereCtx.fill()
      }
    }

    drawCoastlines(sphereCtx)

    sphereCtx.beginPath(); sphereCtx.arc(CX,CY,R_OUT,0,2*Math.PI)
    sphereCtx.strokeStyle='#4b5563'; sphereCtx.lineWidth=1; sphereCtx.stroke()
  }

  // The cheap part: blit the cached sphere, then draw the pins/arc and text fresh on top.
  // This is what runs on every mousemove while dragging a pin or drawing a new section.
  function drawGlobe(){
    ctx.clearRect(0,0,SEC,SEC)
    ctx.drawImage(sphereCache, 0, 0, SEC, SEC)
    drawMarkersAndArc()

    const p = PARAMS[param]
    const modelPrefix = p.lateral ? p.modelName+' ' : ''
    par.querySelector('#tgGlobeInfo').textContent = modelPrefix+p.label+' at '+Math.round(p.depths[depthIdx])+' km'

    updateLegend()
  }

  const FAN_R = SEC*0.85
  const FAN_APEX = [SEC/2, SEC*0.92]
  function fanPoint(sectorAngle, r){
    return [FAN_APEX[0] + r*Math.sin(sectorAngle), FAN_APEX[1] - r*Math.cos(sectorAngle)]
  }

  function drawFan(){
    fctx.clearRect(0,0,SEC,SEC)
    const p = PARAMS[param]
    const range = profileRange(p)
    const vA = llToXYZ(lonA, latA), vB = llToXYZ(lonB, latB)
    const ang = greatCircleAngle(vA, vB)
    const halfAng = Math.max(ang, 0.01)/2
    const maxDepth = p.lateral ? p.depths[p.depths.length-1] : R_EARTH_KM

    const NK = 90, NR = 60
    for(let k=0;k<NK;k++){
      const f0 = k/NK, f1 = (k+1)/NK
      const a0 = -halfAng + f0*2*halfAng, a1 = -halfAng + f1*2*halfAng
      const v = slerp(vA, vB, ang, (f0+f1)/2)
      const [theta, phi] = xyzToThetaPhi(v)
      for(let i=0;i<NR;i++){
        const rf0 = i/NR, rf1 = (i+1)/NR
        const depth0 = maxDepth*rf0, depth1 = maxDepth*rf1
        const depthM = (depth0+depth1)/2
        const r0 = FAN_R*(1-depth0/R_EARTH_KM), r1 = FAN_R*(1-depth1/R_EARTH_KM)
        const val = sampleForFan(p, depthM, theta, phi)
        const c0 = fanPoint(a0,r0), c1 = fanPoint(a1,r0), c2 = fanPoint(a1,r1), c3 = fanPoint(a0,r1)
        fctx.beginPath()
        fctx.moveTo(c0[0],c0[1]); fctx.lineTo(c1[0],c1[1]); fctx.lineTo(c2[0],c2[1]); fctx.lineTo(c3[0],c3[1])
        fctx.closePath()
        fctx.fillStyle = colorForFan(p, val, range)
        fctx.fill()
      }
    }

    // reference depth lines: transition zone, CMB, ICB (whichever fall within this fan's reach)
    const refs = [[660,'660 km'], [CMB_KM,'CMB'], [ICB_KM,'ICB']]
    fctx.font = '11px sans-serif'
    for(const [dk, lab] of refs){
      if(dk > maxDepth) continue
      const rr = FAN_R*(1-dk/R_EARTH_KM)
      fctx.beginPath()
      fctx.setLineDash([4,3])
      const s0 = fanPoint(-halfAng, rr); fctx.moveTo(s0[0], s0[1])
      for(let k=1;k<=40;k++){ const a=-halfAng+2*halfAng*k/40; const s=fanPoint(a,rr); fctx.lineTo(s[0],s[1]) }
      fctx.strokeStyle = 'rgba(156,163,175,0.55)'; fctx.lineWidth = 1; fctx.stroke()
      fctx.setLineDash([])
      const lbl = fanPoint(halfAng, rr)
      fctx.fillStyle = 'rgba(156,163,175,0.85)'
      fctx.fillText(lab, lbl[0]+4, lbl[1]+3)
    }

    // ring marking the depth the globe's own slider is currently showing
    if(PARAMS[param].depths[depthIdx] <= maxDepth){
      const rSel = FAN_R*(1 - PARAMS[param].depths[depthIdx]/R_EARTH_KM)
      fctx.beginPath()
      const s0 = fanPoint(-halfAng, rSel); fctx.moveTo(s0[0], s0[1])
      for(let k=1;k<=40;k++){ const a=-halfAng+2*halfAng*k/40; const s=fanPoint(a,rSel); fctx.lineTo(s[0],s[1]) }
      fctx.strokeStyle = '#f5f3ef'; fctx.lineWidth = 2; fctx.stroke()
    }

    // outline: surface arc + the two radial edges down to the data limit
    const rSurf = FAN_R, rInner = FAN_R*(1-maxDepth/R_EARTH_KM)
    fctx.beginPath()
    let s0 = fanPoint(-halfAng, rSurf); fctx.moveTo(s0[0], s0[1])
    for(let k=1;k<=60;k++){ const a=-halfAng+2*halfAng*k/60; const s=fanPoint(a,rSurf); fctx.lineTo(s[0],s[1]) }
    fctx.strokeStyle = '#4b5563'; fctx.lineWidth = 1.5; fctx.stroke()
    fctx.beginPath()
    const eA0 = fanPoint(-halfAng,rSurf), eA1 = fanPoint(-halfAng,rInner)
    const eB0 = fanPoint(halfAng,rSurf), eB1 = fanPoint(halfAng,rInner)
    fctx.moveTo(eA0[0],eA0[1]); fctx.lineTo(eA1[0],eA1[1])
    fctx.moveTo(eB0[0],eB0[1]); fctx.lineTo(eB1[0],eB1[1])
    fctx.strokeStyle = '#374151'; fctx.lineWidth = 1; fctx.stroke()

    const sA = fanPoint(-halfAng, rSurf), sB = fanPoint(halfAng, rSurf)
    fctx.fillStyle = '#facc15'; fctx.font = 'bold 13px sans-serif'; fctx.fillText('A', sA[0]-14, sA[1]-4)
    fctx.fillStyle = '#f472b6'; fctx.fillText('B', sB[0]+4, sB[1]-4)

    const label = (currentPreset !== 'custom' && PRESETS[currentPreset]) ? PRESETS[currentPreset].label : 'Custom section'
    const fillNote = p.lateral ? 'fills to the CMB — no '+p.modelName+' data below it' : 'fills all the way to the centre (r = 0)'
    par.querySelector('#tgFanInfo').textContent = label+' · '+Math.round(ang*R_EARTH_KM)+' km across — '+fillNote
  }

  snapView()
  updateReadouts()
  drawSphereLayer()
  drawGlobe()
  drawFan()

  // Left-drag is the primary gesture: starting on an existing pin moves just that pin;
  // starting on empty globe draws a brand-new A->B section from scratch (A snaps to
  // the press point, B live-follows the drag). Rotating the view moved to right-drag
  // (contextmenu suppressed below) so it doesn't collide with "draw a new section" on
  // the far more common empty-space gesture.
  cvs.addEventListener('contextmenu', e => e.preventDefault())
  cvs.addEventListener('mousedown', e=>{
    const mx = e.offsetX, my = e.offsetY
    if(e.button === 2){
      dragMode = 'rotate'; lastX = mx; lastY = my
      return
    }
    const sA = markerScreenPos(lonA, latA), sB = markerScreenPos(lonB, latB)
    const near = s => s && Math.hypot(s[0]-mx, s[1]-my) < 14
    if(near(sA)) dragMode = 'A'
    else if(near(sB)) dragMode = 'B'
    else {
      const ll = screenToLL(mx, my)
      if(ll){
        lonA = ll[0]; latA = ll[1]
        lonB = ll[0]; latB = ll[1]
        dragMode = 'draw'
        setCustomPreset()
        updateReadouts()
        drawGlobe()   // sphere unchanged here too
        drawFan()
      }
    }
  })
  cvs.addEventListener('mousemove', e=>{
    const mx = e.offsetX, my = e.offsetY
    if(dragMode === 'rotate'){
      const dx = mx-lastX, dy = my-lastY
      lastX = mx; lastY = my
      yaw += dx*0.008
      pitch = Math.max(-1.3, Math.min(1.3, pitch + dy*0.008))
      drawSphereLayer()   // orientation changed -- the cache must be redrawn
      drawGlobe()
    } else if(dragMode === 'A' || dragMode === 'B' || dragMode === 'draw'){
      const ll = screenToLL(mx, my)
      if(ll){
        if(dragMode === 'A'){ lonA = ll[0]; latA = ll[1] }
        else if(dragMode === 'B'){ lonB = ll[0]; latB = ll[1] }
        else { lonB = ll[0]; latB = ll[1] }   // 'draw': B follows the drag, A stays put
        setCustomPreset()
        updateReadouts()
        drawGlobe()   // sphere unchanged -- just re-blit the cache and redraw the pins
        drawFan()
      }
    }
  })
  window.addEventListener('mouseup', ()=>{ dragMode = null })

  // Each param carries its own depth grid (three different native tomography grids, plus
  // PREM's), so switching params can't just keep the same depthIdx -- it has to re-target
  // the new grid's nearest depth (in km) and refresh the slider's own range to match.
  function retargetDepthIndex(oldDepths){
    const oldKm = oldDepths[depthIdx]
    const newDepths = PARAMS[param].depths
    let bestI = 0, bestD = Infinity
    for(let i=0;i<newDepths.length;i++){
      const d = Math.abs(newDepths[i]-oldKm)
      if(d < bestD){ bestD = d; bestI = i }
    }
    depthIdx = bestI
    const slider = par.querySelector('#tgdepth')
    slider.max = newDepths.length-1
    slider.value = depthIdx
    par.querySelector('#tgdepthv').textContent = Math.round(newDepths[depthIdx])+' km'
  }

  par.querySelector('#tgparam').addEventListener('change', e=>{
    const oldDepths = PARAMS[param].depths
    param = e.target.value
    retargetDepthIndex(oldDepths)
    drawSphereLayer(); drawGlobe(); drawFan()
  })
  par.querySelector('#tgdepth').addEventListener('input', e=>{
    depthIdx = parseInt(e.target.value)
    par.querySelector('#tgdepthv').textContent = Math.round(PARAMS[param].depths[depthIdx])+' km'
    drawSphereLayer(); drawGlobe(); drawFan()
  })
  par.querySelector('#tgpreset').addEventListener('change', e=>{
    const key = e.target.value
    currentPreset = key
    if(key !== 'custom'){
      const pr = PRESETS[key]
      lonA = pr.lonA; latA = pr.latA; lonB = pr.lonB; latB = pr.latB
      // Only bail out of a PREM param (which has no lateral structure to reveal); leave
      // whichever tomography model was already active alone.
      if(!PARAMS[param].lateral){
        const oldDepths = PARAMS[param].depths
        param = 'gypsum'
        par.querySelector('#tgparam').value = 'gypsum'
        retargetDepthIndex(oldDepths)
      }
      snapView()
    }
    updateReadouts()
    drawSphereLayer(); drawGlobe(); drawFan()
  })
</script>
""")
    end

    # See the matching comment on `_epi_ready`: forces the correct execution order on
    # a fresh restart, since Pluto doesn't otherwise detect that the bind cell below
    # depends on this cell defining `TomographyGlobeInput`.
    const _tgi_ready = true
end

# ╔═╡ 2c64d4c0-1b70-4724-8d71-f068ee9c738a
begin
    _tgi_ready
    @bind _tgi TomographyGlobeInput()
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Interpolations = "~0.16.2"
PlutoUI = "~0.7.75"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "cccabc2be65c5bfa97342a77cdb26caab2444b04"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

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

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

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

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "db8a06ef983af758d285665a0398703eb5bc1d66"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.75"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

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

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
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
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

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

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═81b3647e-d741-11f0-b633-478f0d655f0b
# ╠═98f30e60-f1e5-4102-8539-2c6600c3de90
# ╟─169e2164-b62f-4b76-bf12-bec1a6aba243
# ╟─c8266e03-23be-43ea-9292-97359ab7aa41
# ╟─98788e51-80d7-4a50-adf8-f0731d99f2b0
# ╟─3f6270d2-a8f3-404c-87eb-d544a36cbe29
# ╟─5a36b51d-da51-4ae2-95d8-227650643530
# ╟─070fe588-0e0e-4c02-8888-b724670a7f0a
# ╟─2c64d4c0-1b70-4724-8d71-f068ee9c738a
# ╟─98799827-8c07-4706-9354-faf9ee1a3e15
# ╟─22a88956-3e85-4d61-a3fa-aec2da926fd9
# ╟─e2493968-3488-4a46-b0c1-222010c6dcb9
# ╠═67a3d0fe-0a8c-4163-9d79-5e9a511ddab0
# ╟─0defc265-fe5a-42c6-a876-5e5e16cde7fb
# ╠═f2e9afae-368b-4224-bfc9-7aa681e0d4e7
# ╟─bd4aff4f-a4af-4d01-967f-02d7acf2a098
# ╟─259dab7b-4d62-4a5a-aaee-de0de4fd86fb
# ╟─7e1b745c-5174-4e36-90a0-51db489d4768
# ╠═306d2112-cab7-4fdc-b441-e9228e243cdc
# ╟─0fad939d-cb1e-43bf-80d1-bdd8957515f9
# ╠═1faacc81-f508-4832-9a31-148aa3e3c7cc
# ╟─b8d8c763-ccf3-4c29-af91-3959a0b903eb
# ╟─4bef3d37-3bfd-4c7f-83c6-31b4318f9e5a
# ╠═7cb6d994-2002-44f2-bbc7-e73f146d1704
# ╠═5c8a2f14-6e91-4b3a-9d27-8f41c0a9e6b2
# ╠═d3f7b920-1a44-4e6c-8b19-2c5d9a7f0e33
# ╟─3ac914f3-56e4-4025-81d1-9461a0bfc835
# ╠═9bebe536-3dc8-4e45-b27d-160dc0b4f608
# ╟─8438fcce-2d66-4693-b747-641a4c2264b7
# ╠═a41bd2f5-bce3-4a30-8f31-5b86c19e59fb
# ╠═a916b96c-ad21-4cd1-ba7e-6434ce2e5396
# ╠═81868aa1-b26e-4aea-ab9a-2fc2a9e522b0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
