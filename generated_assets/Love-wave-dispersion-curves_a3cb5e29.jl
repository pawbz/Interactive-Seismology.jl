### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ 8d8f440d-4b79-4835-841b-739b5171f979
using PlutoUI, Printf, Roots, LinearAlgebra

# ╔═╡ 29afa36b-391f-4861-aead-d62918daf3c6
using HypertextLiteral: @htl

# ╔═╡ 8d5d9594-2197-4ccc-a7a4-0f0d54a2370a
using PlutoPlotly

# ╔═╡ ee482657-fd34-4c92-95cd-0bf8e254676c
TableOfContents(include_definitions=true)

# ╔═╡ d3bc6646-b22f-4e22-87a8-da37a9aa1d1a
md"""
# Love Wave Dispersion Curves
This notebook provides an environment to explore dispersion of Love waves, a type of surface seismic wave that plays a key role in seismology and Earth structure studies.  

Love waves are horizontally polarized shear waves that are *trapped* in near-surface low-velocity layers. Their dispersion (variation of phase velocity with period) carries information about the elastic properties and layering of the crust and upper mantle.

Here,
- build a layered Earth model interactively, by adjusting thickness, shear velocity (`Vs`), and density (`ρ`) for each layer;
- compute dispersion curves for Love waves and observe how changes in structure affect the shape of the dispersion curve.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 1af153bc-e471-4608-984b-61f1116dfa16
md"""
Number of layers: $(@bind n_layers Slider(1:5, default=2, show_value=true))
"""

# ╔═╡ d2c48242-92ce-11f0-01e6-cb1b8fe1fecb

struct Layer
    thickness::Float64   # km
    vp::Float64          # km/s
    vs::Float64          # km/s
    rho::Float64         # g/cm³
    Qp::Float64          # P-wave quality factor
    Qs::Float64          # S-wave quality factor
end


# ╔═╡ d058a019-40c2-4546-b578-808381506d1f
struct EarthModel
    layers::Vector{Layer}
end

# ╔═╡ 363d8b6d-ddc9-4048-a02b-16346158112c
"""
    starting_phase_velocity(vp, vs) -> c

Compute a starting solution for phase velocity using Newton iteration.
"""
function starting_phase_velocity(vp::Float64, vs::Float64)
    c = 0.95 * vs
    for _ in 1:5
        γ = vs / vp
        κ = c / vs
        k2 = κ^2
        gk2 = (γ * κ)^2
        fac1 = sqrt(1 - gk2)
        fac2 = sqrt(1 - k2)
        fr = (2 - k2)^2 - 4 * fac1 * fac2
        frp = -4 * (2 - k2) * κ +
              4 * fac2 * γ^2 * κ / fac1 +
              4 * fac1 * κ / fac2
        frp /= vs
        c -= fr / frp
    end
    return c
end

# ╔═╡ 6c766c4b-8a1c-4f8b-9467-fec5a91fcd05
struct DispersionResult
    periods::Vector{Float64}
    phase_velocities::Vector{Float64}
end

# ╔═╡ 756e06e3-4f25-4f55-b242-26a45a1b7dd3
periods = 5.0:5.0:100.0

# ╔═╡ 9c9ad5f6-debf-4376-aaf6-759e9b0c949b
"""
    layer_matrix_SH(layer::Layer, ω, c)

Return the 2×2 complex transfer matrix for an SH layer of thickness d:
    [ u_bottom ]   [ A  B ] [ u_top   ]
    [ t_bottom ] = [ C  D ] [ t_top   ]

Where u is horizontal displacement and t is shear traction.
"""
function layer_matrix_SH(layer::Layer, ω::Float64, c::Float64)
    # angular wavenumber horizontally
    k = ω / c

    # shear velocity and density
    vs = layer.vs
    ρ = layer.rho

    # shear modulus (units consistent if rho in mass/volume and vs in velocity units)
    μ = ρ * vs^2

    # vertical wavenumber q (can be complex)
    # q^2 = (ω/vs)^2 - k^2
    q2 = (ω / vs)^2 - k^2
    q = sqrt(complex(q2))          # complex sqrt handles evanescent/propagating

    d = layer.thickness

    # if thickness is effectively zero (half-space), return identity (no propagation)
    if isapprox(d, 0.0; atol=1e-14)
        return Matrix{ComplexF64}(I, 2, 2)
    end

    # arguments for trig/hyperbolic via complex arithmetic
    α = q * d

    # use cos and sin of complex argument (covers cosh/sinh automatically)
    ca = cos(α)
    sa = sin(α)

    # Layer transfer matrix for SH (standard form)
    # [ cos(qd)          sin(qd)/(μ*q)   ]
    # [ -μ*q*sin(qd)     cos(qd)         ]
    M = zeros(ComplexF64, 2, 2)
    M[1, 1] = ca
    M[1, 2] = sa / (μ * q)
    M[2, 1] = -μ * q * sa
    M[2, 2] = ca

    return M
end

# ╔═╡ ac92ec5d-de77-4235-a701-50de13c502b1
"""
    dispersion_function_SH(model, ω, c)

Return complex-valued F(c) whose zero corresponds to a Love-wave mode:
    F(c) = M21 - R * M11
where M = product of layer matrices from top -> bottom and
R = -μ_b * q_b (bottom half-space traction/displacement ratio for decaying solution).
"""
function dispersion_function_SH(model::EarthModel, ω::Float64, c::Float64)
    layers = model.layers
    N = length(layers)
    if N < 1
        error("Model must contain at least one layer (half-space).")
    end

    # Build total propagator M (top -> bottom)
    M = Matrix{ComplexF64}(I, 2, 2)
    for L in layers
        M = layer_matrix_SH(L, ω, c) * M   # multiply top->down (pre-multiply or post? we choose pre)
    end

    # bottom (half-space) properties (assumed last entry)
    Lb = layers[end]
    vs_b = Lb.vs
    ρ_b = Lb.rho
    μ_b = ρ_b * vs_b^2

    # compute bottom vertical wavenumber q_b (for half-space decaying solution)
    q2_b = (ω / vs_b)^2 - (ω / c)^2
    q_b = sqrt(complex(q2_b))

    # radiation/decay relation tb/ub = R  (we take R = -μ_b * q_b)
    R = -μ_b * q_b

    # M maps [u_top; t_top] -> [u_bot; t_bot]. For t_top=0 (free surface), ub = M11*u0, tb = M21*u0.
    # impose tb/ub = R  => M21/M11 = R  => M21 - R*M11 = 0
    F = M[2, 1] - R * M[1, 1]
    return F
end


# ╔═╡ 03e6b940-ba12-41df-a3a2-2b4e0c44d64b
"""
    find_root_bracketed(f, cmin, cmax; nsteps=200)

Search for a sign change of real(f(c)) in [cmin, cmax] by sampling nsteps intervals.
Return (cL, cR) bracket where sign change occurs, or (nothing, nothing) if none found.
"""
function find_root_bracketed(f, cmin::Float64, cmax::Float64; nsteps::Int=200)
    xs = range(cmin, cmax; length=nsteps + 1)
    prev = real(f(first(xs)))
    for i in 2:length(xs)
        x = xs[i]
        val = real(f(x))
        if !isfinite(val)
            prev = val
            continue
        end
        if prev == 0.0
            return (xs[i-1], xs[i-1])
        end
        if sign(prev) != sign(val)
            return (xs[i-1], x)
        end
        prev = val
    end
    return (nothing, nothing)
end

# ╔═╡ e05989a8-f2f1-49ad-a17f-8ddbf5776ab4
"""
    solve_phase_velocity_love(model, T; c_search_pad=0.1)

Find a phase velocity at period T (seconds) for Love (SH) waves.
The search interval is set relative to model layer S velocities.
"""
function solve_phase_velocity_love(model, T; c_search_pad=0.10)
    ω = 2π / T

    # set search bounds using layer shear velocities
    vs_vals = [L.vs for L in model.layers if L.vs > 0.0]
    if isempty(vs_vals)
        throw(ArgumentError("Model has no positive shear velocities."))
    end
    vmin = minimum(vs_vals)
    vmax = maximum(vs_vals)

    # expand search window by fraction
    cmin = max(1e-5, (1 - c_search_pad) * vmin)
    cmax = (1 + c_search_pad) * vmax * 1.5    # allow some margin above vmax
    @show cmin, cmax
    f(c) = dispersion_function_SH(model, ω, c)

    # find bracket by sampling
    (L, R) = find_root_bracketed(f, cmin, cmax; nsteps=800)
    if L === nothing
        # No sign change found -> return NaN to indicate no root
        return NaN
    end

    # If bracketed exactly at same point (rare), return that point
    if L == R
        return L
    end

    # Use robust bisection on real part of f
    g(c) = real(f(c))
    root = find_zero(g, (L, R), Bisection(); tol=1e-6, maxevals=100)
    return root
end

# ╔═╡ cca474e2-7152-43d5-be66-432847de2897
"""
    solve_love_dispersion(model, periods)

Compute Love-wave dispersion (phase velocity for each period).
Returns DispersionResult.
"""
function solve_love_dispersion(model, periods)
    velocities = Float64[]
    for T in periods
        @printf("Solving T=%.3f s ... ", T)
        c = try
            c = solve_phase_velocity_love(model, T)
        catch err
            @warn "error solving period $T: $err"
            NaN
        end
        if isnan(c)
            println("no root found")
        else
            println(@sprintf("c=%.5f km/s", c))
        end
        push!(velocities, c)
    end
    return DispersionResult(periods, velocities)
end


# ╔═╡ bcc394ad-284d-4299-891a-065ec7a6c6b3
"""
    solve_phase_velocity(model, T; c0_guess=nothing)

Find phase velocity at period T (s).
"""
function solve_phase_velocity(model::EarthModel, T::Float64; c0_guess=nothing)
    ω = 2π / T

    # crude initial bracket around Vs
    vsurf = model.layers[1].vs
    cmin, cmax = 0.9 * vsurf, 1.2 * vsurf

    f(c) = real(dispersion_determinant(model, ω, c))

    # Find root of f(c) = 0
    croot = find_zero(f, (cmin, cmax), Bisection(), verbose=false)
    return croot
end

# ╔═╡ 80086cb5-7e39-444b-a85e-ddf46a67ec42
md"## Appendix"

# ╔═╡ 99e9fc2c-fe13-4f94-b613-86cedb0f3653
function layer_table_input(n_layers::Int; vp_vs_ratio=1.73)
    ui = PlutoUI.combine() do Child
        # header
        header = @htl("""
        <tr style="background:#f0f3f8; text-align:center;">
          <th style="padding:6px;">#</th>
          <th style="padding:6px;">Thickness</th>
          <th style="padding:6px;">Vs</th>
          <th style="padding:6px;">Density</th>
        </tr>
        """)

        # rows
        rows = Any[]
        for i in 1:n_layers
            # no thickness slider for last layer (half-space)
            tw = i < n_layers ? Child("thickness_$i", Slider(0.5:0.5:200, default=35, show_value=true)) : "∞"

            vs = Child("vs_$i", Slider(2.5:0.01:5.5, default=(i == 1) ? 3.4 : 4.6, show_value=true))
            rw = Child("rho_$i", Slider(1.0:0.01:6.0, default=(i == 1) ? 2.6 : 3.4, show_value=true))

            push!(rows, @htl("""
                <tr>
                  <td style="text-align:center; padding:6px;"><b>Layer $i</b></td>
                  <td style="padding:6px; text-align:center;">$tw</td>
                  <td style="padding:6px; text-align:center;">$vs</td>
                  <td style="padding:6px; text-align:center;">$rw</td>
                </tr>
            """))
        end

        tbl = @htl("""
        <table style="border-collapse:collapse; border:1px solid #ddd; width:100%;">
          <thead>$header</thead>
          <tbody>
            $(rows...)
          </tbody>
        </table>
        """)

        md"""
        ### Layer Editor
        Number of layers: **$n_layers**

        $tbl
        _Notes:_ Thickness in km.  
        Constraint enforced: **Vp = $(vp_vs_ratio) × Vs > Vs**  
        The last layer is treated as a half-space (no thickness).
        """
    end

    return PlutoUI.Experimental.transformed_value(ui) do vals
        layers = [
            begin
                if i < n_layers
                    t = vals[3*(i-1)+1]
                    vs = vals[3*(i-1)+2]
                    ρ = vals[3*(i-1)+3]
                else
                    # last layer: no thickness slider
                    t = 0.0
                    vs = vals[3*(i-1)+1]
                    ρ = vals[3*(i-1)+2]
                end

                vp = vp_vs_ratio * vs # (not used, anyway for love waves )

                Layer(t, vp, vs, ρ, 100.0, 100.0)
            end
            for i in 1:n_layers
        ]
        le = layers[end]
        # add artificial thick layer to impose radiation condition
        layer_new = Layer(1000, le.vp, le.vs, le.rho, 100.0, 100.0)
        layers = vcat(layers[1:end-1], layer_new, le)
    end
end

# ╔═╡ 7bce51e0-9191-42eb-8fe3-cc68f9b2681f
(@bind layers layer_table_input(n_layers))

# ╔═╡ a3cd5bb3-9b59-48a5-a505-4199a04941ed
model = EarthModel(layers)

# ╔═╡ 45579d35-df74-454f-91d4-03542f77e1ac
# ╠═╡ show_logs = false
res = solve_love_dispersion(model, periods)

# ╔═╡ 6e696382-0a37-46f1-8209-0bf5d7dbc82f
let
    plt = Plot(
        scatter(
            x=res.periods,
            y=res.phase_velocities,
            mode="lines+markers",
            name="Dispersion Curve"
        )
    )

    plt.layout = Layout(
        title="Love Fundamental-Mode Dispersion Curve",
        xaxis=attr(title="Period (s)", showgrid=true),
        yaxis=attr(
            title="Phase Velocity (km/s)",
            showgrid=true,
            range=[2.5, 6.5]
        ),
        hovermode="closest"
    )

    plot(plt)
end

# ╔═╡ fc1fd90c-020a-4f2d-aedd-66115ac6a287
md"""
### Reference
This notebook is inspired by the classical implementation of surface wave theory in  

Herrmann, R. B. (2013), *Computer Programs in Seismology*.

The Julia code here provides a **clean reimplementation** of some of the ideas behind CPS, designed for **interactive exploration** rather than production‐level modeling.
md"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"

[compat]
HypertextLiteral = "~0.9.5"
PlutoPlotly = "~0.3.4"
PlutoUI = "~0.7.71"
Roots = "~2.2.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "9038c24234f0293a4bad7fb0208acaec306d02e2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

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
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

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

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

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

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

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
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "28278bb0053da0fd73537be94afd1682cc5a0a83"
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
deps = ["AbstractPlutoDingetjes", "Dates", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "PlotlyBase", "PlutoUI", "Reexport"]
git-tree-sha1 = "b470931aa2a8112c8b08e66ea096c6c62c60571e"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.3.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
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

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "442b4353ee8c26756672afb2db81894fc28811f3"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.6"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
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
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═ee482657-fd34-4c92-95cd-0bf8e254676c
# ╟─d3bc6646-b22f-4e22-87a8-da37a9aa1d1a
# ╟─1af153bc-e471-4608-984b-61f1116dfa16
# ╟─7bce51e0-9191-42eb-8fe3-cc68f9b2681f
# ╟─6e696382-0a37-46f1-8209-0bf5d7dbc82f
# ╠═45579d35-df74-454f-91d4-03542f77e1ac
# ╠═d2c48242-92ce-11f0-01e6-cb1b8fe1fecb
# ╠═d058a019-40c2-4546-b578-808381506d1f
# ╠═363d8b6d-ddc9-4048-a02b-16346158112c
# ╠═6c766c4b-8a1c-4f8b-9467-fec5a91fcd05
# ╠═a3cd5bb3-9b59-48a5-a505-4199a04941ed
# ╠═756e06e3-4f25-4f55-b242-26a45a1b7dd3
# ╠═9c9ad5f6-debf-4376-aaf6-759e9b0c949b
# ╠═e05989a8-f2f1-49ad-a17f-8ddbf5776ab4
# ╠═cca474e2-7152-43d5-be66-432847de2897
# ╠═ac92ec5d-de77-4235-a701-50de13c502b1
# ╠═03e6b940-ba12-41df-a3a2-2b4e0c44d64b
# ╠═bcc394ad-284d-4299-891a-065ec7a6c6b3
# ╟─80086cb5-7e39-444b-a85e-ddf46a67ec42
# ╠═8d8f440d-4b79-4835-841b-739b5171f979
# ╠═29afa36b-391f-4861-aead-d62918daf3c6
# ╠═8d5d9594-2197-4ccc-a7a4-0f0d54a2370a
# ╠═99e9fc2c-fe13-4f94-b613-86cedb0f3653
# ╟─fc1fd90c-020a-4f2d-aedd-66115ac6a287
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
