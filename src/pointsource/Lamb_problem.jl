### A Pluto.jl notebook ###
# v0.19.46

#> [frontmatter]
#> title = "Lamb's Problem"
#> layout = "layout.jlhtml"
#> tags = ["pointsource"]
#> description = "Superposition of planewaves"

using Markdown
using InteractiveUtils

# ╔═╡ 897afffa-77e8-11ef-1a54-c73d1df8f6a4
using QuadGK, Bessels, PlutoPlotly, FFTW, Symbolics, SpecialFunctions, PlutoUI, ProgressLogging

# ╔═╡ 4621f804-6a44-4c46-af37-3d364ba74cfe
using Groebner, Nemo, Statistics

# ╔═╡ 6188ff3e-cfd5-4c9c-aa39-619dc280d494
TableOfContents()

# ╔═╡ 1bbedd43-9e69-4ba0-a251-855f1f63dcf8
md"# Lamb's Problem"

# ╔═╡ 7a3bd5df-0f9f-489c-b522-4098c325c0c0
md"## Conical Wave"

# ╔═╡ 29020170-3173-40a4-aac6-9d0fd1212aaf
function get_vertical_wavenumber(c, p)
    c2 = abs2(c)
    ξ = if (inv(c2) >= abs2(p))
        -sqrt(inv(c2) - abs2(p))
    else
        im * sqrt(abs2(p) - inv(c2))
    end
    return ξ
end

# ╔═╡ 05ec38ea-3431-490c-bc38-24f8c1b2d54f
function conical_wave(p, ωgrid, x, y, z, c)
    r = sqrt(abs2(x) + abs2(y))
    ξ = get_vertical_wavenumber(c, p)
    return Bessels.besselj0.(ωgrid * p * r) .* exp.(im * ωgrid * ξ * z) * p / ξ
end

# ╔═╡ e2fea7d5-00a3-4796-ad51-26f2dcffa55b
md"## Medium"

# ╔═╡ d430c8ec-fa65-4e72-951e-987ba9f2f160
media = (; α₁=5.0, α₂=7.0, β₁=3.0, β₂=5.0, ρ₁=3.0, ρ₂=7.0)

# ╔═╡ 360de109-d7db-4402-bca8-5b39c6f17da9
md"## Derive Reflection Coefficients"

# ╔═╡ 3dfbfc66-8f87-4906-84f5-495bf113e4ef
@variables ω::Real ξ₁ ν₁ ξ₂ ν₂ x::Real y::Real z::Real p r::Real λ₁::Real μ₁::Real λ₂::Real μ₂::Real t::Real

# ╔═╡ 0de7835c-eb50-4de0-874b-e08270d988aa
@syms α₁ β₁ ρ₁ α₂ β₂ ρ₂

# ╔═╡ 8ef12535-49a1-46fe-a26c-880358755a16
@variables Aₚ Aₛ Bₚ Bₛ

# ╔═╡ c077d6f7-9fc7-4c0b-bf42-26882a374cbf
@variables ı

# ╔═╡ 68e54426-2a9f-409c-acd9-e0e6e9681487
plane_wave(p, η, A) = A * exp(ı * ω * (t - (p * x + η * z)))

# ╔═╡ 4c8103ff-8928-4c78-ad63-d0175738ef7f
incident_potential = plane_wave(p, ξ₁, 1.0)

# ╔═╡ e2438fc6-78e9-4d06-9e4c-ee4ff22d0062
reflected_P_potential = plane_wave(p, -ξ₁, Aₚ)

# ╔═╡ 0ee489a5-3769-4bdf-9e88-cdf54199d692
reflected_S_potential = plane_wave(p, -ν₁, Aₛ)

# ╔═╡ 75fc91ce-840f-4625-87eb-1870f668f5c0
transmitted_P_potential = plane_wave(p, ξ₂, Bₚ)

# ╔═╡ 03a80695-170a-4a2f-a124-630dc3dbf5b6
transmitted_S_potential = plane_wave(p, ν₂, Bₛ)

# ╔═╡ 99ee6828-0a0d-4048-8ffd-d25b7848368b
begin
    Dx = Differential(x)
    Dz = Differential(z)
end

# ╔═╡ 8fe6d5f7-7ec9-4a67-a557-ffe51124066b
begin
    uxInc = expand_derivatives.(Dx(incident_potential))
    uzInc = expand_derivatives.(Dz(incident_potential))
end

# ╔═╡ 277a6190-3f6e-4008-aa2c-aac84d205f42
begin
    uxRefP = expand_derivatives.(Dx.(reflected_P_potential))
    uzRefP = expand_derivatives.(Dz.(reflected_P_potential))
    uxRefS = expand_derivatives.(-Dz.(reflected_S_potential))
    uzRefS = expand_derivatives.(Dx.(reflected_S_potential))
end

# ╔═╡ 1d58e672-5c50-4000-8842-7c592b22c32e
begin
    uxTransP = expand_derivatives.(Dx.(transmitted_P_potential))
    uzTransP = expand_derivatives.(Dz.(transmitted_P_potential))
    uxTransS = expand_derivatives.(-Dz.(transmitted_S_potential))
    uzTransS = expand_derivatives.(Dx.(transmitted_S_potential))
end

# ╔═╡ 65c2f4a8-92f7-4673-ba7a-ba2367b558aa
σzz1 = (expand_derivatives(λ₁ .* Dx.(uxInc + uxRefP + uxRefS) .+ (λ₁ + 2 * μ₁) .* Dz.(uzInc + uzRefP + uzRefS)))

# ╔═╡ ae0e4d63-8348-46ac-994f-11b1ce1001fb
σzz2 = (expand_derivatives(λ₂ .* Dx.(uxTransP + uxTransS) .+ (λ₂ + 2 * μ₂) .* Dz.(uzTransP + uzTransS)))

# ╔═╡ 7b55d65d-7654-4ac5-9b3c-d0f8b52905cd
σxz1 = (expand_derivatives.(μ₁ .* Dx.(uzInc + uzRefP + uzRefS) .+ μ₁ .* Dz.(uxInc + uxRefP + uxRefS)))

# ╔═╡ 0e3d2a52-e357-4b54-80ce-fb7f44c6334f
σxz2 = (expand_derivatives.(μ₂ .* Dx.(uzTransP + uzTransS) .+ μ₂ .* Dz.(uxTransP + uxTransS)))

# ╔═╡ 060c9127-775a-45c0-b31b-2a6dcacac179
Ux_constraint = uxInc + uxRefP + uxRefS - uxTransP - uxTransS

# ╔═╡ 0f034701-7953-4428-bde9-b2a21fbc5c6e
Uz_constraint = uzInc + uzRefP + uzRefS - uzTransP - uzTransS

# ╔═╡ 52149e10-3a89-412a-bd43-e29e26592551
Ux_z0 = simplify(substitute(Ux_constraint, z => 0) / ω / ı / exp(ı * ω * (t - p * x)))

# ╔═╡ 94fc6b36-7e38-4fc1-af47-9c82d6595e6b
Uz_z0 = simplify(substitute(Uz_constraint, z => 0) / ω / ı / exp(ı * ω * (t - p * x)))

# ╔═╡ cca1bbd7-5455-47f9-bef9-448ece10c2b6
σzz_constraint = σzz1 - σzz2

# ╔═╡ d94908f6-1d0f-46ff-8c8e-ef0d407b539d
σzz_z0 = simplify(substitute(σzz_constraint, [z => 0, ı^2 => -1]) / ω^2 / exp(ı * ω * (t - p * x)))

# ╔═╡ a015f013-7333-44e7-b7bb-a293aefeac0b
σxz_constraint = σxz1 - σxz2

# ╔═╡ 3aaca98b-4b75-4970-97cf-b4649f443017
σxz_z0 = simplify(substitute(σxz_constraint, [ı^2 => -1, z => 0]) / ω^2 / exp(ı * ω * (t - p * x)))

# ╔═╡ 93d2f1ad-06ed-4080-a0ce-5ec1eb642e4e
function subs(ex, media)
    α₁, α₂, β₁, β₂, ρ₁, ρ₂ = media
    ex = simplify(substitute(ex, [ı^2 => -1, λ₁ => ρ₁ * (α₁^2 - 2 * β₁^2), μ₁ => β₁^2 * ρ₁, λ₂ => ρ₂ * (α₂^2 - 2 * β₂^2), μ₂ => β₂^2 * ρ₂,
        ξ₁^2 => (1 / α₁^2 - p^2), ξ₂^2 => (1 / α₂^2 - p^2), ν₁^2 => (1 / β₁^2 - p^2), ν₂^2 => (1 / β₂^2 - p^2)]))
    ex = simplify(substitute(ex, [α₁ => α₁, α₂ => α₂, β₁ => β₁, β₂ => β₂, ρ₁ => ρ₁, ρ₂ => ρ₂]))
    return ex
end

# ╔═╡ 02237fbb-f2d3-4436-ad0e-c60df62befaa
function get_A(f, p1, media)
    α₁, α₂, β₁, β₂, ρ₁, ρ₂ = media
    ξ₁p = get_vertical_wavenumber(α₁, p1)
    ξ₂p = get_vertical_wavenumber(α₂, p1)
    ν₁p = get_vertical_wavenumber(β₁, p1)
    ν₂p = get_vertical_wavenumber(β₂, p1)
    return Complex{Float64}(f(p1, ξ₁p, ξ₂p, ν₁p, ν₂p))
    # ex = simplify(substitute(ex, [p=>]))
end

# ╔═╡ 239a50ff-3970-47ff-b772-fc3e0cb247fe
function solve_linear_AB(media)
    constraints = [σzz_z0, σxz_z0, Ux_z0, Uz_z0]
    # Avec = Symbolics.symbolic_linear_solve(constraints, [Aₚ, Bₚ, Aₛ, Bₛ])
    Avec = Symbolics.symbolic_solve(constraints, [Aₚ, Bₚ, Aₛ, Bₛ])
    Avec_subs = Dict(string(key) => subs(Avec[1][key], media) for key in keys(Avec[1]))
end

# ╔═╡ fa4b8951-569e-4517-9b0a-f77a22a31fca
begin
    Avec = solve_linear_AB(media)
    Afunctions_vec = Dict(string(key) => build_function(Avec[key], p, ξ₁, ξ₂, ν₁, ν₂, expression=Val{false}) for key in keys(Avec))
end;

# ╔═╡ f0a4c439-ee9b-4001-97db-66f6ac5afd5a
md"## Sommerfeld Integral — Reflectivity Method"

# ╔═╡ 12f7733a-edd6-481f-8166-ef967520b35a
begin
    Nt = 256
    tgrid = range(0, 50, length=Nt)
    ωgrid = collect(rfftfreq(Nt, inv(step(tgrid)))) * 2 * pi
    f0 = 0.5
    ω0 = 2 * pi * f0
    ωW = @. exp(-0.1(ωgrid - ω0)^2)
    Ω = (; Nt, tgrid, ωgrid, ω0, f0, ωW)
end

# ╔═╡ d1e7d7c1-e794-4426-87b8-12898c13fc4d
plot(scatter(x=ωgrid, y=ωW))

# ╔═╡ 2103d1d4-8041-455d-b513-02643762bd3a
function get_seismogram(x1, y1, z1, Nt, H1, Ω, media)
    U = im * Ω.ωgrid .* Ω.ωW .* quadgk(0, 5, norm=v -> maximum(abs, v)) do p
        Aₚ = get_A(Afunctions_vec["Aₚ"], p, media)
        Aₛ = get_A(Afunctions_vec["Aₛ"], p, media)

        return (
            # direct arrival
            conical_wave(p, ωgrid, x1, y1, z1, media.α₁ + 1e-3)
            # P reflection
            + Aₚ * (conical_wave(p, ωgrid, x1, y1, H1 - z1, media.α₁ + 1e-3))
            # S reflection (wrong, need to add displacement, not potential)
            + Aₛ * (conical_wave(p, ωgrid, x1, y1, H1 - z1, media.β₁ + 1e-3))
        )
    end[1]
    # end
    # U = Vector{Complex{Float64}}(U)
    return real.(irfft(U, Nt))
end

# ╔═╡ 72aca1bf-3997-4be3-b316-921b45480c1a
md"## Seismograms"

# ╔═╡ 27f3136e-f6ac-43a1-9c95-8807cca36455
x1range = collect(range(20, 100, length=15))

# ╔═╡ 3693ffe9-979b-416a-85de-8438f00a6fc2
H1 = 100.0

# ╔═╡ 461ec509-b07d-45e8-85ef-d838190c1247
Threads.nthreads()

# ╔═╡ c695e4d3-c49d-4587-8adc-cdd4c001325e
md"## Appendix"

# ╔═╡ f711a17d-530d-48e3-9e32-2cc2a519e131
y1 = 0

# ╔═╡ 9c12ba1e-0c79-4982-9157-d95c453b4e19
z1 = 0.01

# ╔═╡ 53a45db3-8571-4fc7-a855-5173030b7cb9
begin
    U = [zeros(length(ωgrid)) for x in x1range]
    Threads.@threads for i in 1:length(x1range)
        U[i] = get_seismogram(x1range[i], y1, z1, Nt, H1, Ω, media)
    end
end

# ╔═╡ 96556fec-4bea-4132-946a-92d17372df54
let
    U = hcat(U...)
    Umaxclip = maximum(abs, U) * 0.1
    plot(heatmap(y=tgrid, x=x1range, z=U, zmin=-Umaxclip, zmax=Umaxclip, colorscale=:seismic,), Layout(width=300, xaxis=attr(title="distance to receiver"), yaxis=attr(title="time (s)")))
end

# ╔═╡ e91921c5-4600-4344-85c7-e35820ba2daa
begin
    N = 25

    zgrid = range(5, stop=18, length=N)
    xgrid = range(-50, stop=50, length=N)
    ygrid = range(0, stop=100, length=N)
end

# ╔═╡ a69c4ef6-3b04-4101-b1d4-bb5cf7d40233
# ╠═╡ disabled = true
#=╠═╡
function get_xz_snapshot(xgrid, y1, zgrid, ω, H1, media)
	snapshot = map(Iterators.product(zgrid, xgrid)) do (z1, x1)

		# map(ωgrid) do ω
	# im*ω*quadgk(p -> conical_wave(p, ω, x, y, z, c), 0, inv(c*c))[1]
	# return im*ω*(quadgk(0, Inf) do p
	# 	conical_wave(p, ω, x1, y1, z1, c1)
	# end)[1]
			# im*ω*quadgk(p -> conical_wave(p, ω, x1, y1, z1, c1), 0, Inf)[1]+
	im*ω*quadgk(0, Inf, rtol=1e-3) do p
		Aₚ = get_A(Afunctions_vec["Aₚ"], p, media)
		Aₚ * 
		(conical_wave(p, ω, x1, y1, H1-z1, media.α₁ + 1e-3))
	end[1]
	end
	# U = Matrix{Complex{Float64}}(hcat(U...))
	# return real.(irfft(U, Nt, 1))
end
  ╠═╡ =#

# ╔═╡ 0b3a10e3-b7a2-4e2b-9e88-b8c7bea7834f
# snapshot = get_xz_snapshot(xgrid, 0.0, zgrid, 2,  H1, media)

# ╔═╡ e7c1cbf4-7939-45a2-a754-b36ae9bdcab4
md"""### To Do
- Parallelize the numerical integration over several receivers
- Snapshots to visualize generation of e.g., headwaves
- Plot displacement field, instead of potentials
- Add free surface
"""

# ╔═╡ 9d10ab94-38d5-406e-8697-7f5bc0c11df2
md"### References"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Bessels = "~0.2.8"
FFTW = "~1.8.0"
Groebner = "~0.7.5"
Nemo = "~0.45.7"
PlutoPlotly = "~0.5.0"
PlutoUI = "~0.7.60"
ProgressLogging = "~0.1.4"
QuadGK = "~2.11.1"
SpecialFunctions = "~2.4.0"
Symbolics = "~6.12.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "4950f4f413fdb3b07d9c4689d2a805f865136db5"

[[deps.ADTypes]]
git-tree-sha1 = "5a5eafb8344b81b8c2237f8a6f6b3602b3f6180e"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.8.1"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractAlgebra]]
deps = ["InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "ca88b21f7430be07fb6e9c0500fe4af65e023857"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.41.11"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown"]
git-tree-sha1 = "b392ede862e506d451fc1616e79aa6f4c673dab8"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.38"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsDatesExt = "Dates"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsTestExt = "Test"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "3640d077b6dafd64ceb8fd5c1ec76f7ca53bcf76"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.16.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "cb25e4b105cc927052c2314f8291854ea59bf70a"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.4"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "3e4b134270b372f2ed4d4d0e936aabaefc1802bc"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "e6c693a0e4394f8fda0e51a5bdf5aef26f8235e9"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.111"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "490392af2c7d63183bfa2c8aaa6ab981c5ba7561"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.14"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "bbf1ace0781d9744cb697fb856bd2c3f6568dadb"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "fc3951d4d398b5515f91d7fe5d45fc31dccb3c9b"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.8.5"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FLINT_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "MPFR_jll", "OpenBLAS32_jll"]
git-tree-sha1 = "4bde447e7bc6de73de36870fc942df02f4dd44b5"
uuid = "e134572f-a0d5-539d-bddf-3cad8db41a82"
version = "300.100.300+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"
version = "6.2.1+6"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Atomix", "Combinatorics", "ExprTools", "HostCPUFeatures", "Logging", "MultivariatePolynomials", "Nemo", "PrecompileTools", "PrettyTables", "Primes", "Printf", "Random", "TimerOutputs"]
git-tree-sha1 = "6d28801daff0c972b9bad25d900d60b2adc36c72"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.7.5"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "10bd689145d2c3b2a9844005d01087cc1194e79e"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "f389674c99bfcde17dc57454011aa44d5a260a40"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

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

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MPFR_jll]]
deps = ["Artifacts", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"
version = "4.2.0+1"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "5c1d1d9361e1417e5a065e1f84dc3686cbdaea21"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "d0a6b1096b584a2b88efb70a92f8cb8c881eb38a"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.6"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Nemo]]
deps = ["AbstractAlgebra", "FLINT_jll", "Libdl", "LinearAlgebra", "Pkg", "Random", "RandomExtensions", "SHA"]
git-tree-sha1 = "bae0f1452c648d84535b581103ab8a59cf10a32a"
uuid = "2edaba10-b0f1-5616-af89-8c11ac63239a"
version = "0.45.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6065c4cff8fee6c6770b277af45d5082baacdba1"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.24+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "BaseDirs", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "Reexport", "TOML"]
git-tree-sha1 = "653b48f9c4170343c43c2ea0267e451b68d69051"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.5.0"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "b8a399e95663485820000f26b6a43c794e166a49"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.4"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "b034171b93aebc81b3e1890a036d13a9c4a9e3e0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.27.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "82584ea03bda16156ceaa6af75f9c8e1287029d3"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.53.2"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "e39c5f217f9aca640c8e27ab21acf557a3967db5"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.10"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "25514a6f200219cd1073e4ff23a6324e4a7efe64"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.5.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "87d51a3ee9a4b0d2fe054bdd3fc2436258db2603"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "988e04b34a4c3b824fb656f542473df99a4f610d"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.30"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "fabf4650afe966a2ba646cabd924c3fd43577fc3"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "3927e02dc7648a45ec6aa592bcd8374094a44740"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.7.1"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "8b48697e7fec6d4b7c4a9fe892857a5ed2bae7e8"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.12.0"

    [deps.Symbolics.extensions]
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

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
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "5a13ae8a41237cff5ecf34f73eb1b8f42fff6531"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.24"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

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
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═6188ff3e-cfd5-4c9c-aa39-619dc280d494
# ╟─1bbedd43-9e69-4ba0-a251-855f1f63dcf8
# ╠═96556fec-4bea-4132-946a-92d17372df54
# ╟─7a3bd5df-0f9f-489c-b522-4098c325c0c0
# ╠═29020170-3173-40a4-aac6-9d0fd1212aaf
# ╠═05ec38ea-3431-490c-bc38-24f8c1b2d54f
# ╟─e2fea7d5-00a3-4796-ad51-26f2dcffa55b
# ╠═d430c8ec-fa65-4e72-951e-987ba9f2f160
# ╟─360de109-d7db-4402-bca8-5b39c6f17da9
# ╠═3dfbfc66-8f87-4906-84f5-495bf113e4ef
# ╠═0de7835c-eb50-4de0-874b-e08270d988aa
# ╠═8ef12535-49a1-46fe-a26c-880358755a16
# ╠═c077d6f7-9fc7-4c0b-bf42-26882a374cbf
# ╠═68e54426-2a9f-409c-acd9-e0e6e9681487
# ╠═4c8103ff-8928-4c78-ad63-d0175738ef7f
# ╠═e2438fc6-78e9-4d06-9e4c-ee4ff22d0062
# ╠═0ee489a5-3769-4bdf-9e88-cdf54199d692
# ╠═75fc91ce-840f-4625-87eb-1870f668f5c0
# ╠═03a80695-170a-4a2f-a124-630dc3dbf5b6
# ╠═99ee6828-0a0d-4048-8ffd-d25b7848368b
# ╠═8fe6d5f7-7ec9-4a67-a557-ffe51124066b
# ╠═277a6190-3f6e-4008-aa2c-aac84d205f42
# ╠═1d58e672-5c50-4000-8842-7c592b22c32e
# ╠═65c2f4a8-92f7-4673-ba7a-ba2367b558aa
# ╠═ae0e4d63-8348-46ac-994f-11b1ce1001fb
# ╠═7b55d65d-7654-4ac5-9b3c-d0f8b52905cd
# ╠═0e3d2a52-e357-4b54-80ce-fb7f44c6334f
# ╠═060c9127-775a-45c0-b31b-2a6dcacac179
# ╠═0f034701-7953-4428-bde9-b2a21fbc5c6e
# ╠═52149e10-3a89-412a-bd43-e29e26592551
# ╠═94fc6b36-7e38-4fc1-af47-9c82d6595e6b
# ╠═cca1bbd7-5455-47f9-bef9-448ece10c2b6
# ╠═d94908f6-1d0f-46ff-8c8e-ef0d407b539d
# ╠═a015f013-7333-44e7-b7bb-a293aefeac0b
# ╠═3aaca98b-4b75-4970-97cf-b4649f443017
# ╠═93d2f1ad-06ed-4080-a0ce-5ec1eb642e4e
# ╠═02237fbb-f2d3-4436-ad0e-c60df62befaa
# ╠═239a50ff-3970-47ff-b772-fc3e0cb247fe
# ╠═fa4b8951-569e-4517-9b0a-f77a22a31fca
# ╟─f0a4c439-ee9b-4001-97db-66f6ac5afd5a
# ╠═12f7733a-edd6-481f-8166-ef967520b35a
# ╠═d1e7d7c1-e794-4426-87b8-12898c13fc4d
# ╠═2103d1d4-8041-455d-b513-02643762bd3a
# ╟─72aca1bf-3997-4be3-b316-921b45480c1a
# ╠═27f3136e-f6ac-43a1-9c95-8807cca36455
# ╠═3693ffe9-979b-416a-85de-8438f00a6fc2
# ╠═461ec509-b07d-45e8-85ef-d838190c1247
# ╠═53a45db3-8571-4fc7-a855-5173030b7cb9
# ╟─c695e4d3-c49d-4587-8adc-cdd4c001325e
# ╠═897afffa-77e8-11ef-1a54-c73d1df8f6a4
# ╠═4621f804-6a44-4c46-af37-3d364ba74cfe
# ╠═f711a17d-530d-48e3-9e32-2cc2a519e131
# ╠═9c12ba1e-0c79-4982-9157-d95c453b4e19
# ╠═e91921c5-4600-4344-85c7-e35820ba2daa
# ╠═a69c4ef6-3b04-4101-b1d4-bb5cf7d40233
# ╠═0b3a10e3-b7a2-4e2b-9e88-b8c7bea7834f
# ╟─e7c1cbf4-7939-45a2-a754-b36ae9bdcab4
# ╟─9d10ab94-38d5-406e-8697-7f5bc0c11df2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
