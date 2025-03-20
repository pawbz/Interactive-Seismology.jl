### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "Waves On A String"
#> tags = ["introduction"]
#> layout = "layout.jlhtml"
#> description = "Simulates wave propagation on a string"

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9185cd5d-2557-4f2d-9eeb-3d6865ca9331
begin
    using FFTW, PlutoPlotly, PlutoUI, LinearAlgebra
    using ParallelStencil
    using ParallelStencil.FiniteDifferences1D
    using Printf, Statistics
end

# ╔═╡ 236740c7-9c82-4d23-8133-9f7be87f6520
TableOfContents()

# ╔═╡ 88fefb75-ec4f-43b1-b6f5-f28ce49e0dd0
md"""# Waves On A String
This notebook simulates wave propagation on a string using parallel computing techniques. The core function, `model_string`, initializes the simulation environment and defines two parallel functions: `compute_vy!` and `compute_σ!`, which update the velocity and stress fields, respectively. The simulation parameters include the shear modulus (`μ0`), density (`ρ0`), spatial grid (`xgrid`), and temporal grid (tgrid). The notebook sets up the numerical grid, allocates arrays for stress (`σ`) and velocity (`vy`), and initializes the velocity field with a Gaussian distribution. The medium's properties are heterogeneous, with boundary marked by the red line.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India

"""

# ╔═╡ ab615f30-f095-4d4f-a92b-ecd18e69a20a
md"## Medium"

# ╔═╡ f4179789-c4df-405e-92f5-188ebda7a953
begin
    nx = 1000
    xgrid = range(-100, 100, length=nx)
end

# ╔═╡ 4c7b8ce9-1caf-4a81-a2dd-4b4138a6fef7
ρ0 = Float32(3.22 * 10^-3 * 10^15) # density in kg/km3

# ╔═╡ d72c0058-9138-40fa-82d7-f1e1ea068d4c
μ0 = Float32(82 * 10^9 * 10^3)

# ╔═╡ 569c1f97-92d2-4ef2-a23c-aeb131f86475
invρ0 = inv(ρ0)

# ╔═╡ fb985c7e-f809-47d2-b490-4c91cafb164a
invμ0 = inv(μ0)

# ╔═╡ 77246c3a-0f49-4f43-9e4b-67cf3e4be1dd
medium_ref_values = (; μ0, invμ0, ρ0, invρ0)

# ╔═╡ 10014376-cdf6-4ae4-95e4-32c8ea4ef803
vs0 = sqrt.(μ0 ./ ρ0)

# ╔═╡ ad18ed7d-f050-455b-9a13-a2fab7de4215
begin
    courant_number = 0.1

    # lets calculate the min distance from the center to the edge of the domain
    r = min(xgrid[end] - xgrid[1]) * 0.5

    # choose time stepping dt to satisfy Courant condition
    dt = courant_number * step(xgrid) * inv(vs0)
    nt = Int(floor(r / (vs0 * dt))) * 2
    tgrid = range(0, length=nt, step=dt)
    nothing
end;

# ╔═╡ 73aa375b-4887-4608-b699-17112c238d5a
md"""Time $(@bind T Slider(tgrid, show_value=true)) 

Interface Position $(@bind X Slider(xgrid, show_value=true, default=0))"""

# ╔═╡ 82a7a1ab-88da-4c22-a075-7640991bb6ac
md"## Governing Equations"

# ╔═╡ 32a3036a-14d4-429b-a9a6-2eab02cfa1c5
md"""
Velocity $v_y$ and Stress $\sigma$

```math
\rho\,\partial_t v_y = \partial_x \sigma \quad (1)
```
and 
```math
	\partial_t \sigma = \mu\partial_x v_y \quad (2)
```
"""

# ╔═╡ 39bfcff2-f3c1-43af-902c-67e63e2fa61f
md"""
## Simulation
"""

# ╔═╡ dab87c4b-ac71-47fb-b9a5-3cce2e8c9621
@views function model_string(μ0, ρ0, xgrid, tgrid, X)
    @init_parallel_stencil(Threads, Float32, 1)

    @parallel function compute_vy!(vy::Data.Array, σ::Data.Array, dt::Data.Number, ρ::Data.Array, dx::Data.Number)
        @inn(vy) = @inn(vy) - dt / @all(ρ) * (@d(σ) / dx)
        return
    end

    @parallel function compute_σ!(σ::Data.Array, vy::Data.Array, dt::Data.Number, μ::Data.Array, dx::Data.Number)
        @all(σ) = @all(σ) - dt * @all(μ) * (@d(vy) / dx)
        return
    end

    # Derived numerics
    nx = length(xgrid)   # numerical grid resolution; should be a mulitple of 32-1 for optimal GPU perf
    nt = length(tgrid)       # number of timesteps


    dx = Float32(step(xgrid)) # cell size
    dt = Float32(step(tgrid))

    # Array allocations
    σ = @zeros(nx - 1)
    vy = @zeros(nx)

    # Initial conditions
    vy_initial = exp.(-0.11 .* (xgrid .+ 50) .^ 2)
    copyto!(vy, vy_initial)

    # Medium
    μ = @zeros(nx - 1)
    fill!(μ, Float32(μ0))
    ρ = @zeros(nx)
    fill!(ρ, Float32(ρ0))

    # perturb density of the string
    iX = argmin(abs.(xgrid .- X))
    ρ1 = view(ρ, 1:iX)
    rmul!(ρ1, 4.0f0)

    vysave = @zeros(nx + 1, nt)
    # Time loop
    for it = 1:nt
        @parallel compute_σ!(σ, vy, dt, μ, dx)
        @parallel compute_vy!(vy, σ, dt, ρ, dx)

        # save vy
        vys = view(vysave, :, it)
        copyto!(vys, vy)
    end
    return vysave
end

# ╔═╡ 6b69f053-567a-43db-972a-bdd00d0b464d
vy_save = model_string(μ0, ρ0, xgrid, tgrid, X);

# ╔═╡ e470ca4d-4e06-4808-a513-25ef47c125c8
md"## Appendix"

# ╔═╡ 77710384-fc16-4ece-a7d5-62cfddb5f0f9
md"### Plot"

# ╔═╡ b9287111-7b1f-49be-840e-8ecea5922294
default_plotly_template(:plotly_dark)

# ╔═╡ e7660617-4962-4033-97a9-e16f10968cc3
function plot_string(vy_save, X, T)
	iT = argmin(abs.(tgrid .- T))
    fig = Plot(scatter(x=xgrid, y=vy_save[:, iT]), Layout(title="String's vertical displacement at $(round(T, digits=2)) s, with interface at $(round(X, digits=2)) km", width=700, height=300, xaxis=attr(title="Distance"), yaxis=attr(title="Amplitude", range=(-1.2, 1.2))))
    add_vline!(fig,
        X,
        line_color="red", opacity=1,
        layer="below", line_width=2,
    )
    add_vline!(fig,
        xgrid[1],
        line_color="white", opacity=1,
        layer="below", line_width=2,
    )
    add_vline!(fig,
        xgrid[end],
        line_color="white", opacity=1,
        layer="below", line_width=2,
    )
    plot(fig)
end

# ╔═╡ 647f6d79-818f-468f-900c-cd00eb165664
plot_string(vy_save, X, T)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
ParallelStencil = "94395366-693c-11ea-3b26-d9b7aac5d958"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
FFTW = "~1.8.0"
ParallelStencil = "~0.13.2"
PlutoPlotly = "~0.5.0"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "ecb99d172d0b09ad7712d0c1dd55a9d8aa001af5"

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
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "cb25e4b105cc927052c2314f8291854ea59bf70a"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.4"

[[deps.CellArrays]]
deps = ["Adapt", "StaticArrays"]
git-tree-sha1 = "d1c919d285a876522113bec13611255547b1f8fa"
uuid = "d35fcfd7-7af4-4c67-b1aa-d78070614af4"
version = "0.2.1"

    [deps.CellArrays.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

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

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

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

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14eb2b542e748570b56446f4c50fbfb2306ebc45"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

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

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.ParallelStencil]]
deps = ["CellArrays", "MacroTools", "Random", "StaticArrays"]
git-tree-sha1 = "9b5027f48bb900154e4d2dd5002379b2984b981c"
uuid = "94395366-693c-11ea-3b26-d9b7aac5d958"
version = "0.13.2"

    [deps.ParallelStencil.extensions]
    ParallelStencil_AMDGPUExt = "AMDGPU"
    ParallelStencil_CUDAExt = "CUDA"
    ParallelStencil_EnzymeExt = "Enzyme"

    [deps.ParallelStencil.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"

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
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

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

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

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
# ╠═236740c7-9c82-4d23-8133-9f7be87f6520
# ╟─88fefb75-ec4f-43b1-b6f5-f28ce49e0dd0
# ╟─73aa375b-4887-4608-b699-17112c238d5a
# ╟─647f6d79-818f-468f-900c-cd00eb165664
# ╟─ab615f30-f095-4d4f-a92b-ecd18e69a20a
# ╠═f4179789-c4df-405e-92f5-188ebda7a953
# ╠═ad18ed7d-f050-455b-9a13-a2fab7de4215
# ╠═4c7b8ce9-1caf-4a81-a2dd-4b4138a6fef7
# ╠═d72c0058-9138-40fa-82d7-f1e1ea068d4c
# ╠═569c1f97-92d2-4ef2-a23c-aeb131f86475
# ╠═fb985c7e-f809-47d2-b490-4c91cafb164a
# ╠═77246c3a-0f49-4f43-9e4b-67cf3e4be1dd
# ╠═10014376-cdf6-4ae4-95e4-32c8ea4ef803
# ╟─82a7a1ab-88da-4c22-a075-7640991bb6ac
# ╟─32a3036a-14d4-429b-a9a6-2eab02cfa1c5
# ╟─39bfcff2-f3c1-43af-902c-67e63e2fa61f
# ╠═dab87c4b-ac71-47fb-b9a5-3cce2e8c9621
# ╠═6b69f053-567a-43db-972a-bdd00d0b464d
# ╟─e470ca4d-4e06-4808-a513-25ef47c125c8
# ╠═9185cd5d-2557-4f2d-9eeb-3d6865ca9331
# ╟─77710384-fc16-4ece-a7d5-62cfddb5f0f9
# ╠═b9287111-7b1f-49be-840e-8ecea5922294
# ╠═e7660617-4962-4033-97a9-e16f10968cc3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
