### A Pluto.jl notebook ###
# v0.20.4

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

# ╔═╡ 348240ee-f35f-11ef-01ab-d11599b51ca1
begin
    using PlutoUI, PlutoPlotly, SpecialFunctions, LinearAlgebra, Printf, AssociatedLegendrePolynomials, LaTeXStrings, Bessels



	
end

# ╔═╡ 88195cf5-894b-4676-99bc-a559e0d5ebd9
using ForwardDiff

# ╔═╡ ca388aaf-f515-4c8e-8b39-a7173641dca0
TableOfContents()

# ╔═╡ 2a7c4b45-0153-4095-a67d-b48f8f343d6a
md"""
# Visualization of Spherical Harmonics

Spherical harmonics are fundamental mathematical functions used to describe **wave phenomena** on a sphere. They appear in seismology while modeling Earth's free oscillations. 

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
    """

# ╔═╡ 0e90753b-b2e6-4b41-a5bb-b2dcdab68bd1
md"""
##### **Adjust Spherical Harmonic & Radial Function Parameters**
Use the sliders below to change the properties of the spherical harmonics and radial function.
- **Degree $l$**: Controls the complexity of the spherical harmonic function.
- **Order $m$**: Determines the azimuthal variation (allowed range: $-l \leq m \leq l$).
- **Radial Frequency $ω$**: Affects the oscillations of the radial function.


$( @bind l_value Slider(0:1:20, default=3, show_value=true) ) **Degree**  

$( @bind m_value Slider(-l_value:1:l_value, default=0, show_value=true) ) **Order**  

$( @bind ω_value Slider(range(0, 0.1, length=1000), default=0.05, show_value=true) ) **Radial Frequency**

"""

# ╔═╡ 54b2c235-e693-4f1f-bd1f-1c651e4dbffc
@bind mode Select(["R", "S", "T"])

# ╔═╡ fff4b588-59c3-4cd5-b333-222026b8d666
c = 5.0 # in km/s

# ╔═╡ 6750df2d-459b-4467-b405-1bd8177bee42
md"## Appendix"

# ╔═╡ 4c892f32-ffa3-4b45-a840-6f1a165293b3
md"### Surface Spherical Harmonics "

# ╔═╡ 186c6279-4573-4fda-8589-5fabcbc1cc0e
function Y_lm(l, m, θ, φ)
	P_lm = λlm(l, abs(m), cos(θ))
    return P_lm * exp(im * m * φ)
end

# ╔═╡ b5be41c0-ff78-4f5d-a795-4c5a79135fb6
function Y_lm(l, m, thetaphi)
	return Y_lm(l, m, thetaphi[1], thetaphi[2])
end

# ╔═╡ 7eb05e0a-12e1-4033-b081-1fc79b38218a
function spherical_harmonics(l, m, θ, φ)
        if abs(m) > l
            return zero(θ)
        end
        Y = Y_lm(l, m, θ, φ)
        return Y  # Return the real part for visualization
end

# ╔═╡ 79fa4c7b-5f93-44f9-bf33-311f5efecf32
md"### Vector Spherical Harmonics"

# ╔═╡ ea0d2403-03ea-4c42-a452-59f830c22a16
function spherical_harmonic_gradient(l::Int, m::Int, θ, φ)
    if abs(m) > l
        return (zero(θ), zero(θ))  # Edge case where |m| > l
    end
	Y = Y_lm(l, m, θ, φ)
	
	dY = ForwardDiff.gradient(x->real(Y_lm(l, m, x)), [θ,φ]) + im * ForwardDiff.gradient(x->imag(Y_lm(l, m, x)), [θ,φ])

    return dY
end

# ╔═╡ ca84a03e-053e-47c3-aefc-e81202f0f9f3
spherical_harmonic_gradient(l_value, m_value, 1.2, -0.1)

# ╔═╡ 38e90f9a-5cba-451f-bae3-341f03988e3c
function vector_spherical_harmonics(l, m, θ, φ, mode_type)
        if abs(m) > l
            return zero(θ)
        end
		Y = Y_lm(l, m, θ, φ)
		dY = spherical_harmonic_gradient(l, m, θ, φ)
		s = inv(sqrt(l*(l+1)))
	
 	if mode_type == "R"
        return [Y, 0, 0]  # Purely outward motion (gravity-like)
    elseif mode_type == "S"
        return s.*[0, dY[1], dY[2] / sin(θ)]  # Tangential movement
    elseif mode_type == "T"
        return s.*[0, dY[2] / sin(θ), - dY[1]] # Rotational motion 
    else
        return [0, 0, 0]
    end
end

# ╔═╡ 433d799f-9243-4140-ab1f-006feeaf8816
md"### Grids"

# ╔═╡ 7100ad22-5762-465c-a6a4-bbb06e11ac2f
begin
	    θ_range = range(0, π, length=50)  # Latitude
	    φ_range = range(0, 2π, length=100)  # Longitude
	    θ_grid = first.(Iterators.product(θ_range, φ_range))
		φ_grid = last.(Iterators.product(θ_range, φ_range))
end;

# ╔═╡ 665f39ad-0764-4166-8872-8179eba23374
# Compute spherical harmonics
Ylm_grid = [real(spherical_harmonics(l_value, m_value, θ, φ)) for (θ, φ) in zip(θ_grid, φ_grid)];

# ╔═╡ eaa2c2dd-f109-4a2b-b8b8-4f15d5576ca8
begin
	    # Compute selected VSH function in spherical coordinates
    U_spherical = [vector_spherical_harmonics(l_value, m_value, θ, φ, mode) for (θ, φ) in zip(θ_grid, φ_grid)]
    U_r = [u[1] for u in U_spherical]  # Radial component
    U_θ = [u[2] for u in U_spherical]  # Theta component
    U_φ = [u[3] for u in U_spherical]  # Phi component


    Ux = U_r .* sin.(θ_grid) .* cos.(φ_grid) .+ U_θ .* cos.(θ_grid) .* cos.(φ_grid) .+ U_φ .* -sin.(φ_grid)
    Uy = U_r .* sin.(θ_grid) .* sin.(φ_grid) + U_θ .* cos.(θ_grid) .* sin.(φ_grid) + U_φ .* cos.(φ_grid)
    Uz = U_r .* cos.(θ_grid) + U_θ .* (-sin.(θ_grid)) + U_φ .* 0
	
end;

# ╔═╡ 32fa5f32-72b4-4ce2-8d6b-0ab12acff038
begin
	 # Convert spherical to Cartesian coordinates
	X = sin.(θ_grid) .* cos.(φ_grid) #.* Ylm_grid
	Y = sin.(θ_grid) .* sin.(φ_grid) #.* Ylm_grid
	Z = cos.(θ_grid) #.* Ylm_grid
end;

# ╔═╡ 790e31fe-4cc3-44cd-b33c-0ab0925ce23a
begin
    plot(
        cone(
            x=vec(X), y=vec(Y), z=vec(Z), u=real(vec(Ux)), v=real(vec(Uy)), w=real(vec(Uz)),
            colorscale="Reds", showscale=false, sizeref=(mode == "R" ? 1.0 : 0.15), name="Vector Field"
        )
    , Layout(
        title="$mode <br> l=$l_value, m=$m_value",
        scene=attr(
            xaxis=attr(title="X"),
            yaxis=attr(title="Y"),
            zaxis=attr(title="Z"),
            aspectmode="cube"
        )
    ))
end

# ╔═╡ 71c529bf-77d8-481c-94d3-86a80dfcd641
radius_max = 6000 # in km

# ╔═╡ 8ed29e6e-583c-4fb3-86d3-f34b8d025c0f
radius_grid = range(0., radius_max, length=1000);

# ╔═╡ 4d8bfdc4-175d-4f7a-ad65-338726ad343e
R = Bessels.sphericalbesselj.(l_value, ω_value .* radius_grid ./c);

# ╔═╡ d7ac1aac-2949-4ca3-a0e9-f233a0a6c1b9
begin
	PlutoUI.ExperimentalLayout.hbox([plot(surface(x=X, y=Y, z=Z, surfacecolor=Ylm_grid, colorscale="Seismic", showscale=true), Layout(width=500, title="Spherical Surface Harmonic <br> l=$l_value, m=$m_value")), plot(
scatter(x=R, y=radius_grid, mode="lines", line=attr(color="blue", width=3), name="Radial Function"), Layout(width=200, yaxis=attr(title="radius (km)"),title="Radial Wave Function<br>period = $(round(2*pi/ω_value/60, digits=2)) min"))])
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AssociatedLegendrePolynomials = "2119f1ac-fb78-50f5-8cc0-dda848ebdb19"
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[compat]
AssociatedLegendrePolynomials = "~1.0.1"
Bessels = "~0.2.8"
ForwardDiff = "~0.10.38"
LaTeXStrings = "~1.4.0"
PlutoPlotly = "~0.6.2"
PlutoUI = "~0.7.61"
SpecialFunctions = "~2.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "14265a9cf5fbed589d03875dd12d7c3bb2a4478c"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AssociatedLegendrePolynomials]]
git-tree-sha1 = "3204d769e06c5678b23cf928d850f2f4ad5ec8a5"
uuid = "2119f1ac-fb78-50f5-8cc0-dda848ebdb19"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

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

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

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

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

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

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

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
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

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

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "cc0a5deefdb12ab3a096f00a6d42133af4560d71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

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
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "90af5c9238c1b3b25421f1fdfffd1e8fca7a7133"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.20"

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
git-tree-sha1 = "9ebe25fc4703d4112cc418834d5e4c9a4b29087d"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.2"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "1147f140b4c8ddab224c94efa9569fc23d63ab44"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

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
# ╠═ca388aaf-f515-4c8e-8b39-a7173641dca0
# ╟─2a7c4b45-0153-4095-a67d-b48f8f343d6a
# ╟─0e90753b-b2e6-4b41-a5bb-b2dcdab68bd1
# ╟─d7ac1aac-2949-4ca3-a0e9-f233a0a6c1b9
# ╟─54b2c235-e693-4f1f-bd1f-1c651e4dbffc
# ╟─790e31fe-4cc3-44cd-b33c-0ab0925ce23a
# ╠═4d8bfdc4-175d-4f7a-ad65-338726ad343e
# ╠═fff4b588-59c3-4cd5-b333-222026b8d666
# ╟─6750df2d-459b-4467-b405-1bd8177bee42
# ╠═348240ee-f35f-11ef-01ab-d11599b51ca1
# ╠═88195cf5-894b-4676-99bc-a559e0d5ebd9
# ╟─4c892f32-ffa3-4b45-a840-6f1a165293b3
# ╠═7eb05e0a-12e1-4033-b081-1fc79b38218a
# ╠═186c6279-4573-4fda-8589-5fabcbc1cc0e
# ╠═b5be41c0-ff78-4f5d-a795-4c5a79135fb6
# ╠═665f39ad-0764-4166-8872-8179eba23374
# ╟─79fa4c7b-5f93-44f9-bf33-311f5efecf32
# ╠═ca84a03e-053e-47c3-aefc-e81202f0f9f3
# ╠═ea0d2403-03ea-4c42-a452-59f830c22a16
# ╠═38e90f9a-5cba-451f-bae3-341f03988e3c
# ╠═eaa2c2dd-f109-4a2b-b8b8-4f15d5576ca8
# ╟─433d799f-9243-4140-ab1f-006feeaf8816
# ╠═7100ad22-5762-465c-a6a4-bbb06e11ac2f
# ╠═32fa5f32-72b4-4ce2-8d6b-0ab12acff038
# ╠═71c529bf-77d8-481c-94d3-86a80dfcd641
# ╠═8ed29e6e-583c-4fb3-86d3-f34b8d025c0f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
