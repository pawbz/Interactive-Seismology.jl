### A Pluto.jl notebook ###
# v0.19.25

#> [frontmatter]
#> title = "Born Approximation"
#> description = "Draw subsurface heterogeneities and model their linearized response using the notion of scattering theory."

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8acbffaf-1811-4592-a2d1-a8f561242d85
begin
    using Symbolics
    using SymbolicUtils
    using Latexify
    using FFTW
    using PlutoPlotly
    using PlutoUI
    using PlutoTeachingTools
    using SpecialFunctions
    using Einsum
    using Tullio
end

# ╔═╡ c0e6e258-7fcb-4bd8-9931-be0203c1adc9
ChooseDisplayMode()

# ╔═╡ 4aa9e374-27a1-4d80-9d7f-9a7c1ee859b2
TableOfContents()

# ╔═╡ 5d3692af-dee6-4adf-9276-82f11a1a9544
md"""
# Born Approximation
In this notebook, we construct an integral equation to model scattering due to subsurface heterogeneities using the notion of scattering theory. We assume that the medium may be decomposed into a known reference wave speed profile plus a perturbation called the scatterer. The wavefield similarly may be decomposed into a reference wavefield plus a scattered/perturbed field.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)
ES218; August 2022

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 3d88f8c4-ba31-4d99-8966-3ddf617e5b5f
@bind slowness_pert_draw_input HTML("""
<div id=parent>
	<canvas id=canvas width=500px height=250px></canvas>
	<button id=clearButton>clear</button>
</div>
	
<script>
	const canvasWidth = 500, canvasHeight = 250, background = "#f1f1f1";
	
	const parentDiv = currentScript.previousElementSibling
	const c = parentDiv.querySelector("canvas")
	const ctx = c.getContext("2d");

	ctx.fillStyle = background;
	ctx.fillRect(0, 0, canvasWidth, canvasHeight);

	ctx.fillStyle = "#010101";
	ctx.fillRect(200, 100, 20, 20);

	parentDiv.value = $(vec([[x,y] for x in 200:220, y in 100:120]));

	let drawing = false;
		
	c.addEventListener('mousedown', () => drawing = true);
	c.addEventListener('mouseup', () => drawing = false);
	c.addEventListener('mousemove', (e) => {
		if(drawing) {
			ctx.beginPath();
			ctx.arc(e.offsetX, e.offsetY, 4, 0, 2 * Math.PI);
			ctx.fillStyle = "#010101";
			ctx.fill();
				
			parentDiv.value.push([e.offsetX, (canvasHeight - e.offsetY)]);
			parentDiv.dispatchEvent(new CustomEvent("input"));
		}
	});
	
	function clearCanvas(e) {
		ctx.fillStyle = background;
		ctx.fillRect(0, 0, canvasWidth, canvasHeight);
		parentDiv.value = [];
		parentDiv.dispatchEvent(new CustomEvent("input"));
	}
	
	parentDiv.querySelector("#clearButton").addEventListener('click', clearCanvas);
</script>
""")

# ╔═╡ 4972f23b-88a8-49d9-acad-75a65bdbe101
md"""
## Draw Diffraction Hyperbolas!
The canvas above has a distance of $x∈[0, 1000]\,$m and depth $z∈[0, 500]\,$m. 
The seismic source is located near the surface at $[0, -100]$, and $100$ receivers are equispaced from $[0, -10]$ to $[1000, -10]$ parallel to the surface. You can draw (using mouse clicks) on the canvas to place the scatterers and observe the scattered wavefield below.
"""

# ╔═╡ 30d354da-fceb-42a3-8c4c-59bc211e0022
md"## Scalar Helmholtz Equation"

# ╔═╡ 802eb9d4-98da-4ac8-88de-36365227a971
@syms x z

# ╔═╡ 90c80298-cc7a-4567-bfb2-383570d9e4b1
@syms t

# ╔═╡ deed0f33-0d14-4d28-ac39-e2a75468f696
Dx = Differential(x)

# ╔═╡ 0801d1fa-a22b-4341-b594-245209682f82
Dz = Differential(z)

# ╔═╡ 3a8c62bb-6ecc-4ab7-8e0f-68b5c421767e
Dt = Differential(t)

# ╔═╡ 58ee01e5-bdd5-42de-9670-e2cf3013efec
@variables ρ

# ╔═╡ d08e57c9-8180-4ec1-8a96-4c44170bad5b
@syms ı

# ╔═╡ 36c73570-8cbd-447d-8d1d-20953903ece8
@syms ω # angular frequency

# ╔═╡ e197b835-cf69-46e5-8d0f-a2aa5d02ad04
L(u, s) = ω^2 * s * u + Dx(Dx(u)) + Dz(Dz(u))

# ╔═╡ c4479321-2f68-4359-aa5b-99441a317efe
@syms s(𝐱) δs(𝐱) # reference and perturbed slowness

# ╔═╡ 22b83e2f-a634-441d-a737-279497f08a05
@syms U(𝐱, 𝐱ₛ) # reference wavefield

# ╔═╡ c27bf1a2-a50e-423d-8a81-00272f3069ca
@syms δU(𝐱, 𝐱ₛ) # scattered wavefield

# ╔═╡ 8762984a-eaba-4b3d-9fed-7fc64a68c11a
@syms 𝐱 # spatial coordinate vector

# ╔═╡ d7b8f7e9-49f6-4a7f-b67d-83f00753ef76
@syms 𝐱ₛ # source coordinate vector

# ╔═╡ 75ebf7f5-2f25-4be2-8da4-bc1a483e71be
# TODO: get Helmholtz eq with the following substitution 
U(𝐱, 𝐱ₛ) * exp(ı * ω * t)

# ╔═╡ f478e161-06b3-4a91-b45b-f1e3434e1832
@syms f(t) # source time function

# ╔═╡ ebb0b8c5-94fa-4016-be1f-badf88c09329
@syms F(ω) # source in the frequency domain

# ╔═╡ 703ed5fa-4e67-442d-8423-041b417b9a9e
@syms δ(𝐱) # Dirac delta function

# ╔═╡ bc542855-8281-4b61-b8d9-d758d15f4dcd
L(U(𝐱, 𝐱ₛ), s(𝐱)) ~ F(ω) * δ(𝐱 - 𝐱ₛ)

# ╔═╡ bc57a10f-9d2d-4284-baf6-0875e169dbb4
md"## Perturbation Theory"

# ╔═╡ d71383ea-9bc7-4c35-bbee-be086168f6d9
TwoColumn(
    md"""
 #### Reference State
 Has a medium with slowness
 ```julia
 s(𝐱)
 ```
 and the wavefield 
 ```julia
 U(𝐱,𝐱ₛ)
 ```
 that is initiated at $t=0$ and satisfies
 ```julia
 L(U, s) ~ F(ω)δ(𝐱-𝐱ₛ).
 ```
 	""",
    md"""
#### Perturbed State
Has a perturbed medium with slowness
```julia
s(𝐱)+δs(𝐱)
```
and the perturbed wavefield 
```julia
U(𝐱)+δU(𝐱)
```
that is initiated at $t=0$ and satisfies
```julia
L(U+δU, s+δs) ~ F(ω)δ(𝐱-𝐱ₛ).
```	
"""
)

# ╔═╡ 401baadb-744b-4536-8f07-0ea0f39bde08
L(U(𝐱, 𝐱ₛ) + δU(𝐱, 𝐱ₛ), s(𝐱) + δs(𝐱)) ~ F(ω) * δ(𝐱 - 𝐱ₛ)

# ╔═╡ 7242025b-767c-49f7-a328-8abc85232d6d
scat_wave_eq = L(δU(𝐱, 𝐱ₛ), s(𝐱)) ~ expand(simplify(expand_derivatives(L(U(𝐱, 𝐱ₛ) + δU(𝐱, 𝐱ₛ), s(𝐱) + δs(𝐱)) - L(U(𝐱, 𝐱ₛ), s(𝐱)) - L(δU(𝐱, 𝐱ₛ), s(𝐱)))))

# ╔═╡ 758c029d-68ac-4f1f-8419-f32f7789ab44
expand(scat_wave_eq.rhs)

# ╔═╡ 64eeabf7-2549-4e1d-ac34-df8ac53ea062
md"## Born Approximation"

# ╔═╡ 8e066122-85a6-4493-8642-a39728f9d59b
born_scat_wave_eq = scat_wave_eq.lhs ~ substitute(scat_wave_eq.rhs, ω^2 * δs(𝐱) * δU(𝐱, 𝐱ₛ) => 0)

# ╔═╡ e27ee2ce-cf3c-4b59-be8f-dff23611b619
TwoColumn(
    md"""
   #### Perturbed State
   Has a perturbed medium with slowness
   ```julia
   s(𝐱)
   ```
   and the perturbed wavefield 
   ```julia
   u(𝐱,𝐱ₛ)+δu(𝐱,𝐱ₛ)
   ```
   that is initiated at $t=0$ and satisfies
   ```julia
   L(u+δu, s+δs) ~ F(ω)δ(𝐱-𝐱ₛ).
   ```
   	""",
    md"""
   #### Impulsive State
   Has a perturbed medium with slowness
   ```julia
   s(𝐱)
   ```
   and the perturbed wavefield 
   ```julia
   G(𝐱,𝐱ᵣ)
   ```
   that is initiated at $t=0$ and satisfies
   ```julia
   L(G, s) ~ F(ω)δ(𝐱-𝐱ᵣ).
   ```
   	"""
)

# ╔═╡ 503d28af-6d8a-4919-9aac-85d5dba917e7
born_scat_wave_eq

# ╔═╡ efdecf55-f05c-402d-a860-9c20f369d1ae
# receiver coordinate
@syms 𝐱ᵣ

# ╔═╡ 896660b8-51a6-48f0-9997-72cc3e65c3c0
# Green's function (evaluated at 𝐱, and source at 𝐱ᵣ)
@syms G(𝐱, 𝐱ᵣ)

# ╔═╡ 97f3e97d-5139-4fdb-a3f7-1af86dddd4f4
Green_wave_eq = L(G(𝐱, 𝐱ᵣ), s(𝐱)) ~ δ(𝐱 - 𝐱ᵣ)

# ╔═╡ 2e2dde7f-8120-41bb-a2f0-98a2bfe536e8
ex1 = simplify(expand(Green_wave_eq.lhs * δU(𝐱, 𝐱ₛ) - born_scat_wave_eq.lhs * G(𝐱, 𝐱ᵣ)))

# ╔═╡ 1067f104-ac59-41af-945f-6837ee196778
ex2 = Green_wave_eq.rhs * δU(𝐱, 𝐱ₛ) - born_scat_wave_eq.rhs * G(𝐱, 𝐱ᵣ)

# ╔═╡ b2fefe0e-5d50-4bd0-80f0-6c096943ec53
# volume integral
@syms ∫ᵥ(expression)

# ╔═╡ 39d7b65b-aefa-441d-a47d-59855f4095b3
md"""
Here, we consider only volume scatterers, i.e., we simply consider an unbounded medium. In this case, using the Sommerfeld radiation condition and Green's theorem, the integral below vanishes.
"""

# ╔═╡ 6b991a44-3f54-4018-a115-772cdd4a6e45
∫ᵥ(ex1) ~ 0

# ╔═╡ fae3d3d8-37ee-4a8a-80ca-baa814f8c874
md"Finally, the integral equation to model the scattered field is given below."

# ╔═╡ 62bcf4c5-fc31-47ef-ac38-b0905f3ab0eb
scattered_wavefield = δU(𝐱ᵣ, 𝐱ₛ) ~ ∫ᵥ(arguments(ex2)[1]) ~ ∫ᵥ(-arguments(ex2)[2])

# ╔═╡ b48cf5ec-a033-46e3-aec8-da5864b2386a
md"## Appendix"

# ╔═╡ 2d82c01f-2be9-4628-b926-bab17d451ed8
δx = 2 # spatial sampling interval

# ╔═╡ 33a221fe-ab5d-4a41-80bd-105746f949fb
begin
    nx = 500
    nz = 250
end

# ╔═╡ 2aa9168a-cadb-4060-9eac-d93a1cf3bf0b
# return (x,z) for the element I of user input, assuming spatial sampling ds
function get_location(I, δ)
    [(I[1] - 1) * δ, (nz - I[2] + 1) * δ]
end

# ╔═╡ baeaecf7-ad9c-42d5-8946-3c336d28fb8b
slocs = [[0, -100.0]] # source location(s)

# ╔═╡ 4389a492-2096-4e2a-b05f-91e9820e15a4
rlocs = [[rx, -10.0] for rx in range(0, stop=(500 - 1) * δx, length=100)] # receiver locations

# ╔═╡ d342ea67-3f1f-421d-8850-09c72b0398a4
# time grid
tgrid = range(0, stop=1.5, step=0.005);

# ╔═╡ 968fb80c-0f44-4198-9aa0-87a07f28c6ca
# ... corresponding frequency grid
freqgrid = collect(rfftfreq(length(tgrid), inv(step(tgrid))));

# ╔═╡ 504f09a0-6d86-459a-9b91-f1c8dd684427
begin
    # a Gaussian spectrum for source with peak at 25 Hz
    Fsource = exp.(-abs2.(freqgrid .- 25.0) * 1e-2)
    # remove zero frequency
    Fsource[1] = 0.0
end;

# ╔═╡ 3ebd625f-49dc-4f5e-bcb0-2633092f5062
plot(scatter(x=freqgrid, y=Fsource), Layout(width=400, height=200, title="Source Spectrum", xaxis_title="Frequency [Hz]"),)

# ╔═╡ 9404b52d-c571-4e23-8dda-d20d10de0457
md"### 2-D Green's Function"

# ╔═╡ f78d9f96-0d4c-445b-a433-06f1a9a9eb12
md"### Modeling"

# ╔═╡ 135efe98-672e-4013-8e4c-15f0a29e9bab
# reference P velocity
vp0 = 2500

# ╔═╡ 094b0ad6-1ec0-407a-884f-7147d8a2ec9f
# reference mass density
rho0 = 2000

# ╔═╡ ab29f4e1-cc96-46a3-b0b1-9592378d5c64
# sloc: source location
# rloc: receiver location
# f: frequency 
# returns Green's function for a 
# 2-D homogeneous medium with velocity c
function G0(rloc, sloc, f, c)
    # ω/c wavenumber
    k = 2 * pi * f * inv(c)
    # r distance between source and receiver
    r = sqrt(sum(abs2.(sloc .- rloc)))
    return -0.25 * rho0 * im * hankelh2(0, k * r)
end

# ╔═╡ e5a6e125-c700-4471-9781-abbd0b1c49c7
function get_reference_wavefield(freqgrid, Fsource, rlocs, slocs,)
    @tullio D[iω, ir, is] := G0(rlocs[ir], slocs[is], freqgrid[iω], vp0) * Fsource[iω]
    # remove zero frequencies
    @tullio D[1, i, j] = complex(0.0)
    # transform to time domain
    return irfft(D, length(tgrid), 1,)
end

# ╔═╡ 9c8ce20e-171a-4909-9f3d-c2fe233dda7b
d = get_reference_wavefield(freqgrid, Fsource, rlocs, slocs,);

# ╔═╡ 99fd676e-d0a5-4975-9a0a-82d6438a3143
# δs is the perturbation in slowness
# δx is the spatial sampling used to scale slowness_pertindices
# rlocs and slocs are receiver and source positions
function get_scattered_wavefield(δs, δx, freqgrid, Fsource, slowness_pertindices, rlocs, slocs)
    if(slowness_pertindices == []) 
		return zeros(length(tgrid), length(rlocs), length(slocs))
	else(ismissing(slowness_pertindices))
		# default box in UI
		slowness_pert_indices = vec([[x,y] for x in 200:220, y in 100:120])
	end
    @tullio D[iω, ir, is] := G0(get_location(slowness_pertindices[i], δx), slocs[is], freqgrid[iω], vp0) * G0(rlocs[ir], get_location(slowness_pertindices[i], δx), freqgrid[iω], vp0) * freqgrid[iω] * freqgrid[iω] * δs * Fsource[iω] * 4.0 * pi * pi
    # remove zero frequencies
    @tullio D[1, i, j] = complex(0.0)
    # transform to time domain
    return irfft(D, length(tgrid), 1,)
end

# ╔═╡ e5ee5c1f-f056-47de-9759-077777c6a157
δd = get_scattered_wavefield(1e-6, δx, freqgrid, Fsource, slowness_pert_draw_input, rlocs, slocs);

# ╔═╡ 3c23a484-0b9a-4044-bcec-fec447509991
md"### Plots"

# ╔═╡ a6f59c95-2271-4f2f-894f-b574401cd545
function plot_data(d, δd)

    fig = Plot(Layout(yaxis_autorange="reversed", height=400, width=500, yaxis_title="time [s]", Subplots(shared_xaxes=true, shared_yaxes=true, rows=1, cols=2, subplot_titles=["Reference Wavefield" "Scattered Wavefield"])))
    add_trace!(fig, heatmap(y=tgrid, z=d, showscale=false, colorscale="jet"), row=1, col=1)
    add_trace!(fig, heatmap(y=tgrid, z=δd, showscale=false, colorscale="jet"), row=1, col=2)

    return PlutoPlotly.plot(fig)

end

# ╔═╡ 5b60b72c-73a4-4a90-b50d-c3bd3c8654d5
plot_data(d[:, :, 1], δd[:, :, 1])

# ╔═╡ 98521f53-fbf4-4ba2-b0b9-629f21478a31
md"## References"

# ╔═╡ 8644be64-5901-47df-8840-b234e4ce4a01
md"""
- [Stanford Notes](http://sepwww.stanford.edu/public/docs/sep131/rgunther1/paper_html/node3.html)
- [SEG Wiki](https://wiki.seg.org/wiki/Born-approximate_modeling_formula)
"""

# ╔═╡ f3b78d8f-9226-47e8-beca-38956260a40b
md"""## TODO
- refine the derivation
- Disconnect b/w symbolic derivation and numerical simulation, is it inevitable here?
- Elastic Wavefield w/ density scatterers
- Radiation of velocity and density scatterers
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Einsum = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
Tullio = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"

[compat]
Einsum = "~0.4.1"
FFTW = "~1.5.0"
Latexify = "~0.15.17"
PlutoPlotly = "~0.3.6"
PlutoTeachingTools = "~0.2.3"
PlutoUI = "~0.7.48"
SpecialFunctions = "~2.1.7"
SymbolicUtils = "~0.19.11"
Symbolics = "~4.13.0"
Tullio = "~0.3.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "36d471d82421f1898c5124672d8c091a5e6036e3"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "df23d15b1090a3332a09a7a51da45bd9f0a07f92"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.8"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0310e08cb19f5da31d08341c6120c047598f5b9c"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.5.0"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "SnoopPrecompile", "Static"]
git-tree-sha1 = "dedc16cbdd1d32bead4617d27572f582216ccf23"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.25"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "f12dc65aef03d0a49650b20b2fdaf184928fd886"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.5"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "93c8ba53d8d26e124a5a8d4ec914c3a16e6a0970"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.3"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "7fe6d92c4f281cf4ca6f2fba0ce7b299742da7ca"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.37"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "0e5c14c3bb8a61b3d53b2c0620570c332c8d0663"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.2.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "74911ad88921455c6afcad1eefa12bd7b1724631"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.80"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "988e2db482abeb69efc76ae8b6eba2e93805ee70"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.15"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "d0fa82f39c2a5cdb3ee385ad52bc05c42cb4b9f0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.5"

[[deps.Einsum]]
deps = ["Compat"]
git-tree-sha1 = "4a6b3eee0161c89700b6c1949feae8b851da5494"
uuid = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
version = "0.4.1"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "d3ba08ab64bdfd27234d3f61956c966266757fe6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.7"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "a5e6e7f12607e90d71b09e6ce2c965e41b337968"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.1"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "57f7cde02d7a53c9d1d28443b9f11ac5fbe7ebc9"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.3"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "47f0f03eddecd7ad59c42b1dd46d5f42916aff63"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.11"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "72ab280d921e8a013a83e64709f66bc3e854b2ed"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.20"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "0a92979c14dfa71adbf892f0cd073e34b7189197"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.13.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "45b288af6956e67e621c5cbb2d75a261ab58300b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.20"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0f39bc7f71abdff12ead4fc4a7d998fb2f3c171f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.5"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "4d5917a26ca33c66c8e5ca3247bd163624d35493"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "393fc4d82a73c6fe0e2963dd7c882b09257be537"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "aa532179d4a643d4bd9f328589ca01fa20a0d197"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.1.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "8175fc2b118a3755113c8e68084dc1a9e63c61ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Colors", "Dates", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "PlotlyBase", "PlutoUI", "Reexport"]
git-tree-sha1 = "dec81dcd52748ffc59ce3582e709414ff78d947f"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.3.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "ea3e4ac2e49e3438815f8946fa7673b658e35bdb"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff"]
git-tree-sha1 = "758f3283aba57c53960c8e1900b4c724bf24ba74"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.8"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "de191bc385072cc6c7ed3ffdc1caeed3f22c74d4"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.7.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "ZygoteRules"]
git-tree-sha1 = "f311e004143b4dc7c5492a2e9b9a1d945058fd8c"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.36.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "fd5dba2f01743555d8435f7c96437b29eae81a17"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "50314d2ef65fce648975a8e80ae6d8409ebbf835"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.5"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "c4e3c4bb8544b2f0bc537baa44925b6b39a8c383"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.81.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "c35b107b61e7f34fa3f124026f2a9be97dea9e1c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "6b764c160547240d868be4e961a5037f47ad7379"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "027b43d312f6d52187bb16c2d4f0588ddb8c4bb2"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.11"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "718328e81b547ef86ebf56fbc8716e6caea17b00"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.13.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "34e6bcf36b9ed5d56489600cf9f3c16843fa2aa2"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.11"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f2fd3f288dfc6f507b0c3a2eb3bac009251e548b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.22"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c42fa452a60f022e9e087823b47e5a5f8adc53d5"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.75"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.Tullio]]
deps = ["ChainRulesCore", "DiffRules", "LinearAlgebra", "Requires"]
git-tree-sha1 = "7871a39eac745697ee512a87eeff06a048a7905b"
uuid = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"
version = "0.3.5"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

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
version = "1.2.12+3"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═c0e6e258-7fcb-4bd8-9931-be0203c1adc9
# ╠═4aa9e374-27a1-4d80-9d7f-9a7c1ee859b2
# ╟─5d3692af-dee6-4adf-9276-82f11a1a9544
# ╟─3d88f8c4-ba31-4d99-8966-3ddf617e5b5f
# ╟─4972f23b-88a8-49d9-acad-75a65bdbe101
# ╠═5b60b72c-73a4-4a90-b50d-c3bd3c8654d5
# ╠═9c8ce20e-171a-4909-9f3d-c2fe233dda7b
# ╠═e5ee5c1f-f056-47de-9759-077777c6a157
# ╟─30d354da-fceb-42a3-8c4c-59bc211e0022
# ╠═75ebf7f5-2f25-4be2-8da4-bc1a483e71be
# ╠═e197b835-cf69-46e5-8d0f-a2aa5d02ad04
# ╠═bc542855-8281-4b61-b8d9-d758d15f4dcd
# ╠═802eb9d4-98da-4ac8-88de-36365227a971
# ╠═90c80298-cc7a-4567-bfb2-383570d9e4b1
# ╠═deed0f33-0d14-4d28-ac39-e2a75468f696
# ╠═0801d1fa-a22b-4341-b594-245209682f82
# ╠═3a8c62bb-6ecc-4ab7-8e0f-68b5c421767e
# ╠═58ee01e5-bdd5-42de-9670-e2cf3013efec
# ╠═d08e57c9-8180-4ec1-8a96-4c44170bad5b
# ╠═36c73570-8cbd-447d-8d1d-20953903ece8
# ╠═c4479321-2f68-4359-aa5b-99441a317efe
# ╠═22b83e2f-a634-441d-a737-279497f08a05
# ╠═c27bf1a2-a50e-423d-8a81-00272f3069ca
# ╠═8762984a-eaba-4b3d-9fed-7fc64a68c11a
# ╠═d7b8f7e9-49f6-4a7f-b67d-83f00753ef76
# ╠═f478e161-06b3-4a91-b45b-f1e3434e1832
# ╠═ebb0b8c5-94fa-4016-be1f-badf88c09329
# ╠═703ed5fa-4e67-442d-8423-041b417b9a9e
# ╟─bc57a10f-9d2d-4284-baf6-0875e169dbb4
# ╟─d71383ea-9bc7-4c35-bbee-be086168f6d9
# ╠═401baadb-744b-4536-8f07-0ea0f39bde08
# ╠═7242025b-767c-49f7-a328-8abc85232d6d
# ╠═758c029d-68ac-4f1f-8419-f32f7789ab44
# ╟─64eeabf7-2549-4e1d-ac34-df8ac53ea062
# ╠═8e066122-85a6-4493-8642-a39728f9d59b
# ╟─e27ee2ce-cf3c-4b59-be8f-dff23611b619
# ╠═503d28af-6d8a-4919-9aac-85d5dba917e7
# ╠═efdecf55-f05c-402d-a860-9c20f369d1ae
# ╠═896660b8-51a6-48f0-9997-72cc3e65c3c0
# ╠═97f3e97d-5139-4fdb-a3f7-1af86dddd4f4
# ╠═2e2dde7f-8120-41bb-a2f0-98a2bfe536e8
# ╠═1067f104-ac59-41af-945f-6837ee196778
# ╠═b2fefe0e-5d50-4bd0-80f0-6c096943ec53
# ╟─39d7b65b-aefa-441d-a47d-59855f4095b3
# ╠═6b991a44-3f54-4018-a115-772cdd4a6e45
# ╟─fae3d3d8-37ee-4a8a-80ca-baa814f8c874
# ╠═62bcf4c5-fc31-47ef-ac38-b0905f3ab0eb
# ╟─b48cf5ec-a033-46e3-aec8-da5864b2386a
# ╠═8acbffaf-1811-4592-a2d1-a8f561242d85
# ╠═2d82c01f-2be9-4628-b926-bab17d451ed8
# ╠═33a221fe-ab5d-4a41-80bd-105746f949fb
# ╠═2aa9168a-cadb-4060-9eac-d93a1cf3bf0b
# ╠═baeaecf7-ad9c-42d5-8946-3c336d28fb8b
# ╠═4389a492-2096-4e2a-b05f-91e9820e15a4
# ╠═d342ea67-3f1f-421d-8850-09c72b0398a4
# ╠═968fb80c-0f44-4198-9aa0-87a07f28c6ca
# ╠═504f09a0-6d86-459a-9b91-f1c8dd684427
# ╠═3ebd625f-49dc-4f5e-bcb0-2633092f5062
# ╟─9404b52d-c571-4e23-8dda-d20d10de0457
# ╠═ab29f4e1-cc96-46a3-b0b1-9592378d5c64
# ╟─f78d9f96-0d4c-445b-a433-06f1a9a9eb12
# ╠═135efe98-672e-4013-8e4c-15f0a29e9bab
# ╠═094b0ad6-1ec0-407a-884f-7147d8a2ec9f
# ╠═e5a6e125-c700-4471-9781-abbd0b1c49c7
# ╠═99fd676e-d0a5-4975-9a0a-82d6438a3143
# ╟─3c23a484-0b9a-4044-bcec-fec447509991
# ╠═a6f59c95-2271-4f2f-894f-b574401cd545
# ╟─98521f53-fbf4-4ba2-b0b9-629f21478a31
# ╟─8644be64-5901-47df-8840-b234e4ce4a01
# ╟─f3b78d8f-9226-47e8-beca-38956260a40b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
