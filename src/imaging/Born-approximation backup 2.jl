### A Pluto.jl notebook ###
# v0.20.19

#> [frontmatter]
#> title = "Born Approximation"
#> tags = ["imaging"]
#> layout = "layout.jlhtml"
#> description = "Draw subsurface heterogeneities and model their linearized response using the notion of scattering theory."

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

# ╔═╡ 8acbffaf-1811-4592-a2d1-a8f561242d85
begin
    using Symbolics
    using SymbolicUtils
    using Latexify
    using FFTW
    using PlutoPlotly
    using PlutoUI
    using PlutoTeachingTools
    using Einsum
    using MLUtils
    using PlutoLinks, PlutoHooks
	using Bessels
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
Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 99d5ad02-2897-4f2c-9a80-a146f90e38dd
md"""
---
Draw on the canvas below to place the scatterers and interact with this notebook.
"""

# ╔═╡ 3d88f8c4-ba31-4d99-8966-3ddf617e5b5f
@bind slowness_pert_draw_input HTML("""
<div id=parent>
	<canvas id=canvas width=200px height=100px></canvas>
	<button id=clearButton>clear</button>
</div>
	
<script>
	const canvasWidth = 200, canvasHeight = 100, background = "#f1f1f1";
	
	const parentDiv = currentScript.previousElementSibling
	const c = parentDiv.querySelector("canvas")
	const ctx = c.getContext("2d");

	ctx.fillStyle = background;
	ctx.fillRect(0, 0, canvasWidth, canvasHeight);

	ctx.fillStyle = "#010101";
	ctx.fillRect(50, 50, 5, 5);

	parentDiv.value = $(vec([[x,y] for x in 50:55, y in 50:55]));

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

# ╔═╡ 9d498284-f702-4f6b-aace-54743c7df206
TwoColumn(md"""
$(@bind it_plot Clock(0.1))
	""", md"""$(@bind reset_anim CounterButton("Reset Animation"))
 """)

# ╔═╡ 30d354da-fceb-42a3-8c4c-59bc211e0022
md"## Scalar Helmholtz Equation"

# ╔═╡ 75ebf7f5-2f25-4be2-8da4-bc1a483e71be
# TODO: get Helmholtz eq with the following substitution 
U(𝐱, 𝐱ₛ) * exp(ı * ω * t)

# ╔═╡ 802eb9d4-98da-4ac8-88de-36365227a971
@syms x z

# ╔═╡ 90c80298-cc7a-4567-bfb2-383570d9e4b1
@syms t

# ╔═╡ deed0f33-0d14-4d28-ac39-e2a75468f696
Dx = Differential(x)

# ╔═╡ 0801d1fa-a22b-4341-b594-245209682f82
Dz = Differential(z)

# ╔═╡ e197b835-cf69-46e5-8d0f-a2aa5d02ad04
L(u, s) = ω^2 * s * u + Dx(Dx(u)) + Dz(Dz(u))

# ╔═╡ bc542855-8281-4b61-b8d9-d758d15f4dcd
L(U(𝐱, 𝐱ₛ), s(𝐱)) ~ F(ω) * δ(𝐱 - 𝐱ₛ)

# ╔═╡ 3a8c62bb-6ecc-4ab7-8e0f-68b5c421767e
Dt = Differential(t)

# ╔═╡ 58ee01e5-bdd5-42de-9670-e2cf3013efec
@variables ρ

# ╔═╡ d08e57c9-8180-4ec1-8a96-4c44170bad5b
@syms ı

# ╔═╡ 36c73570-8cbd-447d-8d1d-20953903ece8
@syms ω # angular frequency

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

# ╔═╡ f478e161-06b3-4a91-b45b-f1e3434e1832
@syms f(t) # source time function

# ╔═╡ ebb0b8c5-94fa-4016-be1f-badf88c09329
@syms F(ω) # source in the frequency domain

# ╔═╡ 703ed5fa-4e67-442d-8423-041b417b9a9e
@syms δ(𝐱) # Dirac delta function

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

# ╔═╡ db0b5eb4-406d-4d8c-b0ac-544c1d6f37b1
md"## Modeling"

# ╔═╡ b48cf5ec-a033-46e3-aec8-da5864b2386a
md"## Appendix"

# ╔═╡ 3affdcaf-8b32-496e-81a3-a898578948d2
import PlutoUIExtra

# ╔═╡ 094b0ad6-1ec0-407a-884f-7147d8a2ec9f
# reference mass density
rho0 = 2000.0f0

# ╔═╡ 135efe98-672e-4013-8e4c-15f0a29e9bab
# reference P velocity
vp0 = 2500.0f0

# ╔═╡ 953e46b3-b507-4ff6-a24b-763250917da3
md"### Source Bandwidth"

# ╔═╡ d342ea67-3f1f-421d-8850-09c72b0398a4
# time grid
tgrid = range(0, stop=1.5, step=0.005);

# ╔═╡ 968fb80c-0f44-4198-9aa0-87a07f28c6ca
# ... corresponding frequency grid
freqgrid = Float32.(collect(rfftfreq(length(tgrid), inv(step(tgrid)))));

# ╔═╡ 8962a7d6-1e29-4d99-b8d3-aecd9d69de76
kgrid = Float32.(2.0f0 * pi * freqgrid * inv(vp0)) .+ 1.0f-6;

# ╔═╡ 504f09a0-6d86-459a-9b91-f1c8dd684427
begin
    # a Gaussian spectrum for source with peak frequency
    fpeak = 50.0f0 # Hz
    Fsource = exp.(-abs2.(freqgrid .- fpeak) * 1.0f-3)
    # remove zero frequency
    Fsource[1] = 0.0f0
end;

# ╔═╡ 3ebd625f-49dc-4f5e-bcb0-2633092f5062
plot(scatter(x=freqgrid, y=Fsource), Layout(width=400, height=200, title="Source Spectrum", xaxis_title="Frequency [Hz]"),)

# ╔═╡ 2d82c01f-2be9-4628-b926-bab17d451ed8
begin
    nxUI = 200
    nzUI = 100
    δxUI = 4f0 # spatial sampling interval
    xgridUI = Float32.(collect(range(0.0, step=δxUI, length=nxUI)))
    zgridUI = Float32.(collect(range(0, step=δxUI, length=nzUI)))
end;

# ╔═╡ acba3f10-975b-423d-9ab8-ead6d9f9774b
begin
    slowness_grid_input = zeros(Float32, nzUI, nxUI)
    map(slowness_pert_draw_input) do I
        slowness_grid_input[nzUI-I[2]+1, I[1]] += 1.0f-8
    end
end;

# ╔═╡ 227a30d4-9cec-440a-8238-87cd0d560e33
paUI = (; tgrid, vp0, rho0, freqgrid, kgrid, Fsource, xgrid=xgridUI, zgrid=zgridUI)

# ╔═╡ 33a221fe-ab5d-4a41-80bd-105746f949fb
begin
    δxplot = 11f0
    xgrid = Float32.(collect(range(first(xgridUI), stop=last(xgridUI), step=δxplot)))
    zgrid = Float32.(collect(range(first(xgridUI), stop=last(zgridUI), step=δxplot)))
    nxgrid = length(xgrid)
    nzgrid = length(zgrid)
end;

# ╔═╡ 07a5f979-f83d-45b2-9e19-0c770577f09e
pagrid = (; tgrid, vp0, rho0, freqgrid, kgrid, Fsource, xgrid=xgridUI .+ 1.0f-3, zgrid=zgridUI .+ 1.0f-3)

# ╔═╡ 855bbfe1-39bd-46b0-bafd-be63cc5f67b2
md"### Source Receiver"

# ╔═╡ 0b4b74d4-2fd1-4a52-89c0-4141d047b55a
rlocs_x_grid = last.(collect(Iterators.product(zgrid, xgrid)));

# ╔═╡ d856645d-39bc-4603-837f-3c39d9df64af
rlocs_z_grid = first.(collect(Iterators.product(zgrid, xgrid)));

# ╔═╡ 2aa9168a-cadb-4060-9eac-d93a1cf3bf0b
# return (x,z) for the element I of user input, assuming spatial sampling ds
function get_xlocation(I, δ)
    Float32.((I[1] - 1) * δ)
end

# ╔═╡ 453fcb4f-1cc6-4c57-9bec-2e1d7c76244a
function get_zlocation(I, δ, nzUI)
    Float32.((nzUI - I[2] + 1) * δ)
end

# ╔═╡ 2ac0ea94-aa5d-43da-a240-db4a57dda204
nr = 50;

# ╔═╡ 6e74479f-96ce-4af7-ae7d-603fe32ba88f
md"Select Receivers $(@bind rUI MultiCheckBox(1:nr, select_all=true, default=collect(1:nr)))"

# ╔═╡ 4972f23b-88a8-49d9-acad-75a65bdbe101
md"""
In the experiment above, the distance (x) range $(extrema(xgridUI))m and the depth (z) range $(extrema(zgridUI))m.
The seismic source is located near the surface at (x, z)=(0,0),
and $(length(rUI)) receivers on the surface (z=0).
"""

# ╔═╡ 97b3eaf7-c528-4693-8b89-bd33bdd9184c
ns = 5;

# ╔═╡ ffd19316-d2ba-4bec-abde-19cf05994ecf
md"Select Sources $(@bind sUI MultiCheckBox(1:ns, select_all=true, default=collect(1:ns)))"

# ╔═╡ 4389a492-2096-4e2a-b05f-91e9820e15a4
rlocs_x = [rx for rx in range(0, stop=last(xgridUI), length=nr)] # receiver locations

# ╔═╡ 912a031d-ebed-44b4-a350-2c72fb7623de
# source x location(s)
slocs_x = [sx for sx in range(0, stop=last(xgridUI), length=ns)] # source locations

# ╔═╡ 55631892-919f-4950-acfe-d08bf0a691fc
rlocs_z = fill(-10.0, length(rlocs_x));

# ╔═╡ 0c74076a-a3a1-4b9b-9592-0e67077d2f83
# source z location(s)
slocs_z = fill(-20.0, length(slocs_x));

# ╔═╡ 5ffccf58-7f01-49fc-a6c8-8e284a05b2dc
acq = map((; rlocs_x, rlocs_z, slocs_x, slocs_z)) do x
    Float32.(x)
end;

# ╔═╡ 3cf87a6e-d244-494b-8675-aaaaac6fbf00
acqgrid = map((; rlocs_x=rlocs_x_grid, rlocs_z=rlocs_z_grid, slocs_x, slocs_z)) do x
    Float32.(x)
end;

# ╔═╡ 9404b52d-c571-4e23-8dda-d20d10de0457
md"### 2-D Green's Function"

# ╔═╡ fc363930-8ab5-40bd-9e2c-62aea74cae57
md"""
Calculate Green's function for a 2-D homogeneous medium.
Arguments;
`sx`: Source location in the x-direction;
`sz`: Source location in the z-direction;
`rx`: Receiver location in the x-direction;
`rz`: Receiver location in the z-direction;
`f`: Frequency of the wave;
`c`: Velocity of the medium.
Returns:
Green's function for the given source-receiver pair in the 2-D homogeneous medium.
Green's function describes the response of a medium with velocity `c` and density `rho` to a point source at `(sx, sz)` and a receiver at `(rx, rz)` at a given frequency `f`.
"""

# ╔═╡ e2f1867b-a786-4a0a-a82f-b8a80c2745d5
md"Method to get distance between source and receiver"

# ╔═╡ dc40bf64-587f-452c-abb1-56fc44a740c7
rad(sx, sz, rx, rz) = sqrt(abs2(sx - rx) + abs2(sz - rz))

# ╔═╡ ab29f4e1-cc96-46a3-b0b1-9592378d5c64
function G0(rx, rz, sx, sz, k, rho)
    # k = ω/c wavenumber
    # r distance between source and receiver
    r = rad(sx, sz, rx, rz)
    return -0.25 * rho * im * hankelh2(0, k * r)
end

# ╔═╡ f78d9f96-0d4c-445b-a433-06f1a9a9eb12
md"### Modeling"

# ╔═╡ 7b46eac4-193a-4c25-99b3-b25009ba1260
function get_forward_operator_with_scatterer_locations(pa, acq, scatterer_locations, δxUI, nzUI)
    if (isempty(scatterer_locations))
        return 0.0
    else
        # select only 10 scatterers for speed
        scatterer_locations = randobs(scatterer_locations, 5)
    end
    (; xgrid, zgrid, freqgrid, kgrid, Fsource, vp0, rho0) = pa
    scat_x = map(scatterer_locations) do l
        get_xlocation(l, δxUI)
    end
    scat_z = map(scatterer_locations) do l
        get_zlocation(l, δxUI, nzUI)
    end
    nr = length(acq.rlocs_x)
    nω = length(pa.freqgrid)
    @tullio F[iω] := freqgrid[iω] * freqgrid[iω] * Fsource[iω] * 4.0 * pi * pi
    @tullio G[iω, ir, iS] := G0(scat_x[iS], scat_z[iS], acq.slocs_x[1], acq.slocs_z[1], kgrid[iω], rho0) * G0(acq.rlocs_x[ir], acq.rlocs_z[ir], scat_x[ix], scat_z[iz], kgrid[iω], rho0) * F[iω]
    # remove zero frequencies
    @tullio G[1, i, j] = complex(0.0)
    return reshape(G, nr * nω, :)
end

# ╔═╡ 8218676f-0d09-4beb-bd17-8ae5b692c56d
forward_grid = @use_memo([reset_anim]) do
    get_forward_operator_with_scatterer_locations(pagrid, acqgrid, slowness_pert_draw_input, δxUI, nzUI)
end;

# ╔═╡ 8eba66af-f61f-4453-a9b1-ee8c779ab059
forward_grid

# ╔═╡ 173adc88-5366-49e7-a1af-6478404e082e
function get_forward_operator(pa, acq, sloc_x, sloc_z)
    (; xgrid, zgrid, freqgrid, kgrid, Fsource, vp0, rho0) = pa

    X = Float32.(collect(xgrid))
    nx = length(X)
    Z = Float32.(collect(zgrid))
    nz = length(Z)
    nr = length(acq.rlocs_x)
    nω = length(pa.freqgrid)
    @tullio F[iω] := freqgrid[iω] * freqgrid[iω] * Fsource[iω] * 4.0 * pi * pi
    @tullio G[iω, ir, iz, ix] := G0(X[ix], Z[iz], sloc_x, sloc_z, kgrid[iω], rho0) * G0(acq.rlocs_x[ir], acq.rlocs_z[ir], X[ix], Z[iz], kgrid[iω], rho0) * F[iω]
    # remove zero frequencies
    @tullio G[1, i, j, k] = complex(0.0)
    return reshape(G, nr * nω, nx * nz)
end

# ╔═╡ 9ce0b88a-4e39-4c95-83fb-1758ee661b45
forward_UI = map(slocs_x, slocs_z) do sx, sz
	get_forward_operator(paUI, acq, sx, sz);
end

# ╔═╡ e5a6e125-c700-4471-9781-abbd0b1c49c7
function get_reference_wavefield(pa, acq, tgrid)
    (; xgrid, zgrid, freqgrid, kgrid, Fsource, vp0, rho0) = pa
    @tullio D[iω, ir, is] := G0(acq.rlocs_x[ir], acq.rlocs_z[ir], acq.slocs_x[is], acq.slocs_z[is], kgrid[iω], rho0) * Fsource[iω]
    # remove zero frequencies
    @tullio D[1, i, j] = complex(0.0)
    # transform to time domain
    return irfft(D, length(tgrid), 1,)
end

# ╔═╡ 887cc848-e97e-451d-8bb4-3f0eebb20723
d = get_reference_wavefield(paUI, acq, tgrid);

# ╔═╡ 9f05fce9-9c03-49b4-abeb-70977dcb9892
dgrid = reshape(get_reference_wavefield(pagrid, acqgrid, tgrid)[:, :, 1], :, nzgrid, nxgrid);

# ╔═╡ 04669a0d-9bb8-4106-ab87-7afd3ac86597
function get_scattered_wavefield(slowness, G, acq, pa)
    nr = length(acq.rlocs_x)
    nω = length(pa.freqgrid)
    s = view(slowness, :)
    d = G * s
    d = reshape(d, nω, nr)
    # transform to time domain
    return irfft(d, length(pa.tgrid), 1,)
end

# ╔═╡ 68d92a09-8fc6-485c-aeac-c392a21951e4
begin
    δd = map(forward_UI) do Gmat
		d = get_scattered_wavefield(slowness_grid_input, Gmat, acq, paUI)
		d[:, filter(x -> x ∉ rUI, 1:nr)] .= 0.0f0
		d
	end
end;

# ╔═╡ 2bf0042b-f9b3-4ea3-b3d8-46d89a1ce2fe
δdgrid = reshape(get_scattered_wavefield(1.0f-8 * ones(5), forward_grid, acqgrid, pagrid), :, nzgrid, nxgrid);

# ╔═╡ 661ccded-ddcc-4606-9788-3a0913755d67
get_scattered_wavefield(1.0f-8 * ones(5), forward_grid, acqgrid, pagrid)

# ╔═╡ 2ebe48c7-e2c2-4002-b57a-2c5b94d90165
# δs is the perturbation in slowness
# δx is the spatial sampling used to scale slowness_pertindices
# rlocs and slocs are receiver and source positions
function get_migration_image(δd, G, acq, pa)
    (; xgrid, zgrid) = pa
    nx = length(xgrid)
    nz = length(zgrid)
    δd = rfft(δd, 1,)
    δd = view(δd, :)
    δs = G' * δd
    δs = reshape(δs, nz, nx)
    # transform to time domain
    return real.(δs)
end

# ╔═╡ d8fc5b82-84a7-4f86-88f0-55b09b355a25
images = map(forward_UI, δd) do Gmat, d
	get_migration_image(d, Gmat, acq, paUI)
end;

# ╔═╡ 98d344f7-3ee2-42b7-aefe-d88604abdc8b
image = sum(images[sUI]);

# ╔═╡ 3c23a484-0b9a-4044-bcec-fec447509991
md"### Plots"

# ╔═╡ a6f59c95-2271-4f2f-894f-b574401cd545
function plot_data(d, δd, d1max)
    fig = Plot(Layout(title="Measured Wavefield (Source #$(first(sUI)))", yaxis_autorange="reversed", height=400, width=500, yaxis_title="time [s]", Subplots(shared_xaxes=true, shared_yaxes=true, rows=1, cols=2, subplot_titles=["Reference" "Scattered"], x_title="# Receiver")))
    

    add_trace!(fig, heatmap(y=tgrid, z=d, zmin=-d1max, zmax=d1max, showscale=false, colorscale="jet"), row=1, col=1)
    add_trace!(fig, heatmap(y=tgrid, z=δd, zmin=-d1max, zmax=d1max, showscale=false, colorscale="jet"), row=1, col=2)

    return PlutoPlotly.plot(fig)

end

# ╔═╡ 56d75f42-6a38-4e79-a0a9-e9f005996ac6
plot_data(d[:, :, first(sUI)], δd[:, :, 1][first(sUI)], maximum(abs,d[:,:,1][first(sUI)])/2)

# ╔═╡ 44dfa401-7f1d-467f-99e9-beab9ec52427
function add_ageom!(fig, ageom, row=1, col=1)
    if (!(ageom === nothing))
        add_trace!(fig, scatter(
                x=ageom.rlocs_x[filter(x -> x ∈ rUI, 1:nr)],
                y=ageom.rlocs_z, mode="markers",
                marker_color="black", marker_symbol="triangle-down", showlegend=false), row=row, col=col)
        add_trace!(fig, scatter(
                x=ageom.slocs_x,
                y=ageom.slocs_z, mode="markers",
                marker_color="black", marker_size=10, marker_symbol="star", showlegend=false), row=row, col=col)
    end

end

# ╔═╡ e80193d9-0d99-4cfd-92d6-20afde48bef7
function plot_image(image)
    fig = Plot(Layout(yaxis_autorange="reversed", xaxis=attr(range=extrema(xgridUI) .+ [-20, 20]),
        yaxis=attr(range=extrema(zgridUI) .+ [-10, 10]),
        height=225, width=350, yaxis_title="Depth", xaxis_title="Distance", Subplots(shared_xaxes=true, shared_yaxes=true, rows=1, cols=1, subplot_titles=["Migration Image" ""])))
    add_trace!(fig, heatmap(x=xgridUI, y=zgridUI, z=image, showscale=false, colorscale="seismic"), row=1, col=1)
    add_ageom!(fig, acq)
    return PlutoPlotly.plot(fig)

end

# ╔═╡ 3786bd54-105a-46d7-8232-b93835b01129
plot_image(image)

# ╔═╡ 1c95eb36-adc8-44f1-867b-eb3521684eb9
function plot_animations(d1, d2, d1max)

    fig = Plot(Layout(title="Wavefield (Source #$(first(sUI)))", yaxis_autorange="reversed", xaxis=attr(range=extrema(xgridUI) .+ [-20, 20]),
        yaxis=attr(range=extrema(zgridUI) .+ [-10, 10]),
        height=225, width=700, yaxis_title="Depth", xaxis_title="Distance", Subplots(shared_xaxes=true, shared_yaxes=true, rows=1, cols=2, subplot_titles=["Reference" "Scattered"])))
    add_trace!(fig, heatmap(x=xgrid, y=zgrid, z=d1, zmin=-d1max, zmax=d1max, showscale=false, colorscale="Greys"), row=1, col=1)
    add_trace!(fig, heatmap(x=xgrid, y=zgrid, z=d2, zmin=-d1max, zmax=d1max, showscale=false, colorscale="Greys"), row=1, col=2)
    add_ageom!(fig, acq, 1, 1)
    add_ageom!(fig, acq, 1, 2)
    return PlutoPlotly.plot(fig)

end

# ╔═╡ cff433fb-e1c4-42ae-a310-db85f48d09e3
begin
	reset_anim
	plot_animations(dgrid[mod(it_plot, div(length(tgrid), 2))+1, :, :], δdgrid[mod(it_plot, div(length(tgrid), 2))+1, :, :], maximum(abs, dgrid) / 1.0)
end

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

# ╔═╡ Cell order:
# ╠═c0e6e258-7fcb-4bd8-9931-be0203c1adc9
# ╠═4aa9e374-27a1-4d80-9d7f-9a7c1ee859b2
# ╟─5d3692af-dee6-4adf-9276-82f11a1a9544
# ╟─99d5ad02-2897-4f2c-9a80-a146f90e38dd
# ╟─3d88f8c4-ba31-4d99-8966-3ddf617e5b5f
# ╟─3786bd54-105a-46d7-8232-b93835b01129
# ╟─ffd19316-d2ba-4bec-abde-19cf05994ecf
# ╟─6e74479f-96ce-4af7-ae7d-603fe32ba88f
# ╟─56d75f42-6a38-4e79-a0a9-e9f005996ac6
# ╟─4972f23b-88a8-49d9-acad-75a65bdbe101
# ╟─9d498284-f702-4f6b-aace-54743c7df206
# ╟─cff433fb-e1c4-42ae-a310-db85f48d09e3
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
# ╟─db0b5eb4-406d-4d8c-b0ac-544c1d6f37b1
# ╠═887cc848-e97e-451d-8bb4-3f0eebb20723
# ╠═9f05fce9-9c03-49b4-abeb-70977dcb9892
# ╠═9ce0b88a-4e39-4c95-83fb-1758ee661b45
# ╠═8218676f-0d09-4beb-bd17-8ae5b692c56d
# ╠═acba3f10-975b-423d-9ab8-ead6d9f9774b
# ╠═68d92a09-8fc6-485c-aeac-c392a21951e4
# ╠═2bf0042b-f9b3-4ea3-b3d8-46d89a1ce2fe
# ╠═8eba66af-f61f-4453-a9b1-ee8c779ab059
# ╠═661ccded-ddcc-4606-9788-3a0913755d67
# ╠═d8fc5b82-84a7-4f86-88f0-55b09b355a25
# ╠═98d344f7-3ee2-42b7-aefe-d88604abdc8b
# ╟─b48cf5ec-a033-46e3-aec8-da5864b2386a
# ╠═8acbffaf-1811-4592-a2d1-a8f561242d85
# ╠═3affdcaf-8b32-496e-81a3-a898578948d2
# ╠═094b0ad6-1ec0-407a-884f-7147d8a2ec9f
# ╠═135efe98-672e-4013-8e4c-15f0a29e9bab
# ╟─953e46b3-b507-4ff6-a24b-763250917da3
# ╠═d342ea67-3f1f-421d-8850-09c72b0398a4
# ╠═968fb80c-0f44-4198-9aa0-87a07f28c6ca
# ╠═8962a7d6-1e29-4d99-b8d3-aecd9d69de76
# ╠═504f09a0-6d86-459a-9b91-f1c8dd684427
# ╠═3ebd625f-49dc-4f5e-bcb0-2633092f5062
# ╠═2d82c01f-2be9-4628-b926-bab17d451ed8
# ╠═227a30d4-9cec-440a-8238-87cd0d560e33
# ╠═33a221fe-ab5d-4a41-80bd-105746f949fb
# ╠═07a5f979-f83d-45b2-9e19-0c770577f09e
# ╟─855bbfe1-39bd-46b0-bafd-be63cc5f67b2
# ╠═0b4b74d4-2fd1-4a52-89c0-4141d047b55a
# ╠═d856645d-39bc-4603-837f-3c39d9df64af
# ╠═2aa9168a-cadb-4060-9eac-d93a1cf3bf0b
# ╠═453fcb4f-1cc6-4c57-9bec-2e1d7c76244a
# ╠═2ac0ea94-aa5d-43da-a240-db4a57dda204
# ╠═97b3eaf7-c528-4693-8b89-bd33bdd9184c
# ╠═4389a492-2096-4e2a-b05f-91e9820e15a4
# ╠═912a031d-ebed-44b4-a350-2c72fb7623de
# ╠═55631892-919f-4950-acfe-d08bf0a691fc
# ╠═0c74076a-a3a1-4b9b-9592-0e67077d2f83
# ╠═5ffccf58-7f01-49fc-a6c8-8e284a05b2dc
# ╠═3cf87a6e-d244-494b-8675-aaaaac6fbf00
# ╟─9404b52d-c571-4e23-8dda-d20d10de0457
# ╟─fc363930-8ab5-40bd-9e2c-62aea74cae57
# ╠═ab29f4e1-cc96-46a3-b0b1-9592378d5c64
# ╟─e2f1867b-a786-4a0a-a82f-b8a80c2745d5
# ╠═dc40bf64-587f-452c-abb1-56fc44a740c7
# ╟─f78d9f96-0d4c-445b-a433-06f1a9a9eb12
# ╠═7b46eac4-193a-4c25-99b3-b25009ba1260
# ╠═173adc88-5366-49e7-a1af-6478404e082e
# ╠═e5a6e125-c700-4471-9781-abbd0b1c49c7
# ╠═04669a0d-9bb8-4106-ab87-7afd3ac86597
# ╠═2ebe48c7-e2c2-4002-b57a-2c5b94d90165
# ╟─3c23a484-0b9a-4044-bcec-fec447509991
# ╠═a6f59c95-2271-4f2f-894f-b574401cd545
# ╠═e80193d9-0d99-4cfd-92d6-20afde48bef7
# ╠═1c95eb36-adc8-44f1-867b-eb3521684eb9
# ╠═44dfa401-7f1d-467f-99e9-beab9ec52427
# ╟─98521f53-fbf4-4ba2-b0b9-629f21478a31
# ╟─8644be64-5901-47df-8840-b234e4ce4a01
# ╟─f3b78d8f-9226-47e8-beca-38956260a40b
