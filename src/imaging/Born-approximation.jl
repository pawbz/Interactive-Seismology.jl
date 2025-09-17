### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
Einsum = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
MLUtils = "f1d291b0-491e-4a28-83b9-f70985020b54"
PlutoHooks = "0ff47ea0-7a50-410d-8455-4348d5de0774"
PlutoLinks = "0ff47ea0-7a50-410d-8455-4348d5de0420"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlutoUIExtra = "a011ac08-54e6-4ec3-ad1c-4165f16ac4ce"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
Tullio = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"

[compat]
Bessels = "~0.2.8"
Einsum = "~0.4.1"
FFTW = "~1.7.2"
Latexify = "~0.16.1"
MLUtils = "~0.4.3"
PlutoHooks = "~0.0.5"
PlutoLinks = "~0.1.6"
PlutoPlotly = "~0.4.4"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.54"
PlutoUIExtra = "~0.1.8"
SymbolicUtils = "~1.5.0"
Symbolics = "~5.13.0"
Tullio = "~0.3.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "29271830d51ccb833d6e2ba70c5d0d1c03800bb8"

[[deps.ADTypes]]
git-tree-sha1 = "fb97701c117c8162e84dfcf80215caa904aef44f"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.13.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "d7832de8cf7af26abac741f10372080ac6cb73df"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.34.7"

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
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "0ba8f4c1f06707985ffb4804fdad1bf97b233897"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.41"

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
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50c3c56a52972d78e8be9fd135bfb91c9574c140"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgCheck]]
git-tree-sha1 = "680b3b8759bd4c54052ada14e52355ab69e07876"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "017fcb757f8e921fb44ee063a7aafe5f89b86dd1"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.18.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
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
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

[[deps.BangBang]]
deps = ["Accessors", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires"]
git-tree-sha1 = "e2144b631226d9eeab2d746ca8880b7ccff504ae"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.4.3"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTablesExt = "Tables"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "cb25e4b105cc927052c2314f8291854ea59bf70a"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.4"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

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

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

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
deps = ["LinearAlgebra"]
git-tree-sha1 = "d8a9c0b6ac2d9081bf76324b39c78ca3ce4f0c98"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.6"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataPipes]]
git-tree-sha1 = "29077a8d5c093f4e0988e92c0d76f56c4c581900"
uuid = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
version = "0.3.18"

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
version = "1.11.0"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

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

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "03aa5d44647eaec98e1920635cdfed5d5560a8b9"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.117"

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
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "51b4b84d33ec5e0955b55ff4b748b99ce2c3faa9"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.7"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "30a1848c4f4fc35d1d4bbbd125650f6a11b5bc6c"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.7"

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
git-tree-sha1 = "ec22cbbcd01cba8f41eecd7d44aac1f23ee985e3"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.2"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+3"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "0a2e5873e9a5f54abb06418d57a8df689336a660"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.2"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

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

[[deps.FlexiMaps]]
deps = ["Accessors", "DataPipes", "InverseFunctions"]
git-tree-sha1 = "88fb6ab75454c21be1d75a0a430a0ed95f0d3f1e"
uuid = "6394faf6-06db-4fa8-b750-35ccc60383f7"
version = "0.1.28"

    [deps.FlexiMaps.extensions]
    AxisKeysExt = "AxisKeys"
    DictionariesExt = "Dictionaries"
    IntervalSetsExt = "IntervalSets"
    StructArraysExt = "StructArrays"
    UnitfulExt = "Unitful"

    [deps.FlexiMaps.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dictionaries = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

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
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "ExprTools", "Logging", "MultivariatePolynomials", "PrecompileTools", "PrettyTables", "Primes", "Printf", "Random", "SIMD", "TimerOutputs"]
git-tree-sha1 = "6b505ef15e55bdc5bb3ddbcfebdff1c9e67081e8"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.5.1"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6df9cd6ee79fc59feab33f63a1b3c9e95e2461d5"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.2"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "2bd56245074fab4015b9174f24ceba8293209053"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.27"

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

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

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
version = "1.11.0"

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
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "4bf4b400a8234cff0f177da4a160a90296159ce9"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.41"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "MacroTools", "PrecompileTools", "Requires", "StaticArrays", "UUIDs"]
git-tree-sha1 = "80d268b2f4e396edc5ea004d1e0f569231c71e9e"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.34"

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"
    LinearAlgebraExt = "LinearAlgebra"
    SparseArraysExt = "SparseArrays"

    [deps.KernelAbstractions.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "e459fda6b68ea8684b3fcd513d2fd1e5130c4402"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.16.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd714447457c660382fe634710fb56eb255ee42e"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.6"

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
version = "1.11.0"

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

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "688d6d9e098109051ae33d126fcfc88c4ce4a021"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.1.0"

[[deps.MIMEs]]
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MLCore]]
deps = ["DataAPI", "SimpleTraits", "Tables"]
git-tree-sha1 = "73907695f35bc7ffd9f11f6c4f2ee8c1302084be"
uuid = "c2834f40-e789-41da-a90e-33b280584a8c"
version = "1.0.0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MLUtils]]
deps = ["ChainRulesCore", "Compat", "DataAPI", "DelimitedFiles", "FLoops", "MLCore", "NNlib", "Random", "ShowCases", "SimpleTraits", "Statistics", "StatsBase", "Tables", "Transducers"]
git-tree-sha1 = "6963295133aaa789f5fb18a6dd276c420793cf43"
uuid = "f1d291b0-491e-4a28-83b9-f70985020b54"
version = "0.4.7"

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

[[deps.MicroCollections]]
deps = ["Accessors", "BangBang", "InitialValues"]
git-tree-sha1 = "44d32db644e84c75dab479f1bc15ee76a1a3618f"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.2.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "9c0bc309df575c85422232eedfb74d5a9c155401"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.3"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Pkg", "Random", "Requires", "Statistics"]
git-tree-sha1 = "78de319bce99d1d8c1d4fe5401f7cfc2627df396"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.9.18"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"
    NNlibCUDACUDNNExt = ["CUDA", "cuDNN"]
    NNlibCUDAExt = "CUDA"
    NNlibEnzymeCoreExt = "EnzymeCore"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "cc0a5deefdb12ab3a096f00a6d42133af4560d71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.2"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

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
version = "0.8.5+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "966b85253e959ea89c53a9abebbf2e964fbf593b"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.32"

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
deps = ["AbstractPlutoDingetjes", "BaseDirs", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "Reexport", "TOML"]
git-tree-sha1 = "1ae939782a5ce9a004484eab5416411c7190d3ce"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.4.6"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.PlutoUIExtra]]
deps = ["AbstractPlutoDingetjes", "ConstructionBase", "FlexiMaps", "HypertextLiteral", "InteractiveUtils", "IntervalSets", "Markdown", "PlutoUI", "Random", "Reexport"]
git-tree-sha1 = "b4ff5d24e2dc8fbf319cd44f9f81b5356e27bafb"
uuid = "a011ac08-54e6-4ec3-ad1c-4165f16ac4ce"
version = "0.1.8"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "8765738bc5a6f1554cb61c5ddfae5bf279e8b110"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.25"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

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

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

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
git-tree-sha1 = "fe9d37a17ab4d41a98951332ee8067f8dca8c4c2"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.29.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
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

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "9bb80533cb9769933954ea4ffbecb3025a783198"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.7.2"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

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

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "Expronicon", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "ffed2507209da5b42c6881944ef41a340ab5449b"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.74.1"

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
git-tree-sha1 = "6149620767866d4b0f0f7028639b6e661b6a1e44"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.12"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "0444a37a25fab98adbd90baa806ee492a3af133a"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.6.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShowCases]]
git-tree-sha1 = "7f534ad62ab2bd48591bdeac81994ea8c445e4a5"
uuid = "605ecd9f-84a6-4c9e-81e2-4798472b76a3"
version = "0.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "e3be13f448a43610f978d29b7adf78c76022467a"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.12"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

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
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "fd2d4f0499f6bb4a0d9f5030f5c7d61eed385e03"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.37"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "669e43e90df46fcee4aa859b587da7a7948272ac"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.5.1"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils"]
git-tree-sha1 = "80e67b76699b1ec130ac60946b3d06017f58141f"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.13.0"

    [deps.Symbolics.extensions]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "d7298ebdfa1654583468a487e8e83fae9d72dac3"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.26"

[[deps.Transducers]]
deps = ["Accessors", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "SplittablesBase", "Tables"]
git-tree-sha1 = "7deeab4ff96b85c5f72c824cae53a1398da3d1cb"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.84"

    [deps.Transducers.extensions]
    TransducersAdaptExt = "Adapt"
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.Tullio]]
deps = ["DiffRules", "LinearAlgebra", "Requires"]
git-tree-sha1 = "972698b132b9df8791ae74aa547268e977b55f68"
uuid = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"
version = "0.3.8"

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
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
