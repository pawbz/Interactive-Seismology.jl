### A Pluto.jl notebook ###
# v0.19.25

#> [frontmatter]
#> title = "Seismic Full Waveform Inversion"
#> description = "A notebook that demonstrates how interactivity makes seismic full waveform inversion more accessible and intuitive."

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

# ╔═╡ c5815f5e-9164-11ec-10e1-691834761dff
begin
  using FFTW
  using LinearAlgebra
  using LaTeXStrings
  using PlutoUI
  using Statistics
  using ProgressLogging
  using SparseArrays
  using DSP
  using PlutoPlotly
  using LossFunctions
  using MLUtils
  using FiniteDifferences
  using ForwardDiff
  using PlutoTeachingTools
end

# ╔═╡ 3c540889-49dc-415c-acbc-3494897b260c
TableOfContents()

# ╔═╡ 9c32f5bc-f6d1-4048-903a-27224aaa1f40
md"""
# Full Waveform Inversion

Seismic full waveform inversion (FWI) is a powerful technique for imaging the Earth's subsurface using seismic waves. However, FWI can be computationally intensive and require a deep understanding of numerical methods and algorithms. This notebook demonstrates how interactivity makes FWI more accessible and intuitive.
Pseudo-spectral methods are used to solve the 2-D SH stress-velocity system.
The particle velocity and stress fields are staggered in time and we use a simple explicit time-stepping approach.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 9ab80963-a575-4a1a-a529-f0f2e64af4cc


# ╔═╡ 881c7368-57df-469a-96ec-821512cf98e0
md"""## Earth Medium
Consider a mantle-like medium that has density of $3.22$ $g/cc$ and shear wave velocity of $5$ $km$ $s^{-1}$. Suppose it has a reflector such that the medium has a higher density below the reflector compared to above. You can specify the location and density of the reflector below.
"""

# ╔═╡ fa78af13-e3c6-4d6f-8a3e-1187fe9ae159
md"""## A Priori Information
$(@bind apriori_flag Select([:a=>"Accurate Background Velocity; Unknown Reflector", :b=>"Inaccurate Background Velocity; Known Reflector", :c=>"Inaccurate Background Velocity; Unknown Reflector"]))
"""

# ╔═╡ 37157e04-273a-4be2-9687-2016480a0421
apriori_flag

# ╔═╡ e233afec-6049-4277-8be9-95687c4589b5
md"""
## Measured Displacement Field
Pseudo spectral method is applied for solving the seismic equation given the true medium parameters. Wavefields at the receivers are recorded and observed data is generated.
"""


# ╔═╡ b572f855-db01-4c4f-8922-968bc0ef5fdf
md"## Gradients"

# ╔═╡ b3015aff-bdab-4c94-87da-8744c174263a
md"## Wave Equation"

# ╔═╡ 4adbb7f0-7927-470d-8c00-07d3b3c0cd78
md"""
We shall first consider the equation of motion and simplify it to derive the 2-D SH system

```math
\rho{\ddot{{u}}}= \mathbf{\nabla} \cdot \mathbf{\sigma} + {f}. \

```
Lets assume that the \
1) particles are constrained to move only along y direction, i.e., $u_x=u_z=0$
```math
\rho{\ddot{{u_{y}}}}= \partial_x{\sigma_{yx}} +\partial_z{\sigma_{yz}} + {f_y},
```\
2) stresses don't change along $y$ direction to result in
```math
\partial_{y}{\sigma}_{(\cdot)}=0.
```
The constitutive equations are:
```math
\sigma_{yx}=2\mu\epsilon_{yx}=\mu\partial_{x}u_x,

```
```math
\sigma_{yz}=2\mu\epsilon_{yz}=\mu\partial_{z}u_y.

```
"""


# ╔═╡ 7fd6fe40-0065-4074-b8cf-3d8cd8e68319
md"""
## Forward Experiment
"""

# ╔═╡ 122c2aa4-baec-4288-ab53-afa1d977c486
md"""### Spatial Grids
Let us set up the spatial parameters of the simulation. We define a 2-D spatial grid of extents $x∈[-40,40]$ km and $z∈[0,40]$ km.
"""

# ╔═╡ 7913b296-3010-410d-942f-834a47d599c7
DSP.nextfastfft(106)

# ╔═╡ 77fa76f3-ffda-4d95-8c12-7ccce6a7e52e
# input number of wavelengths, then roughly get the length, then use nextfast fft

# ╔═╡ 7928a88c-f217-4338-a6a5-50ab2d422480
begin
  # domain extends
  zgrid = range(0, stop=40, length=107)
  xgrid = range(-40, stop=40, length=221)

  # Numerics
  nx, nz = length(xgrid), length(zgrid)

  # Derived numerics
  dx, dz = step(xgrid), step(zgrid) # cell sizes
end;

# ╔═╡ cb83dbd1-c423-4b27-b29e-7dc8051f43d5
md"""
### Acquistion Geometry
Choose the number of sources and receivers to be used in the simulation.
"""

# ╔═╡ 2ec95e6e-7c2c-41e5-87b7-84583564f079
md"""### Medium
This will also serve as an initial model during inversion. Medium without the reflector is considered.
"""

# ╔═╡ 5b271e5f-879c-4c43-825a-9660f322febd
md"""### Time Grid
In order to study the propagation of the wave-front, we need to choose a time step that satisfies the Courant condition. Here, we have chosen the Courant number to be $0.2$.
"""

# ╔═╡ 2855c8cf-8364-4c6c-a122-781b99440e89
md"### Body Forces"

# ╔═╡ 8b5372ec-0742-48ea-81c4-d303b96f56c7
md"## Inversion"

# ╔═╡ 14c0cd30-a52c-4f93-9f7b-cae6328c0655
md"### Generate Data"

# ╔═╡ 439df138-5198-4737-915d-7bd2c157aa4b
md"### Bundle Parameters "

# ╔═╡ 0454ce0a-d6de-427f-bbdc-3bcec21327f2
md"### Loss Function"

# ╔═╡ f2fb92bb-33d6-4e15-8caf-b245e000ad69
loss = L2DistLoss()

# ╔═╡ 4a4f1300-94ba-4c1b-9d10-9e955d194ff6
function Jdata(dobs, d, loss=loss)
  return mapreduce(+, dobs, d) do d1, d2
    value(loss, d1, d2)
  end
end

# ╔═╡ a3a9deea-e2d6-4d58-90d7-5a54be176289
md"### Adjoint Simulation"

# ╔═╡ 1dc2ca10-5ba3-4efa-b6d9-d203cf91598b
md"Deviation between the observed and modelled data is referred to as the data error and it acts as the forcing in the case of adjoint simulation."

# ╔═╡ ddb37082-cae0-4a68-ab55-19563d8727ed
function get_adj_source(dobs, d)
  adj_source = deriv.(loss, dobs, d)
  reverse!(adj_source, dims=1)
end

# ╔═╡ 66f9c698-61e3-4b61-aff3-dfc67eb2f6af
md"""
### Gradients
"""

# ╔═╡ 8d161f09-8339-4277-8739-ff76607f7abf
md"""
## Finite-Difference Tests
Tick to perform these tests: $(@bind do_fd_tests CheckBox())
"""

# ╔═╡ 006739fb-24a1-49b0-9619-fe8e2d3c8fca
do_fd_tests && (xs = zeros(Float32, 3))

# ╔═╡ fc49a6d7-a1b1-458a-a9ad-e120282bbabc
md"""
## Appendix
"""

# ╔═╡ a62839d5-837b-4c37-996f-33659c34911c
md"### UI"

# ╔═╡ 3fc0e673-2fa3-489f-a56e-a867ea37cbce
md"""
Function to choose the number of sources and receivers
"""

# ╔═╡ 2ea24e92-d66e-4c60-ad0b-f671d894fef2
function src_rec_ip()
  return PlutoUI.combine() do Child
    src = [md"""Number of sources = $(Child("ns", Slider(range(start=1,stop=10,step=1), default=1, show_value=true)))
    """,]

    rec = [md"""Number of receivers = $(Child("nr", Slider(range(start=5, stop=15, step=1), default=50, show_value=true)))
    """,]

    md"""
    $(src)
    $(rec)
    """
  end
end;

# ╔═╡ 86ed93a1-c7b9-4b3a-8cb5-ca4405cff3df
@bind acq src_rec_ip()

# ╔═╡ 7f797571-055e-4975-9c26-fc968bbc0094
md"""
Function to choose the parameters of the true medium. Z-location of the reflector as the well as the density of the medium below the reflector can be chosen.
"""

# ╔═╡ d7b37c59-e0b3-4e47-86d3-7f1df7400f09
function choose_param_truemed()
  return PlutoUI.combine() do Child
    zloc = [md"""Z location (km) = $(Child("z", Slider(zgrid[floor(Int,0.3*nz):end], default=zgrid[floor(Int,0.5*nz)], show_value=true)))
    """,]

    dens = [md"""Density (g/cc) = $(Child("ρ", Slider(range(start=4, stop=6, step=0.1), default=5, show_value=true)))
    """,]

    md"""
    $(zloc)
    $(dens)
    """
  end
end;

# ╔═╡ 65efba3b-16b0-4113-a59e-809365e7bdd6
@bind param_true_med choose_param_truemed()

# ╔═╡ 6be2f4c2-e9ed-43c2-b66c-ef3176bb9000
md"""
### Fourier Derivatives
We now define methods to compute a 2-D Fourier transform.
```math
\hat{u}(\mathbf{k}, t) = \frac{1}{2\pi}\sum_{\mathbf{x}} u(\mathbf{x}, t) \, e^{-i \mathbf{k}\cdot \mathbf{x}} \ ,
```
*The derivative property of Fourier Transform*: First order differentiation of a function along the $x$ dimension is equivalent to multiplying its Fourier Transform by $\imath k_x$ in the wavenumber domain.\

```math
\partial_xu(\mathbf{x}, t)\leftrightarrow{ik_x}\hat{u}(\mathbf{k}, t).

```
In the following cell, the derivatives are calculated using the functions `Dx!` and `Dz!`. Simply put, these functions compute the Fourier transform, apply the derivative property, and take the inverse Fourier Transform to generate the spatial derivative.
"""

# ╔═╡ aa19e992-2735-4324-8fd7-15eacadf0faa
begin
  Fz = plan_rfft(zeros(Float32, nz, nx), (1))
  Fx = plan_rfft(zeros(Float32, nz, nx), (2))
  kx = reshape(collect(rfftfreq(nx, inv(step(xgrid)))), 1, :) * 2 * pi
  kz = reshape(collect(rfftfreq(nz, inv(step(zgrid)))), :, 1) * 2 * pi
  storagex = zero(Fx * zeros(Float32, nz, nx))
  storagez = zero(Fz * zeros(Float32, nz, nx))
  fp = (; Fx, Fz, kx, kz, storagex, storagez)
  function Dx!(dPdx, P, fp)
    mul!(fp.storagex, fp.Fx, P)
    broadcast!(*, fp.storagex, fp.storagex, fp.kx)
    rmul!(fp.storagex, im)
    ldiv!(dPdx, fp.Fx, fp.storagex)
  end
  Dx!(dP, P) = Dx!(dP, P, fp)
  Dx(P) = (dPdx = zero(P); Dx!(dPdx, P, fp); dPdx)
  function Dz!(dPdz, P, fp)
    mul!(fp.storagez, fp.Fz, P)
    broadcast!(*, fp.storagez, fp.storagez, fp.kz)
    rmul!(fp.storagez, im)
    ldiv!(dPdz, fp.Fz, fp.storagez)
  end
  Dz!(dP, P) = Dz!(dP, P, fp)
  Dz(P) = (dPdz = zero(P); Dz!(dPdz, P, fp); dPdz)
  nothing
end

# ╔═╡ e2127d9b-f2a4-4970-a36e-5fa70c304ca7
# test transpose of Dx
begin
  x1 = rand(Float32, nz, nx)
  y1 = rand(Float32, nz, nx)
  dot(Dx(x1), y1), dot(x1, -Dx(y1))
end

# ╔═╡ e8333b23-53c3-445e-9ca3-6b278359f8ab
md"### Acquisition"

# ╔═╡ 9248af7f-dc1a-4bf6-8f3f-304db73fc604
function get_ageom(xgrid, zgrid, ns, nr; zs=quantile(zgrid, 0.25), zr=quantile(zgrid, 0.25))
  A = (; ns, nr, zs=fill(zs, ns), zr=fill(zr, nr),
    xr=(nr == 1) ? [quantile(xgrid, 0.25)] : range(quantile(xgrid, 0.25), stop=quantile(xgrid, 0.75), length=nr),
    xs=(ns == 1) ? [quantile(xgrid, 0.25)] : range(quantile(xgrid, 0.25), stop=quantile(xgrid, 0.75), length=ns)
  )
  return A
end;

# ╔═╡ d39753e2-5986-4394-9293-9e394f2807f0
ageom = get_ageom(xgrid, zgrid, acq.ns, acq.nr);

# ╔═╡ 00a637a6-ddc4-4830-be65-1891d3cb18bc
get_adj_ageom(ageom) = (; nr=ageom.ns, ns=ageom.nr, xs=ageom.xr, zs=ageom.zr, xr=ageom.xs, zr=ageom.zs)

# ╔═╡ e7e72f61-d79e-4a82-ad9c-129191c8a8c2
adj_ageom = get_adj_ageom(ageom);

# ╔═╡ e08bf013-00c7-4870-82d8-19b899e7208d
# m are the medium properties that will be used 
function get_projection_matrix(xpos, zpos, xgrid, zgrid; transpose_flag=false, m=ones(Float32, length(zgrid), length(xgrid)))
  l = LinearIndices((length(zgrid), length(xgrid)))
  @assert length(xpos) == length(zpos)
  n = length(xpos)
  N = length(xgrid) * length(zgrid)
  I = broadcast(zpos, xpos) do z, x
    iz = argmin(abs.(zgrid .- z))[1]
    ix = argmin(abs.(xgrid .- x))[1]
    return l[iz, ix]
  end
  J = collect(1:n)
  V = m[I]
  return transpose_flag ? sparse(J, I, V, n, N) : sparse(I, J, V, N, n)
end


# ╔═╡ ab8b1a22-ca7a-409e-832e-8d5d08a29a1e
md"### Data"

# ╔═╡ f4d91971-f806-4c5c-8548-b58a20acfb2c
function initialize_data(grid_param, ageom)
  zeros(Float32, length(grid_param.tgrid), length(ageom.xr))
end

# ╔═╡ 9bc38d55-285b-4b83-98d9-d7f9e03405d1
md"### Medium"

# ╔═╡ 27844886-0b54-4b08-a592-a1a38e4b0be2
function bundle_medium(μ, ρ)
  return (; μ=μ, ρ=ρ, invρ=inv.(ρ))
end

# ╔═╡ 489dcf10-b7f2-4544-b80d-3588ff00ff4a
function update_xsρ!(x, xs)
  xρ, xinvμ = chunk(x, 2)
  N = length(xρ)
  N2 = div(N, 2)
  X = view(xρ, N2:N2+length(xs)-1)
  copyto!(X, xs)
  return xs
end

# ╔═╡ 02844bc9-9faf-413f-82d2-28d6ba01668d
function update_xsinvμ!(x, xs)
  xρ, xinvμ = chunk(x, 2)
  N = length(xρ)
  N2 = div(N, 2)
  X = view(xinvμ, N2:N2+length(xs)-1)
  copyto!(X, xs)
  return xs
end

# ╔═╡ 199908d6-f5fd-4461-be53-39919a103835
function get_xs(x, xs)
  xρ, xinvμ = chunk(x, 2)
  N = length(xρ)
  N2 = div(N, 2)
  return vcat(xρ[N2:N2+length(xs)-1], xinvμ[N2:N2+length(xs)-1])
end

# ╔═╡ ad21da29-f6ff-4a94-be87-4e88640cddbf
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	# e.g., get_vs(mean, medium)
	get_vs(op, medium) = op(sqrt.(medium.μ ./ medium.ρ))
	get_vs(medium) = sqrt.(medium.μ ./ medium.ρ)
end
  ╠═╡ =#

# ╔═╡ 2abfb6d5-a201-405e-99c3-42bf0ad5aad5
χ(m, m0) = log(m / m0)

# ╔═╡ a26a5bd6-f4b8-410a-b437-06587011093e
invχ(x, m0) = exp(x) * m0

# ╔═╡ 78c09864-1fb1-427e-b77d-db82945e850c
ρ0 = Float32(3.22 * 10^-3 * 10^15) # density in kg/km3

# ╔═╡ daef685e-6c07-419f-9c68-a1c2c0ec8973
ρ1 = Float32(5.00 * 10^-3 * 10^15) # density in kg/km3

# ╔═╡ 5648367f-8628-4859-a2a5-ab92bab54e9b
μ0 = Float32(82 * 10^9 * 10^3)

# ╔═╡ c7ef6f8b-0c01-4788-982f-17ed6d0098ff
μ1 = Float32(50 * 10^9 * 10^3)

# ╔═╡ dc95afc5-0232-4c2a-801e-dccac012d2e6
invρ0 = inv(ρ0)

# ╔═╡ 13ffff0c-23cd-4a98-8ce8-31644b9429d3
invμ0 = inv(μ0)

# ╔═╡ 63a177c1-034e-4aa6-9951-367570c49850
medium_ref_values = (; μ0, invμ0, ρ0, invρ0, ρ1, μ1)

# ╔═╡ 0b890ec0-1886-4494-b4ac-46de7639f358
medium_ref_values

# ╔═╡ 1c67cda8-7712-4d5b-a2aa-af47f290f745
begin
  # lets add density perturbation (olivine)
  μtrue = ones(Float32, nz, nx) .* medium_ref_values.μ0
  ρtrue = ones(Float32, nz, nx) .* medium_ref_values.ρ0 
  reflect_index = findall(zgrid .== param_true_med.z)[1]
  ρtrue[reflect_index:end, :] .= medium_ref_values.ρ1 
  medium_true = bundle_medium(μtrue, ρtrue)
end;

# ╔═╡ 6a8139c4-12c1-4d18-bd1e-14334290aec1
#=╠═╡
begin
	courant_number = 0.2

	# lets calculate the min distance from the center to the edge of the domain
	r = min(xgrid[end] - xgrid[1], zgrid[end] - zgrid[1]) * 0.5

	# choose time stepping dt to satisfy Courant condition
    dt = courant_number * step(xgrid) * inv(get_vs(minimum, medium_true))
    nt = Int(floor(r / (mean(sqrt.(medium_true.μ ./ medium_true.ρ)) * dt))) * 2
    tgrid = range(0, length=nt, step=dt)
    nothing
end;
  ╠═╡ =#

# ╔═╡ 12a0f7d8-64e0-411f-85ad-14d48bd9492d
#=╠═╡
md"""
Reduce the time window 
$(@bind t_grad RangeSlider(tgrid; left=first(tgrid), right=last(tgrid), show_value=true))
"""
  ╠═╡ =#

# ╔═╡ 92d05447-6d0e-4ee0-a330-244b9c65c871
#=╠═╡
@bind t_forw Slider(1:length(tgrid)-1, show_value=true)
  ╠═╡ =#

# ╔═╡ 6937b103-9ce2-4189-8129-aae1e7936d4f
#=╠═╡
length(tgrid)
  ╠═╡ =#

# ╔═╡ 8b3776bd-509b-4232-9737-36c9ae003350
begin
  # Unperturbed medium
	μref =  (apriori_flag ≠ :a) ? fill(medium_ref_values.μ1, nz, nx)  : fill(medium_ref_values.μ0, nz, nx)   # shear modulus in kg / km/s^2
  ρref = ones(Float32, nz, nx) * medium_ref_values.ρ0 # density in kg/km3
	
		if(apriori_flag == :b) 
		ρref[reflect_index:end, :] .= medium_ref_values.ρ1
		end
		
  medium_ref = bundle_medium(μref, ρref)
end;

# ╔═╡ 8298ae48-0ddc-49a9-a43f-8434e4cc3758
function get_x(medium, ref=medium_ref_values)
  return vcat(vec(χ.(medium.ρ, ref.ρ0)), vec(χ.(inv.(medium.μ), ref.invμ0)))
end

# ╔═╡ 22ac08a7-f4d7-4809-8ee3-903d96c96cd6
x0 = get_x(medium_ref)

# ╔═╡ ff1ba1f6-f127-4c47-8198-aeff5f051ad9
function update_medium!(medium, x, ref=medium_ref_values)
  xρ, xinvμ = chunk(x, 2)
  map!(medium.μ, xinvμ) do xm
    inv(invχ(xm, ref.invμ0))
  end
  map!(medium.ρ, xρ) do xm
    invχ(xm, ref.ρ0)
  end
  map!(medium.invρ, medium.ρ) do ρ
    inv(ρ)
  end
  return medium
end

# ╔═╡ cb6adea2-9ea2-4857-b969-540a3439e700
function update_x!(x, medium, ref=medium_ref_values)
  xρ, xinvμ = chunk(x, 2)
  map!(xinvμ, medium.μ) do m
    χ(inv(m), ref.invμ0)
  end
  map!(xρ, medium.ρ) do m
    χ(m, ref.ρ0)
  end
  return x
end

# ╔═╡ ae8012be-e7ab-4e85-a27f-febf08b3380b
md"### Absorbing Boundaries"

# ╔═╡ ce89710a-48d8-46d4-83b4-163a7eb0a2e5
function get_taper_array(nx, nz; np=50, tapfact=0.20)
  tarray = ones(Float32, nz, nx)
  for ix in 1:np
    tarray[:, ix] .= tarray[:, ix] * abs2(exp(-tapfact * abs(np - ix + 1) / np))
    tarray[:, nx-ix+1] .= tarray[:, nx-ix+1] * abs2(exp(-tapfact * abs(np - ix + 1) / np))
  end
  for iz in 1:np
    tarray[iz, :] .= tarray[iz, :] * abs2(exp(-tapfact * abs(np - iz + 1) / np))
    tarray[nz-iz+1, :] .= tarray[nz-iz+1, :] * abs2(exp(-tapfact * abs(np - iz + 1) / np))
  end
  return tarray
end

# ╔═╡ b7f4078a-ead0-4d42-8b44-4f471eefc6fc
function clip_edges(m, grid_param)
  (; xgrid, zgrid) = grid_param
  I = findall(x -> isequal(x, 1), grid_param.tarray)
  xs = extrema(unique(getindex.(I, 2)))
  zs = extrema(unique(getindex.(I, 1)))
  return xgrid[range(xs...)], zgrid[range(zs...)], m[range(zs...), range(xs...)]
end

# ╔═╡ 9494e2a6-e2fd-4728-94d4-d68816a00e72
md"""
### Fields
"""

# ╔═╡ 9b744e72-07a2-4f8c-b88c-e0b1131794d0
function initialize_fields(pa, nt; snap_store=false)
  nz, nx = pa.nz, pa.nx
  f = (; vy=zeros(Float32, nz, nx),
    dvydx=zeros(Float32, nz, nx),
    dvydz=zeros(Float32, nz, nx),
    σyx=zeros(Float32, nz, nx),
    dσyxdx=zeros(Float32, nz, nx),
    σyz=zeros(Float32, nz, nx),
    dσyzdz=zeros(Float32, nz, nx))
  if snap_store
    fs = (; vys=[zeros(Float32, nz, nx) for i in 1:nt], σyxs=[zeros(Float32, nz, nx) for i in 1:nt], σyzs=[zeros(Float32, nz, nx) for i in 1:nt])
    return merge(f, fs)
  else
    return f
  end
end

# ╔═╡ 31d742a4-100f-4744-afe5-381b265b6f4c
function reset_fields!(fields)
  fill!(fields.vy, zero(Float32))
  fill!(fields.dvydx, zero(Float32))
  fill!(fields.dvydz, zero(Float32))
  fill!(fields.σyx, zero(Float32))
  fill!(fields.σyz, zero(Float32))
  fill!(fields.dσyxdx, zero(Float32))
  fill!(fields.dσyzdz, zero(Float32))
end

# ╔═╡ bd8b9a0c-6664-4ba8-9795-0a496e932d85
md"### Gradients"

# ╔═╡ 7d00ecf7-80ba-4eee-99f6-b03f7b1f6e52
function initialize_grad(pa, nt)
  nz, nx = pa.nz, pa.nx
  return (; g=zeros(Float32, 2 * nz * nx), gρ=zeros(Float32, nz, nx), ginvμ=zeros(Float32, nz, nx))
end

# ╔═╡ d9ddd356-8276-45e3-8f60-87779fae270e
a=zeros(100)

# ╔═╡ 5d7d9a8d-0c96-4533-862b-98418b84566b
function update_gx!(gx, g, x, ref=medium_ref_values)
  xρ, xinvμ = chunk(x, 2)
  gρ, ginvμ = chunk(g, 2)
  gxρ, gxinvμ = chunk(gx, 2)
  map!(gxinvμ, ginvμ, xinvμ) do g, m
    ref.invμ0 * exp(m) * g
  end
  map!(gxρ, gρ, xρ) do g, m
    ref.ρ0 * exp(m) * g
  end
  return x
end

# ╔═╡ 4e1a0d4b-5f25-4b25-8dbb-4069e38dc5c4
# create a function to compute grad_phi and grad_mu 
# choosing t_grad as tgrid will give an accurate gradient of the least-squares objective function
function reduce_gradients!(grad, x, forwfields, adjfields, pa, t_grad=pa.tgrid)

  (; g, gρ, ginvμ) = grad

  fill!(g, zero(Float32))
	fill!(gρ, zero(Float32))
fill!(ginvμ, zero(Float32))
 
	
  (; nx, nz, tarray, dt, tgrid) = pa

  nt1=argmin(abs.(tgrid .- first(t_grad)))
  nt2=argmin(abs.(tgrid .- last(t_grad)))

  # term u[nt]*v[1]

   @. gρ =  gρ + adjfields.vys[nt2] * forwfields.vys[nt2]


  # term τ[nt-1]*σ[1]
 @. ginvμ = ginvμ + adjfields.σyxs[nt2-1] * forwfields.σyxs[nt2] + adjfields.σyzs[nt2-1]* forwfields.σyzs[nt2]

  # for gradρ
  for it = nt1:nt2-1
   
		@. gρ =  gρ + adjfields.vys[it] * (forwfields.vys[it+1] - forwfields.vys[it])
	
    # @. gρ = gρ * pa.tarray

    # for gradμ
    for it = nt1+1:nt2-1
      @. ginvμ = ginvμ + adjfields.σyxs[it-1] * (forwfields.σyxs[it+1] - forwfields.σyxs[it]) +  adjfields.σyzs[it-1] * (forwfields.σyzs[it+1] - forwfields.σyzs[it])
      # @. gμ = gμ * pa.tarray
    end
  end
  # multiply with dt
  rmul!(ginvμ, inv(dt))
rmul!(gρ, inv(dt))
	
gρ1, ginvμ1 = chunk(g, 2)

  copyto!(gρ1, gρ)
  copyto!(ginvμ1, ginvμ)
	# parameterization
  update_gx!(g, g, x)
end

# ╔═╡ f95a08ce-a38d-4b7f-b478-4dbfa607740e
md"### Wavelets"

# ╔═╡ 90953c64-8a87-4065-8c34-d0ead540b728
md"""
Generate a Ricker Wavelet. Reference:
Frequencies of the Ricker wavelet, Yanghua Wang, GEOPHYSICS, VOL. 80, NO. 2.
Its bandwidth is roughly 1.2 * fpeak.

* `fqdom::Float32`: dominant frequency 
* `tgrid`: time-domain grid
* `tpeak::Float32=tgrid[1]+1.5/fqdom`: the peak of the ricker in time (has a default)
"""

# ╔═╡ b993e3e6-8d4e-4de7-aa4b-6a2e3bd12212
function ricker(fqdom,
  tgrid;
  # tpeak is the location of the peak amplitude
  tpeak=tgrid[1] + 1.5 / fqdom, # using approximate half width of ricker
  trim_tol=0.0,
  maxamp=1.0
)
  (tpeak < tgrid[1] + 1.5 / fqdom) && error("cannot output Ricker for given tgrid and tpeak")
  (tpeak > tgrid[end] - 1.5 / fqdom) && error("cannot output Ricker for given tgrid and tpeak")

  isapprox(fqdom, 0.0) && error("dominant frequency cannot be zero")

  # some constants
  pf = (π * π) * (fqdom^2.0)
  nt = length(tgrid)
  δt = step(tgrid)

  # initialize
  wav = zero(tgrid)

  #! ricker wavelet
  for it = 1:nt
    tsquare = (tgrid[it] - tpeak) * (tgrid[it] - tpeak)
    wav[it] = (1.0 - 2.0 * pf * tsquare) * exp(-1.0e0 * pf * tsquare) * maxamp
  end

  isapprox(maximum(abs.(wav)), 0.0) && warn("wavelet is zeros")
  return Float32.(wav)
end

# ╔═╡ 17dd3d57-d5ca-443c-b003-b3a97b963d57
#=╠═╡
begin
	source_fpeak = 2.0 # in Hz
	source_wavelet = ricker(source_fpeak, tgrid, maxamp=1e15)
	source = repeat(source_wavelet, 1, ageom.ns)
end;
  ╠═╡ =#

# ╔═╡ ebab6005-2ad6-4057-9275-bf7d53d41b0b
#=╠═╡
begin
	# Choosing the extent of taper for absorbing boundaries
	true_medium_lambda = get_vs(mean,medium_true)/source_fpeak
	taper_points = floor(Int, 2*true_medium_lambda/dz)

	#NamedTuple for grid-related parameters
	grid_param = (; xgrid, zgrid, tgrid, dt=step(tgrid), nt=length(tgrid), nx=length(xgrid), nz=length(zgrid), tarray=get_taper_array(nx, nz, np=taper_points, tapfact=0.1))
end;
  ╠═╡ =#

# ╔═╡ 661c49e6-b5e0-4d85-a41f-b433936bc522
#=╠═╡
true_medium_lambda / step(xgrid)
  ╠═╡ =#

# ╔═╡ ab25a039-eba2-48d0-9338-bf304e5cdbc8
#=╠═╡
true_medium_lambda / step(zgrid)
  ╠═╡ =#

# ╔═╡ dd78c460-a446-46f5-af57-7e859383348d
#=╠═╡
taper_points
  ╠═╡ =#

# ╔═╡ 53ebecf7-5ea8-4372-9be2-fa48bd2be130
#=╠═╡
gradient = initialize_grad(grid_param, grid_param.nt);
  ╠═╡ =#

# ╔═╡ e586a423-b66b-455d-a88c-8ea70ad7ee2c
#=╠═╡
do_fd_tests && (g2ρ=get_xs(gradient.g, xs)[1:length(xs)]) # Gradients wrt ρ computed via adjoint-state method
  ╠═╡ =#

# ╔═╡ e3f3b379-4add-4866-8472-7bc7e53a7a28
#=╠═╡
do_fd_tests && (g2invμ=get_xs(gradient.g, xs)[length(xs)+1:end]) # Gradients wrt μ⁻¹ computed via adjoint-state method
  ╠═╡ =#

# ╔═╡ c34b1a5d-5078-4b8f-94d1-a088cbe5ab3e
#=╠═╡
plot(heatmap(x=xgrid, y=zgrid, z=(grid_param.tarray)))
  ╠═╡ =#

# ╔═╡ 15bbf544-34bd-4d38-bac5-0f43b1305df3
#=╠═╡
function propagate!(data, fields, pa, medium, ageom, source)

	# source spray and receiver projection matrix
	Rs=get_projection_matrix(ageom.xs, ageom.zs, grid_param.xgrid, grid_param.zgrid, transpose_flag=false, m=medium.invρ)
	Rr=get_projection_matrix(ageom.xr, ageom.zr, grid_param.xgrid, grid_param.zgrid, transpose_flag=true)

    reset_fields!(fields)

    (; vy, dvydx, dvydz, σyx, σyz, dσyxdx, dσyzdz) = fields
    (; nx, nz, tarray, tgrid, dt, nt) = pa

    (; μ, ρ, invρ) = medium

    # time loop
    # @progress 
	for it = 1:nt
		# ============= update vy ================
        Dx!(dσyxdx, σyx)
        Dz!(dσyzdz, σyz)

        @. vy = vy + invρ * (dσyxdx + dσyzdz) * dt

		# ============== body force ==============
		# need to view vy as a vector for source/recording operations
        vyv = view(vy, :)
        f = view(source, it, :)
	    mul!(vyv, Rs, f, dt, 1.0)

		# ============= update σ ================
        Dx!(dvydx, vy)
        Dz!(dvydz, vy)

        @. σyx = σyx + μ * dvydx * dt
        @. σyx = σyx * pa.tarray
        @. σyz = σyz + μ * dvydz * dt
        @. σyz = σyz * pa.tarray


        # ============= record fields ===========
        d = view(data, it, :)
        mul!(d, Rr, vyv)

        (:vys ∈ keys(fields)) && copyto!(fields.vys[it], vy)
        (:σyxs ∈ keys(fields)) && copyto!(fields.σyxs[it], σyx)
        (:σyzs ∈ keys(fields)) && copyto!(fields.σyzs[it], σyz)

    end

    return nothing
end
  ╠═╡ =#

# ╔═╡ d812711d-d02f-44bb-9e73-accd1623dea1
#=╠═╡
plot(tgrid, source_wavelet, size=(500, 200), w=2, label="Source Wavelet", Layout(width=500, height=250))
  ╠═╡ =#

# ╔═╡ b6001e39-db63-4f8e-a975-1185d090a744
#=╠═╡
source_wavelet |> typeof
  ╠═╡ =#

# ╔═╡ 77c9696c-58c5-40bf-acd0-16d5cf877810
#=╠═╡
begin
	# Initialisation of fields and data
	fields_true = initialize_fields(grid_param, grid_param.nt)
	dobs = initialize_data(grid_param, ageom)

	# Running the simulation to generate observed data
	@time propagate!(dobs, fields_true, grid_param, medium_true, ageom, source)	
end;
  ╠═╡ =#

# ╔═╡ 3be62716-f2d9-434c-a69a-ed272b89c85d
#=╠═╡
begin
	# Initialisation of fields and data
	fields_forw = initialize_fields(grid_param, grid_param.nt, snap_store=true)
	dref = initialize_data(grid_param, ageom)

	# Simulation to compute wavefields
	propagate!(dref, fields_forw, grid_param, medium_ref, ageom, source)

	# reverse the time order of fields stored 
	reverse!(fields_forw.vys); reverse!(fields_forw.σyxs); reverse!(fields_forw.σyzs)
end;
  ╠═╡ =#

# ╔═╡ ee679ef2-3cf8-4b3d-a717-ae2d088b5fe8
#=╠═╡
adj_source = get_adj_source(dobs, dref);
  ╠═╡ =#

# ╔═╡ 3f5f9d8a-3647-4a16-89ba-bd7a31c01064
#=╠═╡
begin
	# Initialisation of fields and data
	fields_adj = initialize_fields(grid_param, nt, snap_store=true)
	dadj = initialize_data(grid_param, adj_ageom)

	# Simulating the adjoint field
	propagate!(dadj, fields_adj, grid_param, medium_ref, adj_ageom, adj_source)
end;
  ╠═╡ =#

# ╔═╡ 2df3bd5b-630e-451a-a207-2c3a1719916e
#=╠═╡
reduce_gradients!(gradient, x0, fields_forw, fields_adj, grid_param, t_grad);
  ╠═╡ =#

# ╔═╡ a42d3b46-ae60-41d1-8b2d-e85af895ec14
#=╠═╡
fwi_param = (; fields_forw = initialize_fields(grid_param, nt, snap_store=true), fields_adj = initialize_fields(grid_param, nt, snap_store=true), medium=deepcopy(medium_ref), dref, ageom, adj_ageom, source, adj_source, dobs, grid_param, xbuffer=get_x(medium_ref), gradient=initialize_grad(grid_param, grid_param.nt))
  ╠═╡ =#

# ╔═╡ a24bfc41-8e50-4a68-b5cf-3973f4003221
#=╠═╡
function J(x; fwi_param=fwi_param)

	(; dref, fields_forw, grid_param, medium, ageom, source) = fwi_param
	update_medium!(medium, x)
	propagate!(dref, fields_forw, grid_param, medium, ageom, source)

	return Jdata(dref, dobs)
end
  ╠═╡ =#

# ╔═╡ 008b99f2-c9f3-4d30-a234-393e1ed69840
#=╠═╡
function Jsρ(xs; fwi_param=fwi_param)
	# Function to check gradients wrt ρ
	update_xsρ!(fwi_param.xbuffer, xs)
	return J(fwi_param.xbuffer, fwi_param=fwi_param)
end
  ╠═╡ =#

# ╔═╡ 12a089f0-23b5-4091-8861-d3ba2d0073a0
#=╠═╡
do_fd_tests && @time Jsρ(xs)
  ╠═╡ =#

# ╔═╡ 50733229-38f1-4ac1-acbc-ebb2c92d3891
#=╠═╡
do_fd_tests && (g1ρ=grad(central_fdm(2, 1), Jsρ, xs)) # Gradients wrt ρ using central difference
  ╠═╡ =#

# ╔═╡ b525408f-0d7d-4333-a789-def42565520c
#=╠═╡
do_fd_tests && (g1ρ[1] ./ g2ρ)
  ╠═╡ =#

# ╔═╡ 9d319561-38c9-46c8-aaf6-06d0a41ed0bf
#=╠═╡
function Jsinvμ(xs; fwi_param=fwi_param)
	# Function to check gradients wrt invμ
	update_xsinvμ!(fwi_param.xbuffer, xs)
	return J(fwi_param.xbuffer, fwi_param=fwi_param)
end
  ╠═╡ =#

# ╔═╡ 9aa2e4ac-b221-4f3e-9072-3d4d762f01c7
#=╠═╡
do_fd_tests && @time Jsinvμ(xs)
  ╠═╡ =#

# ╔═╡ 803ac9ba-93d2-4f66-9018-36232b8a3076
#=╠═╡
do_fd_tests && (g1invμ=grad(central_fdm(2, 1), Jsinvμ, xs)) # Gradients wrt μ⁻¹ using central difference
  ╠═╡ =#

# ╔═╡ 99541b49-caf2-40ab-b299-081111e35675
#=╠═╡
do_fd_tests && (g1invμ[1] ./ g2invμ)
  ╠═╡ =#

# ╔═╡ fbe44944-499a-4881-94b6-07855d1165aa
md"""
### Plots
"""

# ╔═╡ a055de92-9eba-4ca8-a579-f00bfa855e4a
md"This function plots the wavefield after clipping the edges."

# ╔═╡ a72184ba-f642-4bed-8a31-c70d5fdb88e2
function add_ageom!(fig, ageom, row, col)
  if (!(ageom === nothing))
    add_trace!(fig, scatter(
        x=ageom.xr,
        y=ageom.zr, mode="markers",
        marker_color="black", marker_symbol="triangle-down", showlegend=false), row=row, col=col)
    add_trace!(fig, scatter(
        x=ageom.xs,
        y=ageom.zs, mode="markers",
        marker_color="black", marker_size=10, marker_symbol="star", showlegend=false), row=row, col=col)
  end

end

# ╔═╡ f1b9f4a0-5554-48fd-9268-4d7b03d74688
function fieldheat(fields, titles, grid_param, ageom=nothing)

  @assert length(fields) == length(titles) == 4
  fig = Plot(Layout(yaxis_autorange="reversed", title=attr(font_size=12,), font=attr(
      size=10), yaxis=attr(scaleanchor="x", scaleratio=1), Subplots(shared_xaxes=true, shared_yaxes=true, horizontal_spacing=0.04, vertical_spacing=0.04, rows=div(length(fields) - 1, 2) + 1, cols=2, subplot_titles=titles)))
  i = 0
  for f in fields
    i = i + 1
    ic = 2 - mod(i, 2)
    ir = div(i - 1, 2) + 1
    x, y, z = clip_edges(f, grid_param)
    dmax = maximum(abs, f)
    add_trace!(fig, heatmap(
        xmin=-30, xmax=30,
        x=x,
        y=y,
        z=z, zmin=-dmax,
        zmax=dmax, colorscale="seismic", showscale=false), row=ir, col=ic)
    add_ageom!(fig, ageom, ir, ic)
  end
  return PlutoPlotly.plot(fig)

end

# ╔═╡ 77e134d8-bd8b-4303-8c44-a4920cf0ee81
#=╠═╡
let t=mod(t_forw, grid_param.nt-1)
	t_grad
fieldheat([fields_forw.vys[t], fields_adj.vys[nt-t], gradient.gρ, gradient.ginvμ], 
	["Forward Field" "Adjoint Field" "Gradient w.r.t. ρ" "Gradient w.r.t. μ⁻¹"], 
	grid_param, ageom)
end
  ╠═╡ =#

# ╔═╡ 99ba8f6a-551a-432b-abb1-a79be233fa46
#=╠═╡
function mediumheat(medium, ageom=nothing)
 	(; μ, ρ) = medium
	c = get_vs(medium)

    fig = Plot(Layout(yaxis_autorange="reversed", height=300, width=650, title=attr(font_size=12,), font=attr(
            size=10), yaxis=attr(scaleanchor="x"), Subplots(shared_xaxes=true, shared_yaxes=true, horizontal_spacing=0.3, rows=1, cols=2, subplot_titles=["Velocity (km/s)" "Density (g/cc)"])))
    add_trace!(fig, heatmap(
            x=xgrid,
            y=zgrid,
            z=c, colorscale="Portland", colorbar_x=0.35), row=1, col=1)


    add_trace!(fig, heatmap(
            x=xgrid,
            y=zgrid,
            z=ρ * 1e-12, colorscale="Portland", colorbar_x=1.0), row=1, col=2)
	if(!(ageom===nothing))
		for ic in 1:2
			add_ageom!(fig, ageom, 1, ic)
		end
	end

    return PlutoPlotly.plot(fig)

end
  ╠═╡ =#

# ╔═╡ 55c7f981-96a7-40e9-811f-37334622565b
#=╠═╡
mediumheat(medium_true, ageom)
  ╠═╡ =#

# ╔═╡ f26f1b1a-18e0-413c-86a9-351ba5dfaebf
#=╠═╡
mediumheat(medium_ref, ageom)
  ╠═╡ =#

# ╔═╡ 54986cc5-2ea0-4097-bbe0-1ed174ec9ae4
function dataheat(d, tgrid; title="Data")
  data = d
  dmax = maximum(abs, data)
  trace = heatmap(
    z=data,
    x=1:size(data, 2),
    y=tgrid,
    colorscale="RdBu",
    reversescale=true,
    zmin=-dmax,
    zmax=dmax,
  )
  layout = Layout(title=title, width=225, height=225, xaxis_range=[0, size(data, 2) + 1], yaxis_scaleanchor="x", yaxis_scaleratio=2.5, yaxis_autorange="reversed",
    xaxis_title="Receiver #",
    yaxis_title="Time (s)",)
  return plot(trace, layout)
end

# ╔═╡ 4171af00-1d14-45ba-9fd3-a2c30d0b759f
#=╠═╡
ThreeColumn(md"""
$(dataheat(dobs, tgrid, title="Observed Data"))""", 
	md"""
	$(dataheat(dref, tgrid, title="Modelled Data"))
	""",
	md"""
	$(dataheat(adj_source, tgrid, title="Adjoint Sources"))
	""")
  ╠═╡ =#

# ╔═╡ 802d9652-7597-43c4-b13a-3c60682d0a69
md"""
We now know how to calculate the derivatives using  functions `Dx` and `Dz`. Let's now discretize the time dimension and formulate an explicit time-stepping (leap-frog) scheme and alternatively update the velocity and stress fields. We will start with the velocity field
```math
\rho{\dot{v}_y}= \mathrm{Dx}(\sigma_{yx}[it-\tfrac{1}{2}])+ \mathrm{Dz}(\sigma_{yz}[it-\tfrac{1}{2}]),
```
where the finite-difference approximation for the time derivative for timestep $dt$ is:
```math
\frac{\rho({{v_y}[it]-{v_y}[it-1]})}{dt}= \mathrm{Dx}(\sigma_{yx}[it-\tfrac{1}{2}])+ \mathrm{Dz}(\sigma_{yz}[it-\tfrac{1}{2}]),
```
```math
{{v_y}[it]-{v_y}[it-1]}= \frac{dt}{\rho}\left(\mathrm{Dx}(\sigma_{yx}[it-\tfrac{1}{2}])+ \mathrm{Dz}(\sigma_{yz}[it-\tfrac{1}{2}])\right),
```
```math
{v_y}[it]= {v_y}[it-1]+\frac{dt}{\rho}(\mathrm{Dx}(\sigma_{yx}[it-\tfrac{1}{2}])+ \mathrm{Dz}(\sigma_{yz}[it-\tfrac{1}{2}])).
```
Similarly, for updating the stresses,
```math
{\dot{\sigma}_{21}}= \mu\frac{dv_y}{dx},
```
```math
{\dot{\sigma}_{23}}= \mu\frac{dv_y}{dz},
```
we can write

```math
\sigma_{21}[it+\tfrac{1}{2}]=\sigma_{21}[it-\tfrac{1}{2}]+\mu\mathrm{Dx}(v_y)dt,
```
```math
\sigma_{23}[it+\tfrac{1}{2}]=\sigma_{23}[it-\tfrac{1}{2}]+\mu\mathrm{Dz}(v_y)dt.
```
Notice that the temporal grid of stress and velocity are staggered. Finally, let's write a function now that leaps by a given number of steps.
"""

# ╔═╡ Cell order:
# ╠═3c540889-49dc-415c-acbc-3494897b260c
# ╟─9c32f5bc-f6d1-4048-903a-27224aaa1f40
# ╠═9ab80963-a575-4a1a-a529-f0f2e64af4cc
# ╟─881c7368-57df-469a-96ec-821512cf98e0
# ╟─65efba3b-16b0-4113-a59e-809365e7bdd6
# ╟─55c7f981-96a7-40e9-811f-37334622565b
# ╟─fa78af13-e3c6-4d6f-8a3e-1187fe9ae159
# ╠═37157e04-273a-4be2-9687-2016480a0421
# ╟─f26f1b1a-18e0-413c-86a9-351ba5dfaebf
# ╠═e233afec-6049-4277-8be9-95687c4589b5
# ╟─4171af00-1d14-45ba-9fd3-a2c30d0b759f
# ╟─b572f855-db01-4c4f-8922-968bc0ef5fdf
# ╟─12a0f7d8-64e0-411f-85ad-14d48bd9492d
# ╠═92d05447-6d0e-4ee0-a330-244b9c65c871
# ╠═77e134d8-bd8b-4303-8c44-a4920cf0ee81
# ╟─b3015aff-bdab-4c94-87da-8744c174263a
# ╟─4adbb7f0-7927-470d-8c00-07d3b3c0cd78
# ╟─7fd6fe40-0065-4074-b8cf-3d8cd8e68319
# ╟─122c2aa4-baec-4288-ab53-afa1d977c486
# ╠═7913b296-3010-410d-942f-834a47d599c7
# ╠═77fa76f3-ffda-4d95-8c12-7ccce6a7e52e
# ╠═7928a88c-f217-4338-a6a5-50ab2d422480
# ╠═6937b103-9ce2-4189-8129-aae1e7936d4f
# ╠═661c49e6-b5e0-4d85-a41f-b433936bc522
# ╠═ab25a039-eba2-48d0-9338-bf304e5cdbc8
# ╟─cb83dbd1-c423-4b27-b29e-7dc8051f43d5
# ╟─86ed93a1-c7b9-4b3a-8cb5-ca4405cff3df
# ╠═d39753e2-5986-4394-9293-9e394f2807f0
# ╟─2ec95e6e-7c2c-41e5-87b7-84583564f079
# ╠═0b890ec0-1886-4494-b4ac-46de7639f358
# ╠═1c67cda8-7712-4d5b-a2aa-af47f290f745
# ╠═8b3776bd-509b-4232-9737-36c9ae003350
# ╟─5b271e5f-879c-4c43-825a-9660f322febd
# ╠═6a8139c4-12c1-4d18-bd1e-14334290aec1
# ╠═dd78c460-a446-46f5-af57-7e859383348d
# ╠═ebab6005-2ad6-4057-9275-bf7d53d41b0b
# ╟─2855c8cf-8364-4c6c-a122-781b99440e89
# ╠═17dd3d57-d5ca-443c-b003-b3a97b963d57
# ╠═d812711d-d02f-44bb-9e73-accd1623dea1
# ╠═b6001e39-db63-4f8e-a975-1185d090a744
# ╟─8b5372ec-0742-48ea-81c4-d303b96f56c7
# ╟─14c0cd30-a52c-4f93-9f7b-cae6328c0655
# ╠═77c9696c-58c5-40bf-acd0-16d5cf877810
# ╠═3be62716-f2d9-434c-a69a-ed272b89c85d
# ╟─439df138-5198-4737-915d-7bd2c157aa4b
# ╠═a42d3b46-ae60-41d1-8b2d-e85af895ec14
# ╠═0454ce0a-d6de-427f-bbdc-3bcec21327f2
# ╠═f2fb92bb-33d6-4e15-8caf-b245e000ad69
# ╠═4a4f1300-94ba-4c1b-9d10-9e955d194ff6
# ╠═a24bfc41-8e50-4a68-b5cf-3973f4003221
# ╠═a3a9deea-e2d6-4d58-90d7-5a54be176289
# ╠═e7e72f61-d79e-4a82-ad9c-129191c8a8c2
# ╟─1dc2ca10-5ba3-4efa-b6d9-d203cf91598b
# ╠═ee679ef2-3cf8-4b3d-a717-ae2d088b5fe8
# ╠═ddb37082-cae0-4a68-ab55-19563d8727ed
# ╠═3f5f9d8a-3647-4a16-89ba-bd7a31c01064
# ╟─66f9c698-61e3-4b61-aff3-dfc67eb2f6af
# ╠═22ac08a7-f4d7-4809-8ee3-903d96c96cd6
# ╠═53ebecf7-5ea8-4372-9be2-fa48bd2be130
# ╠═2df3bd5b-630e-451a-a207-2c3a1719916e
# ╟─8d161f09-8339-4277-8739-ff76607f7abf
# ╠═008b99f2-c9f3-4d30-a234-393e1ed69840
# ╠═9d319561-38c9-46c8-aaf6-06d0a41ed0bf
# ╠═12a089f0-23b5-4091-8861-d3ba2d0073a0
# ╠═9aa2e4ac-b221-4f3e-9072-3d4d762f01c7
# ╠═006739fb-24a1-49b0-9619-fe8e2d3c8fca
# ╠═50733229-38f1-4ac1-acbc-ebb2c92d3891
# ╠═e586a423-b66b-455d-a88c-8ea70ad7ee2c
# ╠═b525408f-0d7d-4333-a789-def42565520c
# ╠═803ac9ba-93d2-4f66-9018-36232b8a3076
# ╠═e3f3b379-4add-4866-8472-7bc7e53a7a28
# ╠═99541b49-caf2-40ab-b299-081111e35675
# ╟─fc49a6d7-a1b1-458a-a9ad-e120282bbabc
# ╠═c5815f5e-9164-11ec-10e1-691834761dff
# ╟─a62839d5-837b-4c37-996f-33659c34911c
# ╟─3fc0e673-2fa3-489f-a56e-a867ea37cbce
# ╠═2ea24e92-d66e-4c60-ad0b-f671d894fef2
# ╟─7f797571-055e-4975-9c26-fc968bbc0094
# ╠═d7b37c59-e0b3-4e47-86d3-7f1df7400f09
# ╟─6be2f4c2-e9ed-43c2-b66c-ef3176bb9000
# ╠═aa19e992-2735-4324-8fd7-15eacadf0faa
# ╠═e2127d9b-f2a4-4970-a36e-5fa70c304ca7
# ╟─e8333b23-53c3-445e-9ca3-6b278359f8ab
# ╠═9248af7f-dc1a-4bf6-8f3f-304db73fc604
# ╠═00a637a6-ddc4-4830-be65-1891d3cb18bc
# ╠═e08bf013-00c7-4870-82d8-19b899e7208d
# ╟─ab8b1a22-ca7a-409e-832e-8d5d08a29a1e
# ╠═f4d91971-f806-4c5c-8548-b58a20acfb2c
# ╟─9bc38d55-285b-4b83-98d9-d7f9e03405d1
# ╠═27844886-0b54-4b08-a592-a1a38e4b0be2
# ╠═489dcf10-b7f2-4544-b80d-3588ff00ff4a
# ╠═02844bc9-9faf-413f-82d2-28d6ba01668d
# ╠═199908d6-f5fd-4461-be53-39919a103835
# ╠═8298ae48-0ddc-49a9-a43f-8434e4cc3758
# ╠═ff1ba1f6-f127-4c47-8198-aeff5f051ad9
# ╠═cb6adea2-9ea2-4857-b969-540a3439e700
# ╠═ad21da29-f6ff-4a94-be87-4e88640cddbf
# ╠═2abfb6d5-a201-405e-99c3-42bf0ad5aad5
# ╠═a26a5bd6-f4b8-410a-b437-06587011093e
# ╠═78c09864-1fb1-427e-b77d-db82945e850c
# ╠═daef685e-6c07-419f-9c68-a1c2c0ec8973
# ╠═5648367f-8628-4859-a2a5-ab92bab54e9b
# ╠═c7ef6f8b-0c01-4788-982f-17ed6d0098ff
# ╠═dc95afc5-0232-4c2a-801e-dccac012d2e6
# ╠═13ffff0c-23cd-4a98-8ce8-31644b9429d3
# ╠═63a177c1-034e-4aa6-9951-367570c49850
# ╟─ae8012be-e7ab-4e85-a27f-febf08b3380b
# ╠═ce89710a-48d8-46d4-83b4-163a7eb0a2e5
# ╠═b7f4078a-ead0-4d42-8b44-4f471eefc6fc
# ╠═c34b1a5d-5078-4b8f-94d1-a088cbe5ab3e
# ╟─9494e2a6-e2fd-4728-94d4-d68816a00e72
# ╠═9b744e72-07a2-4f8c-b88c-e0b1131794d0
# ╠═31d742a4-100f-4744-afe5-381b265b6f4c
# ╠═15bbf544-34bd-4d38-bac5-0f43b1305df3
# ╟─bd8b9a0c-6664-4ba8-9795-0a496e932d85
# ╠═7d00ecf7-80ba-4eee-99f6-b03f7b1f6e52
# ╠═d9ddd356-8276-45e3-8f60-87779fae270e
# ╠═4e1a0d4b-5f25-4b25-8dbb-4069e38dc5c4
# ╠═5d7d9a8d-0c96-4533-862b-98418b84566b
# ╟─f95a08ce-a38d-4b7f-b478-4dbfa607740e
# ╟─90953c64-8a87-4065-8c34-d0ead540b728
# ╠═b993e3e6-8d4e-4de7-aa4b-6a2e3bd12212
# ╟─fbe44944-499a-4881-94b6-07855d1165aa
# ╟─a055de92-9eba-4ca8-a579-f00bfa855e4a
# ╠═f1b9f4a0-5554-48fd-9268-4d7b03d74688
# ╠═a72184ba-f642-4bed-8a31-c70d5fdb88e2
# ╠═99ba8f6a-551a-432b-abb1-a79be233fa46
# ╠═54986cc5-2ea0-4097-bbe0-1ed174ec9ae4
# ╟─802d9652-7597-43c4-b13a-3c60682d0a69
