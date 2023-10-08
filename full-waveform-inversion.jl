### A Pluto.jl notebook ###
# v0.19.27

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
end;

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

# ╔═╡ 1de42ac5-ddaa-4e51-a093-99423a2993ae
#=╠═╡
dobs
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LossFunctions = "30fc2ffe-d236-52d8-8643-a9d8f7c094a7"
MLUtils = "f1d291b0-491e-4a28-83b9-f70985020b54"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DSP = "~0.7.8"
FFTW = "~1.6.0"
FiniteDifferences = "~0.12.26"
ForwardDiff = "~0.10.35"
LaTeXStrings = "~1.3.0"
LossFunctions = "~0.9.0"
MLUtils = "~0.4.1"
PlutoPlotly = "~0.3.6"
PlutoTeachingTools = "~0.2.10"
PlutoUI = "~0.7.50"
ProgressLogging = "~0.1.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "3923f81a1eb85684a638bcb161bbd3e3f6aad7fa"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "16b6dbc4cf7caee4e1e75c49485ec67b667098a0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Accessors]]
deps = ["Compat", "CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Requires", "Test"]
git-tree-sha1 = "c7dddee3f32ceac12abd9a21cd0c4cb489f230d2"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.29"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

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

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "5084cc1a28976dd1642c9f337b28a3cb03e0f7d2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.7"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"
weakdeps = ["ChainRulesCore"]

    [deps.ChangesOfVariables.extensions]
    ChangesOfVariablesChainRulesCoreExt = "ChainRulesCore"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "d730914ef30a06732bdd9f763f6cc32e92ffbff1"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

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

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"

    [deps.ConstructionBase.extensions]
    IntervalSetsExt = "IntervalSets"
    StaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "da8b06f89fce9996443010ef92572b193f8dca1f"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.8"

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
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

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
git-tree-sha1 = "f9818144ce7c8c41edf5c4c179c684d92aa4d9fe"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.6.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "ffb97765602e3cbe59a0589d237bf07f245a8576"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.1"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FiniteDifferences]]
deps = ["ChainRulesCore", "LinearAlgebra", "Printf", "Random", "Richardson", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "3f605dd6db5640c5278f2551afc9427656439f42"
uuid = "26cc04aa-876d-5657-8c51-4c34ba976000"
version = "0.12.26"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.FoldsThreads]]
deps = ["Accessors", "FunctionWrappers", "InitialValues", "SplittablesBase", "Transducers"]
git-tree-sha1 = "eb8e1989b9028f7e0985b4268dabe94682249025"
uuid = "9c68100b-dfe1-47cf-94c8-95104e173443"
version = "0.1.1"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

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

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0cb9352ef2e01574eeebdb102948a58740dcaf83"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

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
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "6a125e6a4cb391e0b9adbd1afa9e771c2179f8ef"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.23"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "LinearAlgebra", "MacroTools", "PrecompileTools", "SparseArrays", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "1e7e27a144936ed6f1b0a01dbc7b7f86afabeb6e"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.3"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "a8960cae30b42b66dd41808beb76490519f6f9e2"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "5.0.0"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "09b7505cc0b1cee87e5d4a26eea61d2e1b0dcd35"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.21+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "98dc144f1e0b299d49e8d23e56ad68d3e4f340a5"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.20"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    DiffEqBiologicalExt = "DiffEqBiological"
    ParameterizedFunctionsExt = "DiffEqBase"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    DiffEqBiological = "eb300fae-53e8-50a0-950c-e21f52c2b7e0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

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
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"
weakdeps = ["ChainRulesCore", "ChangesOfVariables", "InverseFunctions"]

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LossFunctions]]
deps = ["CategoricalArrays", "Markdown"]
git-tree-sha1 = "d4c7ff8c7281943371e1725000fd538a699024d0"
uuid = "30fc2ffe-d236-52d8-8643-a9d8f7c094a7"
version = "0.9.0"

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

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MLUtils]]
deps = ["ChainRulesCore", "Compat", "DataAPI", "DelimitedFiles", "FLoops", "FoldsThreads", "NNlib", "Random", "ShowCases", "SimpleTraits", "Statistics", "StatsBase", "Tables", "Transducers"]
git-tree-sha1 = "f69cdbb5b7c630c02481d81d50eac43697084fe0"
uuid = "f1d291b0-491e-4a28-83b9-f70985020b54"
version = "0.4.1"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "9926529455a331ed73c19ff06d16906737a876ed"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.3"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Pkg", "Random", "Requires", "Statistics"]
git-tree-sha1 = "99e6dbb50d8a96702dc60954569e9fe7291cc55d"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.8.20"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

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
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

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
git-tree-sha1 = "901509f7ae6abfb20281c53c34ecfa95ec3119df"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.10"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "434f66dfbb15606c49a7a21dc670119fdf729fa9"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.10"
weakdeps = ["ChainRulesCore", "MakieCore", "MutableArithmetics"]

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "bc2bda41d798c2e66e7c44a11007bb329b15941b"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.0.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

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
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "feafdc70b2e6684314e188d95fe66d116de834a7"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.2"

[[deps.Richardson]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e03ca566bec93f8a3aeb059c8ef102f268a38949"
uuid = "708f8203-808e-40c0-ba2d-98a6953ed40d"
version = "1.4.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

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
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "63e84b7fdf5021026d0f17f76af7c57772313d99"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.21"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

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

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c42fa452a60f022e9e087823b47e5a5f8adc53d5"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.75"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[deps.UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "ea37e6066bf194ab78f4e747f5245261f17a7175"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.1.2"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

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
# ╠═3c540889-49dc-415c-acbc-3494897b260c
# ╟─9c32f5bc-f6d1-4048-903a-27224aaa1f40
# ╟─881c7368-57df-469a-96ec-821512cf98e0
# ╟─65efba3b-16b0-4113-a59e-809365e7bdd6
# ╟─55c7f981-96a7-40e9-811f-37334622565b
# ╟─fa78af13-e3c6-4d6f-8a3e-1187fe9ae159
# ╠═37157e04-273a-4be2-9687-2016480a0421
# ╟─f26f1b1a-18e0-413c-86a9-351ba5dfaebf
# ╟─e233afec-6049-4277-8be9-95687c4589b5
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
# ╠═1de42ac5-ddaa-4e51-a093-99423a2993ae
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
