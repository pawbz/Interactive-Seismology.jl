### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "Seismic Full Waveform Inversion"
#> tags = ["imaging"]
#> layout = "layout.jlhtml"
#> description = "A notebook that demonstrates how interactivity makes seismic full waveform inversion more accessible and intuitive."

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
\sigma_{yx}=2\mu\epsilon_{yx}=\mu\partial_{x}u_y,

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

# ╔═╡ 12a089f0-23b5-4091-8861-d3ba2d0073a0
do_fd_tests && @time Jsρ(xs)

# ╔═╡ 9aa2e4ac-b221-4f3e-9072-3d4d762f01c7
do_fd_tests && @time Jsinvμ(xs)

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
DSP = "~0.7.9"
FFTW = "~1.8.0"
FiniteDifferences = "~0.12.31"
ForwardDiff = "~0.10.36"
LaTeXStrings = "~1.3.1"
LossFunctions = "~0.11.1"
MLUtils = "~0.4.4"
PlutoPlotly = "~0.4.6"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.58"
ProgressLogging = "~0.1.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "7875261f869a258ec326d65112f7c23fd47b3371"

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
git-tree-sha1 = "297b6b41b66ac7cbbebb4a740844310db9fd7b8c"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.1"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
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
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "7aa7ad1682f3d5754e3491bb59b8103cae28e3a3"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.40"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "3e93fcd95fe8db4704e98dbda14453a0bfc6f6c3"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.3"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "575cd02e080939a33b6df6c5853d14924c08e35b"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.23.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "c0216e792f518b39b22212127d4a84dc31e4e386"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.5"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

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
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

    [deps.CompositionsBase.weakdeps]
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

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
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "97d79461925cdb635ee32116978fc735b9463a39"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.19"

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
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

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
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

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
git-tree-sha1 = "d77e4697046989f44dce3ed66269aaf1611a3406"
uuid = "26cc04aa-876d-5657-8c51-4c34ba976000"
version = "0.12.31"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

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
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "e9648d90370e2d0317f9518c9c6e0841db54a90b"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.31"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "LinearAlgebra", "MacroTools", "PrecompileTools", "Requires", "SparseArrays", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "ed7167240f40e62d97c1f5f7735dea6de3cc5c49"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.18"

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"

    [deps.KernelAbstractions.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Preferences", "Printf", "Requires", "Unicode"]
git-tree-sha1 = "839c82932db86740ae729779e610f07a1640be9a"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "6.6.3"

    [deps.LLVM.extensions]
    BFloat16sExt = "BFloat16s"

    [deps.LLVM.weakdeps]
    BFloat16s = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "88b916503aac4fb7f701bb625cd84ca5dd1677bc"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.29+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cad560042a7cc108f5a4c24ea1431a9221f22c1b"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.2"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
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
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

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

[[deps.LossFunctions]]
deps = ["Markdown", "Requires", "Statistics"]
git-tree-sha1 = "df9da07efb9b05ca7ef701acec891ee8f73c99e2"
uuid = "30fc2ffe-d236-52d8-8643-a9d8f7c094a7"
version = "0.11.1"

    [deps.LossFunctions.extensions]
    LossFunctionsCategoricalArraysExt = "CategoricalArrays"

    [deps.LossFunctions.weakdeps]
    CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "31e27f0b0bf0df3e3e951bfcc43fe8c730a219f6"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.4.5"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MLUtils]]
deps = ["ChainRulesCore", "Compat", "DataAPI", "DelimitedFiles", "FLoops", "NNlib", "Random", "ShowCases", "SimpleTraits", "Statistics", "StatsBase", "Tables", "Transducers"]
git-tree-sha1 = "b45738c2e3d0d402dffa32b2c1654759a2ac35a4"
uuid = "f1d291b0-491e-4a28-83b9-f70985020b54"
version = "0.4.4"

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

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

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

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Pkg", "Random", "Requires", "Statistics"]
git-tree-sha1 = "1fa1a14766c60e66ab22e242d45c1857c83a3805"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.9.13"

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
git-tree-sha1 = "89f57f710cc121a7f32473791af3d6beefc59051"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.14"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "81a2a9462003a423fdc59e2a3ff84cde93c4637b"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.7"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

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
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

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
git-tree-sha1 = "12aa2d7593df490c407a3bbd8b86b8b515017f3e"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.14"

[[deps.Richardson]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "48f038bfd83344065434089c2a79417f38715c41"
uuid = "708f8203-808e-40c0-ba2d-98a6953ed40d"
version = "1.4.2"

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
git-tree-sha1 = "bf074c045d3d5ffd956fa0a461da38a44685d6b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.3"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

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

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

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
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "3064e780dbb8a9296ebb3af8f440f787bb5332af"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.80"

    [deps.Transducers.extensions]
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

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

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[deps.UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "323e3d0acf5e78a56dfae7bd8928c989b4f3083e"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.1.3"

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

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═3c540889-49dc-415c-acbc-3494897b260c
# ╟─9c32f5bc-f6d1-4048-903a-27224aaa1f40
# ╟─881c7368-57df-469a-96ec-821512cf98e0
# ╟─65efba3b-16b0-4113-a59e-809365e7bdd6
# ╟─55c7f981-96a7-40e9-811f-37334622565b
# ╟─fa78af13-e3c6-4d6f-8a3e-1187fe9ae159
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
