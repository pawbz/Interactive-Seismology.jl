### A Pluto.jl notebook ###
# v0.19.23

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
end

# ╔═╡ 3c540889-49dc-415c-acbc-3494897b260c
TableOfContents()

# ╔═╡ 9c32f5bc-f6d1-4048-903a-27224aaa1f40
md"""
# Full Waveform Inversion (SH Waves)

Seismic full waveform inversion (FWI) is a powerful technique for imaging the Earth's subsurface using seismic waves. However, FWI can be computationally intensive and require a deep understanding of numerical methods and algorithms. This notebook demonstrates how interactivity makes FWI more accessible and intuitive.
Pseudo-spectral methods are used to solve the 2-D SH stress-velocity system.
The particle velocity and stress fields are staggered in time and we use a simple explicit time-stepping approach.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 729689f3-7b76-4701-a60d-3471e1f20fce
@bind t_forw Clock(0.1)

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


# ╔═╡ 122c2aa4-baec-4288-ab53-afa1d977c486
md"""## Spatial Grids
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
## Acquistion Geometry
Choose the number of sources and receivers to be used in the simulation.
"""

# ╔═╡ 881c7368-57df-469a-96ec-821512cf98e0
md"""## True Medium
Consider a mantle-like medium that has density of $3.22$ $g/cc$ and shear wave velocity of $5$ $km$ $s^{-1}$. Suppose it has a reflector such that the medium has a higher density below the reflector compared to above. You can specify the location and density of the reflector below.
"""

# ╔═╡ 5b271e5f-879c-4c43-825a-9660f322febd
md"""## Time Grid
In order to study the propagation of the wave-front, we need to choose a time step that satisfies the Courant condition. Here, we have chosen the Courant number to be $0.2$.
"""

# ╔═╡ 2855c8cf-8364-4c6c-a122-781b99440e89
md"## Body Forces"

# ╔═╡ e233afec-6049-4277-8be9-95687c4589b5
md"""
## Generate Observed Data
Pseudo spectral method is applied for solving the seismic equation given the true medium parameters. Wavefields at the receivers are recorded and observed data is generated.
"""


# ╔═╡ 22b8db91-73a0-46df-87fd-7cf0b66ee37d
md"""
## Reference Medium
This will also serve as an initial model during inversion. Medium without the reflector is considered.
"""

# ╔═╡ c73f69d0-69e6-47f1-97b0-3e81218776e6
md"""
## Data Error
Deviation between the observed and model data is referred to as the data error and it acts as the forcing in the case of adjoint simulation.
"""

# ╔═╡ a3a9deea-e2d6-4d58-90d7-5a54be176289
md"## Adjoint Simulation"

# ╔═╡ afca5e37-ebc0-45e7-b27b-53602dbf6672
md"""
## Objective Function
"""

# ╔═╡ 66f9c698-61e3-4b61-aff3-dfc67eb2f6af
md"""
## Gradients
"""

# ╔═╡ 8d161f09-8339-4277-8739-ff76607f7abf
md"""
## Finite-Difference Test (Hockey-Stick Test)
"""

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

        rec = [md"""Number of receivers = $(Child("nr", Slider(range(start=5, stop=15, step=1), default=10, show_value=true)))
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
    Fz = plan_rfft(zeros(nz, nx), (1))
    Fx = plan_rfft(zeros(nz, nx), (2))
    kx = reshape(collect(rfftfreq(nx, inv(step(xgrid)))), 1, :) * 2 * pi
    kz = reshape(collect(rfftfreq(nz, inv(step(zgrid)))), :, 1) * 2 * pi
    storagex = zero(Fx * zeros(nz, nx))
    storagez = zero(Fz * zeros(nz, nx))
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
    x1 = rand(nz, nx)
    y1 = rand(nz, nx)
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

# ╔═╡ e08bf013-00c7-4870-82d8-19b899e7208d
# m are the medium properties that will be used 
function get_restriction_matrix(xpos, zpos, xgrid, zgrid, transpose_flag=false; m=ones(length(zgrid), length(xgrid)))
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


# ╔═╡ dd626606-2a4d-494d-ab74-92753072c773
function get_forcing_transform(ageom, medium, xgrid, zgrid, medium_flag=false)
    if medium_flag
        return (; Rs=get_restriction_matrix(ageom.xs, ageom.zs, xgrid, zgrid, false, m=medium.invρ),
            Rr=get_restriction_matrix(ageom.xr, ageom.zr, xgrid, zgrid, true), Rr_t=get_restriction_matrix(ageom.xr, ageom.zr, xgrid, zgrid, false, m=medium.invρ))
    else
        return (; Rs=get_restriction_matrix(ageom.xs, ageom.zs, xgrid, zgrid, false, m=medium.invρ), Rr=get_restriction_matrix(ageom.xr, ageom.zr, xgrid, zgrid, true))
    end
end;

# ╔═╡ ab8b1a22-ca7a-409e-832e-8d5d08a29a1e
md"### Data"

# ╔═╡ f4d91971-f806-4c5c-8548-b58a20acfb2c
function initialize_data(grid_param, ageom)
    [zeros(length(ageom.xr)) for t in grid_param.tgrid]
end

# ╔═╡ 9bc38d55-285b-4b83-98d9-d7f9e03405d1
md"### Medium"

# ╔═╡ 27844886-0b54-4b08-a592-a1a38e4b0be2
function bundle_medium(μ, ρ)
    return (; μ=μ, ρ=ρ, invρ=inv.(ρ))
end

# ╔═╡ 1c67cda8-7712-4d5b-a2aa-af47f290f745
begin
    # lets add density perturbation (olivine)
    μtrue = ones(nz, nx) .* 82 * 10^9 * 10^3
    ρtrue = ones(nz, nx) .* 3.22 * 10^-3 * 10^15 # density in kg/km3
    reflect_index = findall(zgrid .== param_true_med.z)[1]
    ρtrue[reflect_index:end, :] .= param_true_med.ρ * 10^-3 * 10^15  # density in kg/km3
    medium_true = bundle_medium(μtrue, ρtrue)
end;

# ╔═╡ 8b3776bd-509b-4232-9737-36c9ae003350
begin
    # Unperturbed medium
    μref = ones(nz, nx) * 82 * 10^9 * 10^3      # shear modulus in kg / km/s^2
    ρref = ones(nz, nx) * 3.22 * 10^-3 * 10^15 # density in kg/km3
    medium_ref = bundle_medium(μref, ρref)
end;

# ╔═╡ ad21da29-f6ff-4a94-be87-4e88640cddbf
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	# e.g., get_vs(mean, medium)
	get_vs(op, medium) = op(sqrt.(medium.μ ./ medium.ρ))
	get_vs(medium) = sqrt.(medium.μ ./ medium.ρ)
end
  ╠═╡ =#

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

# ╔═╡ 6937b103-9ce2-4189-8129-aae1e7936d4f
#=╠═╡
length(tgrid)
  ╠═╡ =#

# ╔═╡ be94d0e7-e9a8-4021-a94c-a5b70bbcede7
#=╠═╡
length(tgrid)
  ╠═╡ =#

# ╔═╡ ae8012be-e7ab-4e85-a27f-febf08b3380b
md"### Absorbing Boundaries"

# ╔═╡ ce89710a-48d8-46d4-83b4-163a7eb0a2e5
function get_taper_array(nx, nz; np=50, tapfact=0.20)
    tarray = ones(nz, nx)
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
    f = (; vy=zeros(nz, nx),
        dvydx=zeros(nz, nx),
        dvydz=zeros(nz, nx),
        σyx=zeros(nz, nx),
        dσyxdx=zeros(nz, nx),
        σyz=zeros(nz, nx),
        dσyzdz=zeros(nz, nx))
    if snap_store
        fs = (; vys=[zeros(nz, nx) for i in 1:nt], σyxs=[zeros(nz, nx) for i in 1:nt], σyzs=[zeros(nz, nx) for i in 1:nt])
        return merge(f, fs)
    else
        return f
    end
end

# ╔═╡ 31d742a4-100f-4744-afe5-381b265b6f4c
function reset_fields!(fields)
    fill!(fields.vy, zero(Float64))
    fill!(fields.dvydx, zero(Float64))
    fill!(fields.dvydz, zero(Float64))
    fill!(fields.σyx, zero(Float64))
    fill!(fields.σyz, zero(Float64))
    fill!(fields.dσyxdx, zero(Float64))
    fill!(fields.dσyzdz, zero(Float64))
end

# ╔═╡ 15bbf544-34bd-4d38-bac5-0f43b1305df3
function propagate!(data, fields, pa, medium, forcing_transform_mat, forcing, adj=false)
    reset_fields!(fields)

    (; vy, dvydx, dvydz, σyx, σyz, dσyxdx, dσyzdz) = fields
    (; nx, nz, tarray, tgrid, dt, nt) = pa


    (; μ, ρ, invρ) = medium

    if adj
        (; Rs, Rr, Rr_t) = forcing_transform_mat
    else
        (; Rs, Rr) = forcing_transform_mat
    end

    # time loop
    # @progress
    @progress for it = 1:nt
        Dx!(dvydx, vy)
        Dz!(dvydz, vy)


        @. σyx = σyx + μ * dvydx * dt
        @. σyx = σyx * pa.tarray
        @. σyz = σyz + μ * dvydz * dt
        @. σyz = σyz * pa.tarray


        Dx!(dσyxdx, σyx)
        Dz!(dσyzdz, σyz)

        @. vy = vy + invρ * (dσyxdx + dσyzdz) * dt

        # need to view vy as a vector for source/recording operations
        vyv = view(vy, :)

        # add body force
        f = forcing[it]
        if adj
            mul!(vyv, Rr_t, f, dt, 1.0)
        else
            mul!(vyv, Rs, f, dt, 1.0)
        end

        # record data
        d = data[it]
        mul!(d, Rr, vyv)

        (:vys ∈ keys(fields)) && copyto!(fields.vys[it], vy)
        (:σyxs ∈ keys(fields)) && copyto!(fields.σyxs[it], σyx)
        (:σyzs ∈ keys(fields)) && copyto!(fields.σyzs[it], σyz)

    end

    return nothing
end

# ╔═╡ bd8b9a0c-6664-4ba8-9795-0a496e932d85
md"### Gradients"

# ╔═╡ 82116417-7a0e-4363-9b53-2bd5df250f17
function reset_grad!(grad)
    fill!(grad.▽ρ, zero(Float64))
    fill!(grad.▽μ, zero(Float64))
end

# ╔═╡ 7d00ecf7-80ba-4eee-99f6-b03f7b1f6e52
function initialize_grad(pa, nt, snap_store=true)
    nz, nx = pa.nz, pa.nx
    f = (; ▽ρ=zeros(nz, nx),
        ▽μ=zeros(nz, nx),
        ▽ρs=zeros(nz, nx),
        ▽μs=zeros(nz, nx))
    if snap_store
        fs = (; ▽ρs=[zeros(nz, nx) for i in 1:nt], ▽μs=[zeros(nz, nx) for i in 1:nt])
        return merge(f, fs)
    else
        return f
    end
end

# ╔═╡ 4e1a0d4b-5f25-4b25-8dbb-4069e38dc5c4
# create a function to compute grad_phi and grad_mu 
function propagate_gradients(grad, forwfields, adjfields, pa)
    reset_grad!(grad)

    (; ▽ρ, ▽μ) = grad
    (; nx, nz, tarray, tgrid, dt, nt) = pa

    @progress for it = 1:nt-1
        # for gradρ
        v2 = forwfields.vys[it+1]
        v3 = forwfields.vys[it]
        u2 = adjfields.vys[nt-it]
        @. ▽ρ = ▽ρ + (v2 - v3) * u2
        @. ▽ρ = ▽ρ * pa.tarray

        # for gradμ

        σx2 = forwfields.σyxs[it]
        τx2 = adjfields.σyxs[nt-it]
        τx3 = adjfields.σyxs[nt-it+1]

        σz2 = forwfields.σyzs[it]
        τz2 = adjfields.σyzs[nt-it]
        τz3 = adjfields.σyzs[nt-it+1]

        @. ▽μ = ▽μ + (σx2 * (τx2 - τx3)) + (σz2 * (τz2 - τz3))
        @. ▽μ = ▽μ * pa.tarray

        (:▽ρs ∈ keys(grad)) && copyto!(grad.▽ρs[it], ▽ρ)
        (:▽μs ∈ keys(grad)) && copyto!(grad.▽μs[it], ▽μ)
    end
end

# ╔═╡ f95a08ce-a38d-4b7f-b478-4dbfa607740e
md"### Wavelets"

# ╔═╡ 90953c64-8a87-4065-8c34-d0ead540b728
md"""
Generate a Ricker Wavelet. Reference:
Frequencies of the Ricker wavelet, Yanghua Wang, GEOPHYSICS, VOL. 80, NO. 2.
Its bandwidth is roughly 1.2 * fpeak.

* `fqdom::Float64`: dominant frequency 
* `tgrid`: time-domain grid
* `tpeak::Float64=tgrid[1]+1.5/fqdom`: the peak of the ricker in time (has a default)
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
    return wav
end

# ╔═╡ 17dd3d57-d5ca-443c-b003-b3a97b963d57
#=╠═╡
begin
	source_fpeak = 2.0 # in Hz
	source_wavelet = ricker(source_fpeak, tgrid)
	source_forcing = [fill(source_wavelet[it], ageom.ns) for it in 1:nt]
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

# ╔═╡ c34b1a5d-5078-4b8f-94d1-a088cbe5ab3e
#=╠═╡
plot(heatmap(x=xgrid, y=zgrid, z=(grid_param.tarray)))
  ╠═╡ =#

# ╔═╡ d812711d-d02f-44bb-9e73-accd1623dea1
#=╠═╡
plot(tgrid, source_wavelet, size=(500, 200), w=2, label="Source Wavelet", Layout(width=500, height=250))
  ╠═╡ =#

# ╔═╡ 77c9696c-58c5-40bf-acd0-16d5cf877810
#=╠═╡
begin
	# Getting the forcing due to sources
	true_forcing_transform = get_forcing_transform(ageom, medium_true, grid_param.xgrid, grid_param.zgrid)

	# Initialisation of fields and data
	fields_true = initialize_fields(grid_param, grid_param.nt)
	dobs = initialize_data(grid_param, ageom)

	# Running the simulation to generate observed data
	@time propagate!(dobs, fields_true, grid_param, medium_true, true_forcing_transform, source_forcing)	
end;
  ╠═╡ =#

# ╔═╡ 3be62716-f2d9-434c-a69a-ed272b89c85d
#=╠═╡
begin
	# Getting the forcing due to sources
	ref_forcing_transform = get_forcing_transform(ageom, medium_ref, grid_param.xgrid, grid_param.zgrid, true)

	# Initialisation of fields and data
	fields_forw = initialize_fields(grid_param, grid_param.nt, snap_store=true)
	dref = initialize_data(grid_param, ageom)

	# Simulation to compute wavefields
	propagate!(dref, fields_forw, grid_param, medium_ref, ref_forcing_transform, source_forcing)
end;
  ╠═╡ =#

# ╔═╡ 59155ad5-d341-4e16-b7cc-b6a3def51992
#=╠═╡
data_error = dref .- dobs;
  ╠═╡ =#

# ╔═╡ 3f5f9d8a-3647-4a16-89ba-bd7a31c01064
#=╠═╡
begin
	# Forcing at receivers due to data error
	rec_forcing = reverse(data_error)

	# Initialisation of fields and data
	fields_adj = initialize_fields(grid_param, nt, snap_store=true)
	dadj = initialize_data(grid_param, ageom)

	# Simulating the adjoint field
	propagate!(dadj, fields_adj, grid_param, medium_ref, ref_forcing_transform, rec_forcing, true)
end;
  ╠═╡ =#

# ╔═╡ 2df3bd5b-630e-451a-a207-2c3a1719916e
#=╠═╡
begin
	# Initialisation of grad
	fields_grad = initialize_grad(grid_param, grid_param.nt)
	# Simulation to compute grad 
	@time propagate_gradients(fields_grad, fields_forw, fields_adj, grid_param)
end
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
fieldheat([fields_forw.vys[t], fields_adj.vys[nt-t], fields_grad.▽ρs[t], fields_grad.▽μs[t]], 
	["Forward Field" "Adjoint Field" "Gradient w.r.t. ρ" "Gradient w.r.t. μ"], 
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
    data = cat(d..., dims=2)'
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
    layout = Layout(title=title, width=300, height=450, xaxis_range=[0, size(data, 2) + 1], yaxis_scaleanchor="x", yaxis_scaleratio=2.5, yaxis_autorange="reversed",
        xaxis_title="Receiver #",
        yaxis_title="Time (s)",)
    return plot(trace, layout)
end

# ╔═╡ 4171af00-1d14-45ba-9fd3-a2c30d0b759f
#=╠═╡
dataheat(dobs, tgrid, title="Observed Data")
  ╠═╡ =#

# ╔═╡ 270e5d4a-c666-43e7-8c64-02b9fed977e4
#=╠═╡
dataheat(data_error, tgrid, title="Data Error")
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
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DSP = "~0.7.8"
FFTW = "~1.6.0"
LaTeXStrings = "~1.3.0"
PlutoPlotly = "~0.3.6"
PlutoUI = "~0.7.50"
ProgressLogging = "~0.1.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "1c637fc7dc0b02373e9e75bc0030e1f5a3b10ba3"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "16b6dbc4cf7caee4e1e75c49485ec67b667098a0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.3.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

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

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "da8b06f89fce9996443010ef92572b193f8dca1f"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.8"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

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

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

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

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

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
version = "2.28.0+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

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
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Colors", "Dates", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "PlotlyBase", "PlutoUI", "Reexport"]
git-tree-sha1 = "dec81dcd52748ffc59ce3582e709414ff78d947f"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.3.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.Polynomials]]
deps = ["ChainRulesCore", "LinearAlgebra", "MakieCore", "RecipesBase"]
git-tree-sha1 = "86efc6f761df655f8782f50628e45e01a457d5a2"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.8"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

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
# ╠═3c540889-49dc-415c-acbc-3494897b260c
# ╟─9c32f5bc-f6d1-4048-903a-27224aaa1f40
# ╟─729689f3-7b76-4701-a60d-3471e1f20fce
# ╟─77e134d8-bd8b-4303-8c44-a4920cf0ee81
# ╟─4adbb7f0-7927-470d-8c00-07d3b3c0cd78
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
# ╟─881c7368-57df-469a-96ec-821512cf98e0
# ╟─65efba3b-16b0-4113-a59e-809365e7bdd6
# ╠═1c67cda8-7712-4d5b-a2aa-af47f290f745
# ╟─55c7f981-96a7-40e9-811f-37334622565b
# ╟─5b271e5f-879c-4c43-825a-9660f322febd
# ╠═6a8139c4-12c1-4d18-bd1e-14334290aec1
# ╠═dd78c460-a446-46f5-af57-7e859383348d
# ╠═ebab6005-2ad6-4057-9275-bf7d53d41b0b
# ╟─2855c8cf-8364-4c6c-a122-781b99440e89
# ╠═17dd3d57-d5ca-443c-b003-b3a97b963d57
# ╠═be94d0e7-e9a8-4021-a94c-a5b70bbcede7
# ╠═d812711d-d02f-44bb-9e73-accd1623dea1
# ╟─e233afec-6049-4277-8be9-95687c4589b5
# ╠═77c9696c-58c5-40bf-acd0-16d5cf877810
# ╠═4171af00-1d14-45ba-9fd3-a2c30d0b759f
# ╟─22b8db91-73a0-46df-87fd-7cf0b66ee37d
# ╠═8b3776bd-509b-4232-9737-36c9ae003350
# ╠═f26f1b1a-18e0-413c-86a9-351ba5dfaebf
# ╠═3be62716-f2d9-434c-a69a-ed272b89c85d
# ╟─c73f69d0-69e6-47f1-97b0-3e81218776e6
# ╠═59155ad5-d341-4e16-b7cc-b6a3def51992
# ╠═270e5d4a-c666-43e7-8c64-02b9fed977e4
# ╟─a3a9deea-e2d6-4d58-90d7-5a54be176289
# ╠═3f5f9d8a-3647-4a16-89ba-bd7a31c01064
# ╟─afca5e37-ebc0-45e7-b27b-53602dbf6672
# ╟─66f9c698-61e3-4b61-aff3-dfc67eb2f6af
# ╠═2df3bd5b-630e-451a-a207-2c3a1719916e
# ╟─8d161f09-8339-4277-8739-ff76607f7abf
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
# ╠═e08bf013-00c7-4870-82d8-19b899e7208d
# ╠═dd626606-2a4d-494d-ab74-92753072c773
# ╟─ab8b1a22-ca7a-409e-832e-8d5d08a29a1e
# ╠═f4d91971-f806-4c5c-8548-b58a20acfb2c
# ╟─9bc38d55-285b-4b83-98d9-d7f9e03405d1
# ╠═27844886-0b54-4b08-a592-a1a38e4b0be2
# ╠═ad21da29-f6ff-4a94-be87-4e88640cddbf
# ╟─ae8012be-e7ab-4e85-a27f-febf08b3380b
# ╠═ce89710a-48d8-46d4-83b4-163a7eb0a2e5
# ╠═b7f4078a-ead0-4d42-8b44-4f471eefc6fc
# ╠═c34b1a5d-5078-4b8f-94d1-a088cbe5ab3e
# ╟─9494e2a6-e2fd-4728-94d4-d68816a00e72
# ╠═9b744e72-07a2-4f8c-b88c-e0b1131794d0
# ╠═31d742a4-100f-4744-afe5-381b265b6f4c
# ╠═15bbf544-34bd-4d38-bac5-0f43b1305df3
# ╟─bd8b9a0c-6664-4ba8-9795-0a496e932d85
# ╠═82116417-7a0e-4363-9b53-2bd5df250f17
# ╠═7d00ecf7-80ba-4eee-99f6-b03f7b1f6e52
# ╠═4e1a0d4b-5f25-4b25-8dbb-4069e38dc5c4
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
