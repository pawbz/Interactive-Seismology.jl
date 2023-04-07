### A Pluto.jl notebook ###
# v0.19.14

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
    using Plots
    using PlutoUI
    using Statistics
    using ProgressLogging
    using SparseArrays
end

# ╔═╡ 3c540889-49dc-415c-acbc-3494897b260c
TableOfContents()

# ╔═╡ 9c32f5bc-f6d1-4048-903a-27224aaa1f40
md"""
# Full Waveform Inversion (SH Waves)

Pseudo-spectral methods are used in scientific computing for the solution of partial differential equations. The advantage of these methods is that the evaluation of spatial-derivative operators is considerably faster as we use algorithms such as the [Fast Fourier transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform).

The purpose of this notebook is to introduce these methods for solving the seismic-wave equation, specifically, a 2-D SH stress-velocity system. As discussed in the class, the particle velocity and stress fields are staggered in time and we use a simple explicit time-stepping approach.
We shall assume a 2-D periodic domain and define methods to compute FFT, spatial derivative and iFFT. More efficient pseudo-spectral implementations (especially using GPUs and multiple CPU threads) are in the following packages:
[Fourier Flows](https://fourierflows.github.io/FourierFlowsDocumentation/stable/),
[Geophysical Flows](https://github.com/FourierFlows/GeophysicalFlows.jl).
Importantly, equations are oftentimes time-stepped forward in Fourier space, which is faster, however, we avoid that in this notebook to facilitate storage of snapshots at arbit values of the time.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

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


# ╔═╡ fe821a3c-f528-4773-bcf9-f71513a1eace
md"""
Let us set up the parameters of the simulation. We define a 2-D spatial grid of extents $x∈[-40,40]$ and $z∈[-40,40]$. Let the density (ρ) and shear modulus (μ) be equal to 1 at all the grid points. Let us also add a density perturbation i.e a high density region in the bottom portion of the grid. We then apply an initial impulsive velocity ($v_y$) that has a 2-D asymmetric gaussian profile. In order to study the propagation of the wave-front, we need to choose a time step that satisfies the Courant condition. Here, we have chosen the Courant number to be $0.2$.
"""

# ╔═╡ 122c2aa4-baec-4288-ab53-afa1d977c486
md"## Spatial Grids"

# ╔═╡ 7928a88c-f217-4338-a6a5-50ab2d422480
begin
    # domain extends
    zgrid = range(0, stop=40, length=256)
    xgrid = range(-40, stop=40, length=512)

    # Numerics
    nx, nz = length(xgrid), length(zgrid)
    nout = 1  # plotting frequency

    # Derived numerics
    dx, dz = step(xgrid), step(zgrid) # cell sizes
end;

# ╔═╡ 881c7368-57df-469a-96ec-821512cf98e0
md"## True Medium"

# ╔═╡ 5b271e5f-879c-4c43-825a-9660f322febd
md"## Time Grid"

# ╔═╡ a5db1d99-24e5-41d1-a804-c5dd1d6db94c
# lets calculate the min distance from the center to the edge of the domain
r = min(xgrid[end] - xgrid[1], zgrid[end] - zgrid[1]) * 0.5

# ╔═╡ 45bd4e43-b703-4bbc-8dfa-07bb835e4117
courant_number = 0.2

# ╔═╡ cb83dbd1-c423-4b27-b29e-7dc8051f43d5
md"""
## Acquistion Geometry
"""

# ╔═╡ 2855c8cf-8364-4c6c-a122-781b99440e89
md"## Body Forces"

# ╔═╡ 12db5316-90e4-4582-9303-7ff362c32cea
source_fpeak = 1.0 # in Hz

# ╔═╡ e233afec-6049-4277-8be9-95687c4589b5
md"""
## Generate Observed Data
"""

# ╔═╡ 24525b74-e0d9-484c-9493-d3da03a3ffe4


# ╔═╡ 22b8db91-73a0-46df-87fd-7cf0b66ee37d
md"""
## Reference Medium
This will also serve as an initial model during inversion.
"""

# ╔═╡ 6ba91b28-ee35-4d6f-ad10-dccf991dea83


# ╔═╡ c73f69d0-69e6-47f1-97b0-3e81218776e6
md"""
## Data Error
"""

# ╔═╡ a3a9deea-e2d6-4d58-90d7-5a54be176289
md"## Adjoint Simulation"

# ╔═╡ e41d1801-99f9-4c63-9bb2-96381f114c70
md"
## Gradρ
"

# ╔═╡ 54242145-429f-4480-8e5b-e3b51e347097
function gradρ(fields_forw, fields_adj,pa)
	(; nx, nz, tarray, tgrid, dt, nt)=pa
	g=zeros(nz,nx)
	for it in 1:nt-1
		v1=fields_forw.vys[it+1]
		v2=fields_forw.vys[it]
		v3=fields_adj.vys[nt-it]
		@. g = g + (v2 - v1)* v3
	end
	return g
end


# ╔═╡ b3c89340-6795-492b-9eee-5b06266818f0
md"
## Gradμ
"

# ╔═╡ fc49a6d7-a1b1-458a-a9ad-e120282bbabc
md"""
## Appendix
"""

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
    F = plan_rfft(zeros(nz, nx), (1, 2))
    kx = reshape(collect(fftfreq(nx, inv(step(xgrid)))), 1, :) * 2 * pi
    kz = reshape(collect(rfftfreq(nz, inv(step(zgrid)))), :, 1) * 2 * pi
    storage = zero(F * zeros(nz, nx))
    fp = (; F, kx, kz, storage)
    function Dx!(dPdx, P, fp)
        mul!(fp.storage, fp.F, P)
        broadcast!(*, fp.storage, fp.storage, fp.kx)
        rmul!(fp.storage, im)
        ldiv!(dPdx, fp.F, fp.storage)
    end
    Dx!(dP, P) = Dx!(dP, P, fp)
    Dx(P) = (dPdx = zero(P); Dx!(dPdx, P, fp); dPdx)
    function Dz!(dPdz, P, fp)
        mul!(storage, fp.F, P)
        broadcast!(*, fp.storage, fp.storage, fp.kz)
        rmul!(fp.storage, im)
        ldiv!(dPdz, fp.F, fp.storage)
    end
    Dz!(dP, P) = Dz!(dP, P, fp)
    Dz(P) = (dPdz = zero(P); Dz!(dPdz, P, fp); dPdz)
    nothing
end

# ╔═╡ 93090680-2380-412d-9752-df18251c7dbf
md"""
### Propagate
"""

# ╔═╡ 9bc38d55-285b-4b83-98d9-d7f9e03405d1
md"### Medium"

# ╔═╡ 27844886-0b54-4b08-a592-a1a38e4b0be2
function bundle_medium(μ, ρ)
    return (; μ=μ, ρ=ρ, invρ=inv.(ρ))
end

# ╔═╡ 1c67cda8-7712-4d5b-a2aa-af47f290f745
begin
    # lets add density perturbation
    μtrue = ones(nz, nx) .* 82 * 10^9 * 10^3
    ρtrue = ones(nz, nx) .* 3.22 * 10^-3 * 10^15 # density in kg/km3
    ρtrue[120:256, :] .= 5 * 10^-3 * 10^15  # density in kg/km3
    medium_true = bundle_medium(μtrue, ρtrue)
end;

# ╔═╡ b4085356-62a1-4388-96ae-8b2a8fd7de0f
begin
    # choose time stepping dt to satisfy Courant condition
    dt = courant_number * step(xgrid) * minimum(inv.(sqrt.(medium_true.μ ./ medium_true.ρ)))
    nt = Int(floor(r / (mean(sqrt.(medium_true.μ ./ medium_true.ρ)) * dt)))
    tgrid = range(0, length=nt, step=dt)
    nothing
end

# ╔═╡ 9a77180b-0bfa-4401-af30-d95f14e10d2c
@bind t_forw Slider(range(1,nt),show_value=true)

# ╔═╡ 8b3776bd-509b-4232-9737-36c9ae003350
begin
    # Physics of medium
    # Unperturbed
    μref = ones(nz, nx) * 82 * 10^9 * 10^3      # shear modulus in kg / km/s^2
    ρref = ones(nz, nx) *  3.22 * 10^-3 * 10^15 # density in kg/km3
	initial_medium = bundle_medium(μref, ρref)
end;

# ╔═╡ d880f005-381a-4c72-a9a1-ae3eb84a90eb
vsmean(medium) = mean(sqrt.(medium.μ ./ medium.ρ))

# ╔═╡ 8bc93109-4f3d-43c8-9fed-e30d90e54068
dx, dz, vsmean(medium_true) / source_fpeak

# ╔═╡ e8333b23-53c3-445e-9ca3-6b278359f8ab
md"### Acquisition"

# ╔═╡ 9248af7f-dc1a-4bf6-8f3f-304db73fc604
function get_ageom(xgrid, zgrid, ns, nr; zs=quantile(zgrid, 0.25), zr=quantile(zgrid, 0.25))
    A = (; ns, nr, zs=fill(zs, ns), zr=fill(zr, nr),
        xr=(nr == 1) ? quantile(xgrid, 0.25) : range(quantile(xgrid, 0.25), stop=quantile(xgrid, 0.75), length=nr),
        xs=(ns == 1) ? quantile(xgrid, 0.25) : range(quantile(xgrid, 0.25), stop=quantile(xgrid, 0.75), length=ns)
    )
    return A
end;

# ╔═╡ d39753e2-5986-4394-9293-9e394f2807f0
ageom = get_ageom(xgrid, zgrid, 1, 10);

# ╔═╡ e7b65566-a79e-4101-9471-3656a92e95e6
function get_restriction_matrix(xpos, zpos, xgrid, zgrid, transpose_flag=false)
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
    V = ones(n)
    return transpose_flag ? sparse(J, I, V, n, N) : sparse(I, J, V, N, n)
end

# ╔═╡ e08bf013-00c7-4870-82d8-19b899e7208d
function get_restriction_matrix_density(medium, xpos, zpos, xgrid, zgrid, transpose_flag=false)
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
    V = medium_true.invρ[I]
    return transpose_flag ? sparse(J, I, V, n, N) : sparse(I, J, V, N, n)
end

# ╔═╡ dd626606-2a4d-494d-ab74-92753072c773
function get_forcing_transform(ageom,medium,xgrid,zgrid,medium_flag=false)
	if medium_flag
		return (; Rs=get_restriction_matrix_density(medium,ageom.xs, ageom.zs, xgrid, zgrid),
        Rr=get_restriction_matrix(ageom.xr, ageom.zr, xgrid, zgrid, true), Rr_t=get_restriction_matrix_density(medium,ageom.xr, ageom.zr, xgrid, zgrid))
	else
		return (; Rs=get_restriction_matrix_density(medium,ageom.xs, ageom.zs, xgrid, zgrid), Rr=get_restriction_matrix(ageom.xr, ageom.zr, xgrid, zgrid, true))
	end
end;

# ╔═╡ 3fd67421-1c00-4632-8ca4-181dd1656355
true_forcing_transform = get_forcing_transform(ageom,medium_true,xgrid,zgrid);

# ╔═╡ cfb83df9-1e43-418f-870e-806a57c46407
ref_forcing_transform = get_forcing_transform(ageom,initial_medium,xgrid,zgrid,true);

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
    for k in [:vy, :dvydx, :dvydz, :σyx, :σyz, :dσyxdx, :dσyzdz]
        fill!(getfield(fields, k), 0.0)
    end
end

# ╔═╡ 15bbf544-34bd-4d38-bac5-0f43b1305df3
function propagate!(data, fields, pa, medium, forcing_transform_mat, forcing, adj=false)
    reset_fields!(fields)
    vy, dvydx, dvydz, σyx, σyz, dσyxdx, dσyzdz = fields
    (; nx, nz, tarray, tgrid) = pa
    dt = step(tgrid)
    nt = length(tgrid)
	f_grid = zeros(nz*nx)

    (; μ, ρ, invρ) = medium

	if adj
		(; Rs, Rr, Rr_t) = forcing_transform_mat
	else
		(; Rs, Rr) = forcing_transform_mat
	end

    # time loop
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

# ╔═╡ ab8b1a22-ca7a-409e-832e-8d5d08a29a1e
md"### Data"

# ╔═╡ f4d91971-f806-4c5c-8548-b58a20acfb2c
function initialize_data(grid_param, ageom)
    [zeros(length(ageom.xr)) for t in grid_param.tgrid]
end

# ╔═╡ fbe44944-499a-4881-94b6-07855d1165aa
md"""
### Plots
"""

# ╔═╡ a5cd8e7a-380f-4203-a856-f9e56e04b092
# plot vs and density heatmaps of a given medium
function plot_medium(medium,grid)
	(; μ, ρ) = medium
	vel = sqrt.(μ./ρ)

	fig1 = heatmap(grid.xgrid,grid.zgrid,ρ*1e-12,c=:seismic,xlabel="Distance (x km)",ylabel="Depth (z km)",title="Density (g/cc)",yflip=true)
	fig2 = heatmap(grid.xgrid,grid.zgrid,vel,c=:seismic,xlabel="Distance (x km)",ylabel="Depth (z km)",title="Velocity (km/s)",yflip=true)
	return plot(fig1,fig2,size=(1000,600))
end;

# ╔═╡ 43a77919-d880-40c9-95c9-aaa429a65fb7
begin
    @userplot SeisHeat

    @recipe function f(h::SeisHeat)
        grid := true
        xlabel := "Distance (x)"
        ylabel := "Depth (z)"
        title --> "Wavefield"
        aspect_ratio := 1
        yflip := true
        c := :seismic
        seriestype := :heatmap
        @series begin
            h.args
        end
    end
end

# ╔═╡ 36408698-f85c-4529-972f-6389534fcb88
myheat(x, t="") = heatmap(xgrid, zgrid, x, c=:grays, aspect_ratio=1, title=t, ylim=extrema(zgrid))

# ╔═╡ a62839d5-837b-4c37-996f-33659c34911c
md"### UI"

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
function ricker(fqdom::Float64,
    tgrid::StepRangeLen;
    # tpeak is the location of the peak amplitude
    tpeak::Float64=tgrid[1] + 1.5 / fqdom, # using approximate half width of ricker
    trim_tol::Float64=0.0,
    maxamp::Float64=1.0
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

# ╔═╡ a22ab168-ee9c-4433-b890-b46e9849f127
source_wavelet = ricker(source_fpeak, tgrid)

# ╔═╡ d812711d-d02f-44bb-9e73-accd1623dea1
plot(tgrid, source_wavelet, size=(500, 200), w=2, label="Source Wavelet")

# ╔═╡ c3c0aabd-9df4-41af-9512-bdbb7d4f8700
source_forcing = [fill(source_wavelet[it], 1) for it in 1:nt]

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

# ╔═╡ a120a929-d989-4b88-86af-e735a577db18
# a NamedTuple for grid-related parameters
grid_param = (; xgrid, zgrid, tgrid, nx, nz, tarray=get_taper_array(nx, nz, tapfact=0.3,))

# ╔═╡ 55c7f981-96a7-40e9-811f-37334622565b
plot_medium(medium_true,grid_param)

# ╔═╡ 62a18fa1-b4d6-4fde-b4a7-09fad6b16a22
fields_true = initialize_fields(grid_param, nt);

# ╔═╡ ba4c4fd5-56d0-4004-8858-377139cc3d3c
dobs = initialize_data(grid_param, ageom);

# ╔═╡ 6ece24cb-adc1-4aea-93c2-cc74167b26f6
heatmap(tgrid,range(1,ageom.nr),cat(dobs..., dims=2),xlabel="Time(s)",ylabel="Receiver",title="Observed data",legend=:none)

# ╔═╡ 7b84af8f-939a-47ed-8a0b-5721ce026d79
@time propagate!(dobs, fields_true, grid_param, medium_true, true_forcing_transform, source_forcing)

# ╔═╡ 11ca5d78-1434-41ce-9e3b-88ed63474128
plot_medium(initial_medium,grid_param)

# ╔═╡ af1389bb-b13a-496d-914d-fc99a0b904c9
fields_forw = initialize_fields(grid_param, nt, snap_store=true);

# ╔═╡ 804e23ee-d7bc-4ed9-8b46-ad462442bccc
dref = initialize_data(grid_param, ageom);

# ╔═╡ f2e77de4-f1a7-4cb8-a6b7-f9db0d807720
heatmap(tgrid,range(1,ageom.nr),cat(dref..., dims=2),xlabel="Time(s)",ylabel="Receiver",title="Modelled data",legend=:none)

# ╔═╡ 59155ad5-d341-4e16-b7cc-b6a3def51992
data_error = dref .- dobs;

# ╔═╡ bdd9e0d1-9f99-4957-85fa-ef9e258637d5
heatmap(tgrid,range(1,ageom.nr),cat(data_error..., dims=2),xlabel="Time(s)",ylabel="Receiver",title="Modelled - Observed",legend=:none)

# ╔═╡ 1b6e41ed-75a0-4b12-b390-78538566ce43
rec_forcing = reverse(dref - dobs);

# ╔═╡ b075b29b-0be2-467e-9e27-0280fad35520
@time propagate!(dref, fields_forw, grid_param, initial_medium, ref_forcing_transform, source_forcing)

# ╔═╡ e32d0418-6e18-485b-a6ea-0a982d63e1d2
fields_adj = initialize_fields(grid_param, nt, snap_store=true);

# ╔═╡ 4e5a992d-303a-431c-a2c2-e5f3214a1da6
dadj = initialize_data(grid_param, ageom);

# ╔═╡ f129bfac-b0d5-482b-98e5-856c26c0fc3f
@time propagate!(dadj, fields_adj, grid_param, initial_medium, ref_forcing_transform, rec_forcing, true)

# ╔═╡ d614dea3-e557-4374-9be0-b52321206f10
pa=(; nx, nz, tarray=get_taper_array(nx, nz, tapfact=0.3,) , tgrid, dt, nt)

# ╔═╡ 4b854f3c-e957-4a6e-8fff-8fee52d47691
heatmap(xgrid,zgrid,gradρ(fields_forw, fields_adj, pa ), yflip=true ,xlabel="Distance (x km)", ylabel="Depth (z km)", title="Wavefield" , c=:seismic )

# ╔═╡ 283da509-5437-4405-8aaf-19a56c64b037
function gradμ(fields_forw, fields_adj)
	 (; nx, nz, tarray, tgrid, dt, nt)=pa
	h=zeros(nz,nx)
	for it in 1:nt-1
		v1=fields_forw.dvydx[it]
		v2=fields_adj.σyx[nt-it]
		@. h = h - v1*v2
	end
	return h
end

# ╔═╡ 5e598c38-2fe1-475c-b836-82519af69efe
heatmap(xgrid,zgrid,gradμ(fields_forw, fields_adj), yflip=true ,xlabel="Distance (x km)", ylabel="Depth (z km)", title="" , c=:seismic )

# ╔═╡ c34b1a5d-5078-4b8f-94d1-a088cbe5ab3e
heatmap(xgrid,zgrid,grid_param.tarray,yflip=true)

# ╔═╡ b7f4078a-ead0-4d42-8b44-4f471eefc6fc
function clip_edges(mat,grid_param)
	dum = mat
	dum[grid_param.tarray .!= 1] .= NaN
	return(dum)
end;

# ╔═╡ ae2391be-a7c7-4612-a58b-909c6d5eac0d
begin
	figv = heatmap(xgrid,zgrid,clip_edges(fields_forw.vys[t_forw[1]],grid_param),yflip=true,xlabel="Distance (x km)", ylabel="Depth (z km)", title="Velocity field",legend=:none)
	figadj = heatmap(xgrid,zgrid,clip_edges(fields_adj.vys[nt-t_forw[1]+1],grid_param),yflip=true,xlabel="Distance (x km)", ylabel="Depth (z km)", title="Adjoint field",legend=:none)
	plot(figv,figadj,size=(700,500))
end

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
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
FFTW = "~1.6.0"
LaTeXStrings = "~1.3.0"
Plots = "~1.25.8"
PlutoUI = "~0.7.39"
ProgressLogging = "~0.1.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "dfbed4f959a5021850f10211b42df764c189fbab"

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

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

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
version = "0.5.2+0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

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

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

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

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c98aea696662d09e215ef7cda5296024a9646c75"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.4"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bc9f7725571ddb4ab2c4bc74fa397c1c5ad08943"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.69.1+0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "e07a1b98ed72e3cdd02c6ceaab94b8a606faca40"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.2.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "303202358e38d2b01ba46844b92e48a3c238fd9e"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.6"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

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

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

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

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

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

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

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

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "c95373e73290cf50a8a22c3375e4625ded5c5280"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "d16070abde61120e01b4f30f6f398496582301d6"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.12"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

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

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

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

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

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
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b8d897fe7fa688e93aef573711cb207c08c9e11e"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.19"

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

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

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
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

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
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6edfe154ad7b313c01aceca188c05c835c67360"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.4+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═3c540889-49dc-415c-acbc-3494897b260c
# ╟─9c32f5bc-f6d1-4048-903a-27224aaa1f40
# ╟─4adbb7f0-7927-470d-8c00-07d3b3c0cd78
# ╟─fe821a3c-f528-4773-bcf9-f71513a1eace
# ╟─122c2aa4-baec-4288-ab53-afa1d977c486
# ╠═7928a88c-f217-4338-a6a5-50ab2d422480
# ╠═8bc93109-4f3d-43c8-9fed-e30d90e54068
# ╟─881c7368-57df-469a-96ec-821512cf98e0
# ╠═1c67cda8-7712-4d5b-a2aa-af47f290f745
# ╠═55c7f981-96a7-40e9-811f-37334622565b
# ╟─5b271e5f-879c-4c43-825a-9660f322febd
# ╠═a5db1d99-24e5-41d1-a804-c5dd1d6db94c
# ╠═45bd4e43-b703-4bbc-8dfa-07bb835e4117
# ╠═b4085356-62a1-4388-96ae-8b2a8fd7de0f
# ╠═a120a929-d989-4b88-86af-e735a577db18
# ╠═d614dea3-e557-4374-9be0-b52321206f10
# ╟─cb83dbd1-c423-4b27-b29e-7dc8051f43d5
# ╠═d39753e2-5986-4394-9293-9e394f2807f0
# ╟─2855c8cf-8364-4c6c-a122-781b99440e89
# ╠═12db5316-90e4-4582-9303-7ff362c32cea
# ╠═a22ab168-ee9c-4433-b890-b46e9849f127
# ╠═d812711d-d02f-44bb-9e73-accd1623dea1
# ╠═c3c0aabd-9df4-41af-9512-bdbb7d4f8700
# ╟─e233afec-6049-4277-8be9-95687c4589b5
# ╠═3fd67421-1c00-4632-8ca4-181dd1656355
# ╠═62a18fa1-b4d6-4fde-b4a7-09fad6b16a22
# ╠═ba4c4fd5-56d0-4004-8858-377139cc3d3c
# ╠═7b84af8f-939a-47ed-8a0b-5721ce026d79
# ╠═6ece24cb-adc1-4aea-93c2-cc74167b26f6
# ╠═24525b74-e0d9-484c-9493-d3da03a3ffe4
# ╟─22b8db91-73a0-46df-87fd-7cf0b66ee37d
# ╠═8b3776bd-509b-4232-9737-36c9ae003350
# ╠═11ca5d78-1434-41ce-9e3b-88ed63474128
# ╠═cfb83df9-1e43-418f-870e-806a57c46407
# ╠═af1389bb-b13a-496d-914d-fc99a0b904c9
# ╠═804e23ee-d7bc-4ed9-8b46-ad462442bccc
# ╠═b075b29b-0be2-467e-9e27-0280fad35520
# ╠═6ba91b28-ee35-4d6f-ad10-dccf991dea83
# ╠═f2e77de4-f1a7-4cb8-a6b7-f9db0d807720
# ╠═9a77180b-0bfa-4401-af30-d95f14e10d2c
# ╠═ae2391be-a7c7-4612-a58b-909c6d5eac0d
# ╟─c73f69d0-69e6-47f1-97b0-3e81218776e6
# ╠═59155ad5-d341-4e16-b7cc-b6a3def51992
# ╠═bdd9e0d1-9f99-4957-85fa-ef9e258637d5
# ╟─a3a9deea-e2d6-4d58-90d7-5a54be176289
# ╠═1b6e41ed-75a0-4b12-b390-78538566ce43
# ╠═e32d0418-6e18-485b-a6ea-0a982d63e1d2
# ╠═4e5a992d-303a-431c-a2c2-e5f3214a1da6
# ╠═f129bfac-b0d5-482b-98e5-856c26c0fc3f
# ╟─e41d1801-99f9-4c63-9bb2-96381f114c70
# ╠═54242145-429f-4480-8e5b-e3b51e347097
# ╠═4b854f3c-e957-4a6e-8fff-8fee52d47691
# ╟─b3c89340-6795-492b-9eee-5b06266818f0
# ╠═283da509-5437-4405-8aaf-19a56c64b037
# ╠═5e598c38-2fe1-475c-b836-82519af69efe
# ╟─fc49a6d7-a1b1-458a-a9ad-e120282bbabc
# ╠═c5815f5e-9164-11ec-10e1-691834761dff
# ╟─6be2f4c2-e9ed-43c2-b66c-ef3176bb9000
# ╠═aa19e992-2735-4324-8fd7-15eacadf0faa
# ╟─93090680-2380-412d-9752-df18251c7dbf
# ╠═15bbf544-34bd-4d38-bac5-0f43b1305df3
# ╟─9bc38d55-285b-4b83-98d9-d7f9e03405d1
# ╠═27844886-0b54-4b08-a592-a1a38e4b0be2
# ╠═d880f005-381a-4c72-a9a1-ae3eb84a90eb
# ╟─e8333b23-53c3-445e-9ca3-6b278359f8ab
# ╠═9248af7f-dc1a-4bf6-8f3f-304db73fc604
# ╠═dd626606-2a4d-494d-ab74-92753072c773
# ╠═e7b65566-a79e-4101-9471-3656a92e95e6
# ╠═e08bf013-00c7-4870-82d8-19b899e7208d
# ╟─9494e2a6-e2fd-4728-94d4-d68816a00e72
# ╠═9b744e72-07a2-4f8c-b88c-e0b1131794d0
# ╠═31d742a4-100f-4744-afe5-381b265b6f4c
# ╟─ab8b1a22-ca7a-409e-832e-8d5d08a29a1e
# ╠═f4d91971-f806-4c5c-8548-b58a20acfb2c
# ╟─fbe44944-499a-4881-94b6-07855d1165aa
# ╠═a5cd8e7a-380f-4203-a856-f9e56e04b092
# ╠═43a77919-d880-40c9-95c9-aaa429a65fb7
# ╠═36408698-f85c-4529-972f-6389534fcb88
# ╟─a62839d5-837b-4c37-996f-33659c34911c
# ╟─f95a08ce-a38d-4b7f-b478-4dbfa607740e
# ╟─90953c64-8a87-4065-8c34-d0ead540b728
# ╠═b993e3e6-8d4e-4de7-aa4b-6a2e3bd12212
# ╟─ae8012be-e7ab-4e85-a27f-febf08b3380b
# ╠═ce89710a-48d8-46d4-83b4-163a7eb0a2e5
# ╠═c34b1a5d-5078-4b8f-94d1-a088cbe5ab3e
# ╠═b7f4078a-ead0-4d42-8b44-4f471eefc6fc
# ╟─802d9652-7597-43c4-b13a-3c60682d0a69
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
