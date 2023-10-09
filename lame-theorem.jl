### A Pluto.jl notebook ###
# v0.19.29

#> [frontmatter]
#> title = "Lamé's Theorem"
#> description = "Are you in the near-field with P- and S-wave radiation? Don't worry, the response at your station has a higher-order discontinuity than in the far-field."

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

# ╔═╡ 086e74ac-30c6-11ed-3260-cb96ffde4e5c
begin
    using Symbolics
    using SpecialFunctions
    using SymbolicUtils
    using LinearAlgebra
    using PlutoTeachingTools
    using PlutoPlotly
    using PlutoTest
    using Latexify
    using PlutoUI
    using StatsBase
	using Distributions
end

# ╔═╡ 4bdd2493-70dd-4cf7-bcd4-b2a32aaff474
ChooseDisplayMode()

# ╔═╡ bccd3c27-2f9f-4dea-9d81-83cd6ba0ed9e
TableOfContents()

# ╔═╡ d6c09c3c-7fb6-4598-8897-a8f93b3c725e
md"""
# Lamé's Theorem
The elastic wave equation turns into simpler equations for potentials in the case of an unbounded, isotropic and homogeneous medium.
This interactive notebook delves into the fundamental concepts of near-field and far-field elastic wave radiation that emerge from seismic sources.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 280472ec-50db-49b6-9c96-5f9d5c80b3bb
Markdown.MD(Markdown.Admonition("warning", "Theorem",
    [TwoColumnWideLeft(md"""
If the displacement field $\bf{u} = \bf{u}(\bf{x}, t)$ satisfies,
$\rho\ddot{\bf{u}} = {\bf{f}} + (\lambda+2\mu)\nabla(\nabla\cdot\bf{u}) - \mu\nabla \times (\nabla \times \bf{u})$, then
∃ potentials $\phi$ and $\psi$ for $\bf{u}$ with properties
* ${\bf{u}} = \nabla\phi + \nabla \times \psi$
* $\nabla\cdot\psi = 0$
* $\rho\ddot{\phi} = f_1 + (\lambda+2\mu)\nabla^2\phi$
* $\rho\ddot{\psi} = \mathbf{f}_2 + (\mu)\nabla^2\psi$
""", md"""
###### Helmholtz potentials 
the body force and initial values of $\bf{u}$ and $\dot{\bf{u}}$ satisfy
${\bf{f}} = ∇f_1 + ∇ × \bf{f}_2$, $\nabla\cdot\bf{f}_2=0$
${\bf{u}}({\bf{x}}, 0) = ∇ u_1 + ∇ × \bf{u}_2$
${\bf{\dot{u}}}({\bf{x}}, 0) = ∇ \dot{u}_1 + ∇ × \bf{\dot{u}}_2$
$\nabla\cdot\bf{f}_2=0$,
$\nabla\cdot\bf{u}_2=0$,
$\nabla\cdot\bf{\dot{u}}_2=0$
		
	""")]))

# ╔═╡ 44c900ad-e9bf-4bc3-b456-898899ccd36b
md"The wavefield's analytical expression (which dynamically adjusts based on the chosen parameters above) is"

# ╔═╡ 24d2fecb-fe04-41f5-ab59-e5681e8cdca8
md"""
## Notation
To begin, we will introduce the variables that will be utilized throughout this notebook, and we will construct spatial and temporal differential operators.
"""

# ╔═╡ 25f04c16-a107-4f72-bf64-40ee7db3e24f
@syms x::Real y::Real z::Real  # spatial coordinates

# ╔═╡ 0c748592-b599-4ac2-9370-f3175c09c23c
𝐱 = [x, y, z] # spatial coordinate vector

# ╔═╡ a768aa28-7486-4d2f-9dee-a2f9be58a0e0
@syms t::Real ω::Real # time and angular frequency

# ╔═╡ ed35749b-88fd-428e-889c-95d20b9edd36
@syms α β ρ # P and S wave velocities, mass density

# ╔═╡ 5b649a90-fd7b-4cf3-82c2-b9168273e7f2
r = sqrt(sum(𝐱 .^ 2)) # distance from the origin, where the source is placed

# ╔═╡ 93cd7c3a-e1bf-45f0-9ab8-47980cbbc9c5
begin
    Dx = Differential(x)
    Dy = Differential(y)
    Dz = Differential(z)
    Dt = Differential(t)
    Dr = Differential(r)
end

# ╔═╡ e5addb67-5710-464b-b88f-9aa432df1829
∇ = [Dx, Dy, Dz]

# ╔═╡ b0a98128-4b78-471d-85c6-218c304b9578
md"""
## P Radiation
We shall now study P far-field radiation due to body force, e.g., an impulsive force in the x direction.
We now recall that one of the components of the displacement field is $\nabla \phi$,
$\phi$ is the P potential.
"""

# ╔═╡ 58dae178-1e8f-4c90-baac-ed8ca53aee99
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
* The curl of the P component of the displacement field is zero only when $\frac{1}{r^2}$ terms can be ignored, i.e., far away from the source.
	"""]))

# ╔═╡ cc80bf4e-ba95-4201-a077-e66f2bf79746
md"""
 Now, let's check if this component of the displacement field is like "P waves" with particle motion along the radial direction.
When the source is at the origin, 𝐱 is a vector in the direction from the source.
"""

# ╔═╡ be8f0cc1-83ed-458b-9a0b-75fabd15c701
Markdown.MD(Markdown.Admonition("warning", "Longitudinal Waves",
    [md"""
* The far-field P displacement field at any given location 𝐱 has a direction parallel to the direction from the source.
	"""]))

# ╔═╡ d4c329e6-52e7-4ba7-8941-a5724a038d05
md"""
## S Radiation
We shall now study S far-field radiation due to body force.
"""

# ╔═╡ fbbc5151-7726-4793-9ee2-05662631b3e1
Markdown.MD(Markdown.Admonition("warning", "",
    [md"""
* The divergence of the S component of the displacement field is zero only when $\frac{1}{r^2}$ terms can be ignored, i.e., far away from the source.
	"""]))

# ╔═╡ ec7a4506-698d-421d-a7d1-ae61136883f7
Markdown.MD(Markdown.Admonition("warning", "Transverse Waves",
    [md"""
* The far-field S displacement field at any given location 𝐱 has a direction perpendicular to the direction from the source.
	"""]))

# ╔═╡ aef669b3-4492-4487-a587-d6ffa20bd4f5
function curl(v)
    vx = v[1]
    vy = v[2]
    vz = v[3]
    curl_field = [
        Dy(vz) - Dz(vy),
        Dz(vx) - Dx(vz),
        Dx(vy) - Dy(vx)
    ]
    return simplify.(expand_derivatives.(curl_field))
end

# ╔═╡ 3bdae723-c32d-4264-bd57-184b73cddb3d
function divergence(v)
    vx = v[1]
    vy = v[2]
    vz = v[3]
    div_field = Dx(vx) + Dz(vz) + Dy(vy)
    return simplify.(expand_derivatives.(div_field))
end

# ╔═╡ c570f122-2b06-478d-a66e-6ea0b9871446
md"""
## Understanding Nearfield Radiation
We shall now look at the derivative w.r.t. the radial distance $r$ of the monochromatic potential field.
"""

# ╔═╡ 8c79d25c-589e-4198-b2b0-b2fcaab91b0c
ϕ = 1/r*sin(ω*(t-r/α))

# ╔═╡ 23395f6e-4a12-4550-b536-bb85881dfc81
Dx(ϕ) |> expand_derivatives

# ╔═╡ 9e00e076-1306-451b-a9d2-a80c48a7fb13
md"""
Notice that the far-field term contains a time-derivative of the forcing term. Lets plot the time series.
"""

# ╔═╡ bcc3a6c2-ece8-464d-985a-1e3dabae778f
md"""
## Appendix
"""

# ╔═╡ 8d7fbc8f-2b78-4708-8d3d-be9d4324c8ad
md"Need some structs for multiple dispatch"

# ╔═╡ a8eebe3a-0a95-4811-94de-615a04f3f3bf
begin
    struct X end
    struct Y end
    struct Z end
    struct Pfar end
    struct Sfar end
    struct Near end
    struct Gaussian end
    struct Monochromatic end
    SourceType = Union{Gaussian,Monochromatic}
end;

# ╔═╡ e98402cf-3d2a-4471-8628-61cb83e4aa26
begin
    derivative(::X, ex) = expand_derivatives(Dx(ex))
    derivative(::Y, ex) = expand_derivatives(Dy(ex))
    derivative(::Z, ex) = expand_derivatives(Dz(ex))
end

# ╔═╡ 3975cbd0-2e48-4dcc-a0a2-57db172d14a6
md"""
### Displacement Field
The displacement due to, for example, a body force $X_{0}(t)$ applied in the $x_{1}$ direction at the origin is composed of 3 terms: the first term is the near field term which is composed of both P and S-wave motions, the second term is the P-wave far field term and the third term is the S-wave far-field term. The $i$th component of the field is given by:

```math
{u_i}(\mathbf{x},t)=\frac{1}{4\pi\rho}\left(\frac{\partial^2}{\partial x_i \partial x_1}\frac{1}{r}\right)\int_{r/\alpha}^{r/\beta} \tau X_0(t-\tau) \,d\tau+\frac{1}{4\pi\rho\alpha^2r}\left(\frac{\partial r}{\partial x_i}\frac{\partial r}{\partial x_1}\right)X_0\left(t-\frac{r}{\alpha}\right)+
```
```math
\frac{1}{4\pi\rho\beta^2r}\left(\delta_{i1}-\frac{\partial r}{\partial x_i}\frac{\partial r}{\partial x_1}\right)X_0\left(t-\frac{r}{\beta}\right),
```
where, ${\rho}$ is the density of the medium, ${\alpha}$ is the P-wave velocity, ${\beta}$ is the S-wave velocity, ${r}$ is the distance between the source and receiver, ${\tau}$ is the convolution variable for time and ${\delta}_{i1}$ is the Kronecker Delta function.
"""

# ╔═╡ 30cc4084-b410-4f3d-b9e7-1ee6783295d9
begin
    # compute distance r, given receiver at [xr, yr, zr] and source at [xs, ys, zs]
    rad(xr, yr, zr, xs, ys, zs) = sqrt((xr - xs)^2 + (yr - ys)^2 + (zr - zs)^2)
    invrad(xr, yr, zr, xs, ys, zs) = inv(rad(xr, yr, zr, xs, ys, zs))

    # derivative of r w.r.t. x
    drad(::X, xr, yr, zr, xs, ys, zs,) = (xr - xs) / rad(xr, yr, zr, xs, ys, zs,)
    # derivative of r w.r.t. y
    drad(::Y, xr, yr, zr, xs, ys, zs,) = (yr - ys) / rad(xr, yr, zr, xs, ys, zs,)
    # derivative of r w.r.t. z
    drad(::Z, xr, yr, zr, xs, ys, zs,) = (zr - zs) / rad(xr, yr, zr, xs, ys, zs,)

    # far-field P-wave
    function get_displacement(::Nothing, ::Pfar, xr, yr, zr, xs, ys, zs, rc::Union{X,Y,Z}, sc::Union{X,Y,Z}, ρ, strength=1.0)
        return drad(rc, xr, yr, zr, xs, ys, zs,) * drad(sc, xr, yr, zr, xs, ys, zs,) * invrad(xr, yr, zr, xs, ys, zs,) * strength * inv(4 * pi * ρ * α * α)
    end
    # far-field S-wave
    function get_displacement(::Nothing, ::Sfar, xr, yr, zr, xs, ys, zs, rc::Union{X,Y,Z}, sc::Union{X,Y,Z}, ρ, strength=1.0)
        return (isequal(rc, sc) - drad(rc, xr, yr, zr, xs, ys, zs,) * drad(sc, xr, yr, zr, xs, ys, zs,)) * invrad(xr, yr, zr, xs, ys, zs,) * strength * inv(4.0 * pi * ρ * β * β)
    end
    # dc is derivative component of Pfar and Sfar (useful for moment-tensor solutions)
    # get_displacement(dc::Union{X,Y,Z}, ps::Union{Pfar,Sfar}, xr, yr, zr, xs, ys, zs, rc, sc, ρ, c, strength) = get_displacement(nothing, ps, xr, yr, zr, xs, ys, zs, rc, sc, ρ, c, strength) * drad(dc, xr, yr, zr, xs, ys, zs,)

    # near-field 
    function get_displacement(::Nothing, ::Near, xr, yr, zr, xs, ys, zs, rc::Union{X,Y,Z}, sc::Union{X,Y,Z}, ρ, strength=1.0)
        return derivative(rc, derivative(sc, invrad(xr, yr, zr, xs, ys, zs))) * strength * inv(4 * pi * ρ)
    end
end

# ╔═╡ 307416f3-54d6-4f9f-b29b-5a424d1f8451
md"""
### Forcing
We consider two types of sources in this notebook 1) monochromatic source 2) Gaussian pulse.
"""

# ╔═╡ 96753808-7580-4502-8084-a57642809089
begin
	σ2 = 5 # control std for Gaussian source
    f(::Gaussian, t) = exp(-t^2 / σ2) # work with retarded time later
    f(::Monochromatic, t) = sin(ω * t)
    f1(::Monochromatic, t) = cos(ω * t)
    function f(st::Gaussian, ::Near, t, xr, yr, zr, xs, ys, zs) # computes retarded times for P and S internally, using α and β
        tp = rad(xr, yr, zr, xs, ys, zs) / α
        ts = rad(xr, yr, zr, xs, ys, zs) / β
        return 0.5 * sqrt(σ2 * pi) * (erf((t - tp) / sqrt(σ2)) - erf((t - tp) / sqrt(σ2))) + 0.5 * σ2 * (f(st, t - tp) - f(st, t - ts))
    end
    function f(st::Monochromatic, ::Near, t, xr, yr, zr, xs, ys, zs) # computes retarded times for P and S internally, using α and β
        tp = rad(xr, yr, zr, xs, ys, zs) / α
        ts = rad(xr, yr, zr, xs, ys, zs) / β
        return (-f(st, t - tp) - tp * ω * f1(st, t - tp) + f(st, t - ts) + ts * ω * f1(st, t - ts)) * inv(ω) * inv(ω)
    end
    function f(st::SourceType, ::Pfar, t, xr, yr, zr, xs, ys, zs)
        tp = rad(xr, yr, zr, xs, ys, zs) / α
        f(st, t - tp)
    end
    function f(st::SourceType, ::Sfar, t, xr, yr, zr, xs, ys, zs)
        ts = rad(xr, yr, zr, xs, ys, zs) / β
        f(st, t - ts)
    end
	# include the sourcing in g 
    function get_displacement(sourcetype, dc, ps::Union{Pfar,Sfar,Near}, xr, yr, zr, xs, ys, zs, rc::Union{X,Y,Z}, sc::Union{X,Y,Z}, ρ, strength)
        return get_displacement(dc, ps, xr, yr, zr, xs, ys, zs, rc, sc, ρ, strength) * f(sourcetype, ps, t, xr, yr, zr, xs, ys, zs)
    end
end

# ╔═╡ 12e0afd1-559e-4f22-b76c-c1aa6469f693
upx = get_displacement(nothing, Pfar(), x, y, z, 0, 0, 0, X(), Y(), ρ, 1.0)

# ╔═╡ d98b20d0-853a-4f24-b089-101517d3cc78
upy = get_displacement(nothing, Pfar(), x, y, z, 0, 0, 0, Y(), Y(), ρ, 1.0)

# ╔═╡ f4ad875f-cffd-41cd-8f1f-02b62df5950f
upz = get_displacement(nothing, Pfar(), x, y, z, 0, 0, 0, Z(), Y(), ρ, 1.0)

# ╔═╡ a5a9db82-18df-4bef-9378-7cdcb1af95a1
up = [upx, upy, upz]

# ╔═╡ 91108dff-5a8b-42d2-85b7-5cda6a946406
curl(up)

# ╔═╡ c465a6dd-c103-4b97-ba9b-b1b22bb23c78
cross(𝐱, up)

# ╔═╡ 257d6706-bd8b-4dae-8bab-20900c1684d3
pradiation = sqrt(sum(abs2.([upx, upy, upz])))

# ╔═╡ f4aee41d-d63a-4a90-bfb1-a77d721a576b
usx = get_displacement(nothing, Sfar(), x, y, z, 0, 0, 0, X(), X(), ρ, 1.0)

# ╔═╡ 10535763-8613-49fc-b59e-997236109ba7
Dx(usx) |> expand_derivatives

# ╔═╡ 94adb048-f5c2-485e-833c-921c9a87814e
usy = get_displacement(nothing, Sfar(), x, y, z, 0, 0, 0, Y(), X(), ρ, 1.0)

# ╔═╡ f41ce900-88b5-479b-b7f0-b3cab5ab7bba
usz = get_displacement(nothing, Sfar(), x, y, z, 0, 0, 0, Z(), X(), ρ, 1.0)

# ╔═╡ 87e12020-0e28-46a3-814c-7660fc8c5745
sradiation = sqrt(sum(abs2.([usx, usy, usz])))

# ╔═╡ ccc20f79-1b54-43c7-a759-488a6fc53b0d
us = [usx, usy, usz]

# ╔═╡ e01bfb56-031c-4f6b-9835-45111e09f043
divergence(us) |> simplify

# ╔═╡ 0f88f2b0-294a-47e7-8753-f9450c7e722d
dot(𝐱, us) |> simplify

# ╔═╡ 99e4676f-a469-4243-afbd-7be1b587914a
substitute(dot(𝐱, us), sqrt(x^2+y^2+z^2)^2=>x^2+y^2+z^2) |> simplify

# ╔═╡ 8b707d16-c6b5-4a86-9156-f847ed265334
md"### Seismogram"

# ╔═╡ 68d4c48c-3bf5-49f5-ac32-0012b372eec4
md"### UI"

# ╔═╡ 7ba6bb60-4b0c-450a-af25-76a7b100f8d6
md"Some constants which are not in the main UI"

# ╔═╡ c9a6acbb-bbfa-42d9-ad9f-f511375df4d9
α1 = 4 # in km/s

# ╔═╡ b7e46a19-7273-4717-a737-b94b968e7250
β1 = 2 # in km/s

# ╔═╡ f021e0c8-0938-474f-a418-bb29c0a92052
ρ1 = 5e12 # kg/km3 

# ╔═╡ f2551112-c619-42ff-8b09-ec1271971781
md"### Discretize"

# ╔═╡ a80b0e95-d028-483e-a733-88d31422bd6b
md"Before plotting, we should discretize the space/time and ensure a small length to render 3D plots faster, especially when making changes."

# ╔═╡ 0f4defd4-f0d6-4ff1-8d02-138d5b4b8a99
begin
    # distances are in km
    # Don't choose an odd number of grid points here, due to nan values at (0,0,0)
    xgrid = range(-100, stop=100, length=26)
    zgrid = xgrid
    ygrid = xgrid
end;

# ╔═╡ 2957f6f0-7314-4340-9062-b6da2f1a7089
ThreeColumn(md"""
$(@bind sample_receiver Button("Sample a receiver"))
within radius
""",
md"""$(@bind rmax_plot Slider(round.(range(0, stop=maximum(abs, xgrid), length=100), digits=2), show_value=true, default=maximum(xgrid)/2))
and plot 
""",
	md"""
$(@bind plot_seismogram_field_type MultiCheckBox([Pfar()=>"P", Sfar()=>"S", Near()=>"Near"], default=[Pfar(), Sfar(), Near()], select_all=true))
""")

# ╔═╡ 1d0d7290-e561-4cb1-80d2-fd72c43d774a
begin
	sample_receiver
	rx_plot = rand(Uniform(-rmax_plot, rmax_plot))
	ry_plot = rand(Uniform(-rmax_plot, rmax_plot))
	rz_plot = rand(Uniform(-rmax_plot, rmax_plot))
end;

# ╔═╡ 2b31fe24-347d-4c0e-85e9-4fd8b9cb0e16
rz_plot

# ╔═╡ b90575c9-0e86-48a1-b656-954a0dd71968
plot_seismogram_field_type

# ╔═╡ 24c0331a-f479-4f29-b87d-cfabb75528cb
begin
	
    Xgrid3D = first.(Iterators.product(xgrid, ygrid, zgrid))
    Ygrid3D = map(x -> x[2], (Iterators.product(xgrid, ygrid, zgrid)))
    Zgrid3D = last.(Iterators.product(xgrid, ygrid, zgrid))

	xgrid2 = range(-100, stop=100, length=100)
	Xgrid2D = first.(Iterators.product(xgrid2, xgrid2))
	Zgrid2D = last.(Iterators.product(xgrid2, xgrid2))
end;

# ╔═╡ 6edf6598-8e51-4cb8-9c29-7759c17a52ae
tgrid = range(0, stop=100.0 * inv(α1), length=10)

# ╔═╡ d502d260-c562-476f-9e05-5dc1a2ce26b9
md"""
**Customize Radiation Visualization:** Adjust the visualization by modifying these parameters: Choose the component of the displacement field with $(@bind plot_rc Select([X() => "X", Y() => "Y", Z() => "Z"])), select the force density component using $(@bind plot_sc Select([X() => "X", Y() => "Y", Z() => "Z"])), set the forcing (source) type via $(@bind plot_source_type Select([Monochromatic()=>"Monochromatic", Gaussian()=>"Gaussian"], default=Gaussian())), pick the field type with $(@bind plot_field_type Select([Pfar() => "P", Sfar() => "S", Near() => "Near Field"])), adjust angular frequency for monochromatic sources with 
$(@bind ω1 Slider(range(0.2, stop=1, length=9), show_value=true, default=0.5)), and modify the time instance using $(@bind t_forw Slider(5:length(tgrid), default=7, show_value=true)). These controls allow you to explore different radiation patterns below.
"""

# ╔═╡ 25f3761a-3f94-4aa8-9fb8-0fbbe32fe095
wavefield = get_displacement(plot_source_type, nothing, plot_field_type, x, y, z, 0, 0, 0, plot_rc, plot_sc, ρ, 1e20)

# ╔═╡ 24901509-9dcc-4af7-b708-dcc277565ebb
get_displacement(plot_source_type, nothing, Near(), x, y, z, 0, 0, 0, X(), plot_sc, ρ, 1) |> simplify

# ╔═╡ 1df1eea7-b5cb-477c-991a-14ac21aea154
function get_seismogram(tgrid, sourcetype, dc, xr, yr, zr, xs, ys, zs,  sc::Union{X,Y,Z}, strength)
	ex = map([X(), Y(), Z()]) do rc  
		mapreduce(+, plot_seismogram_field_type) do ps
			get_displacement(sourcetype, nothing, ps, x, y, z, xs, ys, zs, rc, sc, ρ, strength)
		end
	end

	return map(ex) do exp
	map(tgrid) do t1
		substitute(exp,[t=>t1,ρ=>ρ1,α=>α1,β=>β1, ω=>ω1, x=>xr, y=>yr, z=>zr])
	end
	end
end

# ╔═╡ ac7a947f-283f-4743-a42b-082a0eb0c42e
md"""
**Customize Radiation Visualization:** Adjust the visualization by modifying these parameters: Choose the forcing (source) type via $(@bind plot_source_type2 Select([Monochromatic()=>"Monochromatic", Gaussian()=>"Gaussian"])), adjust angular frequency for monochromatic sources with 
$(@bind ω2 Slider(range(0.2, stop=1, length=9), show_value=true, default=0.5)), and modify the time instance using $(@bind t_forw2 Slider(5:length(tgrid), default=7, show_value=true)). These controls allow you to explore the differences between near-field and far-field terms.
"""

# ╔═╡ 65e2b15c-c646-4148-82c8-707350e0e112
example_field = expand_derivatives(Dr(f(plot_source_type2, t - r / α) / r))

# ╔═╡ b0c536eb-9ec9-4745-82fe-50b94330a5c9
md"""
Notice that `example_field`, has two terms: one that decays as $\frac{1}{r^2}$ (near-field term) and the far-field term that decays as $\frac{1}{r}$.

* ...at a distance of $r=100$, the far-field term is stronger by a factor of **$(round((ω2/α1*100), digits=3))**, which is *2π × number of wavelengths between the source and the receiver*. Weak enough? Let's plot the near-field and far-field terms separately after building the necessary functions.
"""

# ╔═╡ c497d1b0-4b29-4938-bdad-3de63e6cd022
tgrid_seismogram = range(-50.0 * inv(β1), stop=300.0 * inv(β1), length=1000)

# ╔═╡ eeff59bb-b32f-423f-ba4a-0c193958a79a
seismogram = get_seismogram(tgrid_seismogram, plot_source_type, nothing, rx_plot, ry_plot, rz_plot, 0, 0, 0,  plot_sc, 1e20) # choosing ry=1, rz=1 to avoid source location (0,0,0)

# ╔═╡ df67edac-e750-4ea8-a131-796f1ee9f86f
md"These functions build and discretize a given expression over spatial and temporal grids.
Arguments: `ex`: A symbolic expression representing the function to be discretized.
Returns: A 3D array representing the discretized field of the provided expression over spatial and temporal grids. P and S wave velocities and angular frequency will be substituted before discretizing."

# ╔═╡ 6dcaa88e-7191-4d6d-bcbf-8195c3f863cf
function build_and_discretize3D(ex)
    test_plot = build_function(ex, x, y, z, t, α, β, ρ, ω, expression=Val{false})
    field = map(tgrid) do t1
        map(Xgrid3D, Ygrid3D, Zgrid3D) do x, y, z
            test_plot(x, y, z, t1, α1, β1, ρ1, ω1) # ω1 is used for 3D
        end
    end
    return field
end

# ╔═╡ dbe87c18-347a-4373-890a-8a0c51fc1ff0
function build_and_discretize2D(ex) # fix y=0
    test_plot = build_function(ex, x, y, z, t, α, β, ρ, ω, expression=Val{false})
    field = map(tgrid) do t1
        map(Xgrid2D, Zgrid2D) do x, z
            test_plot(x, 0.0, z, t1, α1, β1, ρ1, ω2) # ω2 is used for 2D
        end
    end
    return field
end

# ╔═╡ 767543ee-152d-499a-8690-7dba8e28ae88
md"### Plots"

# ╔═╡ 10897515-6242-47b5-8407-22ddbf15f0d6
md"This is a function that creates a 3D volume plot after discretizing the wavefield expression at a given time instance global variable `t_forw`."

# ╔═╡ e0cf9350-325c-4edd-9f77-dcf18c4c8b04
function add_field!(fig, field, t_forw)
    # colorbar axis lims are chosen keeping all the time steps in mind, so that 
    # the attenuation can be observed
    cmax = maximum(map(field) do f
        return maximum(abs, f)
	end) * 0.1

    # select time instance
    values = field[t_forw]


    add_trace!(fig, volume(
        x=Xgrid3D[:],
        y=Ygrid3D[:],
        z=Zgrid3D[:],
        value=values[:],
        isomin=-cmax,
        isomax=cmax,
        colorscale="RdBu",
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=50,# needs to be a large number for good volume rendering
        caps=attr(x_show=false, y_show=false, z_show=false, x_fill=1),
        opacityscale="extremes"
    ))
end;

# ╔═╡ 5dde0ad0-926d-4326-9564-12efc7141edc
function plot_ex_wavefronts3D(ex, title="Elastic-wave Radiation")
    field = build_and_discretize3D(ex)
    layout = Layout(
		uirevision = 1,
        scene=attr(
            xaxis=attr(
                nticks=3,
                range=[-100, 100]
            ),
            yaxis=attr(
                nticks=3,
                range=[-100, 100]
            ),
            zaxis=attr(
                nticks=3,
                range=[-100, 100]
            ),
        ),
        width=650,
        margin=attr(
            r=20,
            l=10,
            b=10,
            t=10
        ),
        title=attr(
            text=title,
            y=0.9,
            x=0.15,)
    )
    fig = Plot(layout)
    add_field!(fig, field, t_forw)

N=10
	add_trace!(fig, scatter3d(x = [rx_plot], y = [ry_plot], z = [rz_plot], mode="markers+text", marker_color="black", marker_size=3, text="receiver"))
	add_trace!(fig, scatter3d(x = [0], y = [0], z = [0], mode="markers+text", marker_color="black", marker_size=3, text="source"))
	
    return fig
end



# ╔═╡ 7917e7eb-ecef-4efe-9cca-606895199a7d
plot(plot_ex_wavefronts3D(wavefield))

# ╔═╡ f2f6bbc8-6c05-4b09-82fd-b4f1ed76da12
function plot_ex_wavefronts2D(nearex, farex)
	nearfield = build_and_discretize2D(nearex)
	farfield = build_and_discretize2D(farex)

	cmax_near = maximum(map(nearfield) do f
        return maximum(abs, f)
	end) * 0.006
	cmax_far = maximum(map(farfield) do f
        return maximum(abs, f)
	end) * 0.1

    fig = Plot(Layout(yaxis_autorange="reversed", height=300, width=650, 
		title=attr(font_size=12,), font=attr(
            size=10), yaxis=attr(scaleanchor="x"), Subplots(horizontal_spacing=0.3, rows=1, cols=2, subplot_titles=["Near Field" "Far Field"])))
    add_trace!(fig, heatmap(
            x=xgrid2,
            y=xgrid2,
            z=nearfield[t_forw2], zmin=-cmax_near,
        zmax=cmax_near, colorscale="seismic", showscale=true, colorbar_x=0.35), row=1, col=1)
    add_trace!(fig, heatmap(
            x=xgrid2,
            y=xgrid2,
            z=farfield[t_forw2], zmin=-cmax_far,
        zmax=cmax_far, colorscale="seismic", showscale=true, colorbar_x=1.0), row=1, col=2)

    return PlutoPlotly.plot(fig)

end

# ╔═╡ 4019938e-defe-4024-a7cd-28f34660e46d
plot_ex_wavefronts2D(first(arguments(example_field)), last(arguments(example_field)))

# ╔═╡ 2f140a6e-4775-441c-9f53-a1f33b897251
function plot_seismogram(seismogram)
fig=PlutoPlotly.Plot(Layout(height=200, width=500, 
		title="Displacement at (x, y, z)=($(round(rx_plot, digits=3)), $(round(ry_plot, digits=3)),$(round(rz_plot, digits=3)))", show_legend=true, font=attr(
            size=10)))
	tgrid=tgrid_seismogram
add_trace!(fig,PlutoPlotly.scatter(x=tgrid,y=seismogram[1],line_width=2,mode="lines",line_color="red",name="X"))
	add_trace!(fig,PlutoPlotly.scatter(x=tgrid,y=seismogram[2],line_width=2,mode="lines",line_color="blue",name="Y"))
	add_trace!(fig,PlutoPlotly.scatter(x=tgrid,y=seismogram[3],line_width=2,mode="lines",line_color="black",name="Z"))
PlutoPlotly.plot(fig)
end

# ╔═╡ 3d43de3c-d60f-4997-aff8-f8027319485c
plot_seismogram(seismogram)

# ╔═╡ 3d86326e-735f-4bcd-86b8-d31b132ee036
md"""
Inhomogeneous wave equation
We can now define the wave operator of a scalar wave equation in three dimensions.
"""

# ╔═╡ 01ac68da-a533-4fba-ad5d-8d12c2014ae1
L(ϕ, α) = Dt(Dt(ϕ)) / α^2 - (Dx(Dx(ϕ)) + Dy(Dy(ϕ)) + Dz(Dz(ϕ)))

# ╔═╡ 854072b0-68c0-4492-b56a-d630fbb4d82a
# print the scalar wave equation 
# latexify(L(u, α) ~ δ(x)δ(y)δ(z) * f(t))

# ╔═╡ 3c2300c5-770b-42ad-a182-d274158b9923
md""" Spherical waves
We will now consider a monochromatic source sin function of space and time with amplitude $A$ and parameters $\omega$ and $K$.
"""

# ╔═╡ 4eb8d905-0ac1-45eb-9543-9dbfc4bad709
Dt(sin(ω * t)) |> expand_derivatives

# ╔═╡ 3c6835b0-bd1c-4e99-b8e1-156e700033fa
sphwav = build_function(ϕ, r, t, α, ω, expression=Val{false})

# ╔═╡ 5629a75b-a436-4127-a880-7bc3fce811a5
gif_Ppotential = plot_ex_wavefronts3D(ϕ, "ϕ");

# ╔═╡ 316460ee-fafc-4057-a0a8-3885335569cf
md"""
Similarly, along $y$ and $z$.
"""

# ╔═╡ 2077aa71-9a0f-4b46-97b2-3157e857fcda


# ╔═╡ 89722347-66ba-4b11-ac95-92ba735e58a0
@bind body_forces confirm(PlutoUI.combine() do Child
    components = ["x", "y", "z"]
    s2 = [md"""
       "x" $(
       	Child(NumberField(range(-1, step=0.1, stop=1), default=1))
       ) "y" $( 
       	Child(NumberField(range(-1, step=0.1, stop=1), default=0))
       ) "z" $(
       	Child(NumberField(range(-1, step=0.1, stop=1), default=0))
       )""" for is in 1:2]

    md"""
    	
    Input Body Force
    Until now, we have considered a spherically symmetric field $ϕ$ as P-potential.
    According to Lame's theorem, the P-potential should be generated using a curl-free component of the forcing field. We refer the reader to Eq. 4.23 Aki Richards for a detailed derivation. Here we will just write down the final expressions.
    	
    $(s2[1])
    """
end)

# ╔═╡ Cell order:
# ╠═4bdd2493-70dd-4cf7-bcd4-b2a32aaff474
# ╠═bccd3c27-2f9f-4dea-9d81-83cd6ba0ed9e
# ╟─d6c09c3c-7fb6-4598-8897-a8f93b3c725e
# ╟─280472ec-50db-49b6-9c96-5f9d5c80b3bb
# ╟─d502d260-c562-476f-9e05-5dc1a2ce26b9
# ╟─7917e7eb-ecef-4efe-9cca-606895199a7d
# ╟─44c900ad-e9bf-4bc3-b456-898899ccd36b
# ╟─25f3761a-3f94-4aa8-9fb8-0fbbe32fe095
# ╟─2957f6f0-7314-4340-9062-b6da2f1a7089
# ╟─3d43de3c-d60f-4997-aff8-f8027319485c
# ╟─24d2fecb-fe04-41f5-ab59-e5681e8cdca8
# ╠═25f04c16-a107-4f72-bf64-40ee7db3e24f
# ╠═0c748592-b599-4ac2-9370-f3175c09c23c
# ╠═a768aa28-7486-4d2f-9dee-a2f9be58a0e0
# ╠═ed35749b-88fd-428e-889c-95d20b9edd36
# ╠═5b649a90-fd7b-4cf3-82c2-b9168273e7f2
# ╠═93cd7c3a-e1bf-45f0-9ab8-47980cbbc9c5
# ╠═e5addb67-5710-464b-b88f-9aa432df1829
# ╟─b0a98128-4b78-471d-85c6-218c304b9578
# ╠═12e0afd1-559e-4f22-b76c-c1aa6469f693
# ╠═d98b20d0-853a-4f24-b089-101517d3cc78
# ╠═f4ad875f-cffd-41cd-8f1f-02b62df5950f
# ╠═a5a9db82-18df-4bef-9378-7cdcb1af95a1
# ╠═91108dff-5a8b-42d2-85b7-5cda6a946406
# ╟─58dae178-1e8f-4c90-baac-ed8ca53aee99
# ╠═257d6706-bd8b-4dae-8bab-20900c1684d3
# ╟─cc80bf4e-ba95-4201-a077-e66f2bf79746
# ╠═c465a6dd-c103-4b97-ba9b-b1b22bb23c78
# ╟─be8f0cc1-83ed-458b-9a0b-75fabd15c701
# ╟─d4c329e6-52e7-4ba7-8941-a5724a038d05
# ╠═f4aee41d-d63a-4a90-bfb1-a77d721a576b
# ╠═94adb048-f5c2-485e-833c-921c9a87814e
# ╠═f41ce900-88b5-479b-b7f0-b3cab5ab7bba
# ╠═87e12020-0e28-46a3-814c-7660fc8c5745
# ╠═ccc20f79-1b54-43c7-a759-488a6fc53b0d
# ╠═10535763-8613-49fc-b59e-997236109ba7
# ╠═e01bfb56-031c-4f6b-9835-45111e09f043
# ╟─fbbc5151-7726-4793-9ee2-05662631b3e1
# ╠═0f88f2b0-294a-47e7-8753-f9450c7e722d
# ╠═99e4676f-a469-4243-afbd-7be1b587914a
# ╟─ec7a4506-698d-421d-a7d1-ae61136883f7
# ╠═aef669b3-4492-4487-a587-d6ffa20bd4f5
# ╠═3bdae723-c32d-4264-bd57-184b73cddb3d
# ╟─c570f122-2b06-478d-a66e-6ea0b9871446
# ╟─ac7a947f-283f-4743-a42b-082a0eb0c42e
# ╠═4019938e-defe-4024-a7cd-28f34660e46d
# ╠═8c79d25c-589e-4198-b2b0-b2fcaab91b0c
# ╠═23395f6e-4a12-4550-b536-bb85881dfc81
# ╠═65e2b15c-c646-4148-82c8-707350e0e112
# ╟─b0c536eb-9ec9-4745-82fe-50b94330a5c9
# ╟─9e00e076-1306-451b-a9d2-a80c48a7fb13
# ╟─bcc3a6c2-ece8-464d-985a-1e3dabae778f
# ╠═086e74ac-30c6-11ed-3260-cb96ffde4e5c
# ╟─8d7fbc8f-2b78-4708-8d3d-be9d4324c8ad
# ╠═a8eebe3a-0a95-4811-94de-615a04f3f3bf
# ╠═e98402cf-3d2a-4471-8628-61cb83e4aa26
# ╟─3975cbd0-2e48-4dcc-a0a2-57db172d14a6
# ╠═30cc4084-b410-4f3d-b9e7-1ee6783295d9
# ╟─307416f3-54d6-4f9f-b29b-5a424d1f8451
# ╠═96753808-7580-4502-8084-a57642809089
# ╟─8b707d16-c6b5-4a86-9156-f847ed265334
# ╠═1d0d7290-e561-4cb1-80d2-fd72c43d774a
# ╠═24901509-9dcc-4af7-b708-dcc277565ebb
# ╠═2b31fe24-347d-4c0e-85e9-4fd8b9cb0e16
# ╠═b90575c9-0e86-48a1-b656-954a0dd71968
# ╠═1df1eea7-b5cb-477c-991a-14ac21aea154
# ╠═eeff59bb-b32f-423f-ba4a-0c193958a79a
# ╟─68d4c48c-3bf5-49f5-ac32-0012b372eec4
# ╟─7ba6bb60-4b0c-450a-af25-76a7b100f8d6
# ╠═c9a6acbb-bbfa-42d9-ad9f-f511375df4d9
# ╠═b7e46a19-7273-4717-a737-b94b968e7250
# ╠═f021e0c8-0938-474f-a418-bb29c0a92052
# ╟─f2551112-c619-42ff-8b09-ec1271971781
# ╟─a80b0e95-d028-483e-a733-88d31422bd6b
# ╠═0f4defd4-f0d6-4ff1-8d02-138d5b4b8a99
# ╠═24c0331a-f479-4f29-b87d-cfabb75528cb
# ╠═6edf6598-8e51-4cb8-9c29-7759c17a52ae
# ╠═c497d1b0-4b29-4938-bdad-3de63e6cd022
# ╟─df67edac-e750-4ea8-a131-796f1ee9f86f
# ╠═6dcaa88e-7191-4d6d-bcbf-8195c3f863cf
# ╠═dbe87c18-347a-4373-890a-8a0c51fc1ff0
# ╟─767543ee-152d-499a-8690-7dba8e28ae88
# ╟─10897515-6242-47b5-8407-22ddbf15f0d6
# ╠═5dde0ad0-926d-4326-9564-12efc7141edc
# ╠═e0cf9350-325c-4edd-9f77-dcf18c4c8b04
# ╠═f2f6bbc8-6c05-4b09-82fd-b4f1ed76da12
# ╠═2f140a6e-4775-441c-9f53-a1f33b897251
# ╟─3d86326e-735f-4bcd-86b8-d31b132ee036
# ╠═01ac68da-a533-4fba-ad5d-8d12c2014ae1
# ╠═854072b0-68c0-4492-b56a-d630fbb4d82a
# ╟─3c2300c5-770b-42ad-a182-d274158b9923
# ╠═4eb8d905-0ac1-45eb-9543-9dbfc4bad709
# ╠═3c6835b0-bd1c-4e99-b8e1-156e700033fa
# ╠═5629a75b-a436-4127-a880-7bc3fce811a5
# ╟─316460ee-fafc-4057-a0a8-3885335569cf
# ╠═2077aa71-9a0f-4b46-97b2-3157e857fcda
# ╠═89722347-66ba-4b11-ac95-92ba735e58a0
