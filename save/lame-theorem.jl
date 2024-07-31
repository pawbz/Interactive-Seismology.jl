### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Distributions = "~0.25.104"
Latexify = "~0.16.1"
PlutoPlotly = "~0.4.4"
PlutoTeachingTools = "~0.2.14"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.54"
SpecialFunctions = "~2.3.1"
StatsBase = "~0.34.2"
SymbolicUtils = "~1.5.0"
Symbolics = "~5.13.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "72afc6d83824c162d8bb31b63c3011ed690959ac"

[[deps.ADTypes]]
git-tree-sha1 = "332e5d7baeff8497b923b730b994fa480601efc7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.5"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "e0cf629808a64842d2b133f7b5c492dc6580f16e"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.34.5"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "793501dcd3fa7ce8d375a2c878dca2296232686e"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Test"]
git-tree-sha1 = "a7055b939deae2455aa8a67491e034f735dd08d3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.33"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "bbec08a37f8722786d87bedf84eae19c020c4efa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "4b41ad09c2307d5f24e36cd6f92eb41b218af22c"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.1"

[[deps.Bijections]]
git-tree-sha1 = "c9b163bd832e023571e86d0b90d9de92a9879088"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.6"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "2118cb2765f8197b08e5958cdd17c165427425ee"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.19.0"
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
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

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

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "886826d76ea9e72b35fcd000e535588f7b60f21d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9242eec9b7e2e14f9952e8ea1c7e31a50501d587"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.104"

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

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "fea68c84ba262b121754539e6ea0546146515d4f"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.3"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

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
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
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

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

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
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "3d8866c029dd6b16e69e0d4a939c4dfcb98fac47"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.8"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "68772f49f54b479fa88ace904f6127f0a3bb2e46"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.12"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

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
git-tree-sha1 = "e49bce680c109bc86e3e75ebcb15040d6ad9e1d3"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.27"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "f12f2225c999886b69273f84713d1b9cb66faace"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.15.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LatticeRules]]
deps = ["Random"]
git-tree-sha1 = "7f5b02258a3ca0221a6a9710b0a0a2e8fb4957fe"
uuid = "73f95e8e-ec14-4e6a-8b18-0d2e271c4e55"
version = "0.0.1"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

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
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

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

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "0b8cf121228f7dae022700c1c11ac1f04122f384"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.2"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "6ffb234d6d7c866a75c1879d2099049d3a35a83a"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.3"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "806eea990fb41f9b36f1253e5697aa645bf6a9f8"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

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

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "a935806434c9d4c506ba941871b327b96d41f2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.0"

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
git-tree-sha1 = "58dcb661ba1e58a13c7adce77435c3c6ac530ef9"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.4.4"

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

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "bd7c69c7f7173097e7b5e1be07cee2b8b7447f51"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.54"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "01ac95fca7daabe77a9cb705862bd87016af9ddb"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.13"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "1d05623b5952aed1307bf8b43bec8b8d1ef94b6e"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.5"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9ebcd48c498668c7fa0e97a9cae873fbee7bfee1"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.1"

[[deps.QuasiMonteCarlo]]
deps = ["Accessors", "ConcreteStructs", "LatticeRules", "LinearAlgebra", "Primes", "Random", "Requires", "Sobol", "StatsBase"]
git-tree-sha1 = "cc086f8485bce77b6187141e1413c3b55f9a4341"
uuid = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
version = "0.3.3"
weakdeps = ["Distributions"]

    [deps.QuasiMonteCarlo.extensions]
    QuasiMonteCarloDistributionsExt = "Distributions"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

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
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "f7d4d01d3321b92c2217d599704d29c45e12905e"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.1.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
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
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "6990168abf3fe9a6e34ebb0e05aaaddf6572189e"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.10"

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
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d8911cc125da009051fb35322415641d02d9e37f"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.6"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "QuasiMonteCarlo", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "2c5cb9fcae1db5ce71ded89cb9314ad09f6c89bc"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.11.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "51ae235ff058a64815e0a2c34b1db7578a06813d"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.7"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "5165dfb9fd131cf0c6957a3a7605dede376e7b63"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.0"

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

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "fba11dbe2562eecdfcac49a05246af09ee64d055"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.8.1"
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
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymbolicIndexingInterface]]
git-tree-sha1 = "65f4ed0f9e3125e0836df12c231cea3dd98bb165"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "849b1dfb1680a9e9f2c6023f79a49b694fb6d0da"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.5.0"

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

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

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

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "21c8fc7cd598ef49f11bc9e94871f5d7740e34b9"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.5"

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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
