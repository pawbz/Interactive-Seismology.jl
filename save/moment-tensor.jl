### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ 4a6eaf04-1efe-11ed-3e32-6bdb3afcc5b3
begin
    using PlutoUI
	using ForwardDiff
	using Interpolations
	import ForwardDiff.derivative
	using PlutoTeachingTools
    # using CairoMakie
	using WGLMakie# a better plotting lib to vis. vector fields
	WGLMakie.activate!()
	using GeoMakie
    using LinearAlgebra
	import PlutoUI: combine
	using Meshes
	using CoordinateTransformations
	using PlutoTest
	using MeshViz
	set_theme!(theme_ggplot2())
end

# ╔═╡ 1169b57f-60e0-47cf-9e87-e987be9d72fd
md"""
# Far-field P & S Radiation Patterns

The purpose of this notebook is to visualize the far-field particle displacement, for both P and S waves, in a homogeneous medium due to an input seismic moment tensor. In the notebook, for a given input moment tensor, we shall first visualize the equivalent body-force couples and dipoles, then compute the far-field displacement for the equivalent forcing terms (see Aki & Richards; eq. 4.23).

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)
ES218; August 2022

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 220ff7e0-2afe-4a29-8781-27f9ba7c8ff7
@bind moment_tensor confirm(PlutoUI.combine() do Child
		components=["x", "y", "z"]
	# components=[1, 2, 3]
		s1=[md"""
		   $(c,"x") $(
		Child(NumberField(range(-1, step=0.1, stop=1), default=0))
		)  $(c,"y") $(
		Child(NumberField(range(-1, step=0.1,stop=1), default=1))
		)  $(c,"z") $(
		Child(NumberField(range(-1, step=0.1, stop=1), default=0))
		)""" for c in components]

		md"""
		
		
		## Input Moment Tensor 
		Each element of the moment tensor is equivalent to a body-force couple.
		$(s1...)

		
		"""
end)

# ╔═╡ 8f53975e-6d81-42f8-8936-c6b952c7913f
warning_box(md"Here, we are **not** forcing the moment tensor to be symmetric!")

# ╔═╡ 01f3e07d-de2d-4efe-82b0-fd1363f32158
@bind body_forces confirm(PlutoUI.combine() do Child
	components=["x", "y", "z"]
	s2=[md"""
		"x" $(
			Child(NumberField(range(-1, step=0.1, stop=1), default=1))
		) "y" $( 
			Child(NumberField(range(-1, step=0.1, stop=1), default=0))
		) "z" $(
			Child(NumberField(range(-1, step=0.1, stop=1), default=0))
		)""" for is in 1:2]
	
md"""
	
## Input Body Forces
	
Add any body forces independent to that of the moment tensor. At the moment, we have 2 body forces.
#### Force I
$(s2[1])
#### Force II
$(s2[2])
"""
end)

# ╔═╡ a220fe12-8611-4630-9fd8-8aaa2a4818bc
md"""
## Force Distribution
"""

# ╔═╡ af057ffd-d750-4f04-8e70-e51b8c2b2a08
md"""
## Displacement Field (2-D projections)
We shall now plot the displacement field projected on either $xy$, $yz$ or $xz$ planes.
"""

# ╔═╡ c25ed391-bf65-472f-9b44-2888a26c18fe
tip(md"The far-field P displacement is curl free, and the far-field S displacement is divergence free.")

# ╔═╡ dc0ed979-2e37-4fbf-972d-b1f30cd7ce23
md"""
## Displacement Field (3D)
"""

# ╔═╡ 89b56471-6339-44c4-a81b-1c58e2fe72fc
tip(md"The above 3D plots from WGLMakie are interactive, so you can use them to train your intuition.")

# ╔═╡ bacfc397-1c38-4184-bbde-f5f81a9a5df3
md"""
## Appendix
"""

# ╔═╡ 5e687f7d-2f91-412d-ab2a-d1da145a5d47
md"""
### Equivalent Forces
Methods the are used to convert the input moment tensor to an equivalent for distribution. 
"""

# ╔═╡ c2983a3b-005e-4de5-b9dc-2692850ff627
xyz=Iterators.product([:x, :y, :z], [:x, :y, :z]); # need an iterator for elements of moment tensor

# ╔═╡ b9fd379b-df75-4ba0-bd4d-0ddc64561e4a
md"""
### Displacement Field in Homogeneous Medium 
"""

# ╔═╡ 3db5e6c8-812f-49d7-a0c7-b00d8bc7f8ff
md"""
The displacement due to, for example, a body force $X_{0}(t)$ applied in the $x_{1}$ direction at the origin is composed of 3 terms: the first term is the near field term which is composed of both P and S-wave motions, the second term is the P-wave far field term and the third term is the S-wave far-field term. The $i$th component of the field is given by:

```math
{u_i}(\mathbf{x},t)=\frac{1}{4\pi\rho}\left(\frac{\partial^2}{\partial x_i \partial x_1}\frac{1}{r}\right)\int_{r/\alpha}^{r/\beta} \tau X_0(t-\tau) \,d\tau+\frac{1}{4\pi\rho\alpha^2r}\left(\frac{\partial r}{\partial x_i}\frac{\partial r}{\partial x_1}\right)X_0\left(t-\frac{r}{\alpha}\right)+
```
```math
\frac{1}{4\pi\rho\beta^2r}\left(\delta_{i1}-\frac{\partial r}{\partial x_i}\frac{\partial r}{\partial x_1}\right)X_0\left(t-\frac{r}{\beta}\right),
```
where, ${\rho}$ is the density of the medium, ${\alpha}$ is the P-wave velocity, ${\beta}$ is the S-wave velocity, ${r}$ is the distance between the source and receiver, ${\tau}$ is the convolution variable for time and ${\delta}_{i1}$ is the Kronecker Delta function.
"""

# ╔═╡ 4e4d324b-4405-402d-bbbf-d9287f7cc29f
α = 5000.0 # P-wave velocity in m/s

# ╔═╡ 96512612-6d21-4f80-b8f6-a9efee9797c8
β = 3000.0 # S-wave velocity in m/s

# ╔═╡ 3c5daacc-baad-4109-83d6-faba39137cf5
rand_receiver=randn(3)

# ╔═╡ eb072a57-352e-47de-af61-f6f6259ea763
rP = 1; rS= 1;

# ╔═╡ 875182de-0955-43ae-9805-6dbb8f44b057
# These are used for multiple dispatch
begin
	struct x end
	struct y end
	struct z end
	struct xy end
	struct yz end
	struct xz end
	struct P end
	struct S end
end

# ╔═╡ 64f2a518-83aa-4770-b421-1b3bc45bc7a2
# Functions that outputs the positions (x, y, z) of the force couples
begin
	dcouple = .5 # half distance b/w the elements of a force couple
	body_force_locations(::x) = [[dcouple, 0.0, 0.0], [-dcouple, 0.0, 0.0]]
	body_force_locations(::y) = [[0.0, dcouple, 0.0], [0.0, -dcouple, 0.0]]
	body_force_locations(::z) = [[0.0, 0.0, dcouple], [0.0, 0.0, -dcouple]]
end

# ╔═╡ e2ed2cef-1767-4586-92ad-fb28947cc4a9
# collect forces (just for plotting)
srclocs=reshape(cat([body_force_locations(eval(x[2])()) for x in xyz]..., dims=1), 1, 1, 9*2);

# ╔═╡ 1839abea-3eef-40cd-8e79-5647bb169364
# Functions that output the components of the forces
begin
	body_force_strengths(f, ::x) = reshape([(1, f), (2,0), (3,0), (1, -f), (2,0), (3,0)], 3, 2)
	body_force_strengths(f, ::y) = reshape([(1, 0), (2,f), (3,0), (1, 0), (2,-f), (3,0)], 3, 2)
	body_force_strengths(f, ::z) = reshape([(1, 0), (2,0), (3,f), (1, 0), (2,0), (3,-f)], 3,2)
end

# ╔═╡ c4f34e9b-8554-409a-aab4-df53e9a1e041
# collect strengths of all elements in the moment tensor X 2 forces (just for plotting)
srcstrength=reshape(cat([body_force_strengths(moment_tensor[i], eval(x[1])()) for (i, x) in enumerate(xyz)]..., dims=2), 1, 3, 9*2);

# ╔═╡ 5bd9de9d-1cbf-4e89-8ab0-ade3181ff4b6
begin
	r(xr, yr, zr, xs=zero(xr), ys=zero(yr), zs=zero(zr)) = sqrt(abs2(xr - xs)+abs2(yr-ys)+abs2(zr-zs))
	dr(::x, xr, yr, zr)= xr / r(xr, yr, zr)
	dr(::y, xr, yr, zr)= yr / r(xr, yr, zr)
	dr(::z, xr, yr, zr)= zr / r(xr, yr, zr)
	# far-field P-wave; ignoring overall amplitude factor
	function g(::Nothing, ::P, xr, yr, zr, rc::Union{x, y, z}, sc::Union{x, y, z}, strength=1.0)		return dr(rc, xr, yr, zr) * dr(sc, xr,yr,zr) / r(xr, yr, zr) * strength
	end
	# far-field S-wave; ignoring overall amplitude factor;
	function g(::Nothing, ::S, xr, yr, zr, rc::Union{x, y, z}, sc::Union{x, y, z}, strength=1.0)
		return (isequal(rc,sc) - dr(rc, xr,yr,zr) * dr(sc, xr,yr,zr)) / r(xr, yr, zr) * strength
	end
	g(dc::Union{x, y, z}, ps::Union{P, S}, xr, yr, zr, rc, sc, strength) = g(nothing, ps, xr, yr, zr, rc, sc, strength)*dr(dc, xr, yr, zr)
	
	# Commented (will include near-field terms as well)
	# g(::x, ps::Union{P, S}, xr, yr, zr, rc, sc, strength) = derivative(x->g(nothing, ps, x, yr, zr, rc, sc, strength), xr)[1]
	# g(::y, ps::Union{P, S}, xr, yr, zr, rc, sc, strength) = derivative(y->g(nothing, ps, xr, y, zr, rc, sc, strength), yr)[1]
	# g(::z, ps::Union{P, S}, xr, yr, zr, rc, sc, strength) = derivative(z->g(nothing, ps, xr, yr, z, rc, sc, strength), zr)[1]
end

# ╔═╡ 559928f4-bac3-4ecd-918a-32c4c3c54119
# combine elements of the moment tensor and body forces together into a simple vector that we can later iterate over
m=vcat(broadcast(moment_tensor, vec(collect(xyz))) do m, c
	(m, eval(c[1])(), eval(c[2])())
end, 
	broadcast(tuple, body_forces, repeat([x(), y(), z()], div(length(body_forces), 3)), nothing)
);

# ╔═╡ a61a40b2-59f7-4aaa-85c2-a159984dfc25
function displacement(ps::Union{P,S}, rx, ry, rz, rc, m=m)
mapreduce(+, m) do c
	g(c[3], ps, rx, ry, rz, rc, c[2], c[1])
end
end

# ╔═╡ 2fabd638-b308-4927-99fe-8bccce404396
# method that computes displacement for all three components
displacement(ps::Union{P,S}, rx, ry, rz)=map(rc->displacement(ps, rx, ry, rz, rc, m), [x(), y(), z()])

# ╔═╡ eb75f980-4d19-438b-83d6-3ea5e1a0a1e0
# test if the P-wave displacement is parallal to r̂
@test all(abs.(cross(displacement(P(), rand_receiver...), rand_receiver)) .< eps())

# ╔═╡ 2d7ca9f8-7575-4b2f-8c2b-80d728a79635
# test if the S-wave displacement is always perpendicular to r̂
@test abs(dot(displacement(S(), rand_receiver...), rand_receiver)) < eps()

# ╔═╡ 2f29f440-4cc6-4ab0-95fe-6baf87c0cd06
begin
	ddisplacement(::x, ps, xr, yr, zr, rc, m=m)=derivative(x->displacement(ps, x, yr, zr, rc, m), xr)[1]
	ddisplacement(::y, ps, xr, yr, zr, rc)=derivative(y->displacement(ps, xr, y, zr, rc, m), yr)[1]
	ddisplacement(::z, ps, xr, yr, zr, rc)=derivative(z->displacement(ps, xr, yr, z, rc, m), zr)[1]
end

# ╔═╡ 4660b617-3d5f-422a-9747-f2405b1a4c29
md"""
We shall now choose the number of points, where the displacement vector is eventually evaluated. For 2D, we uniformly sample the azimuth, and for 3D, we will uniformly sample the surface of a unit sphere.
"""

# ╔═╡ 383bdfca-c621-49d9-9f8d-49221283754f
begin
	# define a cube in R^3
	points_cube=Meshes.Point3.(coordinates.(vertices(boundary(Meshes.Box(Meshes.Point(-1,-1,-1), Meshes.Point(1,1,1))))))

	connec_cube = connect.([(1,4,3,2),(5,6,7,8),(1,2,6,5),(3,4,8,7),(1,5,8,4),(2,3,7,6)])
	sphmesh = refine(refine(SimpleMesh(points_cube, connec_cube), CatmullClark()), CatmullClark());
	sphpts = broadcast(coordinates ∘ centroid, collect(elements(sphmesh)))


    n3 = length(sphpts);    n2 = 50
	# a grid for theta (azimuth)
    thgrid = range(-pi, stop = pi, length = n2)
	nothing
end

# ╔═╡ 9b735440-d831-4249-8083-d1469d8176a0
begin
	cart3D(r, sphpt)=r.*sphpt;
	cart2D(r, th)=[r*sin(th), r*cos(th)];
	cart2D(r, th, ::xy)=[r*sin(th), r*cos(th), 0];
	cart2D(r, th, ::yz)=[0, r*sin(th), r*cos(th)];
	cart2D(r, th, ::xz)=[r*sin(th), 0,  r*cos(th)];
	proj2D(ux, uy, uz, ::xy)=ux, uy
	proj2D(ux, uy, uz, ::yz)=uy, uz
	proj2D(ux, uy, uz, ::xz)=ux, uz
end

# ╔═╡ 8f1cc2eb-5c79-47a5-a35e-995d1e5d8775
function displacement2D(ps::Union{P,S}, plane::Union{xy, yz, xz})
	return map(c->broadcast(th->displacement(ps, cart2D(1.0, th, plane)..., c), thgrid), [x(), y(), z()])
end

# ╔═╡ e9450f23-8879-4b49-9689-c50c755463f0
function displacement3D(ps::Union{P,S})
	return map(c->broadcast(sphpt->displacement(ps, cart3D(1.0, sphpt)..., c), sphpts), [x(), y(), z()])
end

# ╔═╡ 2110ae2c-26ac-4f35-b392-897581ff4f5b
md"""
### Packages
"""

# ╔═╡ a09edfd9-64d6-413a-8766-e20d0b1ef8ba
TableOfContents()

# ╔═╡ 4f3e43cd-4f35-46b2-a972-0cc8f2ac7b3f
md"""
### Plotting Methods
"""

# ╔═╡ e2177159-23c1-42bd-a5f8-7fb46d68a31f
md"""
We are only interested in vis. of the direction of displacement. So lets define some function to normalize the fields and remove geometrical spreading.
"""

# ╔═╡ e453aca7-6d3c-421e-83bb-eb772cb001b2
function unormalize(args...)
	# fact=inv(norm(cat(args..., dims=3)))
	# return map(u->u .* fact, [args...])
	return args
end

# ╔═╡ 9c94f9ef-7241-4dd9-ab59-47bcc148a9bc
function get_arrow_sizes(args...)
	a = vec(sqrt.(sum(map(x->abs2.(x),  [args...]))))
    a = a * 10
end

# ╔═╡ b295b338-f97a-4c75-af23-f5ae291a12f2
md"""
Now lets work on plotting 2D/ 3D arrow heads using `Makie.jl`.
"""

# ╔═╡ 11d8323a-66e9-48e3-83f9-5b94c776dc09
begin
	function plotsources()
	ps = Figure(resolution = (800, 800))
    Axis3(
        ps[1, 1],
        xlabel = "x (m)",
        ylabel = "y (m)",
		zlabel = "z (m)",
        title = "Equivalent Force Distribution",
		    )
	limits!(-5, 5, -5, 5, -5, 5)
	arrows!(
		    vec(getindex.(srclocs, 1)), vec(getindex.(srclocs, 2)), vec(getindex.(srclocs, 3)), vec(getindex.(srcstrength[1,1,:],2)),
		vec(getindex.(srcstrength[1,2,:],2)), 
		vec(getindex.(srcstrength[1,3,:],2)), 
				  linecolor=:black, 
		# arrowsize=vec(ap3),
		arrowsize=0.5*Makie.Vec3f(1, 1, 1),
		lengthscale=1,
				arrowcolor=:black,
	    linewidth = 0.05,
	)
	return ps
	end
end

# ╔═╡ b2c260ad-2577-45cc-83ed-105d3d42f419
plotsources()

# ╔═╡ 873277ce-a1c0-456c-8964-a3b973e9f41e
function plot_displacement_2D()
	f2 = Figure(resolution = (800, 1200))
	for (i,dim) in enumerate([["x","y"], ["y","z"], ["x","z"]])
		plane=eval(Symbol(dim...))()
	up1, up2, up3 = displacement2D(P(), plane);
	us1, us2, us3 = displacement2D(S(), plane);
	
	# prepare for 2-D vis.
	up1_plot, up2_plot, up3_plot, us1_plot, us2_plot, us3_plot= unormalize(up1, up2, up3, us1, us2, us3)
	up1_plot, up2_plot = proj2D(up1_plot, up2_plot, up3_plot, plane)
	us1_plot, us2_plot = proj2D(us1_plot, us2_plot, us3_plot, plane)

	pp2 = vec(broadcast(x->Makie.Point2f(cart2D(rP,x)), thgrid))
	ps2 = vec(broadcast(x->Makie.Point2f(cart2D(rS,x)), thgrid))
    ns2 = broadcast((x, y) -> Makie.Vec2f(x, y), us1_plot, us2_plot)
    np2 = broadcast((x, y) -> Makie.Vec2f(x, y), up1_plot, up2_plot)
    ap2 = get_arrow_sizes(up1_plot, up2_plot)
    as2 = get_arrow_sizes(us1_plot, us2_plot)
	
    Axis(
        f2[i, 1],
        xlabel = "$(dim[1]) (m)",
        ylabel = "$(dim[2]) (m)",
        title = "Displacement field ($(dim...)) slice",
		aspect=1
    )

    arrows!(pp2, np2, lengthscale = 0.1, arrowcolor = :blue, linecolor=:blue,arrowsize = ap2, xlabel = "x", align=:center)
			
    Axis(
        f2[i, 2],
        xlabel = "$(dim[1]) (m)",
        ylabel = "$(dim[2]) (m)",
        title = "Displacement field ($(dim...)) slice",
    )
    arrows!(ps2, ns2, lengthscale = 0.1, arrowcolor = :red, linecolor=:red, arrowsize = as2, align=:center)
	
    
	end
	f2
	
end

# ╔═╡ 3a39ae6e-d4eb-4d1f-99dc-210bb6a22772
plot_displacement_2D()

# ╔═╡ 32aec472-4002-42ff-b67c-2c08fe3fa0bc
function plot_displacement_3D()

	uxp3, uyp3, uzp3 = displacement3D(P());
	uxs3, uys3, uzs3 = displacement3D(S());
	
	# prepare for 3-D vis.
	uxp3_plot, uyp3_plot, uzp3_plot, uxs3_plot, uys3_plot, uzs3_plot = unormalize(uxp3, uyp3, uzp3, uxs3, uys3, uzs3)
	
	pp3 = vec(broadcast(x->Makie.Point3f(cart3D(rP, x)), sphpts))
	ps3 = vec(broadcast(x->Makie.Point3f(cart3D(rS, x)), sphpts))
	ns3 = broadcast((x, y, z) -> Makie.Vec3f(x, y, z), uxs3_plot, uys3_plot, uzs3_plot)
    np3 = broadcast((x, y, z) -> Makie.Vec3f(x, y, z), uxp3_plot, uyp3_plot, uzp3_plot)
    ap3 = get_arrow_sizes(uxp3_plot, uyp3_plot, uzp3_plot)
    as3 = get_arrow_sizes(uxs3_plot, uys3_plot, uzs3_plot)
	f3 = Figure(resolution = (800, 800))
    Axis3(
        f3[1, 1],
        xlabel = "x (m)",
        ylabel = "y (m)",
		zlabel = "z (m)",
        title = "Displacement field"
    )
	global Pcolors=broadcast((ux, uy, uz, s)->dot(s, [ux, uy, uz]), uxp3, uyp3, uzp3, sphpts)
	viz!(sphmesh, color=Pcolors, alpha=0.5)
	arrows!(
		    getindex.(pp3, 1), getindex.(pp3, 2), getindex.(pp3, 3), vec(uxp3_plot), vec(uyp3_plot), vec(uzp3_plot);
		  linecolor=:black, 
		# arrowsize=vec(ap3),
		arrowsize = 0.1*Makie.Vec3f(1, 1, 1),
		lengthscale=0.1,
				arrowcolor=:blue,
	    linewidth = 0.02, align=:center
	)
	 Axis3(
        f3[2, 1],
        xlabel = "x (m)",
        ylabel = "y (m)",
		zlabel = "z (m)",
        title = "Displacement field",
    )
	# limits!(-rP, rP, -rP, rP, -rP, rP)	
viz!(sphmesh, color = :gray, alpha=0.1)
	arrows!(#ps3, ns3,
	    getindex.(ps3, 1), getindex.(ps3, 2), getindex.(ps3, 3), vec(uxs3_plot), vec(uys3_plot), vec(uzs3_plot);
		  linecolor=:black, 
		lengthscale=0.1,
		arrowsize = 0.1*Makie.Vec3f(1, 1, 1),
		# arrowsize=vec(ap3),
		arrowcolor=:red,
	    linewidth = 0.02, align=:center
	)	
	f3
end

# ╔═╡ 98442df9-beac-41c4-aeb8-490b3377ccc8
plot_displacement_3D()

# ╔═╡ 473acfd0-0a69-4426-af47-4ac283965431

TODO()


# ╔═╡ e4721ad5-e407-47ed-89c9-6ca3d099c9ad
md"""
* Have beachballs plotted for a given moment tensor.
* Choose the moment tensors corresponding to strike-slip, normal, reverse, and oblique faults. Explain your choices and plot the corresponding radiation patterns.
"""

# ╔═╡ ada78100-ede9-4b42-a04f-eb87461ce9fc
begin
lons =  rad2deg.(getfield.(SphericalFromCartesian().(sphpts), :θ))
lats =  rad2deg.(getfield.(SphericalFromCartesian().(sphpts), :ϕ))

# 	lons1 = -180:180
# lats1 = -90:90
	I = sortperm(lons)
	lons = lons[I]
	lats = lats[I]
# 
	field=Pcolors[I]
	
	field .-= sum(field)
	field=sign.(field)
	# itp = interpolate((lons,) , field, Gridded(Linear()))
	

	# field = [itp[l, y] for l in lons1, y in lats1]
	
	fig = Figure()
	ax1 = GeoAxis(fig[1, 1], dest = "+proj=vitk1 +lat_1=0 +lat_2=96",
	    coastlines = false, aspect=1)
	heatmap!(ax1, lons, lats, field; shading = false, colormap = :viridis)
	fig
end

# ╔═╡ e10dd2e3-9753-4745-97ee-c40ad1621491
broadcast(coordinates ∘ centroid, collect(elements(sphmesh)))

# ╔═╡ 26390557-abfa-4495-b94f-1dd4e639faa3
# broadcast(coordinates centroid, collect(elements(sphmesh)))

# ╔═╡ d377fe49-08fe-49d4-b85d-990f3f8b4503
 rad2deg.(getfield.(SphericalFromCartesian().(sphpts), :ϕ))

# ╔═╡ 46a5ffa0-7c12-4561-9fb8-8ac6f4fa57f9
unique(lats)

# ╔═╡ 11cf79d3-efd5-494f-abed-c22a4a29ce21
sortperm

# ╔═╡ 51fcd87e-a560-442f-a933-e1c6dc95bba5
Pcolors .- sum(Pcolors)

# ╔═╡ 7b0a5157-ab2c-414e-9f4c-74b041a91cd1
field

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CoordinateTransformations = "150eb455-5306-5404-9cee-2592286d6298"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
GeoMakie = "db073c08-6b98-4ee5-b6a4-5efafb3259c6"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MeshViz = "9ecf9c4f-6e5a-4b5e-83ae-06f2c7a661d8"
Meshes = "eacbb407-ea5a-433e-ab97-5258b1ca43fa"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

[compat]
CoordinateTransformations = "~0.6.3"
ForwardDiff = "~0.10.36"
GeoMakie = "~0.5.1"
Interpolations = "~0.15.1"
MeshViz = "~0.8.8"
Meshes = "~0.32.3"
PlutoTeachingTools = "~0.2.15"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.59"
WGLMakie = "~0.8.16"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "ea33027ba95036f272b0409a252647ec4e91b1e6"

[[deps.ADTypes]]
git-tree-sha1 = "60665b326b75db6517939d0e1875850bc4a54368"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.17.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractLattices]]
git-tree-sha1 = "763b6f3f6bfabd72c7e262cbb5ddfd43fd5c6398"
uuid = "398f06c4-4d28-53ec-89ca-5b2656b7603d"
version = "0.3.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "dbd8c3bbbdbb5c2778f85f4422c39960eac65a42"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.20.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
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
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "1568b28f91293458345dabba6a5ea3f183250a61"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.8"

    [deps.CategoricalArrays.extensions]
    CategoricalArraysJSONExt = "JSON"
    CategoricalArraysRecipesBaseExt = "RecipesBase"
    CategoricalArraysSentinelArraysExt = "SentinelArrays"
    CategoricalArraysStructTypesExt = "StructTypes"

    [deps.CategoricalArrays.weakdeps]
    JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SentinelArrays = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
    StructTypes = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CircularArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "e24a6f390e5563583bb4315c73035b5b3f3e7ab4"
uuid = "7a955b69-7140-5f4e-a0ed-f168c5e2e749"
version = "1.4.0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "5ac098a7c8660e217ffac31dc2af0964a8c3182a"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "2.0.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "e771a63cc8b539eca78c85b0cabd9233d6c8f06f"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

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
git-tree-sha1 = "8010b6bb3388abe68d95743dcbea77650bb2eddf"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "0037835448781bb46feb39866934e243886d756a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a692f5e257d332de1e554e4566a4e5a8a72de2b2"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.4"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["DataStructures", "EnumX", "ExactPredicates", "Random", "SimpleGraphs"]
git-tree-sha1 = "d4e9dc4c6106b8d44e40cd4faf8261a678552c7c"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "0.8.12"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "82adce9b61f01f036934b6fb2b100636df941b37"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "2.2.6+0"

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

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "16946a4d305607c3a4af54ff35d56f0e9444ed0e"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.7"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3e6d038b77f22791b8e3472b7c633acea1ecac06"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.120"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7bb1361afdb33c7f2b085aa49ea8fe1b0fb14e58"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.1+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3a948313e7a41eb1db7a1e733e6335f17b4ab3c4"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "7.1.1+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "797762812ed063b9b94f6cc7742bc8883bb5e69e"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.9.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

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

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "31fd32af86234b6b71add76229d53129aa1b87a9"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.28.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Formatting]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "fb409abab2caf118986fc597ba84b50cbaf00b87"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.3"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "8e233d5167e63d708d41f87597433f59a0f213fe"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.4"

[[deps.GeoInterface]]
deps = ["DataAPI", "Extents", "GeoFormatTypes"]
git-tree-sha1 = "294e99f19869d0b0cb71aef92f19d03649d028d5"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.4.1"

[[deps.GeoInterfaceRecipes]]
deps = ["GeoInterface", "RecipesBase"]
git-tree-sha1 = "fb1156076f24f1dfee45b3feadb31d05730a49ac"
uuid = "0329782f-3d07-4b52-b9f6-d3137cf03c7a"
version = "1.0.2"

[[deps.GeoJSON]]
deps = ["Extents", "GeoFormatTypes", "GeoInterface", "GeoInterfaceRecipes", "JSON3", "StructTypes", "Tables"]
git-tree-sha1 = "29bb9a35cac72ccbeb92e37dbfb2feedd5b5aecd"
uuid = "61d90e0f-e114-555e-ac52-39dfb47a3ef9"
version = "0.7.3"

[[deps.GeoMakie]]
deps = ["Colors", "Downloads", "GeoInterface", "GeoJSON", "GeometryBasics", "ImageIO", "LinearAlgebra", "Makie", "Proj", "Reexport", "Statistics", "StructArrays"]
git-tree-sha1 = "06753d663f23fd344b70b972c98d2f4fe65252c5"
uuid = "db073c08-6b98-4ee5-b6a4-5efafb3259c6"
version = "0.5.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "5694b56ccf9d15addedc35e9a4ba9c317721b788"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.10"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "35fbd0cefb04a516104b8e183ce0df11b70a3f1a"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "f57a64794b336d4990d90f80b147474b869b1bc4"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

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
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Random", "RoundingEmulator"]
git-tree-sha1 = "79342df41c3c24664e5bf29395cfdf2f2a599412"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.36"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

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
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "411eccfe8aba0814ffa0fdf4860913ed09c34975"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.3"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JSServe]]
deps = ["Base64", "CodecZlib", "Colors", "Dates", "Deno_jll", "HTTP", "Hyperscript", "LinearAlgebra", "Markdown", "MsgPack", "Observables", "RelocatableFolders", "SHA", "Sockets", "Tables", "ThreadPools", "URIs", "UUIDs", "WidgetsBase"]
git-tree-sha1 = "4bcf2a78f7c80c6f3d594267bb4e7ec03ac9c172"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f9"
version = "2.3.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "d8337622fe53c05d16f031df24daf0270e53bc64"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.5"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "ba51324b894edaf1df3ab16e2cc6bc3280a2f1a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.10"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "706dfd3c0dd56ca090e86884db6eda70fa7dd4af"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d3c8af829abaeba27181db4acb485b18d15d89c6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.1+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "d5d2e3abfb30ea9c2cff81d243e7235b51315ec2"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.2"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "4adee99b7262ad2a1a4bbbc59d993d24e55ea96f"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.4.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "Primes", "SimplePolynomials"]
git-tree-sha1 = "054fe197a75adcf517933f0c807051d79881f3f3"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.2.10"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "73b98709ad811a6f81d84e105f4f695c229385ba"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.4.3"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Setfield", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "35fa3c150cd96fd77417a23965b7037b90d6ffc9"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.12"

[[deps.MakieCore]]
deps = ["Observables", "REPL"]
git-tree-sha1 = "9b11acd07f21c4d035bd4156e789532e8ee2cc70"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.9"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "96ca8a313eb6437db5ffe946c457a401bbb8ce1d"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.7"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MeshViz]]
deps = ["CategoricalArrays", "ColorSchemes", "Colors", "Dates", "Makie", "Meshes", "ScientificTypes", "Tables", "Unitful"]
git-tree-sha1 = "f92cf44d22f993f7911fe1dfdc3e80d3901b93e4"
uuid = "9ecf9c4f-6e5a-4b5e-83ae-06f2c7a661d8"
version = "0.8.8"

[[deps.Meshes]]
deps = ["Bessels", "CircularArrays", "Distances", "IterTools", "LinearAlgebra", "NearestNeighbors", "Random", "Rotations", "SparseArrays", "StaticArrays", "StatsBase", "Tables", "TransformsBase"]
git-tree-sha1 = "8e2877e5a4596e2d6af69764852e7e0c23147ed7"
uuid = "eacbb407-ea5a-433e-ab97-5258b1ca43fa"
version = "0.32.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Mods]]
git-tree-sha1 = "0a3c164857cfc4defe0ac09e6b0123d61320abb1"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "2.2.6"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[deps.Multisets]]
git-tree-sha1 = "8ef67bad0a3def6d77aebf1d3645119e7c8e5a54"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.5"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2ae7d4ddec2e13ad3bddf5c0796f7547cf682391"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.2+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "61942645c38dd2b5b78e2082c9b51ab315315d10"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.2"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.PROJ_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Libtiff_jll", "SQLite_jll"]
git-tree-sha1 = "29d4224ba271e3b679cca03a8fd72e271525c113"
uuid = "58948b4f-47e0-5654-a9ad-f609743f8632"
version = "902.600.200+0"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Permutations]]
deps = ["Combinatorics", "LinearAlgebra", "Random"]
git-tree-sha1 = "b1f03a4943c62552a12c7f95965b76c3f91cf5b7"
uuid = "2ae35dd2-176d-5d53-8349-f30d82d94d4f"
version = "0.4.23"

[[deps.PikaParser]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "d6ff87de27ff3082131f31a714d25ab6d0a88abf"
uuid = "3bbf5609-3e7b-44cd-8549-7c69f321e792"
version = "0.6.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

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

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "2660ba4e0e9354cf0f3deeefb61531bbdac5495d"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.21"

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

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.Proj]]
deps = ["CEnum", "CoordinateTransformations", "GeoFormatTypes", "GeoInterface", "NetworkOptions", "PROJ_jll"]
git-tree-sha1 = "61188669db4f5b400173e4ec60da8bcb72d6e749"
uuid = "c94c279d-25a6-4763-9509-64d165bea63e"
version = "1.9.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "994cc27cdacca10e68feb291673ec3a76aa2fae9"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "d852eba0cc08181083a58d5eb9dccaec3129cb03"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.9.0"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.RingLists]]
deps = ["Random"]
git-tree-sha1 = "70d8d52c6d4238cace6ae1a7dd2c0c4d43acdf0b"
uuid = "286e9d63-9694-5540-9e3c-4e6708fa07b2"
version = "0.2.9"

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

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.SQLite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "9a325057cdb9b066f1f96dc77218df60fe3007cb"
uuid = "76ed43ae-9a5d-5a62-8c75-30186b810ce8"
version = "3.48.0+0"

[[deps.ScientificTypes]]
deps = ["CategoricalArrays", "ColorTypes", "Dates", "Distributions", "PrettyTables", "Reexport", "ScientificTypesBase", "StatisticalTraits", "Tables"]
git-tree-sha1 = "3ad7f09ae97806e86b3ef28cd50c52ca18bd2d80"
uuid = "321657f4-b219-11e9-178b-2701a2544e81"
version = "3.1.1"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "79123bc60c5507f035e6d1d9e563bb2971954ec8"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleGraphs]]
deps = ["AbstractLattices", "Combinatorics", "DataStructures", "IterTools", "LightXML", "LinearAlgebra", "LinearAlgebraX", "Optim", "Primes", "Random", "RingLists", "SimplePartitions", "SimplePolynomials", "SimpleRandom", "SparseArrays", "Statistics"]
git-tree-sha1 = "5415585dfac38ace09c0641595155b82ff579d9a"
uuid = "55797a34-41de-5266-9ec1-32ac4eb504d3"
version = "0.8.7"

[[deps.SimplePartitions]]
deps = ["AbstractLattices", "DataStructures", "Permutations"]
git-tree-sha1 = "76cdce94e07849d1172ccb8f3a60a0bbbebe58b3"
uuid = "ec83eff0-a5b5-5643-ae32-5cbf6eedec9d"
version = "0.3.3"

[[deps.SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "77b2c5e731f7e50bbd088ccc2810232e06c1afa3"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.18"

[[deps.SimpleRandom]]
deps = ["Distributions", "LinearAlgebra", "Random"]
git-tree-sha1 = "bb4f42b25b87f124478207a82f5b02dfafdb3e63"
uuid = "a6525b86-64cd-54fa-8f65-62fc48bdc0e8"
version = "0.3.2"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["Compat", "PikaParser", "SHA", "StructTypes", "Tables", "TupleTools", "WeakKeyIdDicts"]
git-tree-sha1 = "2c44178ac401356247fa0516371e37d49c957897"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "1.3.4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "89f86d9376acd18a1a4fbef66a56335a3a7633b8"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.5.0"

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
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2c962245732371acd51700dbb268af311bddd719"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.6"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "8e45cecc66f3b42633b8ce14d431e8e57a3e242e"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "9537ef82c42cdd8c5d443cbc359110cbb36bae10"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.21"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

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
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

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

[[deps.ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TransformsBase]]
deps = ["AbstractTrees"]
git-tree-sha1 = "c209ddc5ea678a542aabe482713b24c9280919ed"
uuid = "28dd2a49-a57a-4bfb-84ca-1a49db9b96b8"
version = "1.6.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "41e43b9dc950775eac654b9f845c839cd2f1821e"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.6.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "6258d453843c466d84c17a58732dda5deeb8d3af"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.24.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WGLMakie]]
deps = ["Colors", "FileIO", "FreeTypeAbstraction", "GeometryBasics", "Hyperscript", "JSServe", "LinearAlgebra", "Makie", "Observables", "PNGFiles", "PrecompileTools", "RelocatableFolders", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "2328be9950072c7250d340798b25222779ce55bb"
uuid = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
version = "0.8.16"

[[deps.WeakKeyIdDicts]]
git-tree-sha1 = "349165887b0cfdcbb97261e7eb0e7f11d7b39845"
uuid = "ecbed89c-9d59-4137-ae1f-d1084086d01c"
version = "0.1.0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "c217bad8fccb3bbfef7d7902326eacfbd0d702ad"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.14.4+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4bba74fa59ab0755167ad24f98800fe5d727175b"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.12.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

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

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╟─1169b57f-60e0-47cf-9e87-e987be9d72fd
# ╟─220ff7e0-2afe-4a29-8781-27f9ba7c8ff7
# ╟─8f53975e-6d81-42f8-8936-c6b952c7913f
# ╟─01f3e07d-de2d-4efe-82b0-fd1363f32158
# ╠═a220fe12-8611-4630-9fd8-8aaa2a4818bc
# ╠═b2c260ad-2577-45cc-83ed-105d3d42f419
# ╟─af057ffd-d750-4f04-8e70-e51b8c2b2a08
# ╠═3a39ae6e-d4eb-4d1f-99dc-210bb6a22772
# ╟─c25ed391-bf65-472f-9b44-2888a26c18fe
# ╟─dc0ed979-2e37-4fbf-972d-b1f30cd7ce23
# ╠═98442df9-beac-41c4-aeb8-490b3377ccc8
# ╟─89b56471-6339-44c4-a81b-1c58e2fe72fc
# ╟─bacfc397-1c38-4184-bbde-f5f81a9a5df3
# ╟─5e687f7d-2f91-412d-ab2a-d1da145a5d47
# ╠═c2983a3b-005e-4de5-b9dc-2692850ff627
# ╠═c4f34e9b-8554-409a-aab4-df53e9a1e041
# ╠═e2ed2cef-1767-4586-92ad-fb28947cc4a9
# ╠═64f2a518-83aa-4770-b421-1b3bc45bc7a2
# ╠═1839abea-3eef-40cd-8e79-5647bb169364
# ╠═b9fd379b-df75-4ba0-bd4d-0ddc64561e4a
# ╟─3db5e6c8-812f-49d7-a0c7-b00d8bc7f8ff
# ╠═4e4d324b-4405-402d-bbbf-d9287f7cc29f
# ╠═96512612-6d21-4f80-b8f6-a9efee9797c8
# ╠═5bd9de9d-1cbf-4e89-8ab0-ade3181ff4b6
# ╠═559928f4-bac3-4ecd-918a-32c4c3c54119
# ╠═a61a40b2-59f7-4aaa-85c2-a159984dfc25
# ╠═2fabd638-b308-4927-99fe-8bccce404396
# ╠═3c5daacc-baad-4109-83d6-faba39137cf5
# ╠═eb75f980-4d19-438b-83d6-3ea5e1a0a1e0
# ╠═2d7ca9f8-7575-4b2f-8c2b-80d728a79635
# ╠═2f29f440-4cc6-4ab0-95fe-6baf87c0cd06
# ╠═eb072a57-352e-47de-af61-f6f6259ea763
# ╠═8f1cc2eb-5c79-47a5-a35e-995d1e5d8775
# ╠═e9450f23-8879-4b49-9689-c50c755463f0
# ╠═875182de-0955-43ae-9805-6dbb8f44b057
# ╟─4660b617-3d5f-422a-9747-f2405b1a4c29
# ╠═383bdfca-c621-49d9-9f8d-49221283754f
# ╠═9b735440-d831-4249-8083-d1469d8176a0
# ╟─2110ae2c-26ac-4f35-b392-897581ff4f5b
# ╠═4a6eaf04-1efe-11ed-3e32-6bdb3afcc5b3
# ╠═a09edfd9-64d6-413a-8766-e20d0b1ef8ba
# ╟─4f3e43cd-4f35-46b2-a972-0cc8f2ac7b3f
# ╟─e2177159-23c1-42bd-a5f8-7fb46d68a31f
# ╠═e453aca7-6d3c-421e-83bb-eb772cb001b2
# ╠═9c94f9ef-7241-4dd9-ab59-47bcc148a9bc
# ╟─b295b338-f97a-4c75-af23-f5ae291a12f2
# ╠═11d8323a-66e9-48e3-83f9-5b94c776dc09
# ╠═873277ce-a1c0-456c-8964-a3b973e9f41e
# ╠═32aec472-4002-42ff-b67c-2c08fe3fa0bc
# ╟─473acfd0-0a69-4426-af47-4ac283965431
# ╠═e4721ad5-e407-47ed-89c9-6ca3d099c9ad
# ╠═ada78100-ede9-4b42-a04f-eb87461ce9fc
# ╠═e10dd2e3-9753-4745-97ee-c40ad1621491
# ╠═26390557-abfa-4495-b94f-1dd4e639faa3
# ╠═d377fe49-08fe-49d4-b85d-990f3f8b4503
# ╠═46a5ffa0-7c12-4561-9fb8-8ac6f4fa57f9
# ╠═11cf79d3-efd5-494f-abed-c22a4a29ce21
# ╠═51fcd87e-a560-442f-a933-e1c6dc95bba5
# ╠═7b0a5157-ab2c-414e-9f4c-74b041a91cd1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
