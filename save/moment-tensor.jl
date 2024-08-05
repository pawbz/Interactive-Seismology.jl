### A Pluto.jl notebook ###
# v0.19.43

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

# ╔═╡ 4a6eaf04-1efe-11ed-3e32-6bdb3afcc5b3
begin
    using PlutoUI
	using ForwardDiff
	using Interpolations
	import ForwardDiff.derivative
	using PlutoTeachingTools
    using CairoMakie # a better plotting lib to vis. vector fields
	# WGLMakie.activate!()
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

# ╔═╡ eb75f980-4d19-438b-83d6-3ea5e1a0a1e0
# test if the P-wave displacement is parallal to r̂
@test all(abs.(cross(displacement(P(), rand_receiver...), rand_receiver)) .< eps())

# ╔═╡ 2d7ca9f8-7575-4b2f-8c2b-80d728a79635
# test if the S-wave displacement is always perpendicular to r̂
@test abs(dot(displacement(S(), rand_receiver...), rand_receiver)) < eps()

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
	points_cube=Meshes.Point3.(coordinates.(vertices(Meshes.Box(Meshes.Point(-1,-1,-1), Meshes.Point(1,1,1)))))

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
