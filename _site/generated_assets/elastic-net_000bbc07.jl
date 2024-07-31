### A Pluto.jl notebook ###
# v0.19.36

#> [frontmatter]
#> chapter = "1"
#> layout = "layout.jlhtml"
#> tags = ["module1"]
#> title = "Lasso and Elastic Net"

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

# ╔═╡ d5f7245a-ee3e-4390-ba98-417cf941c92d
using HDF5

# ╔═╡ db63cba8-f6f7-49e5-852c-9de9ef5fc8c7
using PlutoPlotly, LinearAlgebra, PlutoUI, GLMNet, Statistics, StatsBase

# ╔═╡ c8057eec-5255-4aad-ab42-610488f9ea9e
TableOfContents()

# ╔═╡ c3661798-b931-11ee-21dc-7f6ee5fa1562
md"""
# Gravity Inversion (Lasso and Elastic Net)

In this exercise, we are interested in the numerical approximation of the gravity of an extended body. The gravitational potential at a point $P$ external to that body is given by Eq. (3.13),

```math
U_g = -G\, \int_x \int_y \int_z \frac{\rho(x,y,z)}{r(x,y,z)}\,dx\,dy\,dz\,,
```

where $G$ is the gravitational constant, $\rho(x,y,z)$ is the density at some location $(x,y,z)$, and $r(x,y,z)$ is the distance between $(x,y,z)$ and the observation point $P$.

To numerically solve the integral, we discretize the spatial domain into evenly-spaced grid points $(x_i, y_i, z_i)$, separated by some small distance $h$. We then obtain an approximation of $U_g$ by summing over all grid points,

```math
U_g \approx -G\,\sum_i \frac{\rho(x_i,y_i,z_i)}{r(x_i,y_i,z_i)}\,V\,,
```

where the small volume $V$ is given by $V=h^3\approx dx\,dy\,dz$. The product $\rho(x_i,y_i,z_i) V$ equals the mass $m_i$ contained in the small volume $V$. Renaming $r_i=r(x_i,y_i,z_i)$, we may rewrite the above equation as

```math
U_g \approx -G\,\sum_i \frac{m_i}{r_i}\,,
```

which is identical to equation (3.12) in the text.

An approximation that we will make throughout this exercise is that the Earth is contained in a small rectangular box. Though this is obviously not realistic, it greatly simplifies the calculations, while still illustrating the basic principles.
"""

# ╔═╡ 993c8d60-3146-438d-b2b6-c38195baf980
md"""
Elastic Net α $(@bind alphaselect Slider(range(0, 1, length=101), show_value=true, default=0.9))
"""

# ╔═╡ c2c1d192-8848-465e-a2d7-7260cb4e46fe
md"""
Constrain on density? 
$(@bind constraint CheckBox())
with maximum $(@bind maximum_density Slider(show_value=true, range(3000, 10000, length=10)))
"""

# ╔═╡ 4ce52e14-0a40-4cd9-82f1-f95ed2e4c938
@bind invG_type Select(["Pseudo Inverse", "Ridge Regression", "Lasso", "Elastic Net"])

# ╔═╡ 8caa697e-3c72-4d75-819f-a8907bef3fb7
md"## Error Bowls"

# ╔═╡ 74671a5e-4aec-4a51-a759-420d0b0cddf1
md"""

The objective function for Elastic Net proposed by Zou and Hastie (2005)
 is given by
```math
L(m, \lambda, \alpha) = \frac{1}{2}\|Gm-d\|_2^2 + P(m, \lambda, \alpha),
```where
```math
P(m, \lambda, \alpha) = \lambda\left((1-\alpha)\frac{1}{2}\|m\|^2_2 + \alpha \|m\|_1\right).
```
The hyperparameter to be tuned in the Elastic Net is the value for $\alpha$ where, $\alpha \in [0, 1]$
"""

# ╔═╡ 3c3751b1-3236-4ae2-ac4b-90e01076d920


# ╔═╡ 21ed6a59-7a4b-434f-9432-4047ca9520f7


# ╔═╡ e79fc893-233f-4e71-aac7-2c259d57cfcd
md"""
In addition, Zou and Hastie (2005) showed that L1–L2 norm regularization encourages a “grouping effect”, which means that the elements of the solution corresponding to the highly correlated columns of X have similar estimated values. By this enforcing of the grouping effect, the second drawback is mitigated. For the details of this point, please refer to section 2.3 of Zou and Hastie

As the result of L1–L2 norm regularization, the effect that the L1 norm regularized solution is overly concen- trated is resolved by introducing the L2 norm constraint, and at the same time, the sparse nature of the model is also provided by the L1 norm constraint
"""

# ╔═╡ 492e7b9a-8156-476f-aba8-a87e4eef9bea
md"## Physics"

# ╔═╡ 6d0655de-3a68-4990-bb3d-ffe2b2eba4fe
gravity_const=6.67508e-11

# ╔═╡ a30a6d2c-6057-4695-8b20-88a4163227d6
begin
	# Dimension of the computational domain [km].
	x_min=-50.0
	x_max=50.0
	z_min=0
	z_max=50
	
	nx=201
	nz=101
	
	# Coordinate axes.
	xgrid=range(x_min,stop=x_max,length=nx)
	dx=step(xgrid)
	zgrid=range(z_min,stop=z_max,length=nz)
	dz=step(zgrid)
	# Grid spacing [m] and cell volume
	A = step(xgrid)*step(zgrid)
	X = first.(Iterators.product(xgrid, zgrid))
	Z = last.(Iterators.product(xgrid, zgrid))
	# xv,zv=meshgrid(xgrid,zgrid)
	
	# # Define some density distribution.
	mtrue=zeros(nz, nx)
	mtrue[50:80,65:70].=10.0
	# mtrue[70:80,160:190].=5500.0
	mtrue[1:20,140:160].=10.0
	# mtrue[40:60,10:30].=5500.0
end;

# ╔═╡ f039c6e6-f5c4-42eb-95fc-586b711822f4
md"""
$(@bind x_bowl1 Slider(xgrid, show_value=true, default=mean(xgrid)))
$(@bind z_bowl1 Slider(zgrid, show_value=true, default=mean(zgrid)))
"""

# ╔═╡ e2da835d-fcdf-4a56-a632-362c4d62eecb
md"""
$(@bind x_bowl2 Slider(xgrid, show_value=true, default=mean(xgrid)))
$(@bind z_bowl2 Slider(zgrid, show_value=true, default=mean(zgrid)))
"""

# ╔═╡ 51f1666b-374b-43bb-9ac7-613837eddc77
mtrue

# ╔═╡ 2060cd03-e036-4acc-b54a-7fb56b068772
# Plot density distribution.


# ╔═╡ 4de006ac-c4a4-4b0a-b40a-ac9639af66cb
begin
	# Define observation points.
	xobs=xgrid
	zobs=-10.0.*ones(length(xobs))
	
	# Initialize gravitational potential.
	# U=zeros(nx)
	
	# # Loop over all observation points.
	# for k in range(length(x_obs))
	    
	#     r=np.sqrt((x_obs[k]-xv)^2 + (z_obs[k]-zv)^2)
	#     U[k]=-G*V*np.sum(rho/r)
	# end
end

# ╔═╡ cd60baa7-275f-46cc-95b9-9027aa10e828
begin
	G=mapreduce(hcat, xobs, zobs) do x, z
	r = @. sqrt(abs2(X-x) + abs2(Z-z)) * 1e3
	@. -gravity_const*A/r
	end
end

# ╔═╡ 5bb647b8-422b-4656-8df6-05bac55be032
function plot_error_bowl(
)

	ix1 = argmin(abs2.(xgrid .- x_bowl1))
	iz1 = argmin(abs2.(zgrid .- z_bowl1))
	ix2 = argmin(abs2.(xgrid .- x_bowl2))
	iz2 = argmin(abs2.(zgrid .- z_bowl2))

	m=zeros(nz, nx)
	m1=zeros(nz, nx)
	m[iz1, ix1] = 500.0
	m[iz2, ix2] = 500.0
	d=G*vec(m)
	rho_values=range(0, 1000, length=50)
	rho_iterator = Iterators.product(rho_values, rho_values)
	bowl=map(first.(rho_iterator), last.(rho_iterator)) do m11, m22
		fill!(m1, 0.0)
		m1[iz1, ix1] = m11
		m1[iz2, ix2] = m22
		abs2(norm(G*vec(m1)-d, 2))
	end
	@show size(bowl)
	
plot(
    surface(
        contours = attr(
            x=attr(show=true, start= 1.5, size=0.04, color="white"),
            x_end=2,
            z=attr(show=true, start= 0.5, size= 0.05),
            z_end=0.8
        ),
        x=rho_values, y=rho_values,
        z = bowl
    ),
    Layout(
        scene=attr(
            xaxis_nticks=20,
            zaxis_nticks=4,
            camera_eye=attr(x=0, y=-1, z=0.5),
            aspectratio=attr(x=1, y=1, z=0.2)
        )
    )
)
	# plot(heatmap(z=bowl))
end

# ╔═╡ 6ac71a81-5935-4e5a-a635-cf05ced40c6d
plot_error_bowl()

# ╔═╡ 16d73bbb-581a-4883-be01-4e0941f4264b
data_observed=G*vec(mtrue)

# ╔═╡ a9a203e4-484a-47c8-beed-eeab49a72c51
h5open("gravity_data.h5", "w") do file
    write(file, "G", G) 
	 write(file, "d", data_observed)
	# alternatively, say "@write file A"
end

# ╔═╡ fb0433b4-83fe-4b2d-95f4-4ac812bf0768
plot(data_observed)

# ╔═╡ 4b227df5-711c-437b-88ed-ff02f060c2a2
Gpinv = pinv(G);

# ╔═╡ e195c809-d359-4a3d-ba0b-4926fbd841e0
@bind λ Slider(0.1:0.9)

# ╔═╡ 27e69011-3df9-4777-86fb-0412bdea76c2
qrG = qr(G)

# ╔═╡ 90cff54a-6967-4d81-a975-d400b3bb169f
G'*G

# ╔═╡ 6202aaa8-09ae-4202-8c4c-fb2841531eb4


# ╔═╡ 0f715e49-2311-438a-a5dd-1dff5310ad4a
F = svd(G)

# ╔═╡ 51193cf3-2464-48aa-b0a0-633dd05dbf72
plot(F.S)

# ╔═╡ f3cb239c-724b-4c4a-b735-fd5ec90cc614


# ╔═╡ 77ea2f3e-481a-4de4-aaf6-a5a8383bc721
md"## Inversion"

# ╔═╡ 4b9f074e-b263-46cb-8d4e-7fa5e7d5eda6
md"### Constraints"

# ╔═╡ 681f2718-c908-4e90-90ae-649568b8f73f
constraints = if(constraint)
	cat(zeros(nx*nz)', maximum_density*ones(nx*nz)', dims=1)
else
	cat(fill(-Inf, nx*nz)', fill(Inf, nx*nz)', dims=1)
end

# ╔═╡ 1e66f28f-5862-4532-a4a5-1abe1c39e44e
alphaused = if(invG_type == "Lasso")
	1.0
elseif(invG_type == "Ridge Regression")
	0.0
else
	alphaselect
end

# ╔═╡ 8d15596b-6beb-41cb-8c40-005d348ee9ca
path = glmnet(G, data_observed, alpha=alphaused, constraints=copy(constraints))

# ╔═╡ 92a0a1b3-31b1-4c1b-a005-44d288b0266a
md"""
Select Lasso Path λ $(@bind lambdaselect Slider(path.lambda, show_value=true, default=last(path.lambda)))
"""

# ╔═╡ f549df45-5ef8-4cd2-a5ea-8e730bcbfeb6
result = if(invG_type == "Pseudo Inverse")
	mhat = pinv(G) * data_observed
	data_predicted = G*mhat
	(; mhat, data_predicted, data_observed)
else
	mhat = path.betas[:,findfirst(x->x==lambdaselect, path.lambda)]
	data_predicted=GLMNet.predict(path, G)[:,findfirst(x->x==lambdaselect, path.lambda)]
	(; mhat, data_predicted, data_observed)
end

# ╔═╡ 1ee66663-8e88-431b-bbf5-32b17e4dac34
md"## Appendix"

# ╔═╡ 13ac643c-6e5c-4c5b-a7ea-e2748a5d953a
md"## Plots"

# ╔═╡ 475058ea-7e24-46a1-9b10-30cdb984b1ef
function plot_densitymap(rho, title="density distribution")
	R=reshape(rho, length(zgrid), length(xgrid))
	plot(heatmap(x=xgrid, y=zgrid, z=R), Layout(width=500, height=300, title=string(title, " [kg/m³]"), xaxis=attr(title="x [km]"), yaxis=attr(autorange="reversed",title="z [km]")))
end

# ╔═╡ 73f9f906-e9aa-4428-8568-44b5180a6e16
PlutoUI.ExperimentalLayout.vbox([plot_densitymap(mtrue, "true mass density"), plot_densitymap(result.mhat, "predicted mass density")])

# ╔═╡ 60928cf4-5eb5-4897-8cba-515114559e8c
function plot_data(data, title="data", names=fill("", size(data,2)))
	D=cat(data,dims=2)
	plot([scatter(x=xgrid, y=d, name=name) for (d, name) in zip(eachslice(D, dims=2), names)], Layout(width=500,height=300, title=title, xaxis=attr(title="x [km]"), yaxis=attr(title="gravitational potential")))
end

# ╔═╡ 746640db-eda8-4638-ad2b-c3b0d8f25f79
plot_data(hcat(result.data_predicted, data_observed), "observed vs predicted data", ["predicted", "observed"])

# ╔═╡ Cell order:
# ╠═c8057eec-5255-4aad-ab42-610488f9ea9e
# ╟─c3661798-b931-11ee-21dc-7f6ee5fa1562
# ╟─92a0a1b3-31b1-4c1b-a005-44d288b0266a
# ╟─993c8d60-3146-438d-b2b6-c38195baf980
# ╟─c2c1d192-8848-465e-a2d7-7260cb4e46fe
# ╟─4ce52e14-0a40-4cd9-82f1-f95ed2e4c938
# ╟─746640db-eda8-4638-ad2b-c3b0d8f25f79
# ╠═73f9f906-e9aa-4428-8568-44b5180a6e16
# ╠═d5f7245a-ee3e-4390-ba98-417cf941c92d
# ╠═a9a203e4-484a-47c8-beed-eeab49a72c51
# ╟─8caa697e-3c72-4d75-819f-a8907bef3fb7
# ╟─74671a5e-4aec-4a51-a759-420d0b0cddf1
# ╠═6ac71a81-5935-4e5a-a635-cf05ced40c6d
# ╠═f039c6e6-f5c4-42eb-95fc-586b711822f4
# ╠═e2da835d-fcdf-4a56-a632-362c4d62eecb
# ╠═51f1666b-374b-43bb-9ac7-613837eddc77
# ╠═3c3751b1-3236-4ae2-ac4b-90e01076d920
# ╠═5bb647b8-422b-4656-8df6-05bac55be032
# ╠═21ed6a59-7a4b-434f-9432-4047ca9520f7
# ╠═e79fc893-233f-4e71-aac7-2c259d57cfcd
# ╟─492e7b9a-8156-476f-aba8-a87e4eef9bea
# ╠═6d0655de-3a68-4990-bb3d-ffe2b2eba4fe
# ╠═a30a6d2c-6057-4695-8b20-88a4163227d6
# ╠═2060cd03-e036-4acc-b54a-7fb56b068772
# ╠═4de006ac-c4a4-4b0a-b40a-ac9639af66cb
# ╠═cd60baa7-275f-46cc-95b9-9027aa10e828
# ╠═16d73bbb-581a-4883-be01-4e0941f4264b
# ╠═fb0433b4-83fe-4b2d-95f4-4ac812bf0768
# ╠═4b227df5-711c-437b-88ed-ff02f060c2a2
# ╠═e195c809-d359-4a3d-ba0b-4926fbd841e0
# ╠═27e69011-3df9-4777-86fb-0412bdea76c2
# ╠═90cff54a-6967-4d81-a975-d400b3bb169f
# ╠═6202aaa8-09ae-4202-8c4c-fb2841531eb4
# ╠═0f715e49-2311-438a-a5dd-1dff5310ad4a
# ╠═51193cf3-2464-48aa-b0a0-633dd05dbf72
# ╠═f3cb239c-724b-4c4a-b735-fd5ec90cc614
# ╟─77ea2f3e-481a-4de4-aaf6-a5a8383bc721
# ╟─4b9f074e-b263-46cb-8d4e-7fa5e7d5eda6
# ╠═681f2718-c908-4e90-90ae-649568b8f73f
# ╠═1e66f28f-5862-4532-a4a5-1abe1c39e44e
# ╠═8d15596b-6beb-41cb-8c40-005d348ee9ca
# ╠═f549df45-5ef8-4cd2-a5ea-8e730bcbfeb6
# ╟─1ee66663-8e88-431b-bbf5-32b17e4dac34
# ╠═db63cba8-f6f7-49e5-852c-9de9ef5fc8c7
# ╟─13ac643c-6e5c-4c5b-a7ea-e2748a5d953a
# ╠═475058ea-7e24-46a1-9b10-30cdb984b1ef
# ╠═60928cf4-5eb5-4897-8cba-515114559e8c
