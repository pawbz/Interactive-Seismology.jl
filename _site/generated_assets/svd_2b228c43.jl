### A Pluto.jl notebook ###
# v0.19.36

#> [frontmatter]
#> chapter = "1"
#> title = "Singular Value Decomposition"
#> tags = ["module1"]
#> layout = "layout.jlhtml"

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

# ╔═╡ 74418960-a9e8-47ab-bf1b-31cdc2e0acb8
using PlutoPlotly, LinearAlgebra, PlutoUI, PlutoTeachingTools

# ╔═╡ 68546785-c79a-4e9b-a717-74120ead9b88
TableOfContents()

# ╔═╡ d2db3192-b5c9-11ee-0d3c-89a9f069d24c
md"""# Singular Value Decomposition (SVD)
```math
G = U\Sigma V^T
```
- SVD produces two sets of singular vectors (columns of matrices `U` and `V`). 
- For a matrix `G` of size $(m,n)$, these singular vectors form orthogonal axes in $\mathbf{R}^m$ and $\mathbf{R}^n$, respectively.
- In this demo, we consider a square matrix `G`: here, the matrix maps each of the columns in `V` to the corresponding columns in `U`.
- Notice that the orthogonal axes of `U` and `V` are similar (except for the sign) in the case when `G` is symmetric.
- If `G` is singular, we can rotate m to be perpendicular to the columns of G to see zero vector d; yes, the null space is orthogonal to the column space of a matrix.
"""

# ╔═╡ 5d11c5f1-3ecd-4264-902f-fa8fbb2d3c6b
ThreeColumn(md"""
$(@bind Atype MultiCheckBox(["Symmetric"=>"G is Symmetric", "Singular"=> "G is Singular"]))""",md"""$(@bind sampleG CounterButton("Resample G"))""", 
	md""" and  Rotate m$(@bind mphi Slider(range(0,2pi,length=1000)))""")

# ╔═╡ 63e61c32-c235-438b-978a-ae1a09fef886
md"Scale the columns of U and V with singular values? $(@bind scaleUV CheckBox())"

# ╔═╡ c712c1ee-ac74-4e3e-b9ed-f73ee505fdf6
md"## Orthogonal Matrices = Rotation in 2D = Invertible"

# ╔═╡ 7d6a1fad-6fb4-4b25-9c8a-9a9f581edfde
md"""
## Symmetric Matrices GᵀG and GGᵀ
`V` contains orthonormal eigenvectors of GᵀG, and `U` contains orthonormal eigenvectors of GGᵀ, where both of them share eigenvalues.
"""

# ╔═╡ 634364b6-05ca-4a54-a0cd-7d0581eaab55
begin
	sampleG; 
	G = randn(2,2)
	if("Symmetric" ∈ Atype)
		G=0.5*(G+G')
	end
	if("Singular" ∈ Atype)
		G[:,2] .= randn() .* G[:,1] 
	end
end

# ╔═╡ f98dab1c-16f1-418c-acda-3fe0b7554c9c
F=svd(G)

# ╔═╡ 6a2e234d-be6c-44a1-8fdd-9e5c6ad8ba01
begin
	m=[cos(mphi), sin(mphi)]
end

# ╔═╡ 75b51cd7-5d3d-4f0c-8750-d169edf0a3d1
n = F.Vt * m

# ╔═╡ cea5f070-70b4-4870-8a41-b0cb4a6ea26b
o = Diagonal(F.S) * n

# ╔═╡ 98338255-a250-4b22-8bc3-b3f11b22583b
m1=F.U*o

# ╔═╡ d94fc2d9-8881-4654-98d3-7fbb78bc1271
md"## Appendix"

# ╔═╡ 1fee9d48-46be-443a-9793-0b84e6de35a1
md"### Plots"

# ╔═╡ 4f9c48af-0714-40f8-a4e6-3228ba0b1393
function quiverplot(p1, p2=zeros(size(p1)); colors=fill("black", size(p1,2)), xylim=3.0, title="", names=string.(collect(1:size(p1,2))))
	vector_scale = 1 #scale factor in (0, 1] for the vector directions to avoid quiver overlapping
    arrow_scale = 0.25
    angle = π/9
    scaleratio = 1.0 #aspect ratio for the 2d plot
    d =1.0 #a scale factor in (0.9, 1]  for the already scaled direction; 
	p1=cat(p1,dims=2)
	p2=cat(p2,dims=2)
	@assert size(p1,2) == size(p2,2)
	plots = mapreduce(vcat, 1:size(p1,2)) do i
    	x = p2[1:1, i]
   	 	y = p2[2:2, i]
   	 	u = p1[1:1, i]
    	v = p1[2:2, i]
    (length(x) == length(y) == length(u) == length(v)) &&
                  vector_scale > 0 && arrow_scale > 0 ||
                  error("the vects x, y, u, v do not have the same length")
    u = vector_scale * scaleratio *u
    v = vector_scale * v
    end_x = x .+ u
    end_y = y .+ v
	function tuple_interleave(tu::Union{NTuple{3, Vector}, NTuple{4, Vector}}) 
   	 #auxilliary function to interleave elements of a NTuple of vectors, N=3 or 4
   	 zipped_data = collect(zip(tu...))
  	  vv_zdata = [collect(elem) for elem in zipped_data]
  	  return reduce(vcat, vv_zdata)
	end
	
    vect_nans = repeat([NaN],  length(x))
    barb_x = tuple_interleave((x, x .+ d*u, vect_nans))
    barb_y = tuple_interleave((y, y .+ d*v, vect_nans))

    barb_length = sqrt.((u/scaleratio) .^2 .+ v .^2)
    arrow_length = arrow_scale #* barb_length
    barb_angle = atan.(v, u/scaleratio)

    ang1 = barb_angle .+ angle
    ang2 = barb_angle .- angle

    seg1_x = arrow_length .* cos.(ang1)
    seg1_y = arrow_length .* sin.(ang1)

    seg2_x = arrow_length .* cos.(ang2)
    seg2_y = arrow_length .* sin.(ang2)

    arrowend1_x = end_x .- seg1_x *scaleratio
    arrowend1_y = end_y .- seg1_y
    arrowend2_x = end_x .- seg2_x *scaleratio
    arrowend2_y = end_y .- seg2_y
    arrow_x =  tuple_interleave((arrowend1_x, end_x, arrowend2_x, vect_nans))
    arrow_y =  tuple_interleave((arrowend1_y, end_y, arrowend2_y, vect_nans))

    barb = scatter(x=barb_x, y=barb_y, mode="lines", line_color=colors[i], name=names[i],showlegend=!(names[i]==""))
    arrow = scatter(x=arrow_x, y=arrow_y, mode="lines", line_color=colors[i],
                     fill="toself", fillcolor=colors[i], hoverinfo="skip", showlegend=false)
	return [barb, arrow]
	end
    layout = Layout(title=title,width=250, height=260, yaxis=attr(range=[-xylim, xylim]),xaxis=attr(range=[-xylim, xylim]), legend=attr(
        x=0.75,
        y=1,
        yanchor="bottom",
        xanchor="right",
        orientation="h"
    ),
                    showlegend=true)
    return plot(Plot(plots, layout, config=PlotConfig(staticPlot=true)))
end

# ╔═╡ 3539521e-d128-4c66-84ed-a85838ea1684
PlutoUI.ExperimentalLayout.hbox([quiverplot(G, title="Columns of G"), quiverplot(m, title="m", colors=["blue"], names=[""])])

# ╔═╡ 0f7e56b7-e71c-4e04-af8c-05dfb424860b
let
	I = scaleUV ? Diagonal(F.S) : Diagonal(ones(2))
	PlutoUI.ExperimentalLayout.hbox([quiverplot(F.U * I, title="Columns of U"), quiverplot(F.V * I, title="Columns of V")])
end

# ╔═╡ af442ee2-d0fc-4754-a1fb-917650237f4b
PlutoUI.ExperimentalLayout.hbox([quiverplot(hcat(m,o,m1), colors=["blue","magenta", "red"], names=["m","o","Gm"]), quiverplot(hcat(n,o), colors=["green", "magenta"], names=["n","o=Σn"]), quiverplot(hcat(m,n), colors=["blue", "green"], names=["m","n=Vᵀm"]), ])

# ╔═╡ 7a04a1a7-8a73-49a0-bb27-4731f8428ef1
quiverplot(G', colors=["blue", "green"], title="Row f G")

# ╔═╡ 82c2ce66-0d86-414d-868a-5d3aff5168dc
quiverplot(hcat(m,F.Vt*m), colors=["blue", "green"], names=["m","Vᵀm"])

# ╔═╡ f1311c00-8a40-4447-9416-8d0f1c512e33
quiverplot(hcat(F.Vt*m, F.Vt*F.V*m), colors=["blue", "green"], names=["m","VᵀVm"])

# ╔═╡ 126aaff1-428e-45e2-a353-f1afc437db7f
md"""### Resources
* [Wolfram Alpha Demo](https://demonstrations.wolfram.com/SingularValueDecomposition/)
"""

# ╔═╡ Cell order:
# ╠═68546785-c79a-4e9b-a717-74120ead9b88
# ╟─d2db3192-b5c9-11ee-0d3c-89a9f069d24c
# ╟─5d11c5f1-3ecd-4264-902f-fa8fbb2d3c6b
# ╟─3539521e-d128-4c66-84ed-a85838ea1684
# ╟─63e61c32-c235-438b-978a-ae1a09fef886
# ╟─0f7e56b7-e71c-4e04-af8c-05dfb424860b
# ╟─af442ee2-d0fc-4754-a1fb-917650237f4b
# ╠═7a04a1a7-8a73-49a0-bb27-4731f8428ef1
# ╟─c712c1ee-ac74-4e3e-b9ed-f73ee505fdf6
# ╠═82c2ce66-0d86-414d-868a-5d3aff5168dc
# ╠═f1311c00-8a40-4447-9416-8d0f1c512e33
# ╟─7d6a1fad-6fb4-4b25-9c8a-9a9f581edfde
# ╠═f98dab1c-16f1-418c-acda-3fe0b7554c9c
# ╠═75b51cd7-5d3d-4f0c-8750-d169edf0a3d1
# ╠═cea5f070-70b4-4870-8a41-b0cb4a6ea26b
# ╠═98338255-a250-4b22-8bc3-b3f11b22583b
# ╠═634364b6-05ca-4a54-a0cd-7d0581eaab55
# ╠═6a2e234d-be6c-44a1-8fdd-9e5c6ad8ba01
# ╟─d94fc2d9-8881-4654-98d3-7fbb78bc1271
# ╠═74418960-a9e8-47ab-bf1b-31cdc2e0acb8
# ╟─1fee9d48-46be-443a-9793-0b84e6de35a1
# ╠═4f9c48af-0714-40f8-a4e6-3228ba0b1393
# ╟─126aaff1-428e-45e2-a353-f1afc437db7f
