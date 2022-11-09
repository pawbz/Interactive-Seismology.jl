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

# ╔═╡ dcaeb6a8-78d1-11ec-24fb-4509de0a7d7b
begin
    using Pluto
    using PlutoPlotly
    using PlutoTest
    using PlutoTeachingTools
    using LinearAlgebra
    using Distances
    using StatsBase
    using PlutoUI
end

# ╔═╡ 7e72f1fc-345a-4a2d-b03b-8a7549ef6efc
ChooseDisplayMode()

# ╔═╡ 35c21158-fc55-45f6-930d-7b82c2c0685d
TableOfContents()

# ╔═╡ d9d53d21-09ee-47cd-b661-8787de32f2c1
md"""
# Traveltime Tomography
Seismic tomography is a technique for imaging the subsurface of the Earth with seismic waves produced by earthquakes or explosions. P-, S-, and surface waves can be used for tomographic models of different resolutions based on seismic wavelength, wave source distance, and source-receiver geometry. Tomography is solved as an inverse problem. We analyze the 2-D inversion of first-arrival traveltimes in this notebook. 

By interacting with this notebook, we aim to

- visualize the tomogram for different source-receiver configurations;
- understand ill-posedness and the need for regularization;
- understand the concept of backpropagation;
- and more.

##### [Introduction of Seismology](https://pawbz.github.io/ES218.jl/)
ES218; August 2022

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ d6787951-065d-4eac-8c73-c30295f570ca
md"""
### Acquiring Traveltimes
Source and receivers are equispaced on $x=-1000\,$m and $x=1000\,$m, respectively. Slide to choose their number.

Number of sources: $(@bind ns Slider(5:2:50, show_value=true, default=25))

Number of receivers: $(@bind nr Slider(5:2:50, show_value=true, default=25))
"""

# ╔═╡ f6902c2c-4489-4ce9-8a19-86df8b16a81d
md"""
### Model Parameters
The true model is used to generate (synthetic) observations, which are then inverted to obtain the maximum likelihood estimate.

#### Resolution
"""

# ╔═╡ 58e651a0-2d41-40d2-8ab4-12f0d50ea09b
md"""
Slide to adjust the spatial resolution (m) of the true model.
$(@bind ds Slider(range(10,stop=200), default=20, show_value=true))
"""

# ╔═╡ 199e4069-dff0-4205-aecc-65f984f043d5
md"""
Slide to adjust the spatial resolution (m) of the inverted model.
$(@bind ds_inv Slider(range(10,stop=200), default=150, show_value=true))
"""

# ╔═╡ 63ed3537-754f-4080-8cf3-bd6e4f6b2a8c
md"We shall now get the source and receiver positons, plot them, based on the user input."

# ╔═╡ f9089736-3744-4382-8bc0-68fc04b3cddb
md"""
## Forward Problem
We need to generate the _observed_ data. Towards that end, we shall construct the (straight-ray) forward operator `G`. It depends on the model grid (which is fine) and source-receiver setup. The forward problem is `G * s`, where s is the model (slowness) vector.
"""

# ╔═╡ 1c18907e-f82c-445d-9383-6ecc12771f73
md"We shall now generate true seismic velocity and slowness arrays based on the use input."

# ╔═╡ 59726374-ff49-42cb-9002-efda78a38afa
md"_Observed_ data."

# ╔═╡ 924555b0-5ddc-4ee5-ad32-4696f9fb47e4
md"""
## Inverse Problem
Here, we shall simply use the Moore-Penrose pseudoinverse to map the observed travel times to slowness estimates. Note that the forward operator will now be generated again, based on the user input on spatial sampling of the seismic slowness, which is usually smaller compared to what we have used for generating the observed data above. 
"""

# ╔═╡ 1bc87b69-a59c-4074-b57d-011bf5e3df73
md"Finally, we will take inverse of slowness and reshape it produce the final 2-D seismic velcoity model."

# ╔═╡ df1f716a-f961-4f6b-821c-e1038f190449
md"## Uncertainity"

# ╔═╡ b6300fcc-8073-4487-80a2-8739123c1469
md"Check if the observations are satisfied by the estimated model."

# ╔═╡ 8fb45733-6b2b-428a-9536-fe6b6e2f2aa3
md"""
As a sanity check, we will now test the operator `G` by computing travel times in a homogeneous medium analytically. This will make sure that raytracing is done correctly.
"""

# ╔═╡ 792e9b54-3438-4338-913c-190565d38029
md"### Backpropagation"

# ╔═╡ 54b95769-af61-46f8-a49b-e288831baec7
md"Gradient is formed by backpropagation of residuals. Let's plot the gradient of the squared euclidean distance w.r.t. to the slowness vector at the homogeneous model."

# ╔═╡ 99fed06a-37bd-452a-b5b6-9cf552da37e3
md"Visualize forward operator."

# ╔═╡ 1547ff9a-0fa8-4295-94dd-73bd3678129f
# get inverse operator using SVD (damped least squares)
# also returns the singular values for viz.
# commented for now
function get_inverse_operator(G)
    U, S, V = svd(G)
    Gi = V * Diagonal(inv.(S .+ 1e-1)) * U'
    return Gi, S
end

# ╔═╡ 27465f25-d6c7-4855-a20f-1142b8cd3e9f
md"[Scree plot](https://en.wikipedia.org/wiki/Scree_plot#:~:text=In%20multivariate%20statistics%2C%20a%20scree,principal%20component%20analysis%20(PCA).)"

# ╔═╡ 70e78f3d-a261-43b1-a590-966c7c96021c
md"Slider to choose the row of the Hessian matrix."

# ╔═╡ ae30bfd8-6b42-4aa7-90f2-e7303b359b94
md"""
We will now vis. rows of data and model resolution matrices.
"""

# ╔═╡ 1807eb3a-ce0b-46fe-8c70-fa4af3d9ebad
@bind irowd Slider(1:nr*ns, show_value=true)

# ╔═╡ 010a12e2-1abc-4471-a81b-005c30578e63
md"## Appendix"

# ╔═╡ 442255bc-4d49-4602-b0d4-a935871a9fe8
# define a for modelling and inversion

# ╔═╡ da873791-517d-4ac3-80f8-ceae5808be24
begin
    xgrid = range(-1000, stop=1000, length=floor(Int, 2000 / ds))
    zgrid = range(-1000, stop=1000, length=floor(Int, 2000 / ds))
end;

# ╔═╡ 322d1562-2197-4131-bd17-93aed063e55c
begin
    xgrid_inv = range(-1000, stop=1000, length=floor(Int, 2000 / ds_inv))
    zgrid_inv = range(-1000, stop=1000, length=floor(Int, 2000 / ds_inv))
end;

# ╔═╡ fe682b63-2a06-4c32-8dc0-2f99ba48a873
@bind irowm Slider(1:(length(xgrid_inv)-1)*(length(zgrid_inv)-1), show_value=true, default=div(length(xgrid_inv) * length(zgrid_inv), 4))

# ╔═╡ cf99206b-8a78-40de-bd0d-20bb37ec0b09
md"""
`get_medium` generates a medium with a background velocity of 2000 m/s,
and a circular perturbation. It outputs the velocity matrix and the corresponding slowness vector.
"""

# ╔═╡ 2bf78ef8-8fc5-4e0b-a7c0-f72757bae6f6
function get_medium(xgrid, zgrid, pert=nothing)
    nx = length(xgrid)
    nz = length(zgrid)
    ctrue = 2000.0 * ones(nz - 1, nx - 1)
    if (!(pert === nothing)) # check if perturbation is required
        for iz in 1:nz, ix in 1:nx
            if (sqrt(sum(abs2.([xgrid[ix] - pert.x, zgrid[iz] - pert.z]))) < pert.r)
                ctrue[iz, ix] += 100
            end
        end
    end
    strue = vec(inv.(ctrue))
    return ctrue, strue
end

# ╔═╡ 3d1cd3b5-66f5-44a8-805a-496e801be858
md"""
### Methods for straight-ray tracing
"""

# ╔═╡ 1040f2b8-3999-46fb-9e3d-c163264a4f8a
@warn "The straight-ray implementation in this notebook requires sources and receivers on the outer edge of the 2-D grid, i.e., arbitrary positions are not allowed."

# ╔═╡ f71dbf67-c6c2-444b-acf4-0569ee85bc6b
# careful before changing these, only place sources and receivers on the outer edge of the grid
function get_source_receivers_outer_edge(xgrid, zgrid)
    srcz = (ns == 1) ? [sample(zgrid)] : range(zgrid[1], stop=zgrid[end-1], length=ns)
    srcx = fill(xgrid[1], ns)
    recz = (nr == 1) ? [sample(zgrid)] : range(zgrid[1], stop=zgrid[end-1], length=nr)
    recx = fill(xgrid[end], nr)
    return srcz, srcx, recx, recz
end

# ╔═╡ aa9782b5-88be-43a2-b1e1-d68f289a8fec
srcz, srcx, recx, recz = get_source_receivers_outer_edge(xgrid, zgrid);

# ╔═╡ 565af43c-8b85-4ab4-b72d-ac9560efd4fc
# find intersect b/w two segments
# returns nothing if 
# * segments are parallel 
# * segments do not intersect
# else returns 
# * the point of intersection (including C and D i.e., h==0 or h==1)
function find_intersect(A, B, C, D)
    E = B - A
    F = D - C
    P = [-1 * E[2], E[1]]
    if ((dot(F, P) == 0)) # when parallel
        return nothing
    end
    h = dot(A - C, P) * inv(dot(F, P))
    if (((h * (h - 1)) <= 0) | (h ≈ 0) | (h ≈ 1)) # h should be in (0, 1) or h==0 or h==1
        return C + h * F
    else
        return nothing
    end
end

# ╔═╡ 3e2460ec-102a-4d3c-a4c5-5d6c6e2193ec
function find_length_in_cell(xmin, xmax, zmin, zmax, A, B)
    P1 = find_intersect(A, B, [xmin, zmin], [xmin, zmax])
    P2 = find_intersect(A, B, [xmin, zmin], [xmax, zmin])
    P3 = find_intersect(A, B, [xmin, zmax], [xmax, zmax])
    P4 = find_intersect(A, B, [xmax, zmin], [xmax, zmax])

    P = filter(x -> !(x === nothing), [P1, P2, P3, P4])
    P = unique(P)
    if (length(P) >= 2)
        # remove sides that will be otherwise counted twice by adjacent cells
        if (([xmin, zmax] ∈ P) & ([xmax, zmax] ∈ P))
            return 0
        elseif (([xmax, zmin] ∈ P) & ([xmax, zmax] ∈ P))
            return 0
        else
            return maximum(broadcast(x -> sqrt(sum(abs2.(x))), diff(P)))
        end
    else
        return 0
    end
end


# ╔═╡ 82d3a20f-ea2b-47e6-96df-45fa568da9f8
function get_forw_operator(xgrid, zgrid, srcx, srcz, recx, recz)
    nx = length(xgrid)
    nz = length(zgrid)
    ns = length(srcx)
    nr = length(recx)
    G = zeros(ns * nr, (nx - 1) * (nz - 1))
    for is in 1:ns, ir in 1:nr, iz in 1:nz-1, ix in 1:nx-1
        G[ir+(is-1)*nr, iz+(ix-1)*(nz-1)] = find_length_in_cell(xgrid[ix], xgrid[ix+1], zgrid[iz], zgrid[iz+1], [srcx[is], srcz[is]], [recx[ir], recz[ir]])
    end
    return G
end

# ╔═╡ 7df1cd87-40fa-45c1-9d85-1d491c414a18
Gtrue = get_forw_operator(xgrid, zgrid, srcx, srcz, recx, recz);

# ╔═╡ dce75e41-274b-4e6a-8949-5caaeef7238a
G = get_forw_operator(xgrid_inv, zgrid_inv, srcx, srcz, recx, recz);

# ╔═╡ 1f66502a-bfbe-407f-b869-142f446dfdf6
Gi = pinv(G); # compute Moore-Penrose pseudoinverse

# ╔═╡ 58867c8d-af21-48e9-ab0c-4472711e8eb0
begin
    tt_analytic = vec([(sqrt(sum(abs2.([srcx[is] - recx[ir], srcz[is] - recz[ir]])))) * inv(2000) for ir in 1:nr, is in 1:ns])
    tt_G = G * inv.(fill(2000, (length(zgrid_inv) - 1) * (length(xgrid_inv) - 1)))
    @test tt_analytic ≈ tt_G
end

# ╔═╡ ad441089-505f-4da2-a345-548e8c4dd7d2
plot(heatmap(z=G), Layout(xaxis_title="model vector index", yaxis_title="data vector index", width=450))

# ╔═╡ acc8d4d0-a332-478f-8630-b22a10e7063b
function plot_G_scree()
    s = svd(G)
    plot(s.S, Layout(title="Singular values of G"))
end

# ╔═╡ aa36d8d7-a7d4-4aa4-b0ee-b07d36cc453b
plot_G_scree()

# ╔═╡ 3c7cddb2-72c6-45d1-a902-f66cb67d2835
plot(heatmap(z=transpose(G) * G), Layout(title="the Hessian matrix", xaxis_title="model index", yaxis_title="model index", width=450))

# ╔═╡ f408a310-fce3-4876-819e-3457037bd48f
plot(heatmap(x=xgrid_inv, y=zgrid_inv, z=reshape((transpose(G)*G)[irowm, :], length(zgrid_inv) - 1, length(xgrid_inv) - 1)), Layout(title="row of the Hessian matrix", width=450))

# ╔═╡ 5bdd9d83-3911-4ab5-aefb-5ead429ac5a5
plot(heatmap(x=xgrid_inv, y=zgrid_inv, z=reshape((Gi*G)[irowm, :], length(zgrid_inv) - 1, length(xgrid_inv) - 1)),
    Layout(title="Row of Model Resolution Matrix", width=450))

# ╔═╡ a3bf8549-fff8-4423-8a87-a81cb21f9eb1
plot((G*Gi)[irowd, :], Layout(title="Row of data resolution matrix", xaxis_title="raypath index"))

# ╔═╡ 4dd5df1f-f0bc-49ec-a533-4498ed17d223
md"""
### UI
"""

# ╔═╡ 208932c4-a57b-487d-9d3b-f165b4a4a4ed
function perturbation_input(xgrid, zgrid)
    nx = length(xgrid)
    nz = length(zgrid)
    return PlutoUI.combine() do Child
        p = [
            md"""
            x: $(Child("x", Slider(xgrid, default=xgrid[div(nx,2)], show_value=true)))
            """,
            md"""
                     z: $(Child("z", Slider(zgrid, default=zgrid[div(nz,2)], show_value=true)))
                     """
        ]
        r = [
            md"""
            $(Child("r", Slider(range(10, stop=500, length=10), default=330, show_value=true)))
            """,]

        md"""
 #### Perturbation Position
 Slide to adjust the position ∈ [-1000, 1000] m of a circular perturbation.
 $(p)
 #### Perturbation Size
 Slide to adjust the radius ∈ [10, 500] m of the circular perturbation.
 $(r)	   
            """
    end
end

# ╔═╡ 78da0510-bf64-460e-82e7-59ba51c7c7f5
@bind pert perturbation_input(xgrid, zgrid)

# ╔═╡ 092331da-31a4-403d-b0cb-6e7705c6d81b
ctrue, strue = get_medium(xgrid, zgrid, pert);

# ╔═╡ 3d201819-788b-4d90-b3a6-5483fc16ca81
dobs = Gtrue * strue;

# ╔═╡ 329b24c6-fd22-4565-86f6-0881cb11942f
sest = Gi * dobs; # inverse map

# ╔═╡ 68b24556-4731-4b8f-b8a6-7aa7ffa38a92
cest = reshape((inv.(sest)), length(zgrid_inv) - 1, length(xgrid_inv) - 1);

# ╔═╡ a26fd943-dd6c-4e42-b6b0-5c9a17f31b19
plot([scatter(y=dobs, name="observed"), scatter(y=G * sest, name="predicted")], Layout(title="Data Error (Observed Vs. Predicted Traveltimes", xaxis_title="# raypath", yaxis_title="traveltime"))

# ╔═╡ 6152161d-1d6d-4f2c-8ad7-c870458fd34b
plot(heatmap(x=xgrid_inv, y=zgrid_inv, z=reshape(G' * (tt_G - dobs), length(zgrid_inv) - 1, length(xgrid_inv) - 1),), Layout(width=450, title="Gradient w.r.t. slowness"))

# ╔═╡ b338591d-cc11-4e9e-827e-7fcee5b2d38b
md"### Plots"

# ╔═╡ 993e5848-77d4-488d-946d-bfd28b744bcb
function heatmap_model(c, xgrid, zgrid,)
    return heatmap(
        x=xgrid,
        y=zgrid,
        z=c, coloraxis="coloraxis")
end

# ╔═╡ b1350eb1-059e-4f83-a539-2e2befc3dabb
function plot_models()

    fig = Plot(Layout(height=350, Subplots(shared_xaxes=true, shared_yaxes=true, rows=1, cols=2, subplot_titles=["True Earth Model (seismic velocity)" "Estimated Earth Model"])))
    add_trace!(fig, heatmap_model(ctrue, xgrid, zgrid), row=1, col=1)
    add_trace!(fig, heatmap_model(cest, xgrid_inv, zgrid_inv), row=1, col=2)

    return PlutoPlotly.plot(fig)

end

# ╔═╡ 842ec98f-505a-4873-9c64-725e2f92cbb9
plot_models()

# ╔═╡ 1161cc82-6b2d-42e0-99fe-2dfb1a2d00b6
ray_setup = broadcast(reshape(1:ns, 1, ns), reshape(1:nr, nr, 1)) do is, ir
    scatter(x=[srcx[is], recx[ir]], y=[srcz[is], recz[ir]], mode="lines+markers",)
end;

# ╔═╡ 4509e5b8-8d54-47e1-9ba7-b4929fd2d2fc
plot(vec(ray_setup), Layout(showlegend=false, width=450, title="Ray Geometry"))

# ╔═╡ d7dd7859-9489-4bfc-bc77-edec76fe96f2
function plot_ray_setup()
    p = broadcast(reshape(1:ns, 1, ns), reshape(1:nr, nr, 1)) do is, ir
        scatter(x=[srcx[is], recx[ir]], y=[srcz[is], recz[ir]], mode="lines+markers",
            color=:black)
    end
    return plot(vec(p))
end


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distances = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Pluto = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Distances = "~0.10.7"
Pluto = "~0.19.14"
PlutoPlotly = "~0.3.6"
PlutoTeachingTools = "~0.2.3"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.30"
StatsBase = "~0.33.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "27c9c104ad9ea9df3d16bae7f3b2af60ccf00ef0"

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

[[deps.BitFlags]]
git-tree-sha1 = "84259bb6172806304b9101094a7cc4bc6f56dbc6"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.5"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "6e39c91fb4b84dcb870813c91674bdebb9145895"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.5"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "cc4bd91eba9cdbbb4df4746124c22c0832a460d6"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.1.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Configurations]]
deps = ["ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "62a7c76dbad02fdfdaa53608104edf760938c4ca"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.17.4"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

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

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExproniconLite]]
git-tree-sha1 = "2321c9c5a07c2658484dacf8e68e3cd8e2470d5d"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.7.6"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

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

[[deps.FuzzyCompletions]]
deps = ["REPL"]
git-tree-sha1 = "e16dd964b4dfaebcded16b2af32f05e235b354be"
uuid = "fb4132e2-a121-4a70-b8a1-d5b831dcdcc2"
version = "0.5.1"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "8556f4b387fcd1d9b3013d798eecbcfa0d985e66"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.5.2"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

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
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "0f960b1404abb0b244c1ece579a0ec78d056a5d1"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.15"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "eecfbe1bd3f377b7e6caa378392eeed1616c6820"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.2.0"

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
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

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

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "a8cbf066b54d793b9a48c5daa5d586cf2b5bd43d"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "5628f092c6186a80484bfefdf89ff64efdaec552"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.Pluto]]
deps = ["Base64", "Configurations", "Dates", "Distributed", "FileWatching", "FuzzyCompletions", "HTTP", "HypertextLiteral", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "MsgPack", "Pkg", "PrecompileSignatures", "REPL", "RegistryInstances", "RelocatableFolders", "Sockets", "TOML", "Tables", "URIs", "UUIDs"]
git-tree-sha1 = "c3f344a915bc1d67455ecc5e38f4a184ffc4ad96"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.19.14"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "0e8bcc235ec8367a8e9648d48325ff00e4b0a545"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.5"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Colors", "Dates", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "PlotlyBase", "PlutoUI", "Reexport"]
git-tree-sha1 = "dec81dcd52748ffc59ce3582e709414ff78d947f"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.3.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "d8be3432505c2febcea02f44e5f4396fae017503"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.3"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5c0eb9099596090bb3215260ceca687b888a1575"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.30"

[[deps.PrecompileSignatures]]
git-tree-sha1 = "18ef344185f25ee9d51d80e179f8dad33dc48eb1"
uuid = "91cefc8d-f054-46dc-8f8c-26e11d7c5411"
version = "3.0.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "dad726963ecea2d8a81e26286f625aee09a91b7c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.4.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

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
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

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

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

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
# ╠═7e72f1fc-345a-4a2d-b03b-8a7549ef6efc
# ╠═35c21158-fc55-45f6-930d-7b82c2c0685d
# ╠═d9d53d21-09ee-47cd-b661-8787de32f2c1
# ╟─842ec98f-505a-4873-9c64-725e2f92cbb9
# ╟─d6787951-065d-4eac-8c73-c30295f570ca
# ╟─f6902c2c-4489-4ce9-8a19-86df8b16a81d
# ╟─58e651a0-2d41-40d2-8ab4-12f0d50ea09b
# ╟─199e4069-dff0-4205-aecc-65f984f043d5
# ╟─78da0510-bf64-460e-82e7-59ba51c7c7f5
# ╟─63ed3537-754f-4080-8cf3-bd6e4f6b2a8c
# ╠═aa9782b5-88be-43a2-b1e1-d68f289a8fec
# ╠═4509e5b8-8d54-47e1-9ba7-b4929fd2d2fc
# ╟─f9089736-3744-4382-8bc0-68fc04b3cddb
# ╠═7df1cd87-40fa-45c1-9d85-1d491c414a18
# ╟─1c18907e-f82c-445d-9383-6ecc12771f73
# ╠═092331da-31a4-403d-b0cb-6e7705c6d81b
# ╟─59726374-ff49-42cb-9002-efda78a38afa
# ╠═3d201819-788b-4d90-b3a6-5483fc16ca81
# ╟─924555b0-5ddc-4ee5-ad32-4696f9fb47e4
# ╠═dce75e41-274b-4e6a-8949-5caaeef7238a
# ╠═1f66502a-bfbe-407f-b869-142f446dfdf6
# ╠═329b24c6-fd22-4565-86f6-0881cb11942f
# ╟─1bc87b69-a59c-4074-b57d-011bf5e3df73
# ╠═68b24556-4731-4b8f-b8a6-7aa7ffa38a92
# ╟─df1f716a-f961-4f6b-821c-e1038f190449
# ╟─b6300fcc-8073-4487-80a2-8739123c1469
# ╟─a26fd943-dd6c-4e42-b6b0-5c9a17f31b19
# ╟─8fb45733-6b2b-428a-9536-fe6b6e2f2aa3
# ╠═58867c8d-af21-48e9-ab0c-4472711e8eb0
# ╟─792e9b54-3438-4338-913c-190565d38029
# ╟─54b95769-af61-46f8-a49b-e288831baec7
# ╠═6152161d-1d6d-4f2c-8ad7-c870458fd34b
# ╟─99fed06a-37bd-452a-b5b6-9cf552da37e3
# ╠═ad441089-505f-4da2-a345-548e8c4dd7d2
# ╠═1547ff9a-0fa8-4295-94dd-73bd3678129f
# ╟─27465f25-d6c7-4855-a20f-1142b8cd3e9f
# ╠═aa36d8d7-a7d4-4aa4-b0ee-b07d36cc453b
# ╠═acc8d4d0-a332-478f-8630-b22a10e7063b
# ╠═3c7cddb2-72c6-45d1-a902-f66cb67d2835
# ╠═70e78f3d-a261-43b1-a590-966c7c96021c
# ╟─fe682b63-2a06-4c32-8dc0-2f99ba48a873
# ╠═f408a310-fce3-4876-819e-3457037bd48f
# ╟─ae30bfd8-6b42-4aa7-90f2-e7303b359b94
# ╠═5bdd9d83-3911-4ab5-aefb-5ead429ac5a5
# ╠═1807eb3a-ce0b-46fe-8c70-fa4af3d9ebad
# ╠═a3bf8549-fff8-4423-8a87-a81cb21f9eb1
# ╟─010a12e2-1abc-4471-a81b-005c30578e63
# ╠═dcaeb6a8-78d1-11ec-24fb-4509de0a7d7b
# ╠═442255bc-4d49-4602-b0d4-a935871a9fe8
# ╠═da873791-517d-4ac3-80f8-ceae5808be24
# ╠═322d1562-2197-4131-bd17-93aed063e55c
# ╟─cf99206b-8a78-40de-bd0d-20bb37ec0b09
# ╠═2bf78ef8-8fc5-4e0b-a7c0-f72757bae6f6
# ╟─3d1cd3b5-66f5-44a8-805a-496e801be858
# ╟─1040f2b8-3999-46fb-9e3d-c163264a4f8a
# ╠═f71dbf67-c6c2-444b-acf4-0569ee85bc6b
# ╠═565af43c-8b85-4ab4-b72d-ac9560efd4fc
# ╠═3e2460ec-102a-4d3c-a4c5-5d6c6e2193ec
# ╠═82d3a20f-ea2b-47e6-96df-45fa568da9f8
# ╟─4dd5df1f-f0bc-49ec-a533-4498ed17d223
# ╠═208932c4-a57b-487d-9d3b-f165b4a4a4ed
# ╟─b338591d-cc11-4e9e-827e-7fcee5b2d38b
# ╠═993e5848-77d4-488d-946d-bfd28b744bcb
# ╠═b1350eb1-059e-4f83-a539-2e2befc3dabb
# ╠═1161cc82-6b2d-42e0-99fe-2dfb1a2d00b6
# ╠═d7dd7859-9489-4bfc-bc77-edec76fe96f2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
