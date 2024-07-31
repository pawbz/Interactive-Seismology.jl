### A Pluto.jl notebook ###
# v0.19.36

#> [frontmatter]
#> order = "2"
#> title = "Linear Algebra Primer"
#> date = "2024-01-04"
#> tags = ["welcome"]
#> description = "Refresh your concepts of linear algebra. Thank you Gilbert Strang."
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

# ╔═╡ 04f8c8df-f2c3-430c-b4a8-281edc72fd62
begin
    using LinearAlgebra
    using PlutoUI
    using PlutoTeachingTools
    using Symbolics
    using PlutoPlotly
end

# ╔═╡ 539b7f56-91ef-4416-8f7c-d1f636fc2d92
ChooseDisplayMode()

# ╔═╡ 7533151c-ad17-4bbe-83a2-b779c410ee22
TableOfContents()

# ╔═╡ 92735ad6-8b5f-11ed-103a-b1eeeee21a3f
md"""
# Linear Algebra (Primer)
It is one of the mathematical pillars of our course:
- the measured data is often organized to form a vector $d$; 
- in many problems, the physics is often simplified (linearized) to form a matrix $G$ that maps the model vector $m$ to produce measurements $d$:
```math
d = G * m.
```
In this notebook, we shall try to revisit some highlights of linear algebra. Hopefully, as most of you learned these basic concepts, this notebook will serve as a necessary primer! One needs to develop a good geometrical picture of the concepts. For example, we would like to clearly visualize the least-squares solution
```math
\hat{m} = [G^TG]^{-1}G^T d
```
that you will encounter soon.
"""

# ╔═╡ bee44787-2b0e-4d1f-b952-4b3714263cf2
md"""
## Vectors
"""

# ╔═╡ d1999b89-026e-4fab-bbb2-7556af48644d
data = [0.4, 2.4, 5.6, 9]

# ╔═╡ bce812f5-01e7-4ebb-a244-af10eda4e45e
typeof(data)

# ╔═╡ e9850bb2-b024-487c-a7de-2a998306a852
danger(md"`x=(1,2,3)` will result in a simple list called `Tuple` in Julia, it is different from a `Vector{Float64}`!")

# ╔═╡ 562b6d9d-52f6-4234-915a-8db69d829acd
typeof((0.4, 2.4, 5.6, 9))

# ╔═╡ ac5ef8aa-f58a-49bc-b7b7-7f0784f4a86b
md"It is easy to check the length of arrays."

# ╔═╡ 1dc74682-9bec-4f32-ab15-24bf3df371f0
length(data)

# ╔═╡ 737aef17-3c28-4f5a-9e20-6df43a5973a6
md"""
## Scalar-vector multiplication
"""

# ╔═╡ c348b675-fb25-45ea-ade3-de1a340b7fc9
3.0 * data

# ╔═╡ 3e4c2e52-b1f8-48bf-b589-6f375d76ade8
md"""
## Indexing
"""

# ╔═╡ 60659376-b11e-4a50-8776-0b57d4d74115
data[1:2]

# ╔═╡ c7fd9244-a62d-4692-b1cd-7579f7ebd642
data[end]

# ╔═╡ b3764bcd-1543-464a-acbe-44e498a9942d
first(data)

# ╔═╡ 32f67661-a921-42c7-9533-9989bbb6b414
last(data)

# ╔═╡ b4d99d84-bb07-4a6c-a73d-b43759f2d3d8
md"""
## Concatenate
"""

# ╔═╡ 6908b1c4-c4fb-4771-b793-ebacef676858
data2 = [1.2, 3.4, 4.5]

# ╔═╡ 1c695805-f05f-47a5-aac9-5bc0c85b4104
data_full = vcat(data, data2)

# ╔═╡ 7685c3c2-67b7-491f-9b89-c9101a90d9f5
cat(data, data2, dims=1)

# ╔═╡ 588f453a-a155-472e-b4c0-3fbec130f5e0
data3 = [1, 2, 3]

# ╔═╡ 81926f06-523a-4b4e-aab6-f43499e0c486
vcat(data, data3)

# ╔═╡ c1b5d93e-5f66-4516-9a83-bf051ce7d197
aaa = Float32.(data_full)

# ╔═╡ d943d77b-549b-4b25-b8e0-6caf6f11dc57
typeof(aaa)

# ╔═╡ e7890fb8-15ce-42d5-a0ca-36ed53d24a5d
md"""
## Special arrays
"""

# ╔═╡ 7c90d63d-1a27-4438-a180-f3a882dda2fd
typeof(zeros(5))

# ╔═╡ 198bc38d-802f-42b0-8a0a-44d445cc1542
typeof(zeros(Float32, 5))

# ╔═╡ 0ce334b2-0198-4524-8ebf-da136a0a774f
randn(Float32, 10)

# ╔═╡ 2e290c1c-167e-4a2a-aeb9-5f78b0bcdc3a
md"""
## Linear combination
Two arrays of same size can be added together.
"""

# ╔═╡ 00c23763-59c5-4cbd-b7d3-7ee4e962b680
0.5 .* randn(5) .+ 0.3 * randn(5)

# ╔═╡ 31b9a1e0-b097-45c9-b6e6-10da870addf6


# ╔═╡ 7c52fd8a-d94b-4b37-9e70-13699636abf8
md"""
## Inner product
"""

# ╔═╡ bda26e18-6ec8-4b2d-abc2-0fd17310f96c
v1 = [2.0, 3.4, 4.5]

# ╔═╡ 5842ae5e-e575-47a1-8501-9b74967fb769
v2 = [1.0, 3.4, 4.5]

# ╔═╡ f296a560-3bf9-4b3a-868c-8f102f458179
dot(v1, v2)

# ╔═╡ ae4264a3-1b09-499a-8039-6eaa4c176fab
sum(v1 .* v2) # alternatively

# ╔═╡ a4f2e578-5dcc-4906-bffc-1fd8b7b08710
md"## Matrix-vector multiplication"

# ╔═╡ 36d1292a-f402-4eca-b5a8-367156fbc183
@syms m₁ m₂

# ╔═╡ a074ce6d-f387-457c-b5bf-ca97201c24f1
m = [m₁, m₂]

# ╔═╡ 721ef3b0-bf3e-4844-8a5a-e1825c8066c3
g₁ = [0.2, 0.1, 3]

# ╔═╡ 76a0aa87-edf0-4d7b-adfb-9291ab5e843a
g₂ = [1, 0.2, 0.1]

# ╔═╡ 36781305-9c23-4ac4-ab51-996bd408a90f
G = hcat(g₁, g₂) # simply a collection of vectors (columns)

# ╔═╡ 8d4aa317-a78f-4da9-b9d7-310dd566e89c
d = G * m

# ╔═╡ be6cc32a-b105-4fc5-9f0c-31d79d1b4302
Markdown.MD(Markdown.Admonition("warning", "Intuition",
    [md"""
    The multiplication of $G$ times $m$ can be interpreted in two ways:
    1) linear combination of columns $g_1$ and $g_2$ -- which means the vector $G\,m$ belongs to the column space of $G$;
    2) an inner product between each row of $G$ with the vector $m$ -- roughly, what is the correlation of each row of $G$ with the vector $m$?
    """]))

# ╔═╡ cdd66332-8c01-4f73-ad4e-14f008be6f25
# linear combination of columns
m₁ * G[:, 1] + m₂ * G[:, 2]

# ╔═╡ f25e1410-3ce1-4d44-97c5-ac5c4b40c334
aside(Markdown.MD(Markdown.Admonition("therom", "Think",
    [md"""
    Is there any linear combination, other than $m_1=0$ and $m_2=0$, of $g_1$ and $g_2$ that produces a zero vector? Why not?
    """])))

# ╔═╡ dddfda51-ddf6-4f18-9cae-de4535199255
md"""
This means $G\,m$ lies in the space of $g_1$ and $g_2$. You can get convinced by using the sliders below.

Choose:  m₁ $(@bind x₁p Slider(range(-1, 1, length=10), default=0.5, show_value=true)) 
m₂: $(@bind x₂p Slider(range(-1, 1, length=10), default=0.5, show_value=true))
"""

# ╔═╡ 024a1e91-ead5-4c83-8f7b-b89377a2c08b
# alternatively, inner product with rows
[sum(G[1, :] .* m), sum(G[2, :] .* m), sum(G[3, :] .* m)]

# ╔═╡ 770b008b-60a1-410a-b806-e1d497d457fe
md"""
## Outer product
Results in a matrix with dependent columns on $u$
```math
A = u\,v^{T}.
```
"""

# ╔═╡ fec04cd7-6ff9-4780-9972-2369f2a09b22
@bind resample1 Button("Resample")

# ╔═╡ 7444223d-10fe-42f5-b47b-d5d9ee3669e3
begin
    resample1
    A1 = randn(3) * transpose(randn(4))
end

# ╔═╡ af9d6040-a1fc-425a-8cd6-b63116b40093
md"""## Rank; building blocks of a matrix
The rank of a matrix is the dimension of its column space.
Let's check the rank of the matrix that we created using the outer product. We know that it has only one independent column.
"""

# ╔═╡ ffb4c6e8-5d03-4398-adca-2950ae83cf61
rank(A1)

# ╔═╡ 53a89eb0-2b27-4619-a4ab-4779de73f6e6
md"Now, lets add two rank-1 peices together."

# ╔═╡ 3702860e-3a89-45c4-a95e-7b6f56637d3b
A2 = randn(3) * randn(4)' .+ randn(3) * randn(4)'

# ╔═╡ b5d19dbb-16ab-48b7-8bf5-a51921747d95
rank(A2)

# ╔═╡ 12771e27-ebfe-4f88-a650-f5f9fe1b2a71
md"Obviously, a matrix with 3 rows has full rank if we add at least 3 rank-1 (outer products) together. To show this, we will create a matrix using `for` loop in Julia."

# ╔═╡ fb6b5b19-4c70-467e-bc48-3f5a9f89591c
rank(sum([randn(3) * randn(7)' for _ in 1:10]))

# ╔═╡ 324a4aab-16af-462a-bc89-e6345e6b3147
md"""
Conversely, any matrix can be decomposed into rank-1 pieces! We will talk more about this when working with singular-value decomposition (SVD).
"""

# ╔═╡ b3d467ef-439a-4cdd-a631-55ef93432035
md"""## Matrix-matrix multiplication 
Product of two matrices $C$ and $R$
```math
CR = \text{columns of } C \text{ times rows of } R.
```

There are many ways of factoring $A$. Even though $A$ is a fat matrix, we can simply collect the basis for the column space of $A$ to form $C$. Other left-out dependent columns can be generated later using $R$ to form (first great theorem)
```math
A = C\,R.
```
Two views again:
1) linear combination of columns of $C$ using each column of $R$ (earlier, we discussed this when $R$ has only one column, i.e., it is a vector);
2) or linear combination of rows of $R$ using rows of $C$.
Finally, we can show that the number of independent columns of $C$ equals the number of rows of $R$.
"""

# ╔═╡ 806c20ed-0ed2-4087-b6b3-7bce0c40bd53
md"""
## Singular-value decomposition
"""

# ╔═╡ d097c4fd-6589-4821-8f46-2140f7f451f2
md"""
You are now given a fat matrix $X$, with a bunch of columns, much more than 3 (column size). Obviously, there are only 3 independent columns in X. What are they? Can we order them by importance?
"""

# ╔═╡ 580f263c-f1aa-45bf-946c-dc31d24aa5d9
@bind resample2 Button("Resample")

# ╔═╡ f32df348-4c61-4af8-ba1c-e2f9ad8bfa3e
begin
    resample2
    X = randn(3, 3) * randn(3, 20) # create a fat matrix
end

# ╔═╡ 39300e57-97c6-4925-a577-9c02d05ddbec
md"""
Choose the number of singular vectors for plotting
$(@bind ns Slider(1:3, show_value=true, default=1))
"""

# ╔═╡ 736ee1cf-d774-4bda-a87c-4231661732bb
sX = svd(X)

# ╔═╡ 4309e053-291a-4616-8596-b8c2cbd0c22d
md"Let's only pick the first rank-1 piece"

# ╔═╡ 4e7d9c26-c2bf-41c8-8be1-f8c4211606cf
norm(X .- (sX.U[:, 1] * sX.V[:, 1]') .* sX.S[1])

# ╔═╡ 068a5dd5-f13c-4c32-bd46-5e133a554761
md"If we pick the second piece, then the distance is higher."

# ╔═╡ 3065f6ec-b140-4e90-a6e0-2e0623be80a4
norm(X .- (sX.U[:, 2] * sX.V[:, 2]') .* sX.S[2])

# ╔═╡ 83bbe70d-9da4-4bc8-a9ad-2876ac56c61e
md"and so on"

# ╔═╡ b436bf7e-7305-41a9-8434-0bfe60e4b1ab
norm(X .- (sX.U[:, 3] * sX.V[:, 3]') .* sX.S[3])

# ╔═╡ 5c8a7f1a-412e-4987-bac1-3f1db2cc5741
Markdown.MD(Markdown.Admonition("note", "Eckart-Young Theorem", [md"""
The Frobenius Norm of $X$ is equal to the sum of squares of singular values.
```math
||X||^2 = \sigma_1^2 + \sigma_2^2 + \sigma_3^2 + \cdots 
```
"""
]))

# ╔═╡ 894c539b-ce4a-4951-a090-168686b08e5e
md"Let's test this out..."

# ╔═╡ 58a3ad1a-53dc-4c3d-beca-925cf9847ed4
norm(X)^2 # Frobenius norm

# ╔═╡ 79498db3-cbd8-4e0e-9bf9-759951ba02d3
sum(sX.S .^ 2) # sum of squared singular values

# ╔═╡ f56f8acc-fabf-48e7-b2b3-f4e3685224d8
md"Finally, we understand how to arrange singular vectors by importance, at least in the L2 sense."

# ╔═╡ 6f7e15fd-2641-4d58-8937-f60968d76c29
md"""## Nullspace
Now we are interested in a linear combination of the columns of a matrix that generates a zero vector. Other than the trivial case of multiplying all the columns with zero, when is such a non-zero combination possible? Obviously, we cannot combine independent vectors to produce zero. 
"""

# ╔═╡ f5952f79-e5bd-4e8c-b168-7b8a8bcd308a
Markdown.MD(Markdown.Admonition("note", "Counting Law", [md"""
A matrix $G$ has rank $r$; then the equation $G\,m=0$ has $n-r$ independent solutions, 
where $n$ is the number of columns of $G$.
"""
]))

# ╔═╡ 8f36d740-b185-4f10-aa93-1ad514ac0205
md"""
## Appendix
"""

# ╔═╡ 5b3cd4da-82ff-41a6-9121-53574c50ab16
md"### Plots"

# ╔═╡ 0a5717f0-c414-4f7b-9333-1249bb436159
function myquiver(v, names=fill("", length(v)), colors=fill("black", length(v)), title="")
    layout = Layout(
        title=title,
        scene=attr(
			uirevision=1,
            aspectmode="manual", aspectratio=attr(x=1, y=1, z=1),
            xaxis=attr(
                nticks=4,
                range=[-5, 5]
            ),
            yaxis=attr(
                nticks=4,
                range=[-5, 5]
            ),
            zaxis=attr(
                nticks=4,
                range=[-5, 5]
            ),
        ),
    )
    p = [scatter(
        x=[0, vv[1]],
        y=[0, vv[2]],
        z=[0, vv[3]],
        mode="markers+lines",
        marker=attr(
            size=2, symbol=["", "x"],
            color=color,                # set color to an array/list of desired values
            opacity=0.8
        ),
        name=name,
        type="scatter3d"
    ) for (name, color, vv) in zip(names, colors, v)]
    return plot(p, layout)
end

# ╔═╡ b03d542b-64b8-4a99-88e5-a0204277309c
myquiver([g₁, g₂, G * [x₁p, x₂p]], ["g₁", "g₂", "G*m"], ["black", "black", "red"], "Column space of G")

# ╔═╡ 8bbf6130-2ca9-42f1-bcee-3ef5891ffe58
myquiver([A1[:, 1], A1[:, 2], A1[:, 3], A1[:, 4]], ["a₁", "a₂", "a₃", "a₄"], fill("black", 4), "Columns of outerproduct(u, v)")

# ╔═╡ 7c477947-88d5-4146-8663-24a154018096
myquiver(vcat([X[:, i] for i in 1:size(X, 2)], [sX.U[:, i] .* sX.S[i] for i in 1:ns]), fill(nothing, 20 + ns), vcat(fill("black", 20), fill(nothing, ns)), "$ns Orthonormal Singular Vector(s)")

# ╔═╡ 316043ba-b950-486a-a0bb-d84f8a60fb6e
md"""
## Resources
[^Book]:  Linear algebra and learning from data, Gilbert Strang
[^Youtube]: [Essence of linear algebra](https://www.youtube.com/playlist?list=PLZHQObOWTQDPD3MizzM2xVFitgF8hE_ab), 3Blue1Brown
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
PlutoPlotly = "~0.4.4"
PlutoTeachingTools = "~0.2.14"
PlutoUI = "~0.7.54"
Symbolics = "~5.14.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "50ccc8988c18be9c06c5f087a6b12ba2069e6a80"

[[deps.ADTypes]]
git-tree-sha1 = "41c37aa88889c171f1300ceac1313c06e891d245"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.6"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "d7832de8cf7af26abac741f10372080ac6cb73df"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.34.7"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

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
version = "1.0.5+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

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
git-tree-sha1 = "ac67408d9ddf207de5cfa9a97e114352430f01ed"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.16"

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
git-tree-sha1 = "a4532d110ce91bd744b99280193a317310960c46"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.106"

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

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

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
git-tree-sha1 = "38756922d32476c8f41f73560b910fc805a5a103"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.4.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "b211c553c199c111d998ecdaf7623d1b89b69f93"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.12"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

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
version = "0.3.21+4"

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
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
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
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "SparseArrays", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "27ee1c03e732c488ecce1a25f0d7da9b5d936574"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.3.3"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
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
git-tree-sha1 = "116d71e489abc472efa460cfa2bc0ac7cd0bab54"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.12"

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
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "09324a0ae70c52a45b91b236c62065f78b099c37"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.15.2"

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
git-tree-sha1 = "4e17a790909b17f7bf1496e3aec138cf01b60b3b"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.0"
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
version = "1.9.0"

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

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
git-tree-sha1 = "be414bfd80c2c91197823890c66ef4b74f5bf5fe"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "849b1dfb1680a9e9f2c6023f79a49b694fb6d0da"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.5.0"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils"]
git-tree-sha1 = "8d28ebc206dec9e250e21b9502a2662265897650"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.14.1"

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
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

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
# ╠═539b7f56-91ef-4416-8f7c-d1f636fc2d92
# ╠═7533151c-ad17-4bbe-83a2-b779c410ee22
# ╟─92735ad6-8b5f-11ed-103a-b1eeeee21a3f
# ╟─bee44787-2b0e-4d1f-b952-4b3714263cf2
# ╠═d1999b89-026e-4fab-bbb2-7556af48644d
# ╠═bce812f5-01e7-4ebb-a244-af10eda4e45e
# ╟─e9850bb2-b024-487c-a7de-2a998306a852
# ╠═562b6d9d-52f6-4234-915a-8db69d829acd
# ╟─ac5ef8aa-f58a-49bc-b7b7-7f0784f4a86b
# ╠═1dc74682-9bec-4f32-ab15-24bf3df371f0
# ╟─737aef17-3c28-4f5a-9e20-6df43a5973a6
# ╠═c348b675-fb25-45ea-ade3-de1a340b7fc9
# ╟─3e4c2e52-b1f8-48bf-b589-6f375d76ade8
# ╠═60659376-b11e-4a50-8776-0b57d4d74115
# ╠═c7fd9244-a62d-4692-b1cd-7579f7ebd642
# ╠═b3764bcd-1543-464a-acbe-44e498a9942d
# ╠═32f67661-a921-42c7-9533-9989bbb6b414
# ╟─b4d99d84-bb07-4a6c-a73d-b43759f2d3d8
# ╠═6908b1c4-c4fb-4771-b793-ebacef676858
# ╠═1c695805-f05f-47a5-aac9-5bc0c85b4104
# ╠═7685c3c2-67b7-491f-9b89-c9101a90d9f5
# ╠═588f453a-a155-472e-b4c0-3fbec130f5e0
# ╠═81926f06-523a-4b4e-aab6-f43499e0c486
# ╠═c1b5d93e-5f66-4516-9a83-bf051ce7d197
# ╠═d943d77b-549b-4b25-b8e0-6caf6f11dc57
# ╟─e7890fb8-15ce-42d5-a0ca-36ed53d24a5d
# ╠═7c90d63d-1a27-4438-a180-f3a882dda2fd
# ╠═198bc38d-802f-42b0-8a0a-44d445cc1542
# ╠═0ce334b2-0198-4524-8ebf-da136a0a774f
# ╟─2e290c1c-167e-4a2a-aeb9-5f78b0bcdc3a
# ╠═00c23763-59c5-4cbd-b7d3-7ee4e962b680
# ╠═31b9a1e0-b097-45c9-b6e6-10da870addf6
# ╟─7c52fd8a-d94b-4b37-9e70-13699636abf8
# ╠═bda26e18-6ec8-4b2d-abc2-0fd17310f96c
# ╠═5842ae5e-e575-47a1-8501-9b74967fb769
# ╠═f296a560-3bf9-4b3a-868c-8f102f458179
# ╠═ae4264a3-1b09-499a-8039-6eaa4c176fab
# ╟─a4f2e578-5dcc-4906-bffc-1fd8b7b08710
# ╠═36d1292a-f402-4eca-b5a8-367156fbc183
# ╠═a074ce6d-f387-457c-b5bf-ca97201c24f1
# ╠═721ef3b0-bf3e-4844-8a5a-e1825c8066c3
# ╠═76a0aa87-edf0-4d7b-adfb-9291ab5e843a
# ╠═36781305-9c23-4ac4-ab51-996bd408a90f
# ╠═8d4aa317-a78f-4da9-b9d7-310dd566e89c
# ╟─be6cc32a-b105-4fc5-9f0c-31d79d1b4302
# ╠═cdd66332-8c01-4f73-ad4e-14f008be6f25
# ╟─f25e1410-3ce1-4d44-97c5-ac5c4b40c334
# ╟─dddfda51-ddf6-4f18-9cae-de4535199255
# ╠═b03d542b-64b8-4a99-88e5-a0204277309c
# ╠═024a1e91-ead5-4c83-8f7b-b89377a2c08b
# ╟─770b008b-60a1-410a-b806-e1d497d457fe
# ╟─fec04cd7-6ff9-4780-9972-2369f2a09b22
# ╠═7444223d-10fe-42f5-b47b-d5d9ee3669e3
# ╟─8bbf6130-2ca9-42f1-bcee-3ef5891ffe58
# ╟─af9d6040-a1fc-425a-8cd6-b63116b40093
# ╠═ffb4c6e8-5d03-4398-adca-2950ae83cf61
# ╟─53a89eb0-2b27-4619-a4ab-4779de73f6e6
# ╠═3702860e-3a89-45c4-a95e-7b6f56637d3b
# ╠═b5d19dbb-16ab-48b7-8bf5-a51921747d95
# ╟─12771e27-ebfe-4f88-a650-f5f9fe1b2a71
# ╠═fb6b5b19-4c70-467e-bc48-3f5a9f89591c
# ╟─324a4aab-16af-462a-bc89-e6345e6b3147
# ╟─b3d467ef-439a-4cdd-a631-55ef93432035
# ╟─806c20ed-0ed2-4087-b6b3-7bce0c40bd53
# ╟─d097c4fd-6589-4821-8f46-2140f7f451f2
# ╟─580f263c-f1aa-45bf-946c-dc31d24aa5d9
# ╠═f32df348-4c61-4af8-ba1c-e2f9ad8bfa3e
# ╟─39300e57-97c6-4925-a577-9c02d05ddbec
# ╟─7c477947-88d5-4146-8663-24a154018096
# ╠═736ee1cf-d774-4bda-a87c-4231661732bb
# ╟─4309e053-291a-4616-8596-b8c2cbd0c22d
# ╠═4e7d9c26-c2bf-41c8-8be1-f8c4211606cf
# ╟─068a5dd5-f13c-4c32-bd46-5e133a554761
# ╠═3065f6ec-b140-4e90-a6e0-2e0623be80a4
# ╟─83bbe70d-9da4-4bc8-a9ad-2876ac56c61e
# ╠═b436bf7e-7305-41a9-8434-0bfe60e4b1ab
# ╟─5c8a7f1a-412e-4987-bac1-3f1db2cc5741
# ╟─894c539b-ce4a-4951-a090-168686b08e5e
# ╠═58a3ad1a-53dc-4c3d-beca-925cf9847ed4
# ╠═79498db3-cbd8-4e0e-9bf9-759951ba02d3
# ╟─f56f8acc-fabf-48e7-b2b3-f4e3685224d8
# ╟─6f7e15fd-2641-4d58-8937-f60968d76c29
# ╟─f5952f79-e5bd-4e8c-b168-7b8a8bcd308a
# ╟─8f36d740-b185-4f10-aa93-1ad514ac0205
# ╠═04f8c8df-f2c3-430c-b4a8-281edc72fd62
# ╟─5b3cd4da-82ff-41a6-9121-53574c50ab16
# ╠═0a5717f0-c414-4f7b-9333-1249bb436159
# ╟─316043ba-b950-486a-a0bb-d84f8a60fb6e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
