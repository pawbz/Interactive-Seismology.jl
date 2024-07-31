### A Pluto.jl notebook ###
# v0.19.36

#> [frontmatter]
#> chapter = "1"
#> title = "Pseudoinverse"
#> tags = ["module1"]
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ 18b66ba4-7e17-4106-a7f3-ad5731cdedca
md"# Pseudoinverse"

# ╔═╡ 764e3ce8-aa61-4e77-abb0-30ead55763bd
md"""
## $G$ has independent columns
then  $G^+=(G^TG)^{-1}G^T$ is the pseudoinverse
"""

# ╔═╡ 20693b41-90c2-402b-9aea-d94e585e5e8a
md"""
## $G$ has independent rows
then $G^+=G^T(GG^T)^{-1}$ is the pseudo inverse
"""

# ╔═╡ 5cefbdff-b384-4292-803f-a47483f1e47c
md"""
## Diagonal Matrix
For a diagonal matrix $\Sigma$, it is inverted where possible, otherwise $\Sigma^+$ has zeros
"""

# ╔═╡ 640fcda7-8e0b-4558-80f3-3cbc9ee60005
md"""
## SVD
SVD of a matrix $G$ leads to its pseudoinverse $G^+$
* The pseudoinverse of $G=U\Sigma V^T$ is $G^+=V\Sigma^+ U^T$
"""

# ╔═╡ Cell order:
# ╟─18b66ba4-7e17-4106-a7f3-ad5731cdedca
# ╟─764e3ce8-aa61-4e77-abb0-30ead55763bd
# ╟─20693b41-90c2-402b-9aea-d94e585e5e8a
# ╟─5cefbdff-b384-4292-803f-a47483f1e47c
# ╟─640fcda7-8e0b-4558-80f3-3cbc9ee60005
