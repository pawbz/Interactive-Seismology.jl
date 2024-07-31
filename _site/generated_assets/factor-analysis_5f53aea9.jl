### A Pluto.jl notebook ###
# v0.19.36

#> [frontmatter]
#> chapter = "1"
#> title = "Factor Analysis"
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

# ╔═╡ abe35895-c888-4441-9a4e-94c16c91dc27
using CSV, DataFrames, LinearAlgebra, PlutoPlotly, PlutoUI, Distributions

# ╔═╡ ab1f729f-f7a9-4cf1-b1ef-7286e41bd37a
TableOfContents()

# ╔═╡ 1a7b09b7-e614-4dff-b252-dc726d87af81
md"""# Factor Analysis
The purpose of this notebook is to show the importance of factoring the observed data matrix.
"""

# ╔═╡ 8ac94996-30af-41a9-b2fe-2b9982b07d78
md"## Geomorphology: Mountain Profiles"

# ╔═╡ 520c6e47-6df9-4ec9-8f48-f85ced6cb30b
md"## Atmospheric Pressure Fluctuations"

# ╔═╡ 63bdb39b-19da-4f0f-98f4-9e9506e02927
md"""
Perturbation of atmospheric pressure measured in `(x,y)` plane at time `t`= $(@bind it_atm Slider(1:25, show_value=true))
"""

# ╔═╡ 5fd30c0b-eb53-4a54-a17e-6f12be989d4f
@bind it_atm2 Slider(1:25, show_value=true)

# ╔═╡ 2df4d945-2470-45b1-993b-131decdb2eae
md"## Petrologic Database"

# ╔═╡ 1cc2d822-df8f-4fac-b8ca-770b6d1e09dd
md"We apply factor analysis to rock chemistry data taken from a petrologic database (PetDB at www.petdb.org). This database contains chemical information on igneous and metamorphic rocks collected from the floor of all the world’s oceans, but we analyze here N 1⁄4 6356 samples from the Atlantic Ocean that have the following chemical species: SiO2, TiO2, Al2O3, FeOtotal, MgO, CaO, Na2O, and K2O (units of weight percent)."

# ╔═╡ 53f2cf9e-c1e6-4861-8b66-fcbe260dd6c9
@bind resample_data CounterButton("Resample Data")

# ╔═╡ cc80419d-f255-4ba9-a2a9-d52d2d6fd014
element_names = ["SiO₂", "TiO₂", "Al₂O₃","FeO", "MgO", "CaO", "Na₂O", "K₂O"]

# ╔═╡ fd2938e0-a187-41ef-9fe2-f5f0c262de1e
md"## Appendix"

# ╔═╡ e9950ad0-73d1-4fd1-8c37-85fc486d6e21
md"### Data"

# ╔═╡ 2c79c088-b36d-11ee-085e-99f73de45b76
rocks = CSV.read("rocks.txt", DataFrame, header=false)

# ╔═╡ f27ca787-11ed-4870-945a-e1e4fd28f855
nsamples = size(rocks, 1)

# ╔═╡ 8cde1452-e939-44e0-832f-3c08608972f0
nelements = size(rocks, 2)

# ╔═╡ 4fbdb49b-65ca-4b53-999d-5b98472aede1
s = svd(Array(rocks)')

# ╔═╡ a02a64d3-b09e-4cf7-ad08-fec927a5fc53
plot(s.S, Layout(title="Singular Values"))

# ╔═╡ 1e156962-8b97-4fec-9935-afe45c9ae076
s.U

# ╔═╡ 01114f95-f1b7-4410-ae2b-32af961bfee7
begin
	mountains=zeros(15,11);
	mountains[1,:] = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0];
	mountains[2,:] = [0, 2, 3, 4, 5, 6, 5, 4, 3, 1, 0];
	mountains[3,:] = [0, 1, 2, 4, 5, 6, 5, 4, 2, 1, 0];
	mountains[4,:] = [0, 1, 3, 5, 5, 4.5, 4, 3.5, 3, 2.5, 0];
	mountains[5,:] = reverse(mountains[4,:]);
	mountains[6,:] = [0, 6, 6, 6, 5, 5, 5, 3.5, 2, 1, 0];
	mountains[7,:] = reverse(mountains[6,:]);
	mountains[8,:] = [0, 0.25, 0.5, 1, 3, 6, 3, 1, 0.5, 0.25, 0];
	mountains[9,:] = [0, 1, 2.5, 3, 4, 6, 4, 3, 2.5, 1, 0];
	mountains[10,:] = [0, 1, 2, 4, 6, 4, 2, 1, 1, 0.5, 0];
	mountains[1,:] = reverse(mountains[10,:]);
	mountains[12,:] = [0, 0, 0.5, 0, 1, 1, 2, 4, 6, 3, 0];
	mountains[13,:] = reverse(mountains[12,:]);
	mountains[14,:] = [0, 0, 0.5, 1, 3, 6, 5.5, 5.5, 4.5, 3.5, 0];
	mountains[15,:] = reverse(mountains[14,:]);
end;

# ╔═╡ fa431790-c3ce-43e8-9a76-1f8d07153a9c
function get_atm_pressure_data()
	Nx=20
	Nx2=Nx*Nx
	dx=1.0
	L = dx*(Nx-1)
	x = (dx*collect(0:Nx-1))'
	m0 = zeros(Nx,Nx);
	m1 = zeros(Nx,Nx);
	m2 = zeros(Nx,Nx);
	for p = 1:Nx
	for q = 1:Nx
     	m0[p,q] = sin(pi*x[p]/L).*sin(pi*x[q]/L);
     	m1[p,q] = sin(pi*x[p]/L).*sin(2*pi*x[q]/L);
     	m2[p,q] = sin(2*pi*x[p]/L).*sin(3*pi*x[q]/L);
	end
	end
		
	# build data with modes with these amplitudes at different times
	c0 = [1, -1, 1, -1, 1,   -1, 1, -1, 1, -1,     1, -1, 1, -1, 1,     -1, 1, -1, 1, -1,    1, -1, 1, -1, 1 ]';
	c1 = [1, 2, 3, 4, 5,     5, 4, 3, 2, 1,        1, 2, 3, 4, 5,        5, 4, 3, 2, 1,      1, 2, 3, 4, 5 ]';
	c2 = [0, 0, 0, 1, 2,     3, 2, 1, 0, 0,        0, 0, 1, 2, 1,        0, 0, 0, 0, 0,      1, 2, 3, 2, 1 ]';

	Nt = 25;
	D = zeros(Nx, Nx, Nt)
	for it in 1:Nt
		  D[:,:,it] .= c0[it]*m0 + c1[it]*m1 + c2[it]*m2 + 0.7 * randn(Nx,Nx);
	end
	return D
end

# ╔═╡ 88345a4f-aab7-4615-9795-020838f62292
atm_pressure = get_atm_pressure_data()

# ╔═╡ 126cacd6-8508-44e5-ba04-06a0704b4d17
plot(heatmap(z=atm_pressure[:,:, it_atm], colorscale="jet", showscale=false), Layout(title="t=$it_atm", width=200, height=200,xaxis=attr(showticklabels=false), yaxis=attr(showticklabels=false)))

# ╔═╡ 82fd02c2-26db-476b-b989-9d22c9853314
atm_pressure;

# ╔═╡ 0c19ea23-5e47-4724-8fa8-fd6023aa1aa3
size(atm_pressure)

# ╔═╡ ae553ef5-5e41-4b6c-9d7e-43195f86c46e
data_matrix_atm_pressure = reshape(atm_pressure, :, 25)

# ╔═╡ 1d5aaf6e-dad5-431f-bc98-488d6910f6fa
s_atm = svd(data_matrix_atm_pressure);

# ╔═╡ 2eddc7d2-d132-4abd-8d1c-6831ed933925
s_atm.Vt

# ╔═╡ c6d8b56e-50e7-44df-b78c-0b654cb72630
data_matrix_atm_pressure_hat = s_atm.U * Diagonal(s_atm.S) * s_atm.Vt

# ╔═╡ fbaf05be-8e3d-4b0b-8722-1c4838855728
C = s_atm.U; R = Diagonal(s_atm.S) * s_atm.Vt;

# ╔═╡ 573fb34d-98aa-4ed7-a1a8-bd3a2f6448e4
size(C)

# ╔═╡ 6f81aec1-0752-4b88-b844-2468056e301a
plot(heatmap(z=reshape(C[:, it_atm2], 20, 20), colorscale="jet", showscale=false), Layout(title="Singular Vector $it_atm2", width=200, height=200,xaxis=attr(showticklabels=false), yaxis=attr(showticklabels=false)))

# ╔═╡ 91bfbd85-74fb-419b-82e9-d81d42058bdd
plot(R[1, :])

# ╔═╡ ab75c417-232c-4e12-93f6-126a17dc0198
plot(R[2, :])

# ╔═╡ 12a2613e-b34b-472a-a1d2-4a375d24d3cf
plot(R[3, :])

# ╔═╡ be788b71-db91-45c1-95c4-2aa6bbdcdb50
plot(R[4, :])

# ╔═╡ 84ee7809-b9f8-4a22-9ed3-abbd1bdd81b8
plot(R[5, :])

# ╔═╡ 28c524db-f495-4284-b9a9-7198431b9106
C*R ≈ data_matrix_atm_pressure

# ╔═╡ 3e211600-9924-4405-bc53-15a3accb400a
data_matrix_atm_pressure ≈ data_matrix_atm_pressure_hat

# ╔═╡ 09a71c9a-804d-49e7-9d9f-cdb6188b9959
md"### Plots"

# ╔═╡ cccd6ea4-5040-46eb-b67a-106846ad3f72
function sample_plot(z, title)
	plot(heatmap(y=element_names, z=z, showscale=false, zmin=-1,zmax=1,colorscale="Blackbody"), Layout(title=string(title), width=125, xaxis=attr(showticklabels=false)))
end

# ╔═╡ 03f08f45-690b-4bd4-976d-e1d7b1360431
let 
	resample_data
	R = transpose(Array(rocks))
	map(rand(1:nsamples, 5)) do i
		sample_plot(log.(R[:, i:i]), i)
	end
end

# ╔═╡ 1c5f4e96-3ed1-4797-97e2-68d3399c41a4
map(1:5) do i
	sample_plot(s.U[:,i:i], i)
end

# ╔═╡ Cell order:
# ╠═ab1f729f-f7a9-4cf1-b1ef-7286e41bd37a
# ╟─1a7b09b7-e614-4dff-b252-dc726d87af81
# ╟─8ac94996-30af-41a9-b2fe-2b9982b07d78
# ╟─520c6e47-6df9-4ec9-8f48-f85ced6cb30b
# ╟─63bdb39b-19da-4f0f-98f4-9e9506e02927
# ╟─126cacd6-8508-44e5-ba04-06a0704b4d17
# ╠═82fd02c2-26db-476b-b989-9d22c9853314
# ╠═0c19ea23-5e47-4724-8fa8-fd6023aa1aa3
# ╠═ae553ef5-5e41-4b6c-9d7e-43195f86c46e
# ╠═1d5aaf6e-dad5-431f-bc98-488d6910f6fa
# ╠═2eddc7d2-d132-4abd-8d1c-6831ed933925
# ╠═c6d8b56e-50e7-44df-b78c-0b654cb72630
# ╠═fbaf05be-8e3d-4b0b-8722-1c4838855728
# ╠═28c524db-f495-4284-b9a9-7198431b9106
# ╠═5fd30c0b-eb53-4a54-a17e-6f12be989d4f
# ╠═573fb34d-98aa-4ed7-a1a8-bd3a2f6448e4
# ╟─6f81aec1-0752-4b88-b844-2468056e301a
# ╠═3e211600-9924-4405-bc53-15a3accb400a
# ╠═91bfbd85-74fb-419b-82e9-d81d42058bdd
# ╠═ab75c417-232c-4e12-93f6-126a17dc0198
# ╠═12a2613e-b34b-472a-a1d2-4a375d24d3cf
# ╠═be788b71-db91-45c1-95c4-2aa6bbdcdb50
# ╠═84ee7809-b9f8-4a22-9ed3-abbd1bdd81b8
# ╟─2df4d945-2470-45b1-993b-131decdb2eae
# ╟─1cc2d822-df8f-4fac-b8ca-770b6d1e09dd
# ╠═f27ca787-11ed-4870-945a-e1e4fd28f855
# ╠═8cde1452-e939-44e0-832f-3c08608972f0
# ╟─53f2cf9e-c1e6-4861-8b66-fcbe260dd6c9
# ╟─03f08f45-690b-4bd4-976d-e1d7b1360431
# ╟─a02a64d3-b09e-4cf7-ad08-fec927a5fc53
# ╠═1c5f4e96-3ed1-4797-97e2-68d3399c41a4
# ╠═4fbdb49b-65ca-4b53-999d-5b98472aede1
# ╠═1e156962-8b97-4fec-9935-afe45c9ae076
# ╠═cc80419d-f255-4ba9-a2a9-d52d2d6fd014
# ╟─fd2938e0-a187-41ef-9fe2-f5f0c262de1e
# ╠═abe35895-c888-4441-9a4e-94c16c91dc27
# ╟─e9950ad0-73d1-4fd1-8c37-85fc486d6e21
# ╠═2c79c088-b36d-11ee-085e-99f73de45b76
# ╠═01114f95-f1b7-4410-ae2b-32af961bfee7
# ╠═fa431790-c3ce-43e8-9a76-1f8d07153a9c
# ╠═88345a4f-aab7-4615-9795-020838f62292
# ╟─09a71c9a-804d-49e7-9d9f-cdb6188b9959
# ╠═cccd6ea4-5040-46eb-b67a-106846ad3f72
