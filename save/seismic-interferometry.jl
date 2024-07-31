### A Pluto.jl notebook ###
# v0.19.32

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

# ╔═╡ 290b50eb-0faf-4a5f-86bf-36628e61ff06
using PlutoUI,DSP,Statistics,PlutoPlotly,FourierTools,FFTW,StatsBase 

# ╔═╡ f986cb42-7cbe-11ee-22de-cbc3714ed81d
TableOfContents()

# ╔═╡ bd04a8af-7666-47b7-87cc-a212ec8adcbd
md"""
## Numerical Simulation of Seismic Interferometry
We create random sources around the 2 receivers and send single wavelet from all the sources to the receivers. The wavelet response due to a source at the 2 receivers are cross-correlated and the cross-correlations of all the sources are stacked to get the approximate Green's function between the 2 receivers.
"""

# ╔═╡ c648890f-ae42-4658-8e87-ba47c6300d69
md"""
### All responses at Receiver 1
"""

# ╔═╡ 57730f00-d0d6-4f54-8868-eb03c9482a30
md"""
### All response at Receiver 2
"""

# ╔═╡ 1bc0ce25-4de0-4628-a308-8106544fc0ec
md"""
## Cross-Correlation
"""

# ╔═╡ e0ae191c-b674-46c7-bf05-bfecade313a5
md"""
## Generate Wavelets
"""

# ╔═╡ 0279d32c-d884-4990-a1a6-68b43fae2f92
begin
	responsetype = Lowpass(0.5; fs=10)
   designmethod = Chebyshev1(4,1)


function get_wav()

      
    win=(rand(Float32,200) .* rand((-1,1),200)) .*        gaussian(200,rand(0.005:0.001:0.1),padding=0,zerophase=false)
	
    
    hk=filt(digitalfilter(responsetype, designmethod),win);
    hk = hk ./ std(hk)
	# m1=maximum(abs.(hk))
	# hk = hk ./ m1
    
return hk
end
end

# ╔═╡ 577d411c-88b5-495d-a49a-d2379336f533
md"""
## Random Sources
"""

# ╔═╡ 8d45d812-fb70-45c5-a123-4f93844fef6d
md"""
## Receivers
"""

# ╔═╡ 3b144853-0d82-46e2-a4ee-301516f3ce10
vel=0.035

# ╔═╡ 3cf54d17-5931-461c-9e35-bad1fc8d9ac3
function trval(rec,src)
	dis=sqrt((rec[1]-src[1])^2+(rec[2]-src[2])^2)
	tm= dis/vel
	return tm
end

# ╔═╡ 9ac02f1f-25b2-404b-bd3c-1c6faf9b8cab
md"""
### Traveltime between the receivers.
"""

# ╔═╡ 9c06aa0f-6674-477f-945d-f1a081182ec3
md"""
Calculate the maximum possible travel time.
"""

# ╔═╡ 6397b8ae-62ea-404a-8add-08e79a39cd30
import PlutoUI: combine

# ╔═╡ 41870289-3874-43ca-9837-3cead43d0200
function Ricker(;dt::Real=0.002,f0::Real=20.0)

    nw = 2.0/(f0*dt)
    nc = floor(Int, nw/2)
    t = dt*collect(- nc:1:nc)
    b = (pi*f0*t).^2
    w = (1 .- 2 .* b).*exp.(-b)
    w = cat(zeros(Float32,75),w,zeros(Float32,75),dims=1)
end

# ╔═╡ 1083c0c6-beb2-42ba-8c99-e8d8d77ce83a
md"""
### Change the setting:
* Wavelets $(@bind wav Select([Ricker() => "Ricker",get_wav() => "Random"]))
* Number of Sources $(@bind sr Slider(1000:1000:10000,default=1000,show_value=true))
* Inner radius limit of the sources $(@bind inrad Slider(1:2:17,default=1,show_value=true))
* Angular Range $(@bind iang Select(["all" => "Circular","small" => "Sectional","leftarc" => "Right Arc","rightarc" => "Left Arc","new" => "Special case"]))
"""

# ╔═╡ e1f1fa09-87ef-4cdb-94b1-0d2e0e483323
begin
function delaywav(t)
	ts = wav
	# ts = get_wav()
	ts = cat(zeros(Float32,200),ts,zeros(Float32,1000),dims=1)
	td= shift(ts,tuple(t-300))
	return td[1:1000]
	
end
end

# ╔═╡ 09c19be2-867a-4c16-ace5-6e4a9389472c
function get_seis(rec1,srcloc)   # get all the responses from all sources at receiver
	tm=map(1:length(srcloc)) do i 
		trval(rec1,srcloc[i])
	end
	wav=map(1:length(tm)) do h
		delaywav(tm[h])
	end
	wav=stack(wav,dims=2)
	#wav=sum(wav,dims=2)
end

# ╔═╡ e09b72c2-5d7a-4215-ad77-f994f051cc93
ang1=range(-0.1*π,stop=0.1*π,length=div(sr,2))

# ╔═╡ d3542b58-07f4-4016-9d68-75b7978b79a4
ang2=range(π*0.9,stop=(π*1.1),length=div(sr,2))

# ╔═╡ c71bc2d6-cc33-4cae-8195-5fe26ec39d7f
begin

num=sr
rad=20

smallang=cat(ang1,ang2,dims=1)
allang=range(0,stop=2*π,length=num)
	leftarc=range(π/3,stop=-π/3,length=num)
	rightarc=range(2*π/3,stop=4*π/3,length=num)
	R1=range(π/3,stop=2*π/3,length=div(num,2))
	R2=range(4*π/3,stop=5*π/3,length=div(num,2))
    R=cat(R1,R2,dims=1)
	new=R
	#newang=range((cat((0,stop=π/2),(π,stop=3*π/2)),dims=1),length=num)
if iang == "all"
  ang=allang
elseif iang == "small"
	ang=smallang
elseif iang == "leftarc"
	ang=leftarc
elseif iang == "rightarc"
	ang=rightarc
elseif iang== "new"
	ang=new
end

rr=rand(inrad:0.05:rad,num)
ys=map(rr,ang) do r,a
	r*sin(a)
end
xs=map(rr,ang) do r,a
	r*cos(a)
end
 end

# ╔═╡ ddab0b60-818a-4fc0-9aca-05f54c8b8352
srcloc=[[xs[i],ys[i]] for i in 1:num]  # all the source location

# ╔═╡ d9727185-8f69-49a5-a98f-7514562529ab
md"""
## Plots
"""

# ╔═╡ ccf893d9-5b6e-4a79-b09d-3da957bc612c
function plot_line(tr; title="",nt=1999)

	nt=nt
 tgrid = range(-100, stop=100, length=nt)
	
	fig=PlutoPlotly.Plot(Layout(height=300, width=600, 
		title=title, font=attr(
            size=10),))

	#data=[reshape(data1[i],:,2,size(data1[i],2)) for i in 1:length(data1)]
    
	
	
	add_trace!(fig, PlutoPlotly.scatter(x=tgrid, y=tr[1],opacity=1, mode="lines",))
	
	
	PlutoPlotly.plot(fig)
	
end

# ╔═╡ a0658b4c-29ea-4a93-9355-87effa4ea67c
function plot_cross(D; title="",nt=1999, ylims=(-8, 8),scale=1)
spacing=10
	#scale = 9
	nt=nt
 tgrid = range(-100, stop=100, length=nt)
	
	fig=PlutoPlotly.Plot(Layout(yaxis_autorange="reversed", height=600, width=800, 
		title=title, font=attr(
            size=10),))

	#data=[reshape(data1[i],:,2,size(data1[i],2)) for i in 1:length(data1)]
    
	
	for i in 1:size(D,2) 
	add_trace!(fig, PlutoPlotly.scatter(x=(scale*(D[:,i])).+(i-1)*spacing, y=tgrid,opacity=1, mode="lines",), row=1, col=1, )
	
	end
	PlutoPlotly.plot(fig)
	
end


# ╔═╡ eb5cee03-ccad-4b48-8d10-446c8df4c7ee
function plot_trace(D; title="",nt=1000, ylims=(-8, 8),scale=1)
spacing=10
	#scale = 9
	nt=nt
 tgrid = range(0, stop=100, length=nt)
	
	fig=PlutoPlotly.Plot(Layout(yaxis_autorange="reversed", height=600, width=800, 
		title=title, font=attr(
            size=10),))

	#data=[reshape(data1[i],:,2,size(data1[i],2)) for i in 1:length(data1)]
    
	
	for i in 1:size(D,2) 
	add_trace!(fig, PlutoPlotly.scatter(x=(scale*(D[:,i])).+(i-1)*spacing, y=tgrid,opacity=1, mode="lines",), row=1, col=1, )
	
	end
	PlutoPlotly.plot(fig)
	
end

# ╔═╡ 5427df02-9278-41c9-9b3b-b5d1ab5b4c24
function parameters_input(pos::Vector)
	
	return combine() do Child
		
		inputs = [
			md""" $(pos[1]): $(
				Child(pos[1], Slider(-20:1:20,default=-5,show_value=true))
			)"""
		   md""" $(pos[2]): $(
				Child(pos[2], Slider(-20:1:20,default=0,show_value=true))
			)"""
			md""" $(pos[3]): $(
				Child(pos[3], Slider(-20:1:20,default=5,show_value=true))
			)"""
			md""" $(pos[4]): $(
				Child(pos[4], Slider(-20:1:20,default=0,show_value=true))
			)"""
			
			#for pos in directions
		]
		
		md"""
		#### Receiver Position
		$(inputs)
		"""
	end
end

# ╔═╡ c4339754-42eb-4969-9045-010ca1d23ef8
@bind pos confirm(parameters_input(["R1 (x-pos)","R1 (y-pos)","R2 (x-pos)","R2 (y-pos)"]))

# ╔═╡ ca1b63bf-c868-4200-80f0-f8d6e40d9a6c
begin
rec1=[pos[1],pos[2]]
rec2=[pos[3],pos[4]]
end

# ╔═╡ a0f471bd-8139-4f42-bf0e-c019e6d42832
rectime=trval(rec1,rec2)

# ╔═╡ 0458e693-5e33-4372-9951-e3c50b301ae5
begin
real=delaywav(rectime)
plot_line([cat(reverse(real),real,dims=1)],title="Expected Response (Both Causal and Acausal)")
end

# ╔═╡ 280b6269-a022-4846-9071-61481649b008
trval(rec1,[-20,0])

# ╔═╡ 97289274-85f2-46f7-b9bf-916f0a420301
seis1=get_seis(rec1,srcloc);

# ╔═╡ 617fd9b1-c6be-476e-9622-b6482602c7c7
#plot_trace(seis1[:,300:400])
let
step=div(size(seis1,2),200)
rr1=stack(map(1:step:size(seis1,2)) do i
     seis1[:,i]
end,dims=2)
plot(heatmap(z=rr1))
end

# ╔═╡ db6547b2-0939-4d10-b590-4205f8c10ed3
seis2=get_seis(rec2,srcloc);

# ╔═╡ 64a67ad0-2a5f-49c4-8a90-8fe8fb55342f
#plot_trace(seis2[:,1:100])
let
step=div(size(seis2,2),200)
rr2=stack(map(1:step:size(seis2,2)) do i
     seis2[:,i]
end,dims=2)
plot(heatmap(z=rr2))# Layout=attr(yaxis=attr(title="time (s)")))
end

# ╔═╡ 2bad7a74-1b12-499f-b435-7f99c7d64ae9
begin
cross=map(1:size(seis1,2)) do i
	xcorr(seis1[:,i],seis2[:,i],padmode=:longest)
end;
cross=stack(cross,dims=2);
end

# ╔═╡ 46c9b13c-088c-47e7-8fd3-b56aee4a8923
begin
gren=mean(cross,dims=2);
plot_line([gren[:,1]],title="Stacked Cross-Correlation")
end

# ╔═╡ e69db63a-52b4-45c0-8b17-5239468cf5ad
#plot_cross(cross[:,])
let
step=div(size(cross,2),200)
rr1=stack(map(1:step:size(cross,2)) do i
     cross[:,i]
end,dims=2)
plot(heatmap(z=rr1))
end

# ╔═╡ 77d2d0d5-05e8-44e4-910b-b8bc8cf3137f
function plot_point(x1,x2)
    layout=Layout(xaxis_range=(-rad,rad),yaxis_range=(-rad,rad),width=600,height=400)
	fig=Plot(layout)

  add_trace!(fig,scatter(x=x1,y=x2,mode="markers",showlegend=false,marker=attr(size=1,color="black")))
  add_trace!(fig,scatter(x=[rec1[1]],y=[rec1[2]],mode="markers",name="Receiver 1",marker=attr(size=15,color="blue",symbol="triangle-down")))
	add_trace!(fig,scatter(x=[rec2[1]],y=[rec2[2]],mode="markers",name="Receiver 2",marker=attr(size=15,color="green",symbol="triangle-down")))
	return fig
end

# ╔═╡ b933216e-203a-47a1-b4f1-901fa5b38f70
plot_point(xs,ys)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
FourierTools = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
DSP = "~0.7.9"
FFTW = "~1.7.1"
FourierTools = "~0.4.2"
PlutoPlotly = "~0.4.2"
PlutoUI = "~0.7.53"
StatsBase = "~0.34.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "349afb7eda2e83bd4f869ccb554f5b4e567b9e75"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractNFFTs]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "292e21e99dedb8621c15f185b8fdb4260bb3c429"
uuid = "7f219486-4aa7-41d6-80a7-e08ef20ceed7"
version = "0.8.2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "02f731463748db57cc2ebfbd9fbc9ce8280d3433"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.1"

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "e28912ce94077686443433c2800104b061a827ed"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.39"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "1c9b6f39f40dba0ef22244a175e2d4e42c8f6ee7"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.2.0"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BasicInterpolators]]
deps = ["LinearAlgebra", "Memoize", "Random"]
git-tree-sha1 = "3f7be532673fc4a22825e7884e9e0e876236b12a"
uuid = "26cce99e-4866-4b6d-ab74-862489e035e0"
version = "0.7.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e0af648f0692ec1691b5d094b8724ba1346281cf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.18.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

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

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

    [deps.CompositionsBase.weakdeps]
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

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

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "ffb97765602e3cbe59a0589d237bf07f245a8576"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.1"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.FourierTools]]
deps = ["ChainRulesCore", "FFTW", "IndexFunArrays", "LinearAlgebra", "NDTools", "NFFT", "PaddedViews", "Reexport", "ShiftedArrays"]
git-tree-sha1 = "8967a9d259ab1c50e3b3abc6b77d3e3d829d2e6d"
uuid = "b18b359b-aebc-45ac-a139-9c0ccbb2871e"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

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

[[deps.IndexFunArrays]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "e66a2aeb6d5814015004080e5203dfff44d2856f"
uuid = "613c443e-d742-454e-bfc6-1d7f8dd76566"
version = "0.2.6"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad37c091f7d7daf900963171600d7c1c5c3ede32"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

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

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "eb006abbd7041c28e0d16260e50a24f8f9104913"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.2.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

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
version = "2.28.2+0"

[[deps.Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

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

[[deps.NDTools]]
deps = ["LinearAlgebra", "OffsetArrays", "PaddedViews", "Random", "Statistics"]
git-tree-sha1 = "10e35b25261dfd361045e16afa9db5b14a4c1184"
uuid = "98581153-e998-4eef-8d0d-5ec2c052313d"
version = "0.5.2"

[[deps.NFFT]]
deps = ["AbstractNFFTs", "BasicInterpolators", "Distributed", "FFTW", "FLoops", "LinearAlgebra", "Printf", "Random", "Reexport", "SnoopPrecompile", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "93a5f32dd6cf09456b0b81afcb8fc29f06535ffd"
uuid = "efe261a4-0d2b-5849-be55-fc731d526b0d"
version = "0.13.3"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

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
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

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
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "BaseDirs", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "Reexport", "TOML"]
git-tree-sha1 = "d169dae7b340a1b0ca50ef359ed100f07ad677a4"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.4.2"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "db8ec28846dbf846228a32de5a6912c63e2052e3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.53"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "ea78a2764f31715093de7ab495e12c0187f231d1"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.4"

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

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

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

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

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

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "e579d3c991938fecbb225699e8f611fa3fbf2141"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.79"

    [deps.Transducers.extensions]
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

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
# ╠═f986cb42-7cbe-11ee-22de-cbc3714ed81d
# ╟─bd04a8af-7666-47b7-87cc-a212ec8adcbd
# ╟─1083c0c6-beb2-42ba-8c99-e8d8d77ce83a
# ╟─c4339754-42eb-4969-9045-010ca1d23ef8
# ╠═b933216e-203a-47a1-b4f1-901fa5b38f70
# ╠═46c9b13c-088c-47e7-8fd3-b56aee4a8923
# ╠═0458e693-5e33-4372-9951-e3c50b301ae5
# ╠═c648890f-ae42-4658-8e87-ba47c6300d69
# ╠═64a67ad0-2a5f-49c4-8a90-8fe8fb55342f
# ╠═57730f00-d0d6-4f54-8868-eb03c9482a30
# ╠═617fd9b1-c6be-476e-9622-b6482602c7c7
# ╠═1bc0ce25-4de0-4628-a308-8106544fc0ec
# ╠═e69db63a-52b4-45c0-8b17-5239468cf5ad
# ╠═e0ae191c-b674-46c7-bf05-bfecade313a5
# ╠═0279d32c-d884-4990-a1a6-68b43fae2f92
# ╠═e1f1fa09-87ef-4cdb-94b1-0d2e0e483323
# ╠═577d411c-88b5-495d-a49a-d2379336f533
# ╠═e09b72c2-5d7a-4215-ad77-f994f051cc93
# ╠═d3542b58-07f4-4016-9d68-75b7978b79a4
# ╠═c71bc2d6-cc33-4cae-8195-5fe26ec39d7f
# ╠═8d45d812-fb70-45c5-a123-4f93844fef6d
# ╠═ca1b63bf-c868-4200-80f0-f8d6e40d9a6c
# ╠═3b144853-0d82-46e2-a4ee-301516f3ce10
# ╠═3cf54d17-5931-461c-9e35-bad1fc8d9ac3
# ╠═9ac02f1f-25b2-404b-bd3c-1c6faf9b8cab
# ╠═a0f471bd-8139-4f42-bf0e-c019e6d42832
# ╟─9c06aa0f-6674-477f-945d-f1a081182ec3
# ╠═280b6269-a022-4846-9071-61481649b008
# ╠═ddab0b60-818a-4fc0-9aca-05f54c8b8352
# ╠═09c19be2-867a-4c16-ace5-6e4a9389472c
# ╠═97289274-85f2-46f7-b9bf-916f0a420301
# ╠═db6547b2-0939-4d10-b590-4205f8c10ed3
# ╠═2bad7a74-1b12-499f-b435-7f99c7d64ae9
# ╠═290b50eb-0faf-4a5f-86bf-36628e61ff06
# ╠═6397b8ae-62ea-404a-8add-08e79a39cd30
# ╠═41870289-3874-43ca-9837-3cead43d0200
# ╠═d9727185-8f69-49a5-a98f-7514562529ab
# ╠═77d2d0d5-05e8-44e4-910b-b8bc8cf3137f
# ╠═ccf893d9-5b6e-4a79-b09d-3da957bc612c
# ╠═a0658b4c-29ea-4a93-9355-87effa4ea67c
# ╠═eb5cee03-ccad-4b48-8d10-446c8df4c7ee
# ╠═5427df02-9278-41c9-9b3b-b5d1ab5b4c24
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
