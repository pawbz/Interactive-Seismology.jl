### A Pluto.jl notebook ###
# v0.19.11

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

# ╔═╡ 1d229d0c-2ffd-11ed-0cb3-1964ddd3b5af
begin
	using Pkg
	Pkg.add(url="https://github.com/anowacki/Geodesics.jl")
	Pkg.add(url="https://github.com/anowacki/Seis.jl")
	Pkg.add(url="https://github.com/anowacki/SeisRequests.jl")
	Pkg.add("GeoMakie")
	Pkg.add("CairoMakie")
	Pkg.add("PlutoUI")
	using Dates, Geodesics, Seis, SeisRequests
	using GeoMakie, CairoMakie
	using PlutoUI
end

# ╔═╡ f1035bd5-77b8-4a79-8445-4b11d01583e5
md"""
## Plot Earthquake Data
A notebook that quickly downloads sample earthquake data.

##### Introduction of Seismology
ES218; August 2022;

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 45c85015-bf14-4358-9087-8f6a526b95b5
SeisRequests.server_list()

# ╔═╡ d9b24b36-0adc-448c-9571-d08ca1c731c1
md"""
We will now request for events occuring between `starttime` and `endtime` from the IRIS server.
"""

# ╔═╡ e69d617a-f065-43ad-a31a-90ea7122c902
ev=get_events(starttime="2014-03-21T13:40:30.00",endtime="2014-03-21T13:45:40.00",minmagnitude=6,latitude=7.0, longitude=94.0,maxradius=5.0,server="IRIS",verbose=false)

# ╔═╡ 972b1e1f-56af-468a-9140-5562a60f2c61
md"""
We will now query for stations that recorded the first earthquake event, with some constraints on the network and station codes.
"""

# ╔═╡ 50dfe3fa-122b-44b1-850f-4234823228ab
newsta=get_stations(ev[1],network="II",station="UOSS",location="10",channel="BH?",server="IRIS",verbose=false)

# ╔═╡ 38a03c9a-15e6-45ad-936c-6973d5dc9694
md"""
Lets go ahead and download the data.
"""

# ╔═╡ bbf8b9cd-5946-479c-84fb-c47181639355
sta1=get_data(code="II.UOSS.10.BH?",starttime="2014-03-21T13:40:30.00",endtime="2014-03-21T14:40:30.00",verbose=false)

# ╔═╡ bce3ff3e-068d-4972-bf32-bc43213f9205
# need a slider to interact with the time axis
@bind tend PlutoUI.Slider(500:3000, default=3000)

# ╔═╡ d284818e-112c-44ef-82aa-28128fb9f4aa
md"""
Plotting the time series recorded at all the three component sensors.
"""

# ╔═╡ 36d3d32d-efde-4a78-b3ca-4b1534dc45ff
begin
	fig = Figure()
	for (i, s) in enumerate(sta1)
	ax = Axis(fig[i,1], title=string(s), xlabel="time [s]")
	lines!(ax, range(s.b, step=s.delta, length=length(s.t)), s.t)
	xlims!(ax, [300, tend])
	end
	fig
end

# ╔═╡ 329f43e2-c243-40fc-9703-569884ae0411
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╠═f1035bd5-77b8-4a79-8445-4b11d01583e5
# ╠═45c85015-bf14-4358-9087-8f6a526b95b5
# ╟─d9b24b36-0adc-448c-9571-d08ca1c731c1
# ╠═e69d617a-f065-43ad-a31a-90ea7122c902
# ╟─972b1e1f-56af-468a-9140-5562a60f2c61
# ╠═50dfe3fa-122b-44b1-850f-4234823228ab
# ╟─38a03c9a-15e6-45ad-936c-6973d5dc9694
# ╠═bbf8b9cd-5946-479c-84fb-c47181639355
# ╠═bce3ff3e-068d-4972-bf32-bc43213f9205
# ╟─d284818e-112c-44ef-82aa-28128fb9f4aa
# ╠═36d3d32d-efde-4a78-b3ca-4b1534dc45ff
# ╟─329f43e2-c243-40fc-9703-569884ae0411
# ╠═1d229d0c-2ffd-11ed-0cb3-1964ddd3b5af
