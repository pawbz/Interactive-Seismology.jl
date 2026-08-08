### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Global Bodywave Arrivals"
#> date = "2025-08-05"
#> tags = ["raytheory"]
#> description = "What is the arrival you are looking for?"
#> layout = "layout.jlhtml"

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

# ╔═╡ 47b2c09a-2ae8-49f0-ba73-ddb6868417b1
begin
    using CondaPkg
    CondaPkg.add_pip("obspy")
    CondaPkg.add_pip("matplotlib")
    using PythonCall
    using PythonPlot
    using PlutoUI
    pygui(false)
end

# ╔═╡ 025b2827-ed43-45f5-a981-56dd599c72cb
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ 6b3bba88-b693-4e39-8866-8166dfc55c30
md"""
# Global Bodywave Arrivals
This notebook interactively visualizes global seismic bodywave arrivals using the TauP toolkit via Python in Julia. Users can select source depth and receiver distance to explore ray paths, phases, and travel times for various seismic waves.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 66e7dd53-6423-4c4f-ba45-49d0439d51cf
md"""
## Select Source Depth and Receiver Distance

- Receiver Distance (°)
  $(@bind receiver_distance Slider(0.0:1.0:180.0, show_value=true, default=120))

- Source Depth (km) 
  $(@bind source_depth Slider(0.0:10.0:700.0, show_value=true, default=20))
"""

# ╔═╡ 5a439ed2-6333-4e44-bade-67e227975b92
@bind selected_phases MultiCheckBox(phases, default=phases[1:1])

# ╔═╡ 9b0b4e7e-fd8f-4573-af80-ed76ff2848f5
md"## Appendix"

# ╔═╡ dd4cb9d8-8d6e-4ea8-b6bf-545631fecff8
taup = pyimport("obspy.taup")

# ╔═╡ 23f2f44c-144a-4fd1-a429-656bf0af4cca
model = taup.TauPyModel(model="iasp91")

# ╔═╡ 7568eb42-fe1b-44bc-86a1-dda9200bb49b
arrivals = model.get_ray_paths(source_depth, receiver_distance)

# ╔═╡ 5299340f-020c-4040-adc8-d6b308c3738e
phases = unique(pyconvert(String, arrival.name) for arrival in arrivals)

# ╔═╡ 3c203808-3886-4d94-9e05-b893d0ba6c4d
arrivals_filtered = isempty(selected_phases) ? nothing : model.get_ray_paths(
    source_depth, receiver_distance; phase_list=selected_phases,
)

# ╔═╡ 460863de-314c-4196-b72b-eb6382cce6f7
if isnothing(arrivals_filtered)
    md"""
    !!! note "Choose a phase"
        Select at least one arrival phase to draw its ray path and travel time.
    """
else
    fig = figure(figsize=(10, 10))
    arrivals_filtered.plot_rays(plot_type="spherical", fig=fig,
        legend=true, label_arrivals=false,
        plot_all=true, indicate_wave_type=true)
    fig.tight_layout(pad=15)
    fig
end

# ╔═╡ 78b1bb95-14a6-461f-a668-8fbbaa7e5ee0
phases_selected = isnothing(arrivals_filtered) ? String[] :
    [pyconvert(String, arrival.name) for arrival in arrivals_filtered]

# ╔═╡ 6dfc1926-4da7-405b-aca4-7eac0aa814c8
traveltimes = isnothing(arrivals_filtered) ? Float64[] :
    [pyconvert(Float64, arrival.time) for arrival in arrivals_filtered]

# ╔═╡ 10cc4de1-f635-4a4c-b195-19dca1e44d4c
md"""
Traveltime in seconds: $(Dict(pyconvert(String, ph)=>t for (ph, t) in zip(phases_selected, traveltimes)))
"""

# ╔═╡ 5d94e67e-5334-4e1c-9838-749b2318c66d
md"""
## Credits
- [https://www.seis.sc.edu/taup/](https://www.seis.sc.edu/taup/)
- [https://docs.obspy.org/packages/obspy.taup.html](https://docs.obspy.org/packages/obspy.taup.html)
"""

# ╔═╡ Cell order:
# ╠═025b2827-ed43-45f5-a981-56dd599c72cb
# ╟─6b3bba88-b693-4e39-8866-8166dfc55c30
# ╟─66e7dd53-6423-4c4f-ba45-49d0439d51cf
# ╟─5a439ed2-6333-4e44-bade-67e227975b92
# ╟─10cc4de1-f635-4a4c-b195-19dca1e44d4c
# ╟─460863de-314c-4196-b72b-eb6382cce6f7
# ╠═7568eb42-fe1b-44bc-86a1-dda9200bb49b
# ╠═5299340f-020c-4040-adc8-d6b308c3738e
# ╠═78b1bb95-14a6-461f-a668-8fbbaa7e5ee0
# ╠═6dfc1926-4da7-405b-aca4-7eac0aa814c8
# ╠═3c203808-3886-4d94-9e05-b893d0ba6c4d
# ╟─9b0b4e7e-fd8f-4573-af80-ed76ff2848f5
# ╠═47b2c09a-2ae8-49f0-ba73-ddb6868417b1
# ╠═dd4cb9d8-8d6e-4ea8-b6bf-545631fecff8
# ╠═23f2f44c-144a-4fd1-a429-656bf0af4cca
# ╟─5d94e67e-5334-4e1c-9838-749b2318c66d
