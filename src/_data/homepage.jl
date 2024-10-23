Dict(
    "title" => @htl("Interactive Seismology</strong>"),

    # # add a disclaimer to the course webpage. Remove it if you dont want to include it.
    "disclaimer" => md"""
A [collection](https://pawbz.github.io/Interactive-Seismology.jl/) of *reactive notebooks* using the [Julia](https://julialang.org/) notebook system, Pluto. These notebooks are designed to illustrate fundamental concepts in seismology in an interactive manner. They are primarily intended to support our course on Seismology (ES218), and are accessible to students with minimal programming experience. The notebooks contain all the essential mathematical expressions, which are easily manipulated using [Symbolics.jl](https://symbolics.juliasymbolics.org/dev/), a modern Computer Algebra System (CAS). By using these standalone notebooks, students can intuitively explore and gain a better understanding of wave phenomena in seismology.
    """,

    # Highlights the key features of your class to make it more engaging. Remove it if you dont want to include it.
    "highlights" => [
        Dict("name" => "Born Approximation", 
             "text" => md"""
             Students can draw the scatterers on a whiteboard, visualize the scattering response, and immediately image these scatterers.
             **This hands-on approach helps students grasp the concept of wave scattering and its applications in seismology.**
             """,
             "img" => "$(root_url)/assets/movies/born-approximation.gif"
        ),
        Dict("name" => "Seismic Anisotropy",
             "text" => md"""
            In anisotropic media, wavefronts are no longer perfect spheres but can take on ellipsoidal shapes, reflecting the directional dependence of wave speeds. 
            **This exploration helps students visualize and understand the impact of anisotropy on seismic wave behavior.**
             """,
             "img" => "$(root_url)/assets/movies/anisotropy.gif"
             ),
        Dict("name" => "Full Waveform Inversion",
             "text" => md"""
             In Full Waveform Inversion, the aim is to update the medium parameters (like velocity or density) 
             of a subsurface model to better match the observed seismic data. **Grasp why the correlation 
             between the forward propagated source wavefield and the backpropagated residual wavefield is crucial for these updates.**
             """,
             "img" => "$(root_url)/assets/movies/full-waveform-inversion.gif"
             ),
        Dict("name" => "Ray Tomography",
             "text" => md"""
             Seismic tomography is a technique used to image the interior of the Earth (or other planetary bodies) by analyzing the propagation of seismic waves.
                **These notebooks helps students understand the basic principles of tomography.**
             """,
             "img" => "$(root_url)/assets/movies/ray-tomography.gif"
        )
    ]
)