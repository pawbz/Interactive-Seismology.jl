### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ 8d8f440d-4b79-4835-841b-739b5171f979
using PlutoUI, Printf, Roots, LinearAlgebra, FFTW

# ╔═╡ 29afa36b-391f-4861-aead-d62918daf3c6
using HypertextLiteral: @htl

# ╔═╡ 8d5d9594-2197-4ccc-a7a4-0f0d54a2370a
using PlutoPlotly

# ╔═╡ 2e91160e-33fa-4e2f-876e-8cccfebd27a5
using Turing

# ╔═╡ ee482657-fd34-4c92-95cd-0bf8e254676c
TableOfContents(include_definitions=true)

# ╔═╡ d3bc6646-b22f-4e22-87a8-da37a9aa1d1a
md"""
# Love Wave Dispersion Curves
This notebook provides an environment to explore dispersion of Love waves, a type of surface seismic wave that plays a key role in seismology and Earth structure studies.  

Love waves are horizontally polarized shear waves that are *trapped* in near-surface low-velocity layers. Their dispersion (variation of phase velocity with period) carries information about the elastic properties and layering of the crust and upper mantle.

Here,
- build a layered Earth model interactively, by adjusting thickness, shear velocity (`Vs`), and density (`ρ`) for each layer;
- compute dispersion curves for Love waves and observe how changes in structure affect the shape of the dispersion curve.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ d2c48242-92ce-11f0-01e6-cb1b8fe1fecb
struct Layer
    thickness::Float64   # km
    vp::Float64          # km/s
    vs::Float64          # km/s
    rho::Float64         # g/cm³
    Qp::Float64          # P-wave quality factor
    Qs::Float64          # S-wave quality factor
end

# ╔═╡ cee3bbac-9de8-46e7-8ba9-2cc1fe18802c
"""
    create_layers(n_layers::Int, thickness::Vector{Float64}, vs::Vector{Float64}; 
                  vp_vs_ratio::Float64=1.73, rho::Union{Float64, Vector{Float64}}=2.8, 
                  Qp::Float64=1000.0, Qs::Float64=1000.0)

Create a vector of Layer structs from input parameters.

**Arguments:**
- `n_layers`: number of layers
- `thickness`: vector of layer thicknesses (km), last layer ignored (half-space)
- `vs`: vector of shear velocities (km/s)
- `vp_vs_ratio`: Vp/Vs ratio for computing Vp (default: 1.73)
- `rho`: scalar density (g/cm³) or vector of densities per layer
- `Qp`: P-wave quality factor (default: 1000)
- `Qs`: S-wave quality factor (default: 1000)

**Returns:**
- `Vector{Layer}`: array of Layer structs

**Example:**
```julia
layers = create_layers(3, [10.0, 20.0, 0.0], [3.5, 4.0, 4.5])
```
"""
function create_layers(n_layers::Int, thickness::Vector{Float64}, vs::Vector{Float64}; 
                       vp_vs_ratio::Float64=1.73, rho::Union{Float64, Vector{Float64}}=2.8, 
                       Qp::Float64=1000.0, Qs::Float64=1000.0)
    
    # Validate input lengths
    if length(thickness) < n_layers
        error("thickness vector must have at least $n_layers elements")
    end
    if length(vs) < n_layers
        error("vs vector must have at least $n_layers elements")
    end
    
    # Handle density: scalar or vector
    rho_vals = if isa(rho, Vector)
        if length(rho) < n_layers
            error("rho vector must have at least $n_layers elements")
        end
        rho
    else
        fill(rho, n_layers)
    end
    
    # Create layers
    layers = Layer[]
    for i in 1:n_layers
        t = (i < n_layers) ? thickness[i] : 0.0  # Last layer is half-space (no thickness)
        v_s = vs[i]
        v_p = vp_vs_ratio * v_s
        ρ = rho_vals[i]
        
        push!(layers, Layer(t, v_p, v_s, ρ, Qp, Qs))
    end
    
    return layers
end

# ╔═╡ 363d8b6d-ddc9-4048-a02b-16346158112c
"""
    starting_phase_velocity(vp, vs) -> c

Compute a starting solution for phase velocity using Newton iteration.
"""
function starting_phase_velocity(vp::Float64, vs::Float64)
    c = 0.95 * vs
    for _ in 1:5
        γ = vs / vp
        κ = c / vs
        k2 = κ^2
        gk2 = (γ * κ)^2
        fac1 = sqrt(1 - gk2)
        fac2 = sqrt(1 - k2)
        fr = (2 - k2)^2 - 4 * fac1 * fac2
        frp = -4 * (2 - k2) * κ +
              4 * fac2 * γ^2 * κ / fac1 +
              4 * fac1 * κ / fac2
        frp /= vs
        c -= fr / frp
    end
    return c
end

# ╔═╡ 6c766c4b-8a1c-4f8b-9467-fec5a91fcd05
struct Velocities
    periods::Vector{Float64}
    phase_velocities::Vector{Float64}
    group_velocities::Vector{Float64}
end

# ╔═╡ 5510a81e-1ac0-45b9-9657-13f0db32f3e1
struct EigenFunctions
	eigenfunctions_u::Matrix{Float64}  # displacement eigenfunctions [depth, frequency]
    eigenfunctions_t::Matrix{Float64}  # stress eigenfunctions [depth, frequency]
    depths::Vector{Float64}            # depth sampling points
    energy_integrals::Vector{Float64}  # energy normalization factors
end

# ╔═╡ 7e9acd7c-9cf5-45e0-8892-77c77a1fd983
md"""
Maximum period: $(@bind period_max Slider(20.0:5.0:200.0, default=100.0, show_value=true)) s
"""

# ╔═╡ 756e06e3-4f25-4f55-b242-26a45a1b7dd3
periods = collect(range(5.0, stop=period_max, length=64))

# ╔═╡ 9c9ad5f6-debf-4376-aaf6-759e9b0c949b
"""
    layer_matrix_SH(layer::Layer, ω, c)

Return the 2×2 complex transfer matrix for an SH layer of thickness d:
    [ u_bottom ]   [ A  B ] [ u_top   ]
    [ t_bottom ] = [ C  D ] [ t_top   ]

Where u is horizontal displacement and t is shear traction.
"""
function layer_matrix_SH(layer::Layer, ω::Float64, c::Float64)
    # angular wavenumber horizontally
    k = ω / c

    # shear velocity and density
    vs = layer.vs
    ρ = layer.rho

    # shear modulus (units consistent if rho in mass/volume and vs in velocity units)
    μ = ρ * vs^2

    # vertical wavenumber q (can be complex)
    # q^2 = (ω/vs)^2 - k^2
    q2 = (ω / vs)^2 - k^2
    q = sqrt(complex(q2))          # complex sqrt handles evanescent/propagating

    d = layer.thickness

    # if thickness is effectively zero (half-space), return identity (no propagation)
    if isapprox(d, 0.0; atol=1e-14)
        return Matrix{ComplexF64}(I, 2, 2)
    end

    # arguments for trig/hyperbolic via complex arithmetic
    α = q * d

    # use cos and sin of complex argument (covers cosh/sinh automatically)
    ca = cos(α)
    sa = sin(α)

    # Layer transfer matrix for SH (standard form)
    # [ cos(qd)          sin(qd)/(μ*q)   ]
    # [ -μ*q*sin(qd)     cos(qd)         ]
    M = zeros(ComplexF64, 2, 2)
    M[1, 1] = ca
    M[1, 2] = sa / (μ * q)
    M[2, 1] = -μ * q * sa
    M[2, 2] = ca

    return M
end






# ╔═╡ 4f122bf3-0693-4f13-b716-c0dc2863af25
# # complete_eigenfunction_implementation
# """
#     compute_love_eigenfunctions_complete(model::EarthModel, period::Float64, phase_velocity::Float64)

# Complete implementation of Love wave eigenfunction computation using layer matrix propagation.
# This follows the methodology from CPS rsh.f and implements the full matrix approach.
# """
# function compute_love_eigenfunctions(model::EarthModel, period::Float64, phase_velocity::Float64)
#     ω = 2π / period
#     c = phase_velocity
#     layers = model.layers
#     n_layers = length(layers)
    
#     # Create comprehensive depth sampling
#     depths = create_depth_sampling(model)
#     n_depths = length(depths)
    
#     # Initialize eigenfunction arrays
#     u_eigenfunction = zeros(ComplexF64, n_depths)
#     t_eigenfunction = zeros(ComplexF64, n_depths)
    
#     # Compute layer matrices and propagation
#     layer_matrices = Vector{Matrix{ComplexF64}}(undef, n_layers)
#     layer_properties = []
    
#     # Build layer matrices
#     for i in 1:n_layers
#         layer = layers[i]
#         vs = layer.vs
#         ρ = layer.rho
#         μ = ρ * vs^2
        
#         # Vertical wavenumber
#         k = ω / c
#         q2 = (ω / vs)^2 - k^2
#         q = sqrt(complex(q2))
        
#         push!(layer_properties, (vs=vs, ρ=ρ, μ=μ, q=q))
        
#         if i < n_layers  # Not half-space
#             d = layer.thickness
#             α = q * d
            
#             # Layer transfer matrix
#             ca = cos(α)
#             sa = sin(α)
            
#             layer_matrices[i] = [
#                 ca              sa/(μ*q);
#                 -μ*q*sa         ca
#             ]
#         else  # Half-space
#             # Radiation condition matrix
#             layer_matrices[i] = [
#                 1.0     0.0;
#                 0.0    1.0
#             ]
#         end
#     end
    
#     # Forward propagation to establish boundary conditions
#     # Start from free surface: u(0) = 1, t(0) = 0
#     surface_vector = [complex(1.0), complex(0.0)]
    
#     # Propagate through all layers to check radiation condition
#     total_matrix = Matrix{ComplexF64}(I, 2, 2)
#     for i in 1:(n_layers-1)
#         total_matrix = layer_matrices[i] * total_matrix
#     end
    
#     # Apply radiation condition at bottom
#     bottom_vector = total_matrix * surface_vector
#     μ_bottom = layer_properties[end].μ
#     q_bottom = layer_properties[end].q
    
#     # Check dispersion relation: t_bottom + μ*q*u_bottom ≈ 0
#     dispersion_check = bottom_vector[2] + μ_bottom * q_bottom * bottom_vector[1]

# 	@show dispersion_check
#     # Now compute eigenfunctions at all depths
#     current_depth_idx = 1
#     cumulative_depth = 0.0
    
#     for layer_idx in 1:(n_layers-1)  # Exclude half-space for now
#         layer = layers[layer_idx]
#         layer_thickness = layer.thickness
#         vs = layer_properties[layer_idx].vs
#         ρ = layer_properties[layer_idx].ρ
#         μ = layer_properties[layer_idx].μ
#         q = layer_properties[layer_idx].q
        
#         # Find depths within this layer
#         layer_start = cumulative_depth
#         layer_end = cumulative_depth + layer_thickness
        
#         # Get state vector at top of this layer by propagating from surface
#         layer_top_matrix = Matrix{ComplexF64}(I, 2, 2)
#         for j in 1:(layer_idx-1)
#             layer_top_matrix = layer_matrices[j] * layer_top_matrix
#         end
#         state_at_top = layer_top_matrix * surface_vector
        
#         # Compute eigenfunctions within this layer
#         for depth_idx in current_depth_idx:n_depths
#             depth = depths[depth_idx]
            
#             if depth > layer_end
#                 break
#             end
            
#             # Depth within current layer
#             depth_in_layer = depth - layer_start
#             α = q * depth_in_layer
            
#             # Propagation matrix within layer
#             ca = cos(α)
#             sa = sin(α)
#             intra_layer_matrix = [
#                 ca              sa/(μ*q);
#                 -μ*q*sa         ca
#             ]
            
#             # State at this depth
#             state_vector = intra_layer_matrix * state_at_top
            
#             u_eigenfunction[depth_idx] = state_vector[1]
#             t_eigenfunction[depth_idx] = state_vector[2]
            
#             current_depth_idx = depth_idx + 1
#         end
        
#         cumulative_depth += layer_thickness
#     end
    
#     # Handle half-space
#     if current_depth_idx <= n_depths
#         # Get state at top of half-space
#         halfspace_top_matrix = Matrix{ComplexF64}(I, 2, 2)
#         for j in 1:(n_layers-1)
#             halfspace_top_matrix = layer_matrices[j] * halfspace_top_matrix
#         end
#         state_at_halfspace_top = halfspace_top_matrix * surface_vector
        
#         # Half-space properties
#         vs_hs = layer_properties[end].vs
#         ρ_hs = layer_properties[end].ρ
#         μ_hs = layer_properties[end].μ
#         q_hs = layer_properties[end].q
        
#         # Exponential decay in half-space
#         for depth_idx in current_depth_idx:n_depths
#             depth = depths[depth_idx]
#             depth_in_halfspace = depth - cumulative_depth
            
#             # Decaying solution
#             decay_factor = exp(-real(q_hs) * depth_in_halfspace)
#             phase_factor = exp(-1im * imag(q_hs) * depth_in_halfspace)
            
#             u_eigenfunction[depth_idx] = state_at_halfspace_top[1] * decay_factor * phase_factor
#             t_eigenfunction[depth_idx] = state_at_halfspace_top[2] * decay_factor * phase_factor
#         end
#     end
    
#     # Normalize eigenfunctions
#     u_real = real.(u_eigenfunction)
#     t_real = real.(t_eigenfunction)
    
#     # Energy normalization integral
#     energy_integral = compute_energy_integral(depths, u_real, t_real, model)
    
#     if energy_integral > 0
#         norm_factor = sqrt(energy_integral)
#         u_eigenfunction ./= norm_factor
#         t_eigenfunction ./= norm_factor
#     end
    
#     return depths, real.(u_eigenfunction), real.(t_eigenfunction), energy_integral
# end






# ╔═╡ 07792153-a4ec-44b2-8d4f-cd9620436331
# """
#     create_depth_sampling(model::EarthModel) -> Vector{Float64}

# Create optimal depth sampling for eigenfunction computation.
# """
# function create_depth_sampling(model::EarthModel)
#     depths = Float64[]
#     current_depth = 0.0
    
#     # Surface point
#     push!(depths, 0.0)
    
#     # Sample within each layer
#     for (i, layer) in enumerate(model.layers[1:end-1])  # Exclude half-space
#         layer_thickness = layer.thickness
        
#         # Adaptive sampling based on layer properties
#         if i <= 2
#             # Fine sampling in upper layers
#             n_samples = max(20, Int(ceil(layer_thickness / 1.0)))  # 1 km sampling
#         elseif i <= 5
#             # Medium sampling in middle layers  
#             n_samples = max(10, Int(ceil(layer_thickness / 2.0)))  # 2 km sampling
#         else
#             # Coarse sampling in deep layers
#             n_samples = max(5, Int(ceil(layer_thickness / 5.0)))   # 5 km sampling
#         end
        
#         # Create sampling points within layer
#         if layer_thickness > 0
#             layer_depths = range(current_depth, current_depth + layer_thickness, length=n_samples+1)
#             append!(depths, layer_depths[2:end])  # Skip first point to avoid duplication
#         end
        
#         current_depth += layer_thickness
#     end
    
#     # Add points into half-space for visualization
#     half_space_depths = current_depth .+ [2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 150.0]
#     append!(depths, half_space_depths)
    
#     return sort(unique(depths))
# end

# ╔═╡ d17163ee-57d8-43fe-aa96-42c271a7f4b3
# """
#     compute_energy_integral(depths, u_real, t_real, model) -> Float64

# Compute energy integral for eigenfunction normalization.
# This implements the energy velocity calculation from CPS.
# """
# function compute_energy_integral(depths::Vector{Float64}, u_real::Vector{Float64}, 
#                                 t_real::Vector{Float64}, model::EarthModel)
    
#     energy = 0.0
#     current_depth = 0.0
    
#     # Integrate through each layer
#     for (layer_idx, layer) in enumerate(model.layers[1:end-1])
#         layer_thickness = layer.thickness
#         ρ = layer.rho
#         μ = ρ * layer.vs^2
        
#         # Find depth indices for this layer
#         layer_start = current_depth
#         layer_end = current_depth + layer_thickness
        
#         # Integration within layer
#         for i in 2:length(depths)
#             if depths[i-1] >= layer_start && depths[i] <= layer_end
#                 dz = depths[i] - depths[i-1]
                
#                 # Kinetic energy density: (1/2) * ρ * u²
#                 kinetic_density_1 = 0.5 * ρ * u_real[i-1]^2
#                 kinetic_density_2 = 0.5 * ρ * u_real[i]^2
                
#                 # Strain energy density: (1/2) * t²/μ  
#                 strain_density_1 = 0.5 * t_real[i-1]^2 / μ
#                 strain_density_2 = 0.5 * t_real[i]^2 / μ
                
#                 # Trapezoidal integration
#                 layer_energy = 0.5 * (kinetic_density_1 + kinetic_density_2 + 
#                                      strain_density_1 + strain_density_2) * dz
                
#                 energy += layer_energy
#             end
#         end
        
#         current_depth += layer_thickness
#     end
    
#     # Add half-space contribution (if significant)
#     halfspace_start = current_depth
#     ρ_hs = model.layers[end].rho
#     μ_hs = ρ_hs * model.layers[end].vs^2
    
#     for i in 2:length(depths)
#         if depths[i-1] >= halfspace_start
#             dz = depths[i] - depths[i-1]
            
#             # Exponentially decaying contribution
#             kinetic_1 = 0.5 * ρ_hs * u_real[i-1]^2
#             kinetic_2 = 0.5 * ρ_hs * u_real[i]^2
#             strain_1 = 0.5 * t_real[i-1]^2 / μ_hs
#             strain_2 = 0.5 * t_real[i]^2 / μ_hs
            
#             halfspace_energy = 0.5 * (kinetic_1 + kinetic_2 + strain_1 + strain_2) * dz
#             energy += halfspace_energy
#         end
#     end
    
#     return energy
# end


# ╔═╡ c7d305d1-3ca5-4506-8163-cfc8954b0168
# """
#     generate_love_synthetic_seismogram(dispersion_result::DispersionResult, distance::Float64;
#                                       source_depth::Float64=10.0, dt::Float64=0.1, 
#                                       duration::Float64=1000.0, f0::Float64=0.1)

# Generate synthetic Love wave seismogram using modal summation with computed eigenfunctions.
# This implements the modal summation method similar to CPS tpulse96/spulse96.
# """
# function generate_love_synthetic_seismogram(dispersion_result::DispersionResult, distance::Float64;
#                                            source_depth::Float64=10.0, dt::Float64=0.1, 
#                                            duration::Float64=1000.0, f0::Float64=0.1)
    
#     # Time vector
#     npts = Int(duration / dt)
#     time = (0:npts-1) * dt
    
#     # Initialize seismogram  
#     displacement = zeros(Float64, npts)
    
#     # Ricker wavelet source time function
#     function ricker_wavelet(t, f0)
#         a = π * f0
#         return (1.0 - 2.0 * (a * t)^2) * exp(-(a * t)^2)
#     end
    
#     # Modal summation over all computed modes
#     for (i, period) in enumerate(dispersion_result.periods)
#         phase_vel = dispersion_result.phase_velocities[i]
#         group_vel = dispersion_result.group_velocities[i]
        
#         if isnan(phase_vel) || isnan(group_vel)
#             continue
#         end
        
#         frequency = 1.0 / period
#         ω = 2π * frequency
        
#         # Find source excitation from eigenfunction at source depth
#         if !isempty(dispersion_result.depths) && !isempty(dispersion_result.eigenfunctions_u)
#             # Interpolate eigenfunction at source depth
#             source_excitation = 1.0  # Default
#             if source_depth <= maximum(dispersion_result.depths)
#                 source_idx = argmin(abs.(dispersion_result.depths .- source_depth))
#                 u_source = dispersion_result.eigenfunctions_u[source_idx, i]
#                 source_excitation = abs(u_source)  # Magnitude of eigenfunction at source
#             end
#         else
#             source_excitation = exp(-source_depth / 50.0)  # Simple depth decay
#         end
        
#         # Travel times
#         group_arrival = distance / group_vel
#         phase_arrival = distance / phase_vel
        
#         # Amplitude factors
#         # Geometric spreading
#         amplitude = 1.0 / sqrt(max(distance, 1.0))
        
#         # Frequency-dependent amplitude (surface wave scaling)
#         amplitude *= 1.0 / frequency
        
#         # Attenuation (simple Q model)
#         Q = 100.0
#         attenuation = exp(-π * distance / (Q * phase_vel * period))
#         amplitude *= attenuation
        
#         # Source excitation
#         amplitude *= source_excitation
        
#         # Energy normalization
#         if !isnan(dispersion_result.energy_integrals[i]) && dispersion_result.energy_integrals[i] > 0
#             amplitude /= sqrt(dispersion_result.energy_integrals[i])
#         end
        
#         # Add this mode's contribution to seismogram
#         for (j, t) in enumerate(time)
#             if t >= group_arrival
#                 # Source time function delayed by group arrival
#                 source_time = t - group_arrival - 20.0  # 20s delay for source
#                 source_amplitude = ricker_wavelet(source_time, f0)
                
#                 # Phase-shifted sinusoid
#                 phase_shift = ω * phase_arrival
#                 wave_contribution = amplitude * source_amplitude * sin(ω * t - phase_shift)
                
#                 displacement[j] += wave_contribution
#             end
#         end
#     end
    
#     return time, displacement
# end

# ╔═╡ ac92ec5d-de77-4235-a701-50de13c502b1
"""
    dispersion_function_SH(model, ω, c)

Return complex-valued F(c) whose zero corresponds to a Love-wave mode:
    F(c) = M21 - R * M11
where M = product of layer matrices from top -> bottom and
R = -μ_b * q_b (bottom half-space traction/displacement ratio for decaying solution).
"""
function dispersion_function_SH(layers, ω::Float64, c::Float64)
    N = length(layers)
    if N < 1
        error("Model must contain at least one layer (half-space).")
    end

    # Build total propagator M (top -> bottom)
    M = Matrix{ComplexF64}(I, 2, 2)
    for L in layers
        M = layer_matrix_SH(L, ω, c) * M   # multiply top->down (pre-multiply or post? we choose pre)
    end

    # bottom (half-space) properties (assumed last entry)
    Lb = layers[end]
    vs_b = Lb.vs
    ρ_b = Lb.rho
    μ_b = ρ_b * vs_b^2

    # compute bottom vertical wavenumber q_b (for half-space decaying solution)
    q2_b = (ω / vs_b)^2 - (ω / c)^2
    q_b = sqrt(complex(q2_b))

    # radiation/decay relation tb/ub = R  (we take R = -μ_b * q_b)
    R = -μ_b * q_b

    # M maps [u_top; t_top] -> [u_bot; t_bot]. For t_top=0 (free surface), ub = M11*u0, tb = M21*u0.
    # impose tb/ub = R  => M21/M11 = R  => M21 - R*M11 = 0
    F = M[2, 1] - R * M[1, 1]
    return F
end


# ╔═╡ 03e6b940-ba12-41df-a3a2-2b4e0c44d64b
"""
    find_root_bracketed(f, cmin, cmax; nsteps=200)

Search for a sign change of real(f(c)) in [cmin, cmax] by sampling nsteps intervals.
Return (cL, cR) bracket where sign change occurs, or (nothing, nothing) if none found.
"""
function find_root_bracketed(f, cmin::Float64, cmax::Float64; nsteps::Int=200)
    xs = range(cmin, cmax; length=nsteps + 1)
    prev = real(f(first(xs)))
    for i in 2:length(xs)
        x = xs[i]
        val = real(f(x))
        if !isfinite(val)
            prev = val
            continue
        end
        if prev == 0.0
            return (xs[i-1], xs[i-1])
        end
        if sign(prev) != sign(val)
            return (xs[i-1], x)
        end
        prev = val
    end
    return (nothing, nothing)
end

# ╔═╡ e05989a8-f2f1-49ad-a17f-8ddbf5776ab4
"""
    solve_phase_velocity_love(model, T; c_search_pad=0.1)

Find a phase velocity at period T (seconds) for Love (SH) waves.
The search interval is set relative to model layer S velocities.
"""
function solve_phase_velocity_love(layers, T; c_search_pad=0.10)
    ω = 2π / T

    # set search bounds using layer shear velocities
    vs_vals = [L.vs for L in layers if L.vs > 0.0]
    if isempty(vs_vals)
        throw(ArgumentError("Model has no positive shear velocities."))
    end
    vmin = minimum(vs_vals)
    vmax = maximum(vs_vals)

    # expand search window by fraction
    cmin = max(1e-5, (1 - c_search_pad) * vmin)
    cmax = (1 + c_search_pad) * vmax * 1.5    # allow some margin above vmax
    @show cmin, cmax
    f(c) = dispersion_function_SH(layers, ω, c)

    # find bracket by sampling
    (L, R) = find_root_bracketed(f, cmin, cmax; nsteps=800)
    if L === nothing
        # No sign change found -> return NaN to indicate no root
        return NaN
    end

    # If bracketed exactly at same point (rare), return that point
    if L == R
        return L
    end

    # Use robust bisection on real part of f
    g(c) = real(f(c))
    root = find_zero(g, (L, R), Bisection(); tol=1e-6, maxevals=100)
    return root
end

# ╔═╡ b2add770-6eef-41f6-b4f1-ff93e299d408
"""
    compute_group_velocity_love(model::EarthModel, period::Float64, phase_velocity::Float64)

Compute Love wave group velocity using energy method.
"""
function compute_group_velocity_love(layers, period::Float64, phase_velocity::Float64)
    ω = 2π / period
    c = phase_velocity
    
    # Small perturbation for numerical derivative
    δω = ω * 1e-6
    
    # Compute phase velocities at ω ± δω
    c_plus = try
        solve_phase_velocity_love(layers, 2π/(ω + δω))
    catch
        NaN
    end
    
    c_minus = try
        solve_phase_velocity_love(layers, 2π/(ω - δω))
    catch
        NaN
    end
    
    if isnan(c_plus) || isnan(c_minus)
        # Fallback: approximate group velocity
        return 0.9 * c
    end
    
    # Group velocity using U = c - λ(dc/dλ) where λ = 2π/ω
    dc_dω = (c_plus - c_minus) / (2 * δω)
    group_velocity = c + ω * dc_dω
    
    return group_velocity
end

# ╔═╡ cca474e2-7152-43d5-be66-432847de2897
"""
    solve_love_dispersion(model, periods)

Compute Love-wave dispersion (phase velocity for each period) along with eigenfunctions.
Returns DispersionResult with phase velocities, group velocities, and eigenfunctions.
"""
function solve_love_dispersion(layers, periods)
    velocities = Float64[]
    group_velocities = Float64[]
    all_eigenfunctions_u = []
    all_eigenfunctions_t = []
    all_depths = Nothing
    energy_integrals = Float64[]
    
    for T in periods
        @printf("Solving T=%.3f s ... ", T)
        c = try
            c = solve_phase_velocity_love(layers, T)
        catch err
            @warn "error solving period $T: $err"
            NaN
        end
        
        if isnan(c)
            println("no root found")
            push!(velocities, NaN)
            push!(group_velocities, NaN)
            push!(all_eigenfunctions_u, Float64[])
            push!(all_eigenfunctions_t, Float64[])
            push!(energy_integrals, NaN)
        else
            println(@sprintf("c=%.5f km/s", c))
            push!(velocities, c)
            
            # Compute group velocity
            U = compute_group_velocity_love(layers, T, c)
            push!(group_velocities, U)
            
            # Compute eigenfunctions
            # depths, u_eigen, t_eigen, energy = compute_love_eigenfunctions(model, T, c)
            
            if all_depths === nothing
                all_depths = depths
            end
            
            # push!(all_eigenfunctions_u, u_eigen)
            # push!(all_eigenfunctions_t, t_eigen)
            # push!(energy_integrals, energy)
            
            @printf("  U=%.5f km/s\n", U)
        end
    end
	velocities = Velocities(periods, velocities, group_velocities)

	return velocities
		#all_depths, all_eigenfunctions_u, all_eigenfunctions_t, energy_integrals
    # # Convert eigenfunction arrays to matrices
    # if all_depths !== nothing
    #     n_depths = length(all_depths)
    #     n_periods = length(periods)
        
    #     eigenfunctions_u = zeros(Float64, n_depths, n_periods)
    #     eigenfunctions_t = zeros(Float64, n_depths, n_periods)
        
    #     for (i, period) in enumerate(periods)
    #         if length(all_eigenfunctions_u[i]) == n_depths
    #             eigenfunctions_u[:, i] = all_eigenfunctions_u[i]
    #             eigenfunctions_t[:, i] = all_eigenfunctions_t[i]
    #         end
    #     end
    # else
    #     all_depths = Float64[]
    #     eigenfunctions_u = zeros(Float64, 0, length(periods))
    #     eigenfunctions_t = zeros(Float64, 0, length(periods))
    # end
    
    # return DispersionResult(periods, velocities, group_velocities, 
    #                        eigenfunctions_u, eigenfunctions_t, all_depths, energy_integrals)
end


# ╔═╡ bcc394ad-284d-4299-891a-065ec7a6c6b3
"""
    solve_phase_velocity(model, T; c0_guess=nothing)

Find phase velocity at period T (s).
"""
function solve_phase_velocity(layers, T::Float64; c0_guess=nothing)
    ω = 2π / T

    # crude initial bracket around Vs
    vsurf = layers[1].vs
    cmin, cmax = 0.9 * vsurf, 1.2 * vsurf

    f(c) = real(dispersion_determinant(model, ω, c))

    # Find root of f(c) = 0
    croot = find_zero(f, (cmin, cmax), Bisection(), verbose=false)
    return croot
end

# ╔═╡ ab1fe653-c3bc-4da0-97c8-8a4dffb1ef01
md"## Posterior"

# ╔═╡ 84b9bd9b-c214-43d7-bc3b-3a36a5601b91
MvNormal(randn(10), I(10))

# ╔═╡ 814b9f2f-d9d2-4cea-9b75-1940a867944d
pdf(LogNormal(10), 1)

# ╔═╡ 27bc22cd-75f3-4f12-9de0-c9a125f2d548
2I(2)

# ╔═╡ b8bace04-1142-441d-a868-e528c7cd0cfc
rand(Poisson(10), 1)

# ╔═╡ 178dfb7e-b3d9-4336-b883-cef234aadb27
DiscreteUniform(5, 100) |> rand

# ╔═╡ 7f8bec68-bf38-4893-921b-189b8ca738eb
@model function invert_group_velocities(periods)
    # Priors
    N ~ DiscreteUniform(5, 100)
    h ~ MvLogNormal(zeros(N), 3I(N))
	vs ~ MvLogNormal(zeros(N), 3I(N))

	layers = create_layers(N, h, vs)
	
    # Likelihood
    velocities = solve_love_dispersion(layers, periods)
	σ = 1.0
    U ~ MvNormal(velocities.group_velocities, σ * σ * I(length(periods)))
end


# ╔═╡ ff5ab9ba-46dd-4761-98a6-b5d9127a4301
periods

# ╔═╡ 80086cb5-7e39-444b-a85e-ddf46a67ec42
md"## Appendix"

# ╔═╡ 5a43781a-b3e4-4bac-86d1-bd822c169804
const GUTENBERG_MODEL = [
    Layer(19.0, 6.14, 3.55, 2.74, 1000.0, 1000.0),
    Layer(19.0, 6.58, 3.80, 3.00, 1000.0, 1000.0),
    Layer(12.0, 8.20, 4.65, 3.32, 1000.0, 1000.0),
    Layer(10.0, 8.17, 4.62, 3.34, 1000.0, 1000.0),
    Layer(10.0, 8.14, 4.57, 3.35, 1000.0, 1000.0),
    Layer(10.0, 8.10, 4.51, 3.36, 1000.0, 1000.0),
    Layer(10.0, 8.07, 4.46, 3.37, 1000.0, 1000.0),
    Layer(10.0, 8.02, 4.41, 3.38, 1000.0, 1000.0),
    Layer(25.0, 7.93, 4.37, 3.39, 1000.0, 1000.0),
    Layer(25.0, 7.85, 4.35, 3.41, 1000.0, 1000.0),
    Layer(25.0, 7.89, 4.36, 3.43, 1000.0, 1000.0),
    Layer(25.0, 7.98, 4.38, 3.46, 1000.0, 1000.0),
    Layer(25.0, 8.10, 4.42, 3.48, 1000.0, 1000.0),
    Layer(25.0, 8.21, 4.46, 3.50, 1000.0, 1000.0),
    Layer(50.0, 8.38, 4.54, 3.53, 1000.0, 1000.0),
    Layer(50.0, 8.62, 4.68, 3.58, 1000.0, 1000.0),
    Layer(50.0, 8.87, 4.85, 3.62, 1000.0, 1000.0),
    Layer(50.0, 9.15, 5.04, 3.69, 1000.0, 1000.0),
    Layer(50.0, 9.45, 5.21, 3.82, 1000.0, 1000.0),
    Layer(100.0, 9.88, 5.45, 4.01, 1000.0, 1000.0),
    Layer(100.0, 10.30, 5.76, 4.21, 1000.0, 1000.0),
    Layer(100.0, 10.71, 6.03, 4.40, 1000.0, 1000.0),
    Layer(100.0, 11.10, 6.23, 4.56, 1000.0, 1000.0),
    Layer(100.0, 11.35, 6.32, 4.63, 1000.0, 1000.0)
]

# ╔═╡ 05d910d2-527c-416f-9d63-c259fbe8a45d
default_layers = GUTENBERG_MODEL;

# ╔═╡ 1af153bc-e471-4608-984b-61f1116dfa16
md"""
Number of layers: $(@bind n_layers Slider(1:length(default_layers), default=length(default_layers), show_value=true))
"""

# ╔═╡ 99e9fc2c-fe13-4f94-b613-86cedb0f3653
function layer_table_input(n_layers::Int; vp_vs_ratio=1.73)
	function gutenberg_layer_defaults(n_layers::Int)
		n_layers <= 1 && return [default_layers[end]]
		finite_layers = default_layers[1:end-1]
		edges = round.(Int, range(1, length(finite_layers) + 1, length=n_layers))
		grouped = Layer[]
		for i in 1:(n_layers - 1)
			block = finite_layers[edges[i]:(edges[i + 1] - 1)]
			h = sum(layer.thickness for layer in block)
			vp = sum(layer.vp * layer.thickness for layer in block) / h
			vs = sum(layer.vs * layer.thickness for layer in block) / h
			rho = sum(layer.rho * layer.thickness for layer in block) / h
			push!(grouped, Layer(h, vp, vs, rho, 100.0, 100.0))
		end
		return vcat(grouped, default_layers[end])
	end

	gutenberg_defaults = gutenberg_layer_defaults(n_layers)

    ui = PlutoUI.combine() do Child
        # header
        header = @htl("""
        <tr style="text-align:center;">
          <th style="padding:6px;">#</th>
          <th style="padding:6px;">Thickness (km)</th>
          <th style="padding:6px;">Vs (km/s)</th>
          <th style="padding:6px;">Density (gm/cc)</th>
        </tr>
        """)

		 rows = Any[]
        for i in 1:n_layers
            layer = gutenberg_defaults[i]
            tw = i < n_layers ? Child("thickness_$i", Slider(0.5:0.5:200, default=layer.thickness, show_value=true)) : "∞"
            vs = Child("vs_$i", Slider(2.5:0.01:8.0, default=layer.vs, show_value=true))
            rw = Child("rho_$i", Slider(1.0:0.01:6.0, default=layer.rho, show_value=true))

            push!(rows, @htl("""
                <tr>
                  <td style="text-align:center; padding:6px;"><b>Layer $i</b></td>
                  <td style="padding:6px; text-align:center;">$tw</td>
                  <td style="padding:6px; text-align:center;">$vs</td>
                  <td style="padding:6px; text-align:center;">$rw</td>
                </tr>
            """))
        end


        tbl = @htl("""
        <table style="border-collapse:collapse; border:1px solid #ddd; width:100%;">
          <thead>$header</thead>
          <tbody>
            $(rows...)
          </tbody>
        </table>
        """)

        md"""
        ### Layer Editor
        Number of layers: **$n_layers**

        $tbl
        _Notes:_  
        The last layer is treated as a half-space (no thickness).
        """
    end

    return PlutoUI.Experimental.transformed_value(ui) do vals
        layers = [
            begin
                if i < n_layers
                    t = vals[3*(i-1)+1]
                    vs = vals[3*(i-1)+2]
                    ρ = vals[3*(i-1)+3]
                else
                    # last layer: no thickness slider
                    t = 0.0
                    vs = vals[3*(i-1)+1]
                    ρ = vals[3*(i-1)+2]
                end

                vp = vp_vs_ratio * vs # (not used, anyway for love waves )

                Layer(t, vp, vs, ρ, 100.0, 100.0)
            end
            for i in 1:n_layers
        ]
        le = layers[end]
        # add artificial thick layer to impose radiation condition
        layer_new = Layer(1000, le.vp, le.vs, le.rho, 100.0, 100.0)
        layers = vcat(layers[1:end-1], layer_new, le)
    end
end

# ╔═╡ 7bce51e0-9191-42eb-8fe3-cc68f9b2681f
(@bind layers layer_table_input(n_layers))

# ╔═╡ 45579d35-df74-454f-91d4-03542f77e1ac
# ╠═╡ show_logs = false
res = solve_love_dispersion(layers, periods)

# ╔═╡ f3653043-33f5-4155-aec4-f412a4f305eb
begin
	struct ObservedDispersion
	    periods::Vector{Float64}
	    phase_velocities::Vector{Float64}
	    group_velocities::Vector{Float64}
	end
	
	empty_observed_dispersion() = ObservedDispersion(Float64[], Float64[], Float64[])
	load_observed_dispersion(::Missing; kind::Symbol=:auto) = empty_observed_dispersion()
	load_observed_dispersion(::Nothing; kind::Symbol=:auto) = empty_observed_dispersion()
	
	function _parse_observed_number(token::AbstractString)
	    t = strip(token)
	    isempty(t) && return NaN
	    lowercase(t) in ("nan", "missing", "none") && return NaN
	    return tryparse(Float64, t) === nothing ? NaN : parse(Float64, t)
	end
	
	function _numeric_rows_from_string(text::AbstractString)
	    rows = Vector{Vector{Float64}}()
	    for line in split(text, '\n')
	        body = strip(first(split(first(split(line, "#")), "//")))
	        isempty(body) && continue
	        tokens = split(replace(body, ',' => ' '))
	        nums = [_parse_observed_number(tok) for tok in tokens]
	        length(nums) >= 2 && push!(rows, nums)
	    end
	    return rows
	end
	
	function _observed_rows_from_path(path::AbstractString)
		path = expanduser(strip(path))
		isempty(path) && return Vector{Vector{Float64}}()
		isfile(path) || error("Observed dispersion file not found: $path")
		return _numeric_rows_from_string(read(path, String))
	end
	
	function _observed_rows_from_filepicker(file)
		data = if haskey(file, "data")
			file["data"]
		elseif haskey(file, :data)
			file[:data]
		else
			nothing
		end
		data === nothing && return Vector{Vector{Float64}}()
		text = data isa AbstractVector{UInt8} ? String(data) : String(data)
		return _numeric_rows_from_string(text)
	end
	
	function _dispersion_from_rows(rows; kind::Symbol=:auto)
		isempty(rows) && return empty_observed_dispersion()
	    period = Float64[]
	    phase = Float64[]
	    group = Float64[]
	    for row in rows
	        T = row[1]
	        isfinite(T) && T > 0 || continue
	
	        if length(row) >= 6
	            # pDSurfTomo flat file: period lat1 lon1 lat2 lon2 velocity.
	            v = row[6]
	            if kind == :phase
	                push!(period, T); push!(phase, v); push!(group, NaN)
	            elseif kind == :group
	                push!(period, T); push!(phase, NaN); push!(group, v)
	            else
	                push!(period, T); push!(phase, v); push!(group, NaN)
	            end
	        elseif length(row) >= 3
	            # Simple combined table: period phase_velocity group_velocity.
	            push!(period, T); push!(phase, row[2]); push!(group, row[3])
	        else
	            v = row[2]
	            if kind == :group
	                push!(period, T); push!(phase, NaN); push!(group, v)
	            else
	                push!(period, T); push!(phase, v); push!(group, NaN)
	            end
	        end
	    end
	    return ObservedDispersion(period, phase, group)
	end
	
	function load_observed_dispersion(path::AbstractString; kind::Symbol=:auto)
		return _dispersion_from_rows(_observed_rows_from_path(path); kind=kind)
	end
	
	function load_observed_dispersion(file::AbstractDict; kind::Symbol=:auto)
		return _dispersion_from_rows(_observed_rows_from_filepicker(file); kind=kind)
	end
	
	function merge_observed_dispersion(phase_obs::ObservedDispersion, group_obs::ObservedDispersion)
	    return ObservedDispersion(
	        vcat(phase_obs.periods, group_obs.periods),
	        vcat(phase_obs.phase_velocities, group_obs.phase_velocities),
	        vcat(phase_obs.group_velocities, group_obs.group_velocities),
	    )
	end
end

# ╔═╡ 227742dd-3cc0-4615-bde2-09dc7b890095
md"""
### Observed dispersion overlay

Phase file: $(@bind observed_phase_file FilePicker())

Group file: $(@bind observed_group_file FilePicker())
"""

# ╔═╡ 46e9bdbf-8618-4937-b789-c03fbe5396f7
observed_dispersion = merge_observed_dispersion(
	load_observed_dispersion(observed_phase_file; kind=:phase),
	load_observed_dispersion(observed_group_file; kind=:group),
)

# ╔═╡ 6e696382-0a37-46f1-8209-0bf5d7dbc82f
let
    # Filter out NaN values for plotting
    valid_indices = findall(.!isnan.(res.phase_velocities))
    observed_phase_indices = findall(i -> isfinite(observed_dispersion.phase_velocities[i]), eachindex(observed_dispersion.periods))
    observed_group_indices = findall(i -> isfinite(observed_dispersion.group_velocities[i]), eachindex(observed_dispersion.periods))
    
    if !isempty(valid_indices)
		function vs_depth_profile(layers)
			model_layers = length(layers) >= 3 && layers[end - 1].thickness > 500 ? vcat(layers[1:end-2], layers[end]) : layers
			positive_thickness = [layer.thickness for layer in model_layers if layer.thickness > 0]
			halfspace_thickness = isempty(positive_thickness) ? 25.0 : max(25.0, 0.15 * sum(positive_thickness))
			depths = Float64[]
			vs = Float64[]
			depth = 0.0
			for (i, layer) in enumerate(model_layers)
				h = (i == length(model_layers) || layer.thickness <= 0) ? halfspace_thickness : layer.thickness
				push!(depths, depth); push!(vs, layer.vs)
				depth += h
				push!(depths, depth); push!(vs, layer.vs)
			end
			return vs, depths
		end

		model_vs, model_depths = vs_depth_profile(layers)
		traces = typeof(scatter())[]
        
        # Phase velocity
        push!(traces,
			scatter(
                x=res.periods[valid_indices],
                y=res.phase_velocities[valid_indices],
                mode="lines+markers",
                name="Phase Velocity",
				xaxis="x",
				yaxis="y",
                line=attr(color="blue", width=3),
                marker=attr(size=6)
            )
        )

        if !isempty(observed_phase_indices)
            push!(traces,
				scatter(
                    x=observed_dispersion.periods[observed_phase_indices],
                    y=observed_dispersion.phase_velocities[observed_phase_indices],
                    mode="markers",
                    name="Observed Phase Velocity",
					xaxis="x",
					yaxis="y",
                    marker=attr(size=8, color="black", symbol="x")
                )
            )
        end
        
        # Group velocity
        push!(traces,
			scatter(
                x=res.periods[valid_indices],
                y=res.group_velocities[valid_indices],
                mode="lines+markers",
                name="Group Velocity", 
				xaxis="x2",
				yaxis="y2",
                line=attr(color="red", width=3),
                marker=attr(size=6)
            )
        )

        if !isempty(observed_group_indices)
            push!(traces,
				scatter(
                    x=observed_dispersion.periods[observed_group_indices],
                    y=observed_dispersion.group_velocities[observed_group_indices],
                    mode="markers",
                    name="Observed Group Velocity",
					xaxis="x2",
					yaxis="y2",
                    marker=attr(size=8, color="black", symbol="x")
                )
            )
        end

		push!(traces,
			scatter(
				x=model_vs,
				y=model_depths,
				mode="lines",
				name="Model Vs",
				xaxis="x3",
				yaxis="y3",
				line=attr(color="green", width=3, shape="hv")
			)
		)
        
        Plot(traces, Layout(
            title="Love Wave Dispersion Curves",
            xaxis=attr(domain=[0.0, 0.68], title="", showgrid=true, anchor="y"),
            yaxis=attr(domain=[0.56, 1.0], title="Phase Velocity (km/s)", showgrid=true, range=[2.5, 6.5], anchor="x"),
            xaxis2=attr(domain=[0.0, 0.68], title="Period (s)", showgrid=true, anchor="y2"),
            yaxis2=attr(domain=[0.0, 0.44], title="Group Velocity (km/s)", showgrid=true, range=[2.5, 6.5], anchor="x2"),
            xaxis3=attr(domain=[0.78, 1.0], title="Vs (km/s)", showgrid=true, anchor="y3"),
            yaxis3=attr(domain=[0.0, 1.0], title="Depth (km)", autorange="reversed", showgrid=true, anchor="x3", side="right", ticks="outside", showticklabels=true),
			legend=attr(orientation="h", x=0.5, xanchor="center", y=-0.16, yanchor="top"),
			margin=attr(l=70, r=40, t=60, b=100),
            showlegend=true,
            height=650
        ))
    else
        md"_No valid dispersion solutions found._"
    end
end

# ╔═╡ 40ed8766-c3db-4bb6-bcd2-abf54dbba726
res.group_velocities

# ╔═╡ 7c668163-2386-40dd-ad81-7a9db54ddeac
posterior = invert_group_velocities(periods) | (; U = res.group_velocities)

# ╔═╡ b223de7e-6ee4-4ce1-9647-6e024091a63d
chain = sample(posterior, NUTS(), 10)

# ╔═╡ fc1fd90c-020a-4f2d-aedd-66115ac6a287
md"""
### Reference
This notebook is inspired by the classical implementation of surface wave theory in  

Herrmann, R. B. (2013), *Computer Programs in Seismology*.

The Julia code here provides a **clean reimplementation** of some of the ideas behind CPS, designed for **interactive exploration**.
"""

# ╔═╡ a882824a-a4da-11f0-8360-85b7e5471cb3
md"""
### Eigenfunction Display

Select period for eigenfunction plot (s): $(@bind selected_period Slider(5:5:100, default=20, show_value=true))
"""

# ╔═╡ a8828434-a4da-11f0-ae4f-7707dd8da2ef
# let
#     if !isempty(res.depths) && !isempty(res.eigenfunctions_u)
#         # Find closest period to selected
#         period_idx = argmin(abs.(res.periods .- selected_period))
#         actual_period = res.periods[period_idx]
        
#         if !isnan(res.phase_velocities[period_idx])
#             u_eigen = res.eigenfunctions_u[:, period_idx]
#             t_eigen = res.eigenfunctions_t[:, period_idx]
            
#             fig = make_subplots(
#                 rows=1, cols=2,
#                 subplot_titles=["Displacement Eigenfunction U(z)", "Stress Eigenfunction T(z)"],
#                 horizontal_spacing=0.15
#             )
            
#             # Displacement eigenfunction
#             add_trace!(fig,
#                 scatter(
#                     x=u_eigen,
#                     y=res.depths,
#                     mode="lines+markers",
#                     name="U(z)",
#                     line=attr(color="blue", width=3),
#                     marker=attr(size=4)
#                 ),
#                 row=1, col=1
#             )
            
#             # Stress eigenfunction  
#             add_trace!(fig,
#                 scatter(
#                     x=t_eigen,
#                     y=res.depths,
#                     mode="lines+markers",
#                     name="T(z)",
#                     line=attr(color="red", width=3),
#                     marker=attr(size=4)
#                 ),
#                 row=1, col=2
#             )
            
#             # Add layer boundaries
#             cumulative_depth = 0.0
#             for (i, layer) in enumerate(layers[1:end-1])
#                 cumulative_depth += layer.thickness
                
#                 # Add vertical line at layer boundary
#                 add_trace!(fig,
#                     scatter(
#                         x=[minimum(u_eigen), maximum(u_eigen)],
#                         y=[cumulative_depth, cumulative_depth],
#                         mode="lines",
#                         line=attr(color="gray", width=1, dash="dash"),
#                         showlegend=false
#                     ),
#                     row=1, col=1
#                 )
                
#                 add_trace!(fig,
#                     scatter(
#                         x=[minimum(t_eigen), maximum(t_eigen)],
#                         y=[cumulative_depth, cumulative_depth],
#                         mode="lines",
#                         line=attr(color="gray", width=1, dash="dash"),
#                         showlegend=false
#                     ),
#                     row=1, col=2
#                 )
#             end
            
#             relayout!(fig,
#                 title="Love Wave Eigenfunctions (T=$(round(actual_period, digits=1)) s, c=$(round(res.phase_velocities[period_idx], digits=3)) km/s)",
#                 xaxis1=attr(title="Displacement U(z)", showgrid=true),
#                 xaxis2=attr(title="Stress T(z)", showgrid=true),
#                 yaxis1=attr(title="Depth (km)", autorange="reversed", showgrid=true),
#                 yaxis2=attr(title="", autorange="reversed", showgrid=true),
#                 showlegend=false,
#                 height=500
#             )
            
#             plot(fig)
#         else
#             md"_No valid eigenfunction for period $(actual_period) s_"
#         end
#     else
#         md"_No eigenfunction data available_"
#     end
# end

# ╔═╡ a88296c2-a4da-11f0-af15-77b6844f5837
# let
#     # Plot velocity model
#     depths = Float64[0]
#     vs_profile = Float64[]
#     vp_profile = Float64[]
#     rho_profile = Float64[]
    
#     current_depth = 0.0
#     for (i, layer) in enumerate(layers)
#         if i < length(layers)  # Not the half-space
#             push!(depths, current_depth)
#             push!(depths, current_depth + layer.thickness)
#             push!(vs_profile, layer.vs)
#             push!(vs_profile, layer.vs)
#             push!(vp_profile, layer.vp)
#             push!(vp_profile, layer.vp)
#             push!(rho_profile, layer.rho)
#             push!(rho_profile, layer.rho)
#             current_depth += layer.thickness
#         else  # Half-space
#             push!(depths, current_depth)
#             push!(depths, current_depth + 50)  # Arbitrary depth for visualization
#             push!(vs_profile, layer.vs)
#             push!(vs_profile, layer.vs)
#             push!(vp_profile, layer.vp)
#             push!(vp_profile, layer.vp)
#             push!(rho_profile, layer.rho)
#             push!(rho_profile, layer.rho)
#         end
#     end
    
#     fig = make_subplots(
#         rows=1, cols=3,
#         subplot_titles=["Shear Velocity", "P Velocity", "Density"],
#         horizontal_spacing=0.1
#     )
    
#     # Vs profile
#     add_trace!(fig,
#         scatter(x=vs_profile, y=depths, mode="lines",
#                line=attr(width=3, color="blue"), name="Vs"),
#         row=1, col=1)
    
#     # Vp profile
#     add_trace!(fig,
#         scatter(x=vp_profile, y=depths, mode="lines",
#                line=attr(width=3, color="red"), name="Vp"),
#         row=1, col=2)
    
#     # Density profile
#     add_trace!(fig,
#         scatter(x=rho_profile, y=depths, mode="lines", 
#                line=attr(width=3, color="green"), name="ρ"),
#         row=1, col=3)
    
#     relayout!(fig,
#         title="Earth Model - Velocity and Density Profiles",
#         xaxis1=attr(title="Vs (km/s)", showgrid=true),
#         xaxis2=attr(title="Vp (km/s)", showgrid=true),
#         xaxis3=attr(title="Density (g/cm³)", showgrid=true),
#         yaxis1=attr(title="Depth (km)", autorange="reversed", showgrid=true),
#         yaxis2=attr(title="", autorange="reversed", showgrid=true),
#         yaxis3=attr(title="", autorange="reversed", showgrid=true),
#         showlegend=false,
#         height=400
#     )
    
#     plot(fig)
# end

# ╔═╡ 3cce37c8-a4db-11f0-8445-2b16c3968493
md"""
## Synthetic Seismogram Generation

Generate synthetic Love wave seismograms using modal summation with the computed eigenfunctions:
"""

# ╔═╡ 3cce396c-a4db-11f0-bc5a-e107cd4639d6
md"""
**Distance (km):** $(@bind distance Slider(10.0:10.0:500.0, default=100.0, show_value=true))

**Source Depth (km):** $(@bind source_depth Slider(0.0:1.0:20.0, default=5.0, show_value=true))

**Source Frequency (Hz):** $(@bind source_freq Slider(0.05:0.01:0.2, default=0.1, show_value=true))

**Duration (s):** $(@bind duration Slider(500.0:100.0:2000.0, default=1000.0, show_value=true))
"""

# ╔═╡ 3cce3a82-a4db-11f0-9333-3da6f5a8cb24
@bind compute_synthetic Button("Generate Synthetic Seismogram")

# ╔═╡ 3cce3aca-a4db-11f0-aedc-edcc81a83eeb
# begin
#     compute_synthetic  # Trigger computation when button pressed
    
#     # Generate synthetic seismogram
#     time_vec, synthetic_displacement = generate_love_synthetic_seismogram(
#         res, distance;
#         source_depth=source_depth,
#         dt=0.1,
#         duration=duration,
#         f0=source_freq
#     )
    
#     md"✓ Generated synthetic Love wave seismogram for distance = $(distance) km"
# end

# ╔═╡ 3cce3bd8-a4db-11f0-baab-cf3610efa1c8
# # Plot synthetic seismogram
# let
#     fig = plot(
#         Layout(
#             title="Synthetic Love Wave Seismogram (Distance: $(distance) km)",
#             xaxis_title="Time (s)",
#             yaxis_title="Displacement (arbitrary units)",
#             height=400,
#             width=800
#         )
#     )
    
#     add_trace!(fig, scatter(
#         x=time_vec,
#         y=synthetic_displacement,
#         mode="lines",
#         name="Love Waves",
#         line=attr(color="black", width=1)
#     ))
    
#     # Add theoretical group velocity arrival times for reference
#     if any(.!isnan.(res.group_velocities))
#         valid_indices = .!isnan.(res.group_velocities)
#         min_arrival = minimum(distance ./ res.group_velocities[valid_indices])
#         max_arrival = maximum(distance ./ res.group_velocities[valid_indices])
        
#         # Add vertical lines for arrival time window
#         add_vline!(fig, min_arrival, line_color="red", line_dash="dash", annotation_text="Fastest Group Velocity")
#         add_vline!(fig, max_arrival, line_color="blue", line_dash="dash", annotation_text="Slowest Group Velocity")
#     end
    
#     fig
# end

# ╔═╡ 5b5d9fb8-2eb4-4b9a-8cf7-119efa543386
md"### Plotting"

# ╔═╡ 0c794e97-3f78-4f37-b370-2df91149ffd2
# Convert a layered model to stepwise depth/value arrays for plotting
function _step_profile(layers::Vector{Layer}, getval; zmax_cap::Float64=1e3)
    z, v = Float64[], Float64[]
    zcum = 0.0
    for ly in layers[1:end-2]
        h = isfinite(ly.thickness) ? ly.thickness : zmax_cap - zcum
        push!(z, zcum);         push!(v, getval(ly))
        zcum += h
        push!(z, zcum);         push!(v, getval(ly))
        if zcum >= zmax_cap; break; end
    end
    z, v
end

# ╔═╡ 88f6c30a-e522-4611-92dc-9590e037d128
""""
    plot_velocity_models(models; labels=[], zmax=500.0, plot_rho=false)

Plot one or more layered velocity models with Plotly step curves (Vp, Vs, optionally ρ).
- `models` : vector of layer stacks (e.g. `[model1, model2]`), each a `Vector{Layer}`
- `labels` : optional names per model
- `zmax`   : depth cap (km)
- `plot_rho` : also plot density if true
Returns a PlotlyJS `Plot`.
"""
function plot_velocity_models(models::Vector{Vector{Layer}}; labels=[], zmax=100000.0, plot_rho=false, plot_vp=false, title="")
    lbls = isempty(labels) ? ["model $(i)" for i in 1:length(models)] : labels
    traces = [scatter()]

    for (m, name) in zip(models, lbls)
       
        z_vs, vs = _step_profile(m, ly -> ly.vs; zmax_cap=zmax)
          push!(traces,
            scatter(x=vs, y=z_vs, mode="lines", line=attr(shape="hv", color="blue"), name="$name Vs",
                    hovertemplate="Vs=%{x:.2f} km/s<br>z=%{y:.1f} km"),
        )
		if(plot_vp)
			 z_vp, vp = _step_profile(m, ly -> ly.vp; zmax_cap=zmax)
			 push!(traces,
            scatter(x=vp, y=z_vp, mode="lines", line=attr(shape="hv", dash="dash"), name="$name Vp",
                    hovertemplate="Vp=%{x:.2f} km/s<br>z=%{y:.1f} km"))
		end
        if plot_rho
            z_rho, rho = _step_profile(m, ly -> ly.rho; zmax_cap=zmax)
            push!(traces,
                scatter(x=rho, y=z_rho, mode="lines", line=attr(shape="hv", dash="dot"),
                        name="$name ρ", hovertemplate="ρ=%{x:.2f} g/cc<br>z=%{y:.1f} km")
            )
        end
    end

    layout = Layout(
        xaxis_title = "Velocity (km/s)",
		width=250,
		title=title,
		height=500,
        yaxis = attr(title="Depth (km)", autorange="reversed"),
        legend = attr(orientation="h"),
        margin = attr(l=60, r=20, t=30, b=60)
    )
    plot(traces, layout)
end

# ╔═╡ cbff7199-4249-49dd-a7bb-0d2e52d0d611
plot_velocity_models([layers])

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"

[compat]
FFTW = "~1.10.0"
HypertextLiteral = "~0.9.5"
PlutoPlotly = "~0.6.5"
PlutoUI = "~0.7.72"
Roots = "~2.2.10"
Turing = "~0.42.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "bef6b3a1b484c103287794a3328723d6f113337d"

[[deps.ADTypes]]
git-tree-sha1 = "8b2b045b22740e4be20654175cc38291d48539db"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.20.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractMCMC]]
deps = ["BangBang", "ConsoleProgressMonitor", "Distributed", "FillArrays", "LogDensityProblems", "Logging", "LoggingExtras", "ProgressLogging", "Random", "StatsBase", "TerminalLoggers", "Transducers", "UUIDs"]
git-tree-sha1 = "64c15bc9be4a9ab588d8b90c5ebf633ef67afb73"
uuid = "80f14c24-f653-4e6a-9b94-39d6b0f70001"
version = "5.10.0"

[[deps.AbstractPPL]]
deps = ["AbstractMCMC", "Accessors", "DensityInterface", "JSON", "LinearAlgebra", "Random", "StatsBase"]
git-tree-sha1 = "66aed89871b6a8458bb407327fa997d980ea894f"
uuid = "7a57a42e-76ec-4ea3-a279-07e840d6d9cf"
version = "0.13.6"
weakdeps = ["Distributions"]

    [deps.AbstractPPL.extensions]
    AbstractPPLDistributionsExt = ["Distributions"]

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdvancedHMC]]
deps = ["AbstractMCMC", "ArgCheck", "DocStringExtensions", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "ProgressMeter", "Random", "Setfield", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "15e1bafe97eb0eeb77c8dae63cdcfa133986b254"
uuid = "0bf59076-c3b1-5ca4-86bd-e02cd72cde3d"
version = "0.8.3"

    [deps.AdvancedHMC.extensions]
    AdvancedHMCADTypesExt = "ADTypes"
    AdvancedHMCCUDAExt = "CUDA"
    AdvancedHMCComponentArraysExt = "ComponentArrays"
    AdvancedHMCMCMCChainsExt = "MCMCChains"
    AdvancedHMCOrdinaryDiffEqExt = "OrdinaryDiffEq"

    [deps.AdvancedHMC.weakdeps]
    ADTypes = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ComponentArrays = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
    MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
    OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"

[[deps.AdvancedMH]]
deps = ["AbstractMCMC", "Distributions", "DocStringExtensions", "FillArrays", "LinearAlgebra", "LogDensityProblems", "Random", "Requires"]
git-tree-sha1 = "c9be6a2d91cfe681f2cf454749c00ea5ceb8f07e"
uuid = "5b7e9947-ddc0-4b3f-9b55-0d8042f74170"
version = "0.8.9"
weakdeps = ["DiffResults", "ForwardDiff", "MCMCChains", "StructArrays"]

    [deps.AdvancedMH.extensions]
    AdvancedMHForwardDiffExt = ["DiffResults", "ForwardDiff"]
    AdvancedMHMCMCChainsExt = "MCMCChains"
    AdvancedMHStructArraysExt = "StructArrays"

[[deps.AdvancedPS]]
deps = ["AbstractMCMC", "Distributions", "Random", "Random123", "Requires", "SSMProblems", "StatsFuns"]
git-tree-sha1 = "5d34d826ece67ce790d4a7f3f97d837e52aba7f8"
uuid = "576499cb-2369-40b2-a588-c64705576edc"
version = "0.7.0"
weakdeps = ["Libtask"]

    [deps.AdvancedPS.extensions]
    AdvancedPSLibtaskExt = "Libtask"

[[deps.AdvancedVI]]
deps = ["ADTypes", "Accessors", "ChainRulesCore", "DiffResults", "DifferentiationInterface", "Distributions", "DocStringExtensions", "FillArrays", "Functors", "LinearAlgebra", "LogDensityProblems", "Optimisers", "ProgressMeter", "Random", "StatsBase"]
git-tree-sha1 = "9196795a7474c48ec7d73056a52a786e61dac062"
uuid = "b5ca4192-6429-45e5-a2d9-87aec30a685c"
version = "0.6.0"

    [deps.AdvancedVI.extensions]
    AdvancedVIBijectorsExt = ["Bijectors", "Optimisers"]
    AdvancedVIEnzymeExt = ["Enzyme", "ChainRulesCore"]
    AdvancedVIMooncakeExt = ["Mooncake", "ChainRulesCore"]
    AdvancedVIReverseDiffExt = ["ReverseDiff", "ChainRulesCore"]

    [deps.AdvancedVI.weakdeps]
    Bijectors = "76274a88-744f-5084-9051-94815aaf08c4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgCheck]]
git-tree-sha1 = "f9e9a66c9b7be1ad7372bbd9b062d9230c30c5ce"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.BangBang]]
deps = ["Accessors", "ConstructionBase", "InitialValues", "LinearAlgebra"]
git-tree-sha1 = "a49f9342fc60c2a2aaa4e0934f06755464fcf438"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.4.6"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTablesExt = "Tables"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijectors]]
deps = ["ArgCheck", "ChainRulesCore", "ChangesOfVariables", "Distributions", "DocStringExtensions", "Functors", "InverseFunctions", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "MappedArrays", "Random", "Reexport", "Roots", "SparseArrays", "Statistics"]
git-tree-sha1 = "bbdaafb28e31f2948916670a8eb7fd94ec12406f"
uuid = "76274a88-744f-5084-9051-94815aaf08c4"
version = "0.15.14"

    [deps.Bijectors.extensions]
    BijectorsDistributionsADExt = "DistributionsAD"
    BijectorsEnzymeCoreExt = "EnzymeCore"
    BijectorsForwardDiffExt = "ForwardDiff"
    BijectorsLazyArraysExt = "LazyArrays"
    BijectorsMooncakeExt = "Mooncake"
    BijectorsReverseDiffChainRulesExt = ["ChainRules", "ReverseDiff"]
    BijectorsReverseDiffExt = "ReverseDiff"

    [deps.Bijectors.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    DistributionsAD = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "SparseInverseSubset", "Statistics", "StructArrays", "SuiteSparse"]
git-tree-sha1 = "3b704353e517a957323bd3ac70fa7b669b5f48d4"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.72.6"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Chairmarks]]
deps = ["Printf", "Random"]
git-tree-sha1 = "9a49491e67e7a4d6f885c43d00bb101e6e5a434b"
uuid = "0ca39b1e-fe0b-4e98-acfc-b1656634c4de"
version = "1.3.1"
weakdeps = ["Statistics"]

    [deps.Chairmarks.extensions]
    StatisticsChairmarksExt = ["Statistics"]

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "3aa4bf1532aa2e14e0374c4fd72bed9a9d0d0f6c"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.10"
weakdeps = ["InverseFunctions", "Test"]

    [deps.ChangesOfVariables.extensions]
    ChangesOfVariablesInverseFunctionsExt = "InverseFunctions"
    ChangesOfVariablesTestExt = "Test"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "8010b6bb3388abe68d95743dcbea77650bb2eddf"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.3"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

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

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "80bd15222b3e8d0bc70d921d2201aa0084810ce5"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.12"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3bc002af51045ca3b47d2e1787d6ce02e68b943a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.122"
weakdeps = ["ChainRulesCore", "DensityInterface", "Test"]

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

[[deps.DistributionsAD]]
deps = ["Adapt", "ChainRules", "ChainRulesCore", "Compat", "Distributions", "FillArrays", "LinearAlgebra", "PDMats", "Random", "Requires", "SpecialFunctions", "StaticArrays", "StatsFuns", "ZygoteRules"]
git-tree-sha1 = "4acbf909e892ce1f94c39a138541566c1aad5e66"
uuid = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
version = "0.6.58"

    [deps.DistributionsAD.extensions]
    DistributionsADForwardDiffExt = "ForwardDiff"
    DistributionsADLazyArraysExt = "LazyArrays"
    DistributionsADReverseDiffExt = "ReverseDiff"
    DistributionsADTrackerExt = "Tracker"

    [deps.DistributionsAD.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPPL]]
deps = ["ADTypes", "AbstractMCMC", "AbstractPPL", "Accessors", "BangBang", "Bijectors", "Chairmarks", "Compat", "ConstructionBase", "DifferentiationInterface", "Distributions", "DocStringExtensions", "InteractiveUtils", "LinearAlgebra", "LogDensityProblems", "MacroTools", "OrderedCollections", "Printf", "Random", "Statistics", "Test"]
git-tree-sha1 = "b87bc878118011066248821177ba26b62a999248"
uuid = "366bfd00-2699-11ea-058f-f148b4cae6d8"
version = "0.39.1"

    [deps.DynamicPPL.extensions]
    DynamicPPLChainRulesCoreExt = ["ChainRulesCore"]
    DynamicPPLEnzymeCoreExt = ["EnzymeCore"]
    DynamicPPLForwardDiffExt = ["ForwardDiff"]
    DynamicPPLJETExt = ["JET"]
    DynamicPPLMCMCChainsExt = ["MCMCChains"]
    DynamicPPLMarginalLogDensitiesExt = ["MarginalLogDensities"]
    DynamicPPLMooncakeExt = ["Mooncake"]

    [deps.DynamicPPL.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    JET = "c3a54625-cd67-489e-a8e7-0a5a0ff4e31b"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
    MarginalLogDensities = "f0c3360a-fb8d-11e9-1194-5521fd7ee392"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.EllipticalSliceSampling]]
deps = ["AbstractMCMC", "ArrayInterface", "Distributions", "Random", "Statistics"]
git-tree-sha1 = "e611b7fdfbfb5b18d5e98776c30daede41b44542"
uuid = "cad2338a-1db2-11e9-3401-43bc07c9ede2"
version = "2.0.0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "5bfcd42851cf2f1b303f51525a54dc5e98d408a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.15.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "9340ca07ca27093ff68418b7558ca37b05f8aeb1"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.29.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cd33c7538e68650bd0ddbb3f5bd50a4a0fa95b50"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.0"
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

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "ba51324b894edaf1df3ab16e2cc6bc3280a2f1a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.10"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "95ba48564903b43b2462318aa243ee79d81135ff"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libtask]]
deps = ["MistyClosures", "Test"]
git-tree-sha1 = "a88cfb11c45f350bbe51df05a1af7c691d58306d"
uuid = "6f1fad26-d15e-5dc8-ae53-837a1d7b8c9f"
version = "0.9.10"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Printf"]
git-tree-sha1 = "9ea3422d03222c6de679934d1c08f0a99405aa03"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.5.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogDensityProblems]]
deps = ["ArgCheck", "DocStringExtensions", "Random"]
git-tree-sha1 = "d9625f27ded4ad726ceca7819394a4cc77ed25b3"
uuid = "6fdf6af0-433a-55f7-b3ed-c6c6e0b8df7c"
version = "2.2.0"

[[deps.LogDensityProblemsAD]]
deps = ["DocStringExtensions", "LogDensityProblems"]
git-tree-sha1 = "7b83f3ad0a8105f79a067cafbfd124827bb398d0"
uuid = "996a588d-648d-4e1f-a8f0-a84b347e47b1"
version = "1.13.1"

    [deps.LogDensityProblemsAD.extensions]
    LogDensityProblemsADADTypesExt = "ADTypes"
    LogDensityProblemsADDifferentiationInterfaceExt = ["ADTypes", "DifferentiationInterface"]
    LogDensityProblemsADEnzymeExt = "Enzyme"
    LogDensityProblemsADFiniteDifferencesExt = "FiniteDifferences"
    LogDensityProblemsADForwardDiffBenchmarkToolsExt = ["BenchmarkTools", "ForwardDiff"]
    LogDensityProblemsADForwardDiffExt = "ForwardDiff"
    LogDensityProblemsADReverseDiffExt = "ReverseDiff"
    LogDensityProblemsADTrackerExt = "Tracker"
    LogDensityProblemsADZygoteExt = "Zygote"

    [deps.LogDensityProblemsAD.weakdeps]
    ADTypes = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
    BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
    DifferentiationInterface = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"
weakdeps = ["ChainRulesCore", "ChangesOfVariables", "InverseFunctions"]

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MCMCChains]]
deps = ["AbstractMCMC", "AxisArrays", "DataAPI", "Dates", "Distributions", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "4f5b84761bbfd1e99c2568815b1108858b760f4c"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "7.6.0"

[[deps.MCMCDiagnosticTools]]
deps = ["AbstractFFTs", "DataAPI", "DataStructures", "Distributions", "LinearAlgebra", "MLJModelInterface", "Random", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "526c98cd41028da22c01cb8a203246799ad853a8"
uuid = "be115224-59cd-429b-ad48-344e309966f0"
version = "0.3.15"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MLJModelInterface]]
deps = ["InteractiveUtils", "REPL", "Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "ccaa3f7938890ee8042cc970ba275115428bd592"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.12.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MicroCollections]]
deps = ["Accessors", "BangBang", "InitialValues"]
git-tree-sha1 = "44d32db644e84c75dab479f1bc15ee76a1a3618f"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.2.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.MistyClosures]]
git-tree-sha1 = "d1a692e293c2a0dc8fda79c04cad60582f3d4de3"
uuid = "dbe65cb8-6be2-42dd-bbc5-4196aaced4f4"
version = "2.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NamedArrays]]
deps = ["Combinatorics", "DelimitedFiles", "InvertedIndices", "LinearAlgebra", "OrderedCollections", "Random", "Requires", "SparseArrays", "Statistics"]
git-tree-sha1 = "33d258318d9e049d26c02ca31b4843b2c851c0b0"
uuid = "86f7a689-2022-50b4-a561-43c23ac3c673"
version = "0.10.5"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "48968edaf014f67e58fe4c8a4ce72d392aed3294"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.3"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Optimisers]]
deps = ["ChainRulesCore", "ConstructionBase", "Functors", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "131dc319e7c58317e8c6d5170440f6bdaee0a959"
uuid = "3bd65402-5787-11e9-1adc-39752487f4e2"
version = "0.4.6"

    [deps.Optimisers.extensions]
    OptimisersAdaptExt = ["Adapt"]
    OptimisersEnzymeCoreExt = "EnzymeCore"
    OptimisersReactantExt = "Reactant"

    [deps.Optimisers.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Reactant = "3c362404-f566-11ee-1572-e11a4b42c853"

[[deps.Optimization]]
deps = ["ADTypes", "ArrayInterface", "ConsoleProgressMonitor", "DocStringExtensions", "LinearAlgebra", "Logging", "LoggingExtras", "OptimizationBase", "Printf", "Reexport", "SciMLBase", "SparseArrays", "TerminalLoggers"]
git-tree-sha1 = "9c8a9bc453fa1b5a5913ffa77d4b9b8b586f1de8"
uuid = "7f7a1694-90dd-40f0-9382-eb1efda571ba"
version = "5.2.0"

[[deps.OptimizationBase]]
deps = ["ADTypes", "ArrayInterface", "DifferentiationInterface", "DocStringExtensions", "FastClosures", "LinearAlgebra", "PDMats", "Reexport", "SciMLBase", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings"]
git-tree-sha1 = "96e15dffd0499eba632c10a4d864b98071042809"
uuid = "bca83a33-5cc9-4baa-983d-23429ab6bcbb"
version = "4.0.2"

    [deps.OptimizationBase.extensions]
    OptimizationEnzymeExt = "Enzyme"
    OptimizationFiniteDiffExt = "FiniteDiff"
    OptimizationForwardDiffExt = "ForwardDiff"
    OptimizationMLDataDevicesExt = "MLDataDevices"
    OptimizationMLUtilsExt = "MLUtils"
    OptimizationMTKExt = "ModelingToolkit"
    OptimizationReverseDiffExt = "ReverseDiff"
    OptimizationSymbolicAnalysisExt = "SymbolicAnalysis"
    OptimizationZygoteExt = "Zygote"

    [deps.OptimizationBase.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLDataDevices = "7e8f7934-dd98-4c1a-8fe8-92b47a384d40"
    MLUtils = "f1d291b0-491e-4a28-83b9-f70985020b54"
    ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SymbolicAnalysis = "4297ee4d-0239-47d8-ba5d-195ecdf594fe"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.OptimizationOptimJL]]
deps = ["Optim", "OptimizationBase", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays"]
git-tree-sha1 = "ba8f23a14d861a303e829afaf67a7b3f75da945f"
uuid = "36348300-93cb-4f02-beb5-3c3902f8871e"
version = "0.4.8"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d922b4d80d1e12c658da7785e754f4796cc1d60d"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.36"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "28278bb0053da0fd73537be94afd1682cc5a0a83"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.21"

    [deps.PlotlyBase.extensions]
    DataFramesExt = "DataFrames"
    DistributionsExt = "Distributions"
    IJuliaExt = "IJulia"
    JSON3Ext = "JSON3"

    [deps.PlotlyBase.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    JSON3 = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "ColorSchemes", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "PrecompileTools", "Reexport", "ScopedValues", "Scratch", "TOML"]
git-tree-sha1 = "8acd04abc9a636ef57004f4c2e6f3f6ed4611099"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.5"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "c05b4c6325262152483a1ecb6c69846d2e01727b"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.34"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "c5a07210bd060d6a8491b0ccdee2fa0235fc00bf"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.1.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "f0803bc1171e455a04124affa9c21bba5ac4db32"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.6"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "dbe5fd0b334694e905cb9fda73cd8554333c46e2"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.1"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "51bdb23afaaa551f923a0e990f7c44a4451a26f1"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.39.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "8a433b1ede5e9be9a7ba5b1cc6698daa8d718f1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "2f609ec2295c452685d3142bc4df202686e555d2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.16"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SSMProblems]]
deps = ["AbstractMCMC", "Distributions", "Random"]
git-tree-sha1 = "5dd0431563784b468db06335ce8777653092d621"
uuid = "26aad666-b158-4e64-9d35-0e672562fa48"
version = "0.5.2"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLLogging", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "b06382e2f1ff0c9aff4f88f560673c800df6099f"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.128.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDifferentiationInterfaceExt = "DifferentiationInterface"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DifferentiationInterface = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLLogging]]
deps = ["Logging", "LoggingExtras", "Preferences"]
git-tree-sha1 = "1713f77f37f6809cd5b051197ea9adf139101ae5"
uuid = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
version = "1.6.0"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "9deb07bb4e6f95712bda5617a717caed701fbe4b"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.14.0"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLPublic]]
git-tree-sha1 = "ed647f161e8b3f2973f24979ec074e8d084f1bee"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.0"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "c3b2323466378a2ba15bea4b2f73b081e022f473"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "322365aa23098275562cbad6a1c2539ee40d9618"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "1.1.3"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerChainRulesCoreExt = "ChainRulesCore"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseInverseSubset]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "52962839426b75b3021296f7df242e40ecfc0852"
uuid = "dc90abb0-5640-4711-901d-7e5b23a2fada"
version = "0.1.2"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "6ed48d9a3b22417c765dc273ae3e1e4de035e7c8"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.23"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"
    SparseMatrixColoringsJuMPExt = ["JuMP", "MathOptInterface"]

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "b8693004b385c842357406e3af647701fe783f98"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.15"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "89f86d9376acd18a1a4fbef66a56335a3a7633b8"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.5.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "be5733d4a2b03341bdcab91cea6caa7e31ced14b"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.9"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a3c1536470bf8c5e02096ad4853606d7c8f62721"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.2"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"
weakdeps = ["PrettyTables"]

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "f133fab380933d042f6796eda4e130272ba520ca"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.7"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Transducers]]
deps = ["Accessors", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "SplittablesBase", "Tables"]
git-tree-sha1 = "4aa1fdf6c1da74661f6f5d3edfd96648321dade9"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.85"

    [deps.Transducers.extensions]
    TransducersAdaptExt = "Adapt"
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.Turing]]
deps = ["ADTypes", "AbstractMCMC", "AbstractPPL", "Accessors", "AdvancedHMC", "AdvancedMH", "AdvancedPS", "AdvancedVI", "BangBang", "Bijectors", "Compat", "DataStructures", "Distributions", "DistributionsAD", "DocStringExtensions", "DynamicPPL", "EllipticalSliceSampling", "ForwardDiff", "Libtask", "LinearAlgebra", "LogDensityProblems", "MCMCChains", "NamedArrays", "Optimization", "OptimizationOptimJL", "OrderedCollections", "Printf", "Random", "Reexport", "SciMLBase", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "69b42b532b4a92c6093e46be2b1a9d62a155489d"
uuid = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"
version = "0.42.0"

    [deps.Turing.extensions]
    TuringDynamicHMCExt = "DynamicHMC"

    [deps.Turing.weakdeps]
    DynamicHMC = "bbc10e6e-7c05-544b-b16e-64fede858acb"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "434b3de333c75fc446aa0d19fc394edafd07ab08"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.7"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═ee482657-fd34-4c92-95cd-0bf8e254676c
# ╟─d3bc6646-b22f-4e22-87a8-da37a9aa1d1a
# ╟─1af153bc-e471-4608-984b-61f1116dfa16
# ╟─7bce51e0-9191-42eb-8fe3-cc68f9b2681f
# ╠═45579d35-df74-454f-91d4-03542f77e1ac
# ╠═f3653043-33f5-4155-aec4-f412a4f305eb
# ╟─227742dd-3cc0-4615-bde2-09dc7b890095
# ╠═46e9bdbf-8618-4937-b789-c03fbe5396f7
# ╟─6e696382-0a37-46f1-8209-0bf5d7dbc82f
# ╠═cbff7199-4249-49dd-a7bb-0d2e52d0d611
# ╠═d2c48242-92ce-11f0-01e6-cb1b8fe1fecb
# ╠═cee3bbac-9de8-46e7-8ba9-2cc1fe18802c
# ╠═363d8b6d-ddc9-4048-a02b-16346158112c
# ╠═6c766c4b-8a1c-4f8b-9467-fec5a91fcd05
# ╠═5510a81e-1ac0-45b9-9657-13f0db32f3e1
# ╟─7e9acd7c-9cf5-45e0-8892-77c77a1fd983
# ╠═756e06e3-4f25-4f55-b242-26a45a1b7dd3
# ╠═9c9ad5f6-debf-4376-aaf6-759e9b0c949b
# ╠═4f122bf3-0693-4f13-b716-c0dc2863af25
# ╠═07792153-a4ec-44b2-8d4f-cd9620436331
# ╠═d17163ee-57d8-43fe-aa96-42c271a7f4b3
# ╠═b2add770-6eef-41f6-b4f1-ff93e299d408
# ╠═c7d305d1-3ca5-4506-8163-cfc8954b0168
# ╠═e05989a8-f2f1-49ad-a17f-8ddbf5776ab4
# ╠═cca474e2-7152-43d5-be66-432847de2897
# ╠═ac92ec5d-de77-4235-a701-50de13c502b1
# ╠═03e6b940-ba12-41df-a3a2-2b4e0c44d64b
# ╠═bcc394ad-284d-4299-891a-065ec7a6c6b3
# ╟─ab1fe653-c3bc-4da0-97c8-8a4dffb1ef01
# ╠═84b9bd9b-c214-43d7-bc3b-3a36a5601b91
# ╠═814b9f2f-d9d2-4cea-9b75-1940a867944d
# ╠═27bc22cd-75f3-4f12-9de0-c9a125f2d548
# ╠═b8bace04-1142-441d-a868-e528c7cd0cfc
# ╠═178dfb7e-b3d9-4336-b883-cef234aadb27
# ╠═7f8bec68-bf38-4893-921b-189b8ca738eb
# ╠═40ed8766-c3db-4bb6-bcd2-abf54dbba726
# ╠═ff5ab9ba-46dd-4761-98a6-b5d9127a4301
# ╠═7c668163-2386-40dd-ad81-7a9db54ddeac
# ╠═b223de7e-6ee4-4ce1-9647-6e024091a63d
# ╟─80086cb5-7e39-444b-a85e-ddf46a67ec42
# ╠═8d8f440d-4b79-4835-841b-739b5171f979
# ╠═29afa36b-391f-4861-aead-d62918daf3c6
# ╠═8d5d9594-2197-4ccc-a7a4-0f0d54a2370a
# ╠═2e91160e-33fa-4e2f-876e-8cccfebd27a5
# ╠═99e9fc2c-fe13-4f94-b613-86cedb0f3653
# ╠═5a43781a-b3e4-4bac-86d1-bd822c169804
# ╠═05d910d2-527c-416f-9d63-c259fbe8a45d
# ╟─fc1fd90c-020a-4f2d-aedd-66115ac6a287
# ╠═a882824a-a4da-11f0-8360-85b7e5471cb3
# ╠═a8828434-a4da-11f0-ae4f-7707dd8da2ef
# ╠═a88296c2-a4da-11f0-af15-77b6844f5837
# ╠═3cce37c8-a4db-11f0-8445-2b16c3968493
# ╠═3cce396c-a4db-11f0-bc5a-e107cd4639d6
# ╠═3cce3a82-a4db-11f0-9333-3da6f5a8cb24
# ╠═3cce3aca-a4db-11f0-aedc-edcc81a83eeb
# ╠═3cce3bd8-a4db-11f0-baab-cf3610efa1c8
# ╟─5b5d9fb8-2eb4-4b9a-8cf7-119efa543386
# ╠═0c794e97-3f78-4f37-b370-2df91149ffd2
# ╠═88f6c30a-e522-4611-92dc-9590e037d128
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
