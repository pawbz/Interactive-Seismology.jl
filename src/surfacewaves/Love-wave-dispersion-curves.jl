### A Pluto.jl notebook ###
# v0.20.19

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


# ╔═╡ d058a019-40c2-4546-b578-808381506d1f
struct EarthModel
    layers::Vector{Layer}
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
struct PhaseVelocities
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

# ╔═╡ 756e06e3-4f25-4f55-b242-26a45a1b7dd3
periods = collect(range(5.0, stop=100.0, length=64))

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
function dispersion_function_SH(model::EarthModel, ω::Float64, c::Float64)
    layers = model.layers
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
function solve_phase_velocity_love(model, T; c_search_pad=0.10)
    ω = 2π / T

    # set search bounds using layer shear velocities
    vs_vals = [L.vs for L in model.layers if L.vs > 0.0]
    if isempty(vs_vals)
        throw(ArgumentError("Model has no positive shear velocities."))
    end
    vmin = minimum(vs_vals)
    vmax = maximum(vs_vals)

    # expand search window by fraction
    cmin = max(1e-5, (1 - c_search_pad) * vmin)
    cmax = (1 + c_search_pad) * vmax * 1.5    # allow some margin above vmax
    @show cmin, cmax
    f(c) = dispersion_function_SH(model, ω, c)

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
function compute_group_velocity_love(model::EarthModel, period::Float64, phase_velocity::Float64)
    ω = 2π / period
    c = phase_velocity
    
    # Small perturbation for numerical derivative
    δω = ω * 1e-6
    
    # Compute phase velocities at ω ± δω
    c_plus = try
        solve_phase_velocity_love(model, 2π/(ω + δω))
    catch
        NaN
    end
    
    c_minus = try
        solve_phase_velocity_love(model, 2π/(ω - δω))
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
function solve_love_dispersion(model, periods)
    velocities = Float64[]
    group_velocities = Float64[]
    all_eigenfunctions_u = []
    all_eigenfunctions_t = []
    all_depths = Nothing
    energy_integrals = Float64[]
    
    for T in periods
        @printf("Solving T=%.3f s ... ", T)
        c = try
            c = solve_phase_velocity_love(model, T)
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
            U = compute_group_velocity_love(model, T, c)
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
	velocities = PhaseVelocities(periods, velocities, group_velocities)

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
function solve_phase_velocity(model::EarthModel, T::Float64; c0_guess=nothing)
    ω = 2π / T

    # crude initial bracket around Vs
    vsurf = model.layers[1].vs
    cmin, cmax = 0.9 * vsurf, 1.2 * vsurf

    f(c) = real(dispersion_determinant(model, ω, c))

    # Find root of f(c) = 0
    croot = find_zero(f, (cmin, cmax), Bisection(), verbose=false)
    return croot
end

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
            # Use Gutenberg model values if available, else fallback
            layer = i <= length(default_layers) ? default_layers[i] : default_layers[end]
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

# ╔═╡ a3cd5bb3-9b59-48a5-a505-4199a04941ed
model = EarthModel(layers)

# ╔═╡ 45579d35-df74-454f-91d4-03542f77e1ac
# ╠═╡ show_logs = false
res = solve_love_dispersion(model, periods)

# ╔═╡ 6e696382-0a37-46f1-8209-0bf5d7dbc82f
let
    # Filter out NaN values for plotting
    valid_indices = findall(.!isnan.(res.phase_velocities))
    
    if !isempty(valid_indices)
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=["Phase Velocity" "Group Velocity"],
            vertical_spacing=0.15
        )
        
        # Phase velocity
        add_trace!(fig,
            scatter(
                x=res.periods[valid_indices],
                y=res.phase_velocities[valid_indices],
                mode="lines+markers",
                name="Phase Velocity",
                line=attr(color="blue", width=3),
                marker=attr(size=6)
            ),
            row=1, col=1
        )
        
        # Group velocity
        add_trace!(fig,
            scatter(
                x=res.periods[valid_indices],
                y=res.group_velocities[valid_indices],
                mode="lines+markers",
                name="Group Velocity", 
                line=attr(color="red", width=3),
                marker=attr(size=6)
            ),
            row=2, col=1
        )
        
        relayout!(fig,
            title="Love Wave Dispersion Curves",
            xaxis1=attr(title="", showgrid=true),
            xaxis2=attr(title="Period (s)", showgrid=true),
            yaxis=attr(title="Phase Velocity (km/s)", showgrid=true, range=[2.5, 6.5]),
            yaxis2=attr(title="Group Velocity (km/s)", showgrid=true, range=[2.5, 6.5]),
            showlegend=false,
            height=600
        )
        
        fig
    else
        md"_No valid dispersion solutions found._"
    end
end

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

[compat]
FFTW = "~1.9.0"
HypertextLiteral = "~0.9.5"
PlutoPlotly = "~0.6.4"
PlutoUI = "~0.7.71"
Roots = "~2.2.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.0"
manifest_format = "2.0"
project_hash = "5e79032ef0cc11d1fdbce4d6e3738ec22f3c9500"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

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

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

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

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "797762812ed063b9b94f6cc7742bc8883bb5e69e"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.9.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

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

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

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

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.11.1+1"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

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
version = "1.12.0"
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
git-tree-sha1 = "232630fee92e588c11c2b260741b4fa70784b4c5"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.4"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

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

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "442b4353ee8c26756672afb2db81894fc28811f3"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.6"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

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

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.13.1+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.5.0+2"
"""

# ╔═╡ Cell order:
# ╠═ee482657-fd34-4c92-95cd-0bf8e254676c
# ╟─d3bc6646-b22f-4e22-87a8-da37a9aa1d1a
# ╟─1af153bc-e471-4608-984b-61f1116dfa16
# ╟─7bce51e0-9191-42eb-8fe3-cc68f9b2681f
# ╟─6e696382-0a37-46f1-8209-0bf5d7dbc82f
# ╠═45579d35-df74-454f-91d4-03542f77e1ac
# ╠═d2c48242-92ce-11f0-01e6-cb1b8fe1fecb
# ╠═d058a019-40c2-4546-b578-808381506d1f
# ╠═363d8b6d-ddc9-4048-a02b-16346158112c
# ╠═6c766c4b-8a1c-4f8b-9467-fec5a91fcd05
# ╠═5510a81e-1ac0-45b9-9657-13f0db32f3e1
# ╠═a3cd5bb3-9b59-48a5-a505-4199a04941ed
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
# ╟─80086cb5-7e39-444b-a85e-ddf46a67ec42
# ╠═8d8f440d-4b79-4835-841b-739b5171f979
# ╠═29afa36b-391f-4861-aead-d62918daf3c6
# ╠═8d5d9594-2197-4ccc-a7a4-0f0d54a2370a
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
