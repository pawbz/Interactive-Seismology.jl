### A Pluto.jl notebook ###
# v0.20.19

#> [frontmatter]
#> title = "Fourier Series And Dispersion"
#> tags = ["introduction"]
#> layout = "layout.jlhtml"
#> description = "Exploration of Fourier series, wave dispersion"

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

# ╔═╡ 903844b8-a500-11f0-afa6-f187b28c7860
TableOfContents()

# ╔═╡ 90384550-a500-11f0-9976-15113ffba311
md"""
# Fourier Series And Dispersion

This basic interactive notebook will help you understand

- **Fourier Series**: How complex waveforms are built from simple sinusoids and 
- **Dispersion**: How different frequencies travel at different speeds.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,  
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 90384aca-a500-11f0-9af9-738a0287bdbb
let
    # Create time and space grids
    t = 0.0  # Fixed time for now
    x = range(-2π, 2π, length=1000)

	wave_type = "impulse"
    
    # Get Fourier coefficients
    coeffs = get_fourier_coefficients(wave_type, n_harmonics)
    
    # Get true/analytical signal
    true_signal = get_true_signal(wave_type, x, f0)
    
    # Extract phase shifts from UI (handle case where phase_shifts might not be available)
    phases = zeros(n_harmonics)
    # if show_phase_controls && @isdefined(phase_shifts)
        for n in 1:n_harmonics
            phase_key = "phase_$n"
            if haskey(phase_shifts, phase_key)
                phases[n] = phase_shifts[phase_key]
            end
        end
    # end
    
    # Initialize arrays for plotting
    individual_waves = []
    cumulative_sum = zeros(length(x))
    
    # Create plot
    fig = plot(Layout(
        title="Fourier Series Decomposition",
        xaxis_title="Position (radians)",
        yaxis_title="Amplitude",
        height=600,
        showlegend=true
    ))
    
    # Plot individual harmonics (show up to 10 for visual clarity)
    colors = ["red", "blue", "green", "orange", "purple", "brown", "pink", "gray", "olive", "cyan"]
    display_limit = min(n_harmonics, 10)  # Limit display for clarity
    
    for n in 1:display_limit
        if coeffs[n] ≠ 0
            # Apply phase shift: cos(n*f0*x + phase) where phase = -f0*n*time_shift
            phase_shift = -f0 * n * phases[n]  # Convert time shift to phase
            harmonic = coeffs[n] * cos.(n * f0 * x .+ phase_shift)
            cumulative_sum .+= harmonic
            
            add_trace!(fig, scatter(
                x=x, y=harmonic,
                mode="lines",
                name="Harmonic $n (Δt=$(phases[n])s)",
                line=attr(color=colors[mod(n-1, length(colors))+1], width=1, dash="dot"),
                opacity=0.75
            ))
        end
    end
    
    # Add remaining harmonics to cumulative sum (without plotting individually)
    for n in (display_limit+1):n_harmonics
        if coeffs[n] ≠ 0
            phase_shift = -f0 * n * phases[n]
            harmonic = coeffs[n] * cos.(n * f0 * x .+ phase_shift)
            cumulative_sum .+= harmonic
        end
    end
    
    # Plot true signal (only if no phase shifts applied)
   
        add_trace!(fig, scatter(
            x=x, y=true_signal,
            mode="lines",
            name="True Signal",
            line=attr(color="darkred", width=3, dash="dash"),
            opacity=0.8
        ))
    # end
    
    # Plot cumulative sum (Fourier series approximation)
    add_trace!(fig, scatter(
        x=x, y=cumulative_sum,
        mode="lines",
        name="Fourier Series (N=$n_harmonics)",
        line=attr(color="black", width=3, dash="solid")
    ))
    
    fig
end

# ╔═╡ 903847d0-a500-11f0-b347-cb097b21036e
md"""
**Number of harmonics:** $(@bind n_harmonics Slider(1:20, default=5, show_value=true))

"""

# ╔═╡ 93948b6b-d2ce-40c5-976d-0e565ddc0ca4
@bind phase_shifts harmonic_phase_input(n_harmonics; default_shift=0.0, max_shift=2.0)

# ╔═╡ d8f3a2b0-a586-11f0-9c84-3b2f5e1a8d7c
md"""
## Group vs Phase Velocity Demonstration

This section illustrates the fundamental difference between **phase velocity** and **group velocity** using a narrow frequency band - key concepts in wave dispersion.

**Time:** $(@bind time_demo Slider(0:0.1:20, default=0, show_value=true)) seconds

**Center Frequency:** $(@bind fc_demo Slider(0.1:0.1:2.0, default=0.5, show_value=true)) Hz

**Frequency Bandwidth:** $(@bind df_demo Slider(0.01:0.01:0.3, default=0.1, show_value=true)) Hz

**Phase Velocity 1:** $(@bind vp_demo Slider(2:0.1:8, default=4, show_value=true)) m/s

**Group Velocity:** $(@bind vg_demo Slider(1:0.1:6, default=2.5, show_value=true)) m/s
"""

# ╔═╡ e9f4b3c0-a586-11f0-8a75-4c3d6f2e9b8f
let
    # Create spatial grid
    x = range(-50, 50, length=1000)
    
    # Parameters from the mathematical formula
    Δω = 2π * df_demo  # Angular frequency bandwidth
    ω0 = 2π * fc_demo  # Center angular frequency
    k_n_ω0 = ω0 / vp_demo  # Wave number at center frequency
    
    # Calculate the derivative dk_n/dω (inverse of group velocity)
    # For demonstration, we'll use the relation: dk/dω = 1/vg
    dk_dω = 1 / vg_demo
    
    # Calculate Y parameter from the formula
    Y = (Δω / 2) * (time_demo .- dk_dω * x)
    
    # Implement the exact formula from equation (7.11, Aki and Richards):
    # f₀(x,t) ≈ (Δω sin Y)/(π Y) * cos[ω₀t - k_n(ω₀)x]
    
    # Handle sinc function properly (sinc(0) = 1)
    sinc_Y = zeros(length(Y))
    for i in eachindex(Y)
        if abs(Y[i]) < 1e-10
            sinc_Y[i] = 1.0
        else
            sinc_Y[i] = sin(Y[i]) / Y[i]
        end
    end
    
    # Calculate the envelope (sinc function)
    envelope = abs.(Δω * sinc_Y / π)
    
    # Calculate the carrier wave
    carrier = cos.(ω0 * time_demo .- k_n_ω0 * x)
    
    # Total wave using exact formula
    wave_total = envelope .* carrier
    
    # For comparison, also calculate the two-frequency superposition
    f1 = fc_demo - df_demo/2
    f2 = fc_demo + df_demo/2
    k1 = 2π * f1 / vp_demo
    k2 = 2π * f2 / vp_demo
    wave1 = cos.(2π * f1 * time_demo .- k1 * x)
    wave2 = cos.(2π * f2 * time_demo .- k2 * x)
    wave_superposition = wave1 + wave2
    
    # Calculate theoretical positions
    xg_pos = vg_demo * time_demo  # Group velocity position
    xp_pos = vp_demo * time_demo  # Phase velocity position
    
    # Create plot
    fig = plot(Layout(
        title="Group vs Phase Velocity (Exact Sinc Formula)<br>Time = $(round(time_demo, digits=1)) s",
        xaxis_title="Position (m)",
        yaxis_title="Amplitude",
        height=600,
        showlegend=true
    ))
    
    # Plot envelope (sinc function)
    add_trace!(fig, scatter(
        x=x, y=envelope,
        mode="lines",
        name="Envelope: (Δω/π)sinc(Y)",
        line=attr(color="red", width=3, dash="dash"),
        opacity=0.8
    ))
    
    # Plot carrier wave (faded)
    add_trace!(fig, scatter(
        x=x, y=0.5 * carrier,
        mode="lines",
        name="Carrier: cos(ω₀t - k₀x)",
        line=attr(color="blue", width=1, dash="dot"),
        opacity=0.4
    ))
    
    # Plot total wave (exact sinc formula)
    add_trace!(fig, scatter(
        x=x, y=wave_total,
        mode="lines",
        name="f₀(x,t) = Envelope × Carrier",
        line=attr(color="black", width=2)
    ))
    
    # # Plot two-frequency superposition for comparison
    # add_trace!(fig, scatter(
    #     x=x, y=wave_superposition,
    #     mode="lines",
    #     name="Two-frequency sum (comparison)",
    #     line=attr(color="gray", width=1, dash="dashdot"),
    #     opacity=0.6
    # ))
    
    # Mark group velocity position (envelope maximum)
    add_trace!(fig, scatter(
        x=[xg_pos, xg_pos],
        y=[-maximum(envelope), maximum(envelope)],
        mode="lines",
        name="Group Velocity Position",
        line=attr(color="red", width=4),
        opacity=0.7
    ))
    
    # Mark phase velocity position (carrier wave crest)
    phase_marker_y = 0.5 * cos(ω0 * time_demo - k_n_ω0 * xp_pos)
    add_trace!(fig, scatter(
        x=[xp_pos],
        y=[phase_marker_y],
        mode="markers",
        name="Phase Velocity (Wave Crest)",
        marker=attr(color="blue", size=12, symbol="circle")
    ))
    
    # Add vertical line for phase velocity
    add_trace!(fig, scatter(
        x=[xp_pos, xp_pos],
        y=[-maximum(envelope), maximum(envelope)],
        mode="lines",
        name="Phase Velocity Position", 
        line=attr(color="blue", width=2, dash="dashdot"),
        opacity=0.7
    ))
    
    fig
end

# ╔═╡ f1a2c4d0-a586-11f0-9b86-5d4e7f3a9c9d
md"""
### Key Observations:

**Phase Velocity (Blue)**: 
- Tracks individual wave crests
- Speed = $(vp_demo) m/s
- Controls propagation of individual oscillations

**Group Velocity (Red)**: 
- Tracks the envelope 
- Speed = $(vg_demo) m/s  
- Controls propagation of wave energy/information

**Experiment**: Try setting $v_g < v_p$ or $v_g > v_p$ to see different dispersion scenarios!
"""

# ╔═╡ 8411b089-85f5-4d27-850e-9747896eea36
md"---"

# ╔═╡ 903846f6-a500-11f0-9bec-955ede71b788
md"""
Any periodic function can be represented as a sum of sine and cosine waves:

$$f(x) = a_0 + \sum_{n=1}^{\infty} [a_n \cos(nx) + b_n \sin(nx)]$$

Let's build complex waveforms by adding simple harmonics!
"""

# ╔═╡ 2d41bce0-2035-4afa-8846-74061ac94333
md"## Appendix"

# ╔═╡ c31dd2bc-562a-4b04-86fa-99e3d241fdde
md"### UI"

# ╔═╡ 384fa2f0-38fb-464e-b72e-2ec26276f1cf
"""
Create a table-style input for harmonic phase shifts with sliders
"""
function harmonic_phase_input(n_harmonics::Int; default_shift=0.0, max_shift=2.0)
    ui = PlutoUI.combine() do Child
        # header
        header = @htl("""
        <tr style="text-align:center;">
          <th style="padding:8px; border:1px solid #dee2e6;">#</th>
          <th style="padding:8px; border:1px solid #dee2e6;">Frequency (Hz)</th>
          <th style="padding:8px; border:1px solid #dee2e6;">Time Shift (s)</th>
        </tr>
        """)

        rows = Any[]
        for n in 1:n_harmonics
            freq = n * f0
            shift_slider = Child("phase_$n", Slider(-max_shift:0.1:max_shift, default=default_shift, show_value=true))
            

            # Alternate row colors for better readability
            row_style = ""
            
            push!(rows, @htl("""
                <tr style="$row_style">
                  <td style="text-align:center; padding:8px; border:1px solid #dee2e6; font-weight:bold;">$n</td>
                  <td style="text-align:center; padding:8px; border:1px solid #dee2e6;">$(round(freq, digits=2))</td>
                  <td style="padding:8px; text-align:center; border:1px solid #dee2e6;">$shift_slider</td>
                </tr>
            """))
        end

        # Add summary info row
        summary_row = @htl("""
            <tr style="background-color:#e3f2fd; border:2px solid #90caf9;">
              <td colspan="4" style="text-align:center; padding:8px; font-style:italic; color:#0d47a1;">
                📊 Total harmonics: $n_harmonics | Range: ±$(max_shift)s | Step: 0.1s
              </td>
            </tr>
        """)
        push!(rows, summary_row)

        tbl = @htl("""
        <table style="border-collapse:collapse; border:2px solid #dee2e6; width:100%; margin:10px 0;">
          <thead>$header</thead>
          <tbody>
            $(rows...)
          </tbody>
        </table>
        """)

        md"""
        #### Individual Harmonic Phase Control
        
        $tbl
        
        **Instructions:**
        - **Positive values** (+): Harmonic arrives **later** (slower velocity)
        - **Negative values** (−): Harmonic arrives **earlier** (faster velocity)  
        - **Zero values**: No phase shift (original timing)
        
        **Experiment Ideas:**
        - **Normal Dispersion**: Set progressive delays (0, 0.2, 0.4, 0.6...)
        - **Anomalous Dispersion**: Set progressive advances (0, -0.2, -0.4, -0.6...)
        - **Random Medium**: Use random values to simulate complex media
        """
    end

    return PlutoUI.Experimental.transformed_value(ui) do vals
        # Create dictionary of phase shifts for all harmonics
        phase_dict = Dict{String, Float64}()
        
        for n in 1:n_harmonics
            phase_dict["phase_$n"] = vals[n]
        end
        
        return phase_dict
    end
end

# ╔═╡ 7f5b0439-0ad7-4ab6-83ca-093fc2d48915
using HypertextLiteral: @htl

# ╔═╡ 90384046-a500-11f0-8c2f-5ff4e899d2de
using PlutoUI, PlutoPlotly, FFTW, LinearAlgebra

# ╔═╡ c8806df8-6aab-4018-9c5e-7230195e8580
f0 = 0.6

# ╔═╡ 5c99f48d-94f6-46e8-b3a3-6c7a2a5c8b7b
"""Generate Fourier series coefficients for different wave types"""
function get_fourier_coefficients(wave_type, n_harmonics)
    if wave_type == "square"
        # Square wave: only odd harmonics, amplitude = 4/(π*n)
        return [n % 2 == 1 ? 4/(π*n) : 0.0 for n in 1:n_harmonics]
    elseif wave_type == "sawtooth" 
        # Sawtooth wave: all harmonics, amplitude = 2/(π*n) * (-1)^(n+1)
        return [2/(π*n) * (-1)^(n+1) for n in 1:n_harmonics]
    elseif wave_type == "triangle"
        # Triangle wave: odd harmonics only, amplitude = 8/(π²*n²) * (-1)^((n-1)/2)
        return [n % 2 == 1 ? 8/(π^2*n^2) * (-1)^((n-1)÷2) : 0.0 for n in 1:n_harmonics]
    elseif wave_type == "impulse"
        # Impulse (Dirac comb): all harmonics with equal amplitude = 2/T = 2*f0
        # For visualization, we scale it down to make it reasonable
        return [2*f0 * 0.1 for n in 1:n_harmonics]  # Scale factor 0.1 for visibility
    else # custom
        # Custom: decreasing amplitude with some randomness
        return [1/n * (1 + 0.3*sin(n)) for n in 1:n_harmonics]
    end
end

# ╔═╡ 48f67b56-3ff0-4c5f-a694-36ef9d09793f
"""Generate true/analytical signal for comparison"""
function get_true_signal(wave_type, x, f0)
    period = 2π / f0
    x_mod = mod.(x, period)  # Map to one period
    
    if wave_type == "square"
        # Square wave: +1 for first half period, -1 for second half
        return [x_val < period/2 ? 1.0 : -1.0 for x_val in x_mod]
    elseif wave_type == "sawtooth"
        # Sawtooth wave: linear ramp from -1 to +1
        return 2 * (x_mod ./ period) .- 1
    elseif wave_type == "triangle"
        # Triangle wave: ramp up then down
        return [x_val < period/2 ? 4*x_val/period - 1 : 3 - 4*x_val/period for x_val in x_mod]
    elseif wave_type == "impulse"
        # Impulse train: narrow spikes at period intervals
        spike_width = period * 0.05  # 5% of period
        return [abs(x_val) < spike_width || abs(x_val - period) < spike_width ? 2.0 : 0.0 for x_val in x_mod]
    else # custom
        # Custom: a more complex but smooth function
        return sin.(x) + 0.3*sin.(3*x) + 0.1*sin.(5*x)
    end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
FFTW = "~1.10.0"
HypertextLiteral = "~0.9.5"
PlutoPlotly = "~0.6.5"
PlutoUI = "~0.7.72"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "c38ba138082eed4cee4c442553757650e05d8ac7"

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

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

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
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

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
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

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

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

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
version = "1.11.0"
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
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
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
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═903844b8-a500-11f0-afa6-f187b28c7860
# ╟─90384550-a500-11f0-9976-15113ffba311
# ╟─90384aca-a500-11f0-9af9-738a0287bdbb
# ╟─903847d0-a500-11f0-b347-cb097b21036e
# ╟─93948b6b-d2ce-40c5-976d-0e565ddc0ca4
# ╟─d8f3a2b0-a586-11f0-9c84-3b2f5e1a8d7c
# ╟─e9f4b3c0-a586-11f0-8a75-4c3d6f2e9b8f
# ╟─f1a2c4d0-a586-11f0-9b86-5d4e7f3a9c9d
# ╟─8411b089-85f5-4d27-850e-9747896eea36
# ╟─903846f6-a500-11f0-9bec-955ede71b788
# ╟─2d41bce0-2035-4afa-8846-74061ac94333
# ╟─c31dd2bc-562a-4b04-86fa-99e3d241fdde
# ╠═384fa2f0-38fb-464e-b72e-2ec26276f1cf
# ╠═7f5b0439-0ad7-4ab6-83ca-093fc2d48915
# ╠═90384046-a500-11f0-8c2f-5ff4e899d2de
# ╠═c8806df8-6aab-4018-9c5e-7230195e8580
# ╠═5c99f48d-94f6-46e8-b3a3-6c7a2a5c8b7b
# ╠═48f67b56-3ff0-4c5f-a694-36ef9d09793f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
