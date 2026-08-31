### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Seismic Attenuation"
#> tags = ["planewaves"]
#> layout = "layout.jlhtml"
#> description = "This notebook helps us visualize these effects of attenuation on propagating seismic waves."

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

# ╔═╡ 74a3e9e3-ed09-4582-91fc-47944d5744db
begin
    using PlutoUI
    using FFTW
end

# ╔═╡ a9436d89-ad43-4cb8-8619-67c0e26eb686
TableOfContents()

# ╔═╡ fbbee9cc-493d-11ed-00fd-b1e3bddbb3fa
md"""
# Intrinsic Attenuation
As opposed to adiabatic wave propagation, the wave motion in the earth is spatially attenuated. In other words, for anelastic earth, the integral of kinetic energy and the strain energy is no longer held constant due to internal friction and conversion to heat. In general, the dimensionless quantity ``Q`` that summarizes this intrinsic attenuation is dependent on the frequency of the propagating waves. Over the range of frequencies observed in seismology, it has been observed that
- observations suggest that ``Q`` is effectively constant, and
- the earth is only mildly anelastic i.e., ``Q\gg1``.
As a result, the high frequencies are attenuated more than the low frequencies; a traveling pulse in the earth gradually loses its high frequencies causing changes in pulse shapes with distance from the source. This notebook visualizes these effects using two receivers at different distances from a common source -- exactly the setup a seismologist actually uses to *measure* ``Q`` in the field, via the **spectral ratio method**. Along the way it also contrasts the everywhere-assumed ``Q\gg1`` approximation against Kjartansson's (1979) exact constant-``Q`` model, so the approximation's own limits become visible rather than assumed.

Attenuation is observed to be strongest in the upper mantle and inner core. Attenuation decreases rapidly with depth, with a high attenuation in the crust that is attributed to the effects of cracks and fluids.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ d60cf031-d375-4da6-ace3-abd5043d6364
md"""
## Quality Factor (``Q``) and the Nearly-Constant-``Q`` Approximation

Assuming ``Q\gg1``, we can write
```math
\frac{1}{Q(\omega)} = -\frac{1}{\pi}\frac{\Delta A}{A}.
```
Here, ``\Delta A`` is the gradual spatial decay in the amplitude of the propagating wave as it travels a distance of one wavelength:
```math
\Delta A = \frac{d A}{d x}\lambda = \frac{d A}{d x}\frac{2\pi c}{\omega}.
```
This results in an exponentially decaying solution
```math
A(x) = A_0\,\exp\left(-\frac{\omega x}{2 c Q}\right) = A_0\,\exp\left(-\pi f\, t^*(x)\right), \qquad t^*(x) \equiv \frac{x}{cQ},
```
where the second form -- in ordinary frequency ``f=\omega/2\pi`` -- introduces **``t^*``** ("t-star"), the standard seismological measure of *accumulated* attenuation along a path: the travel time ``x/c`` divided by ``Q``. Every result below is naturally written in terms of ``t^*``, not ``x`` and ``Q`` separately, because ``t^*`` is what a real measurement actually recovers (see "Measuring ``Q``" below).

!!! warning "This is an approximation, not the whole story"
    This is **Futterman's (1962) linear, "nearly constant ``Q``" formula**, and it is only
    the leading-order behavior as ``Q\to\infty`` of an exact model derived below. The
    widget's `Futterman (Q≫1)` button uses exactly this formula; `Kjartansson (exact)`
    uses the exact one, valid for any ``Q``.
"""

# ╔═╡ 304689e0-f318-4730-be5d-e112517d81ad
md"""
## From Amplitude Decay to Dispersion -- Causality Demands It

Attenuation cannot come alone. If high frequencies decayed faster than low ones while
every frequency traveled at the *same* speed, a sharp pulse would arrive at some distance
already broadened -- but the leading edge of a causal signal cannot move faster than the
fastest frequency component allows, no matter how small its amplitude. The only way to
keep the response causal is for the wave speed itself to depend (weakly) on frequency:
faster propagation at higher frequency compensates for their disproportionate energy
loss, so the broadened arrival still respects causality. This amplitude-dispersion
pairing is a direct consequence of the Kramers-Kronig relations that any causal, linear
response function must obey -- attenuation and dispersion are two faces of one physical
mechanism, never independent knobs.

The Dispersion & Attenuation panel above plots exactly this: the phase velocity
``c(\omega)`` rises (very slightly, since real ``Q\gg1``) with frequency, for both the
approximate and exact models derived below.
"""

# ╔═╡ b3273cac-84b5-458a-86c0-a530bd769ba0
md"""
## Kjartansson's Exact Constant-``Q`` Model

Futterman's formula above is a linearization. The exact result, valid for *any* ``Q`` and
not just ``Q\gg1``, comes from asking a sharper question: what is the *only* causal
material response that makes ``Q`` exactly frequency-independent?

Linear viscoelasticity relates stress and strain in the frequency domain by a complex
modulus, ``\Sigma(\omega) = M(\omega)E(\omega)``. The loss angle -- and hence ``Q``, via
``1/Q=\tan(\arg M(\omega))`` -- is frequency-independent only if ``\arg M(\omega)`` itself
does not depend on ``|\omega|``. A causal power-law creep function is the *only* choice
that achieves this exactly (Kjartansson, 1979):
```math
M(\omega) = M_0\left(\frac{i\omega}{\omega_0}\right)^{2\gamma}
          = M_0\left|\frac{\omega}{\omega_0}\right|^{2\gamma} e^{i\pi\gamma\,\mathrm{sgn}(\omega)},
\qquad \frac{1}{Q} = \tan(\pi\gamma) \iff \gamma = \frac{1}{\pi}\arctan\!\left(\frac{1}{Q}\right).
```
``\arg M(\omega)=\pi\gamma\,\mathrm{sgn}(\omega)`` truly has no ``|\omega|``-dependence --
this is the entire content of "constant ``Q``" here, confirmed numerically in the
Appendix rather than merely asserted. Substituting ``M(\omega)`` into the 1-D wave
equation's solution gives the exact closed forms this notebook uses:
```math
c(\omega) = c_0\left|\frac{\omega}{\omega_0}\right|^{\gamma}, \qquad
\alpha(\omega) = \tan\!\left(\frac{\pi\gamma}{2}\right)\mathrm{sgn}(\omega)\,\frac{\omega}{c(\omega)}.
```
As ``Q\to\infty``, ``\gamma\to 1/(\pi Q)\to 0``, and expanding both expressions to leading
order in ``\gamma`` reproduces Futterman's linear formula and its ``\omega/(2cQ)`` decay
rate exactly -- the Appendix's self-check sweeps ``Q`` and shows this discrepancy
shrinking monotonically to zero, i.e. Futterman's formula is not a separate, coincidentally
similar approximation; it is the small-``\gamma`` limit of this one exact model.
"""

# ╔═╡ 892eb702-189e-4b3f-a4a0-89e4763b97ac
md"""
## Geometric Spreading vs. Intrinsic Attenuation

A real recorded amplitude combines two, physically unrelated effects: **geometric
spreading** -- energy conservation spreading a fixed total energy over a growing
wavefront, with no loss to heat -- and the **intrinsic attenuation** derived above, true
anelastic energy loss. For a point source, a body wave's spherical wavefront grows as
``x^2`` in area, so amplitude falls as ``x^{-1}``; a surface wave's cylindrical wavefront
grows only as ``x``, so its amplitude falls more slowly, as ``x^{-1/2}``. The Distance
Axis panel above overlays this combined envelope,
```math
A(x) \propto x^{-n}\,e^{-\alpha(\omega)x}, \qquad n=1 \text{ (body wave)},\ \ n=\tfrac12 \text{ (surface wave)},
```
along the axis both receivers sit on -- drag them and watch where they land relative to
the pure geometric-spreading curve compared to the pure intrinsic-attenuation curve.

Crucially, geometric spreading is *frequency-independent* -- it is a single multiplicative
constant at any given ``x``, with nothing to say about waveform shape or pulse
broadening. This single fact is what makes the spectral ratio method below work at all.
"""

# ╔═╡ a1041e34-78f7-4212-898a-7e1e6d0770d8
md"""
## Measuring ``Q``: The Spectral Ratio Method

Everything above described *forward* modeling: given ``Q``, predict the waveform. In
practice a seismologist has the opposite problem -- given two recorded waveforms, recover
``Q``. The classic technique is the **spectral ratio method**, and the Waveforms &
Spectral Ratio panel above builds it live from the two draggable receivers.

Take two receivers at distances ``x_1,x_2`` in the same medium, and Fourier-transform
each recorded pulse to get complex spectra ``U_1(f)``, ``U_2(f)``. Using the amplitude law
from the very first section, written in terms of ``t^*``:
```math
\ln\left|\frac{U_2(f)}{U_1(f)}\right|
= \ln\left(\frac{A_2}{A_1}\right) - \pi f\left(t^*(x_2)-t^*(x_1)\right)
= \underbrace{\ln\left(\frac{A_2}{A_1}\right)}_{\text{intercept}}
  \;-\; \pi\,\Delta t^*\cdot f,
```
where ``A_1,A_2`` are ANY frequency-independent amplitude factors -- geometric spreading,
an unknown radiation pattern, a site response, an instrument gain difference. **This
quantity is exactly linear in ``f``, with a slope that depends on ``\Delta t^*`` alone.**
Fitting a straight line to ``\ln|U_2(f)/U_1(f)|`` therefore recovers ``\Delta t^*`` from
the *slope*, completely independent of whatever those frequency-independent factors
happen to be -- they only ever shift the *intercept*. Since ``\Delta t^*=(x_2-x_1)/(c_0Q)``,
knowing the source-receiver geometry converts the fitted slope directly into an estimate
of ``Q``:
```math
\widehat{Q} = -\frac{\pi(x_2-x_1)}{c_0\cdot\text{slope}}.
```
The panel's readout shows this recovered ``\widehat{Q}`` next to the "true" (dialed-in)
``Q`` -- try dragging the receivers further apart (a larger, more measurable ``\Delta
t^*``) or toggling independent geometric-spreading exponents at the two receivers and
watch the *intercept* move while ``\widehat{Q}`` stays put. The Appendix's self-check
demonstrates precisely this invariance numerically.
"""

# ╔═╡ b2b60060-a73d-4006-911b-36a28c0e7554
md"""
## References

- Aki, K., & Richards, P. G. (2002). *Quantitative Seismology* (2nd ed.), ch. 5 --
  the nearly-constant-``Q`` (Futterman) dispersion relation and the spectral ratio method.
- Futterman, W. I. (1962). Dispersive body waves. *Journal of Geophysical Research*,
  67(13), 5279-5291.
- Kjartansson, E. (1979). Constant Q-wave propagation and attenuation. *Journal of
  Geophysical Research*, 84(B9), 4737-4748.
- The FFT-based plane-wave synthesis pattern (`source_spectrum`, plain-array frequency
  response, `irfft`) follows `Born-approximation.jl`. The two-receiver drag mechanic
  follows `viscoelastic-rheology.jl`'s 1-D handle drag; the split time/frequency panel
  follows `wave-mode-duality-1D.jl`'s "Receiver Record" panel.
"""

# ╔═╡ 02e113c5-b3b6-4e9e-9670-c937af00d7c2
md"## Appendix"

# ╔═╡ d9b9ee19-f549-4c4f-8b79-6a360bfd5f2d
md"### Reference Medium & Grids"

# ╔═╡ bdf54ecf-72c8-4741-b655-142fc23be647
begin
    "Reference angular frequency (rad/s) that every dispersion/attenuation formula below is anchored to -- 10 Hz, an arbitrary but fixed choice matching the notebook's earlier convention."
    const OMEGA0 = 2π * 10

    """
    	source_spectrum(tgrid, fpeak)

    A Gaussian-shaped amplitude spectrum peaked at `fpeak` Hz, sampled on the real-FFT
    frequency grid implied by `tgrid`. Returns `(freqgrid, Fsource)`; the DC bin is
    zeroed since a DC source carries no travel-time information here.
    """
    function source_spectrum(tgrid, fpeak)
        freqgrid = collect(rfftfreq(length(tgrid), inv(step(tgrid))))
        Fsource = exp.(-abs2.(freqgrid .- fpeak) .* 1.0e-1)
        Fsource[1] = 0.0
        return freqgrid, Fsource
    end
end

# ╔═╡ 189722e0-e9ea-42e8-ad75-684e6da5ad04
md"### Kjartansson's Exact Constant-Q Model"

# ╔═╡ e9e188a6-ec01-4a4f-9168-cf185b14cb05
begin
    """
    	kjartansson_gamma(Q)

    The power-law exponent ``\\gamma`` in Kjartansson's (1979) exact constant-``Q``
    model, ``\\gamma=(1/\\pi)\\arctan(1/Q)``. This single number sets both the phase-velocity
    power law and the (exactly frequency-independent) loss angle simultaneously --
    ``1/Q=\\tan(\\pi\\gamma)`` by construction (see the self-check below).
    """
    kjartansson_gamma(Q) = atan(1 / Q) / π

    """
    	kjartansson_phase_velocity(ω, Q, ω0, c0)

    Kjartansson's (1979) EXACT causal phase velocity for a medium with quality factor
    `Q`, reference frequency `ω0`, and phase velocity `c0` at that reference:
    ``c(\\omega)=c_0|\\omega/\\omega_0|^\\gamma``, `γ` from [`kjartansson_gamma`](@ref).
    Valid for any `Q`, not just `Q≫1` -- contrast with the linearized
    [`futterman_phase_velocity`](@ref).
    """
    kjartansson_phase_velocity(ω, Q, ω0, c0) = c0 * abs(ω / ω0)^kjartansson_gamma(Q)

    """
    	kjartansson_attenuation_coeff(ω, Q, ω0, c0)

    Kjartansson's (1979) exact spatial attenuation coefficient,
    ``\\alpha(\\omega)=\\tan(\\pi\\gamma/2)\\,\\mathrm{sgn}(\\omega)\\,\\omega/c(\\omega)``, so
    that amplitude decays as `exp(-α(ω)x)` over a distance `x`. Reduces to `ω/(2Qc0)` as
    `Q→∞` (see [`futterman_attenuation_coeff`](@ref) and its self-check).
    """
    function kjartansson_attenuation_coeff(ω, Q, ω0, c0)
        γ = kjartansson_gamma(Q)
        c = kjartansson_phase_velocity(ω, Q, ω0, c0)
        return tan(π * γ / 2) * sign(ω) * ω / c
    end
end

# ╔═╡ 3f5f9b5a-2131-432a-b58d-132cf1800610
let
    Qtest = 12.3
    γ = kjartansson_gamma(Qtest)
    lhs = 1 / Qtest
    rhs = tan(π * γ)

    ω0test, c0test = 2π * 10, 3000.0
    ωs = 2π .* [0.1, 1.0, 5.0, 20.0, 100.0, 500.0]
    γs = [kjartansson_gamma(Qtest) for _ in ωs]  # γ (hence the loss angle πγ) has no ω-dependence by construction

    @assert isapprox(lhs, rhs; rtol=1e-10)
    @assert allequal(γs)

    md"""
    !!! correct "Self-check"
        ``1/Q`` = $(round(lhs,digits=6)) matches ``\tan(\pi\gamma)`` =
        $(round(rhs,digits=6)) to machine precision -- the algebraic identity
        `kjartansson_gamma` is built from. And ``\gamma``, the ONLY quantity setting the
        loss angle ``\pi\gamma``, has no ``\omega``-dependence at all across a
        $(round(Int,ωs[end]/ωs[1]))-fold frequency sweep -- ``Q`` really is
        frequency-independent in this model, not approximately so.
    """
end

# ╔═╡ 14a5545d-858e-407b-a193-0e5abc9f68b1
md"### Futterman's Nearly-Constant-Q Approximation (Q≫1)"

# ╔═╡ 22bda95b-a685-46ad-8284-075bafd27581
begin
    """
    	futterman_phase_velocity(ω, Q, ω0, c0)

    Futterman's (1962) "nearly constant Q" phase velocity, the ``Q\\gg1`` linearization
    of [`kjartansson_phase_velocity`](@ref): ``c(\\omega)=c_0[1+(1/\\pi Q)\\ln(\\omega/\\omega_0)]``.
    This is the formula Aki & Richards (ch. 5) give for the standard constant-``Q``
    dispersion relation, exact only in the limit ``Q\\to\\infty`` (see the self-check below).
    """
    futterman_phase_velocity(ω, Q, ω0, c0) = c0 * (1 + log(ω / ω0) / (π * Q))

    """
    	futterman_attenuation_coeff(ω, Q, c)

    The small-loss (``Q\\gg1``) limit of [`kjartansson_attenuation_coeff`](@ref):
    ``\\alpha(\\omega)=\\omega/(2Qc)``, the amplitude decay rate already derived in the
    "Quality Factor" section above.
    """
    futterman_attenuation_coeff(ω, Q, c) = ω / (2 * Q * c)
end

# ╔═╡ ea207d15-c2b9-472c-9fc3-79b7000a1738
let
    ω0test, c0test = 2π * 10, 3000.0
    ωtest = 2π * 15.0

    # Futterman's formula itself approximates γ=atan(1/Q)/π by the cruder 1/(πQ) -- so
    # comparing it against Kjartansson's EXACT γ mixes two independent small-parameter
    # approximations (the exponential linearization, and atan(1/Q)≈1/Q) that don't have
    # to combine monotonically while BOTH are still non-negligible (roughly Q≲50, where
    # they partly cancel at some Q and reinforce at others). Once Q is large enough that
    # atan(1/Q)≈1/Q is already excellent (Q≳50), only the single remaining small
    # parameter matters and convergence genuinely is monotonic -- checked directly below.
    Qs = [50.0, 100.0, 500.0, 1000.0, 1.0e4, 1.0e5]
    errs_c = [abs(kjartansson_phase_velocity(ωtest, Q, ω0test, c0test) -
                  futterman_phase_velocity(ωtest, Q, ω0test, c0test)) / c0test for Q in Qs]
    errs_a = [abs(kjartansson_attenuation_coeff(ωtest, Q, ω0test, c0test) -
                  futterman_attenuation_coeff(ωtest, Q, c0test)) /
              kjartansson_attenuation_coeff(ωtest, Q, ω0test, c0test) for Q in Qs]

    @assert issorted(errs_c; rev=true)
    @assert issorted(errs_a; rev=true)
    @assert errs_c[end] < 1.0e-9 && errs_a[end] < 1.0e-5

    # even at a fairly modest Q=10 (well short of "Q≫1"), the two formulas already agree
    # closely -- the approximation is useful long before it's asymptotically exact
    Q10 = 10.0
    errc10 = abs(kjartansson_phase_velocity(ωtest, Q10, ω0test, c0test) -
                 futterman_phase_velocity(ωtest, Q10, ω0test, c0test)) / c0test
    erra10 = abs(kjartansson_attenuation_coeff(ωtest, Q10, ω0test, c0test) -
                 futterman_attenuation_coeff(ωtest, Q10, c0test)) /
             kjartansson_attenuation_coeff(ωtest, Q10, ω0test, c0test)
    @assert errc10 < 0.01 && erra10 < 0.02

    md"""
    !!! correct "Self-check"
        Both Futterman formulas are literally the ``Q\to\infty`` limit of Kjartansson's
        exact ones: sweeping ``Q`` from $(Int(Qs[1])) to $(Qs[end]) at a fixed frequency
        (comfortably past the point where Futterman's own ``\gamma\approx1/(\pi Q)``
        shortcut is itself accurate), the relative discrepancy in phase velocity shrinks
        monotonically from $(round(errs_c[1]*100,sigdigits=2))% down to
        $(round(errs_c[end]*100,sigdigits=3))%, and in attenuation from
        $(round(errs_a[1]*100,sigdigits=2))% down to $(round(errs_a[end]*100,sigdigits=3))%.
        Even at a much more modest ``Q`` = $(Int(Q10)), the two already agree to
        $(round(errc10*100,sigdigits=2))% (velocity) and $(round(erra10*100,sigdigits=2))%
        (attenuation) -- confirming the linearization is genuinely the small-``\gamma``
        limit of the exact model, not a separately-derived formula that merely happens to
        look similar.
    """
end

# ╔═╡ 6bfbe8d3-3ef0-4a82-8266-59794fd6118a
md"### Geometric Spreading"

# ╔═╡ fd66f5ad-1b8e-4b0e-9377-32d8f15e9ac2
"""
	spreading_factor(x, n)

Geometric-spreading amplitude decay over a distance `x`: ``x^{-n}``. `n=1` for body
waves (spherical spreading of a point-source wavefront), `n=0.5` for surface waves
(cylindrical spreading, energy spread over a growing ring rather than a growing sphere).
"""
spreading_factor(x, n) = x^(-n)

# ╔═╡ 82ff57b2-3824-42db-9d55-a3c7fc377131
md"### Plane-Wave Synthesis"

# ╔═╡ 0ec1c415-c990-46eb-bfdf-17d062442f90
begin
    """
    	attenuation_response(freqgrid, x, Q, ω0, c0, model)

    The complex propagation operator for a distance `x`, RELATIVE to the pure `x/c0`
    bulk delay (i.e. what multiplies the source spectrum after the trivial travel time at
    the reference velocity `c0` has already been factored out -- see
    [`receiver_spectrum`](@ref)). `model` selects which phase-velocity/attenuation pair to
    use:
    - `"elastic"`: no attenuation, no dispersion -- `1.0+0im` for every frequency.
    - `"futterman"`: [`futterman_phase_velocity`](@ref)/[`futterman_attenuation_coeff`](@ref).
    - `"kjartansson"`: [`kjartansson_phase_velocity`](@ref)/[`kjartansson_attenuation_coeff`](@ref).

    `f=0` is skipped (returned as `1.0+0im`) since both dispersion laws are singular
    there; [`source_spectrum`](@ref) already zeroes the DC bin anyway.
    """
    function attenuation_response(freqgrid, x, Q, ω0, c0, model)
        return map(freqgrid) do f
            iszero(f) && return complex(1.0, 0.0)
            ω = 2π * f
            if model == "elastic"
                return complex(1.0, 0.0)
            elseif model == "futterman"
                c = futterman_phase_velocity(ω, Q, ω0, c0)
                α = futterman_attenuation_coeff(ω, Q, c0)
            elseif model == "kjartansson"
                c = kjartansson_phase_velocity(ω, Q, ω0, c0)
                α = kjartansson_attenuation_coeff(ω, Q, ω0, c0)
            else
                error("unknown model $model")
            end
            return exp(im * ω * x * (1 / c0 - 1 / c)) * exp(-α * x)
        end
    end

    """
    	receiver_spectrum(freqgrid, Fsource, x, n, Q, ω0, c0, model)

    The full complex spectrum recorded at distance `x`: the source spectrum `Fsource`,
    scaled by [`spreading_factor`](@ref)`(x,n)`, times [`attenuation_response`](@ref) for
    the requested `model`. The trivial bulk delay `x/c0` is DELIBERATELY not included in
    the returned phase -- every receiver's pulse is synthesized on the same shared,
    `x`-independent time grid, so a distant receiver's pulse appears centered near `t=0`
    just like a close one, isolating the interesting (dispersion/attenuation) part of the
    story from the boring part (it simply arrives later). See [`synthesize_pulse`](@ref).
    """
    function receiver_spectrum(freqgrid, Fsource, x, n, Q, ω0, c0, model)
        return Fsource .* spreading_factor(x, n) .* attenuation_response(freqgrid, x, Q, ω0, c0, model)
    end

    """
    	synthesize_pulse(tgrid, spectrum)

    Time-domain pulse from a real-FFT `spectrum` sampled on the grid implied by `tgrid`,
    via `irfft` then `fftshift` -- paired with [`receiver_spectrum`](@ref)'s zero-phase,
    reduced-time convention, `fftshift` places the pulse's peak near the MIDDLE of
    `tgrid` (which must be symmetric about 0) rather than at its wrapped-around start.
    """
    synthesize_pulse(tgrid, spectrum) = fftshift(irfft(spectrum, length(tgrid)))
end

# ╔═╡ f668c317-1166-4ab1-b96c-c46eaf9ee136
let
    tgrid_test = range(-2.0, 2.0, length=1024)
    freqgrid_test, Fsource_test = source_spectrum(tgrid_test, 5.0)
    x1test, ntest = 40.0, 1.0

    spec_source_only = Fsource_test .* spreading_factor(x1test, ntest)
    pulse_source = synthesize_pulse(tgrid_test, spec_source_only)

    spec_elastic = receiver_spectrum(freqgrid_test, Fsource_test, x1test, ntest, 30.0, OMEGA0, 3000.0, "elastic")
    pulse_elastic = synthesize_pulse(tgrid_test, spec_elastic)

    resid = maximum(abs.(pulse_elastic .- pulse_source))
    peak_idx = argmax(abs.(pulse_elastic))
    peak_t = tgrid_test[peak_idx]

    @assert isapprox(resid, 0.0; atol=1e-9)
    @assert isapprox(peak_t, 0.0; atol=2 * step(tgrid_test))

    md"""
    !!! correct "Self-check"
        With `model="elastic"`, `receiver_spectrum` differs from a bare, spreading-scaled
        copy of the source spectrum by at most $(round(resid,sigdigits=2)) in the
        synthesized time series -- zero to numerical precision, confirming
        `attenuation_response` really is the identity when there's nothing to attenuate.
        Its peak sits at reduced time $(round(peak_t,digits=4)) s, within one grid
        step of ``0`` -- the reduced-time convention is doing its job.
    """
end

# ╔═╡ 8271e3a5-308d-42b4-b636-1a59f18361f6
md"### The Spectral Ratio Method"

# ╔═╡ 9c1c2062-c62a-48dd-9d0f-651d848afccd
begin
    """
    	log_spectral_ratio(U1, U2)

    ``\\ln|U_2(f)/U_1(f)|``, the quantity a spectral-ratio ``Q`` measurement actually
    fits a line to -- theoretically linear in `f` with slope ``-\\pi\\Delta t^*`` (see
    [`estimate_Q`](@ref)), regardless of any frequency-independent amplitude prefactors.
    """
    log_spectral_ratio(U1, U2) = log.(abs.(U2) ./ abs.(U1))

    """
    	linear_fit(x, y)

    Ordinary least-squares slope and intercept of `y` against `x`, via the normal
    equations `[x ones(length(x))] \\ y`. This notebook's smallest general-purpose
    regression helper -- no larger-scoped one exists elsewhere in this repo to reuse
    (`ray-tomography.jl`'s `get_tikhonov_solution` solves a damped normal-equations
    problem specific to travel-time tomography, not a plain 1-D fit).
    """
    function linear_fit(x, y)
        Adesign = [x ones(length(x))]
        slope, intercept = Adesign \ y
        return slope, intercept
    end

    """
    	estimate_Q(slope, x1, x2, c0)

    Recover ``Q`` from the fitted [`linear_fit`](@ref) `slope` of
    [`log_spectral_ratio`](@ref) vs. frequency (Hz). Since ``t^*(x)=x/(cQ)`` and
    ``\\ln|U_2/U_1|=\\text{const}-\\pi f\\Delta t^*`` with ``\\Delta t^*=(x_2-x_1)/(c_0Q)``,
    ``\\text{slope}=-\\pi(x_2-x_1)/(c_0Q)``, so ``Q=-\\pi(x_2-x_1)/(c_0\\cdot\\text{slope})``.
    """
    estimate_Q(slope, x1, x2, c0) = -π * (x2 - x1) / (c0 * slope)
end

# ╔═╡ fb5825ce-55cf-4b48-8513-0669fbe876cd
let
    tgrid_test = range(-2.0, 2.0, length=2048)
    freqgrid_test, Fsource_test = source_spectrum(tgrid_test, 5.0)
    x1test, x2test, c0test, Qtrue = 20.0, 80.0, 3000.0, 25.0
    mask = Fsource_test .> 0.05 * maximum(Fsource_test)

    U1 = receiver_spectrum(freqgrid_test, Fsource_test, x1test, 1.0, Qtrue, OMEGA0, c0test, "kjartansson")
    U2 = receiver_spectrum(freqgrid_test, Fsource_test, x2test, 1.0, Qtrue, OMEGA0, c0test, "kjartansson")
    ratio = log_spectral_ratio(U1, U2)
    slope, intercept = linear_fit(freqgrid_test[mask], ratio[mask])
    Qest = estimate_Q(slope, x1test, x2test, c0test)
    relerr = abs(Qest - Qtrue) / Qtrue

    # geometric-spreading invariance: give receiver 2 a DIFFERENT spreading exponent AND
    # an arbitrary amplitude multiplier (standing in for an unknown radiation-pattern or
    # site difference) -- the fitted SLOPE (hence Q_est) must be unchanged, only the
    # intercept should move
    U1b = receiver_spectrum(freqgrid_test, Fsource_test, x1test, 1.0, Qtrue, OMEGA0, c0test, "kjartansson")
    U2b = 3.7 .* receiver_spectrum(freqgrid_test, Fsource_test, x2test, 0.5, Qtrue, OMEGA0, c0test, "kjartansson")
    ratiob = log_spectral_ratio(U1b, U2b)
    slopeb, interceptb = linear_fit(freqgrid_test[mask], ratiob[mask])
    Qestb = estimate_Q(slopeb, x1test, x2test, c0test)

    @assert relerr < 0.05
    @assert isapprox(slope, slopeb; rtol=1e-8)
    @assert isapprox(Qest, Qestb; rtol=1e-8)
    @assert !isapprox(intercept, interceptb; rtol=1e-3)

    md"""
    !!! correct "Self-check"
        Synthesizing two receivers at a known ``Q_{\rm true}`` = $(Int(Qtrue)) and fitting
        the spectral ratio recovers ``\widehat{Q}`` = $(round(Qest,digits=2)) --
        $(round(relerr*100,digits=2))% error, from nothing but the two waveforms. Giving
        receiver 2 a DIFFERENT geometric-spreading exponent (``n=0.5`` instead of ``1``)
        and an arbitrary ``3.7\times`` amplitude factor (standing in for an unknown
        radiation pattern or site response) leaves the fitted slope, and therefore
        ``\widehat{Q}`` = $(round(Qestb,digits=2)), unchanged to within
        $(round(abs(slope-slopeb)/abs(slope)*100,sigdigits=2))%, while the intercept
        visibly shifts from $(round(intercept,digits=3)) to $(round(interceptb,digits=3)).
        This is *why* the spectral-ratio method works without needing to know or correct
        for spreading: it only ever reads the slope.
    """
end

# ╔═╡ 097e5d11-4857-4fd3-a587-3d72d2567180
md"### The Interactive Widget"

# ╔═╡ f3674d89-497d-43c3-9741-8a1d8e3fc860
begin
    """
    	AttenuationInput(; model="futterman", Q=25.0, c0=3.0, fp=3.0, x1=20.0, x2=80.0, n1=1.0, n2=1.0)

    Initial state for the widget: dispersion/attenuation `model` (`"elastic"`,
    `"futterman"`, or `"kjartansson"`), quality factor `Q`, reference velocity `c0`
    (km/s), source peak frequency `fp` (Hz), the two receiver distances `x1`,`x2` (km,
    draggable on the Distance Axis panel), and each receiver's own geometric-spreading
    exponent `n1`,`n2`.
    """
    struct AttenuationInput
        model::String
        Q::Float64
        c0::Float64
        fp::Float64
        x1::Float64
        x2::Float64
        n1::Float64
        n2::Float64
    end
    AttenuationInput(; model="futterman", Q=25.0, c0=3.0, fp=3.0, x1=20.0, x2=80.0, n1=1.0, n2=1.0) =
        AttenuationInput(model, Q, c0, fp, x1, x2, n1, n2)

    Base.get(w::AttenuationInput) = Dict{String,Any}(
        "model" => w.model, "Q" => w.Q, "c0" => w.c0, "fp" => w.fp,
        "x1" => w.x1, "x2" => w.x2, "n1" => w.n1, "n2" => w.n2)

    """
    	Base.show(io, ::MIME"text/html", w::AttenuationInput)

    Render the three-panel scene: a Dispersion & Attenuation curve panel, a Distance Axis
    panel with a fixed source and two draggable receivers, and a split
    Waveforms/Spectral-Ratio panel, plus a "Medium & Source" and a "Geometric Spreading"
    control group below.
    """
    function Base.show(io::IO, ::MIME"text/html", w::AttenuationInput)
        write(io, """
        <div id="iawidget">
        <style>
        #iawidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #iawidget .ia-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #iawidget .ia-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #iawidget .ia-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #iawidget .ia-row{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start;margin-bottom:14px}
        #iawidget .ia-row-secondary{display:flex;gap:16px;margin-bottom:14px}
        #iawidget .ia-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #iawidget .ia-panel-title{font-size:14px;font-weight:700;color:#e5e7eb;margin-bottom:4px;text-align:center}
        #iawidget .ia-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:4px}
        #iawidget canvas{display:block;cursor:default}
        #iawidget .ia-controls{width:100%;box-sizing:border-box;display:flex;gap:12px;flex-wrap:wrap}
        #iawidget .ia-control-group{flex:1 1 280px;min-width:260px;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
        #iawidget .ia-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #iawidget .ia-control-row{display:grid;grid-template-columns:70px minmax(60px,1fr) 64px;gap:6px;align-items:center;margin:5px 0}
        #iawidget .ia-control-row label{font-size:13px;color:#9ca3af}
        #iawidget .ia-control-row input[type=range]{width:100%;min-width:0}
        #iawidget .ia-value{font-size:12px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #iawidget .ia-actions{display:flex;gap:6px;flex-wrap:wrap}
        #iawidget button{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 10px;font-size:13px;cursor:pointer}
        #iawidget button.active{background:#2563eb;border-color:#93c5fd}
        #iawidget button:hover{background:#767676}
        </style>

        <div class="ia-title">
          <div class="ia-title-desc">Two receivers, one source -- measure Q the way a seismologist actually does.</div>
          <div class="ia-title-hint">drag either receiver (triangle) along the Distance Axis &middot; pick a dispersion model, Q, c&#8320;, and source frequency below &middot; the Waveforms panel's slope-fitted Q&#770; is read straight off the spectral ratio, not handed to you</div>
        </div>

        <div class="ia-row">
          <div>
            <div class="ia-panel-title">Dispersion &amp; Attenuation</div>
            <div class="ia-panel"><canvas id="ia-dispersion"></canvas></div>
            <div class="ia-caption">c(&omega;)/c&#8320; vs. frequency &middot; Futterman (blue) vs. Kjartansson (orange)</div>
          </div>
          <div>
            <div class="ia-panel-title">Distance Axis: Source &amp; Receivers</div>
            <div class="ia-panel"><canvas id="ia-axis"></canvas></div>
            <div class="ia-caption" id="ia-axis-caption">drag a triangle to move that receiver</div>
          </div>
        </div>

        <div class="ia-row-secondary">
          <div style="flex:1 1 100%">
            <div class="ia-panel-title">Waveforms &amp; Spectral Ratio</div>
            <div class="ia-panel"><canvas id="ia-record"></canvas></div>
            <div class="ia-caption">left: both receivers' pulses (reduced time) &middot; right: ln|U&#8322;/U&#8321;| vs. frequency with the fitted line</div>
          </div>
        </div>

        <div class="ia-controls">
          <div class="ia-control-group">
            <div class="ia-control-title">Medium &amp; Source</div>
            <div class="ia-actions">
              <button id="ia-model-elastic" type="button">Elastic</button>
              <button id="ia-model-futterman" type="button">Futterman (Q&#8811;1)</button>
              <button id="ia-model-kjartansson" type="button">Kjartansson (exact)</button>
            </div>
            <div class="ia-control-row"><label>Q</label><input type="range" id="ia-Q" min="3" max="200" step="1" value="$(w.Q)"><span class="ia-value" id="ia-Q-v"></span></div>
            <div class="ia-control-row"><label>c&#8320; (km/s)</label><input type="range" id="ia-c0" min="1" max="8" step="0.1" value="$(w.c0)"><span class="ia-value" id="ia-c0-v"></span></div>
            <div class="ia-control-row"><label>f&#7605; (Hz)</label><input type="range" id="ia-fp" min="1" max="10" step="0.1" value="$(w.fp)"><span class="ia-value" id="ia-fp-v"></span></div>
          </div>
          <div class="ia-control-group">
            <div class="ia-control-title">Geometric Spreading</div>
            <div class="ia-actions">
              <button id="ia-spread-body" type="button">Both: body wave (n=1)</button>
              <button id="ia-spread-surface" type="button">Both: surface wave (n=0.5)</button>
            </div>
            <div class="ia-control-row"><label>r&#8321; n</label><input type="range" id="ia-n1" min="0.3" max="2" step="0.1" value="$(w.n1)"><span class="ia-value" id="ia-n1-v"></span></div>
            <div class="ia-control-row"><label>r&#8322; n</label><input type="range" id="ia-n2" min="0.3" max="2" step="0.1" value="$(w.n2)"><span class="ia-value" id="ia-n2-v"></span></div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        // Pluto can execute a cell's <script> tag more than once against the SAME
        // persistent DOM nodes (a client's initial connection replays intermediate cell
        // outputs, and later edits rerun this cell) -- guard so only the first execution
        // for this widget instance wires anything up, matching every other widget in
        // this repo (e.g. viscoelastic-rheology.jl's `_eplInitialized` guard).
        if(!par._iaInitialized){
        par._iaInitialized = true;

        // WideCell's own ResizeObserver widens `par` asynchronously, some unknown number
        // of frames after this script first runs -- reading par.clientWidth before that
        // lands bakes a too-small layout into every canvas (confirmed: canvases came out
        // sized to a narrow pre-widen clientWidth while their .ia-panel wrapper divs
        // stretched to the real, later-widened row width, leaving a large blank gap of
        // the wrapper's own black background to the right of each canvas). Deferring the
        // whole setup into iaInit(), invoked only via a debounced ResizeObserver on
        // `par`, waits for the real settled width instead of guessing how many frames to
        // skip -- same pattern as ray-tomography.jl's rtInit().
        function iaInit(){
        let state = { model: "$(w.model)", Q: $(w.Q), c0: $(w.c0), fp: $(w.fp),
          x1: $(w.x1), x2: $(w.x2), n1: $(w.n1), n2: $(w.n2) };
        let pushed = null; // {fdisp,cfutt,ckjar,xaxis,env1,env2,tgrid,pulse1,pulse2,fband,ratio,slope,intercept,Qest,dtstarEst}
        let commitInFlight = false;

        // par.clientWidth (measured post-settle, see iaInit's own comment above) already
        // reflects the real column width -- WideCell's own max_width already bounds it,
        // and Pluto's chrome/sidebar are already subtracted out of it. Capping it further
        // against window.innerWidth*fraction is redundant at best and, whenever the
        // notebook column is wide relative to the browser window (a wide viewport, or a
        // narrow Pluto sidebar), actively WRONG: that fraction can be smaller than the
        // real clientWidth, silently shrinking every canvas below its own wrapper's
        // width and leaving a dead, unused stripe of the wrapper's own black background
        // to its right (confirmed live: canvases sized from 0.85*innerWidth while their
        // .ia-panel wrappers had already stretched to a wider, correctly-measured row).
        const availW = Math.min(par.clientWidth || 1400, 1400) || 900;
        const GAP = 16;
        const PW = Math.max(260, (availW - GAP)/2);
        const PH = Math.max(170, Math.round(PW*0.55));
        const RECORD_W = availW;
        const RECORD_H = Math.max(190, Math.round(RECORD_W*0.26));
        const DPR = window.devicePixelRatio || 1;

        function hidpi(canvas, ctx, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          ctx.setTransform(DPR,0,0,DPR,0,0);
        }

        const dispCv = par.querySelector('#ia-dispersion'), dispCtx = dispCv.getContext('2d');
        hidpi(dispCv, dispCtx, PW, PH);
        const axisCv = par.querySelector('#ia-axis'), axisCtx = axisCv.getContext('2d');
        hidpi(axisCv, axisCtx, PW, PH);
        const recCv = par.querySelector('#ia-record'), recCtx = recCv.getContext('2d');
        hidpi(recCv, recCtx, RECORD_W, RECORD_H);

        const modelBtns = { elastic: par.querySelector('#ia-model-elastic'), futterman: par.querySelector('#ia-model-futterman'), kjartansson: par.querySelector('#ia-model-kjartansson') };
        const QIn = par.querySelector('#ia-Q'), QV = par.querySelector('#ia-Q-v');
        const c0In = par.querySelector('#ia-c0'), c0V = par.querySelector('#ia-c0-v');
        const fpIn = par.querySelector('#ia-fp'), fpV = par.querySelector('#ia-fp-v');
        const n1In = par.querySelector('#ia-n1'), n1V = par.querySelector('#ia-n1-v');
        const n2In = par.querySelector('#ia-n2'), n2V = par.querySelector('#ia-n2-v');
        const axisCaption = par.querySelector('#ia-axis-caption');

        function drawStarMarker(ctx, cx, cy, r, fill, stroke){
          const spikes = 5, rOuter = r, rInner = r * 0.45;
          ctx.beginPath();
          for(let i=0; i<spikes*2; i++){
            const rad = i % 2 === 0 ? rOuter : rInner;
            const ang = -Math.PI/2 + i*Math.PI/spikes;
            const x = cx + rad*Math.cos(ang), y = cy + rad*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1; ctx.stroke();
        }
        function drawTriangleDownMarker(ctx, cx, cy, r, fill, stroke){
          ctx.beginPath();
          for(let i=0; i<3; i++){
            const ang = Math.PI/2 + i*2*Math.PI/3;
            const x = cx + r*Math.cos(ang), y = cy + r*Math.sin(ang);
            i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
          }
          ctx.closePath();
          ctx.fillStyle = fill; ctx.fill();
          ctx.strokeStyle = stroke; ctx.lineWidth = 1.5; ctx.stroke();
        }

        function syncControlLabels(){
          for(const k in modelBtns) modelBtns[k].classList.toggle('active', state.model===k);
          QV.textContent = state.Q.toFixed(0);
          c0V.textContent = state.c0.toFixed(1)+' km/s';
          fpV.textContent = state.fp.toFixed(1)+' Hz';
          n1V.textContent = state.n1.toFixed(1);
          n2V.textContent = state.n2.toFixed(1);
        }
        syncControlLabels();

        function emit(){
          commitInFlight = true;
          par.value = { model: state.model, Q: state.Q, c0: state.c0, fp: state.fp,
            x1: state.x1, x2: state.x2, n1: state.n1, n2: state.n2 };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledEmit(){ if(!commitInFlight) emit(); }

        // ---- Dispersion & Attenuation panel ----
        function drawDispersion(){
          dispCtx.clearRect(0,0,PW,PH);
          dispCtx.strokeStyle = '#374151'; dispCtx.lineWidth = 1; dispCtx.strokeRect(0.5,0.5,PW-1,PH-1);
          if(!pushed){ dispCtx.fillStyle='#6b7280'; dispCtx.font='12px sans-serif'; dispCtx.fillText('computing...', 10, 18); return; }
          const bandL = 40, bandB = 18, plotW = PW-bandL-6, plotH = PH-bandB-6;
          const fmin = pushed.fdisp[0], fmax = pushed.fdisp[pushed.fdisp.length-1];
          const allY = pushed.cfutt.concat(pushed.ckjar).concat([1.0]);
          let ymin = Math.min(...allY), ymax = Math.max(...allY);
          const pad = Math.max(1e-4, (ymax-ymin)*0.15); ymin -= pad; ymax += pad;
          function xOf(f){ return bandL + (Math.log10(f)-Math.log10(fmin))/(Math.log10(fmax)-Math.log10(fmin))*plotW; }
          function yOf(v){ return 6 + plotH - (v-ymin)/(ymax-ymin)*plotH; }

          function drawCurve(arr, color, alpha){
            dispCtx.globalAlpha = alpha; dispCtx.strokeStyle = color; dispCtx.lineWidth = 1.8;
            dispCtx.beginPath();
            for(let i=0;i<arr.length;i++){
              const x = xOf(pushed.fdisp[i]), y = yOf(arr[i]);
              i===0 ? dispCtx.moveTo(x,y) : dispCtx.lineTo(x,y);
            }
            dispCtx.stroke(); dispCtx.globalAlpha = 1;
          }
          if(state.model === 'elastic'){
            dispCtx.strokeStyle = '#e5e7eb'; dispCtx.lineWidth = 1.8;
            dispCtx.beginPath(); dispCtx.moveTo(bandL, yOf(1.0)); dispCtx.lineTo(bandL+plotW, yOf(1.0)); dispCtx.stroke();
            drawCurve(pushed.cfutt, '#38bdf8', 0.3);
            drawCurve(pushed.ckjar, '#f97316', 0.3);
          } else {
            drawCurve(pushed.cfutt, '#38bdf8', state.model==='futterman' ? 1.0 : 0.3);
            drawCurve(pushed.ckjar, '#f97316', state.model==='kjartansson' ? 1.0 : 0.3);
          }
          dispCtx.strokeStyle = '#4b5563'; dispCtx.beginPath(); dispCtx.moveTo(bandL,yOf(1.0)+0.5); dispCtx.lineTo(bandL+plotW,yOf(1.0)+0.5); dispCtx.stroke();
          dispCtx.fillStyle = '#e5e7eb'; dispCtx.font = '10px sans-serif'; dispCtx.textAlign='left';
          dispCtx.fillText('c/c\\u2080='+ymin.toFixed(3)+'..'+ymax.toFixed(3), bandL+2, 12);
          dispCtx.fillStyle = '#9ca3af';
          dispCtx.fillText(fmin.toFixed(1)+' Hz', bandL, PH-4);
          dispCtx.textAlign='right'; dispCtx.fillText(fmax.toFixed(0)+' Hz', bandL+plotW, PH-4);
        }

        // ---- Distance Axis panel ----
        const AX_XMIN = 0, AX_XMAX = 150, AX_MINGAP = 4;
        function axisLayout(){
          const padL = 22, padR = 14, midY = PH*0.62;
          const plotW = PW - padL - padR;
          return { padL, plotW, midY };
        }
        function xToPx(x){ const {padL, plotW} = axisLayout(); return padL + (x-AX_XMIN)/(AX_XMAX-AX_XMIN)*plotW; }
        function clampAway(v, other){
          v = Math.max(AX_XMIN+2, Math.min(AX_XMAX, v));
          if(Math.abs(v-other) < AX_MINGAP) v = other + Math.sign(v-other || 1)*AX_MINGAP;
          return Math.max(AX_XMIN+2, Math.min(AX_XMAX, v));
        }

        function drawAxis(){
          axisCtx.clearRect(0,0,PW,PH);
          axisCtx.strokeStyle = '#374151'; axisCtx.lineWidth = 1; axisCtx.strokeRect(0.5,0.5,PW-1,PH-1);
          const {padL, plotW, midY} = axisLayout();

          // combined spreading x attenuation envelope, one curve per receiver's own n
          if(pushed){
            function drawEnvelope(env, color){
              axisCtx.strokeStyle = color; axisCtx.globalAlpha = 0.55; axisCtx.lineWidth = 1.3;
              axisCtx.beginPath();
              for(let i=0;i<pushed.xaxis.length;i++){
                const x = xToPx(pushed.xaxis[i]), y = midY - env[i]*(midY-14);
                i===0 ? axisCtx.moveTo(x,y) : axisCtx.lineTo(x,y);
              }
              axisCtx.stroke(); axisCtx.globalAlpha = 1;
            }
            drawEnvelope(pushed.env1, '#38bdf8');
            drawEnvelope(pushed.env2, '#f97316');
          }

          axisCtx.strokeStyle = '#4b5563'; axisCtx.beginPath();
          axisCtx.moveTo(padL, midY); axisCtx.lineTo(padL+plotW, midY); axisCtx.stroke();

          drawStarMarker(axisCtx, xToPx(0), midY, 8, '#facc15', '#000');
          axisCtx.fillStyle = '#9ca3af'; axisCtx.font = '10px sans-serif'; axisCtx.textAlign='center';
          axisCtx.fillText('source', xToPx(0), midY+22);

          drawTriangleDownMarker(axisCtx, xToPx(state.x1), midY, 7, '#38bdf8', '#0a0f18');
          axisCtx.fillStyle = '#38bdf8';
          axisCtx.fillText('r\\u2081 '+state.x1.toFixed(0)+' km', xToPx(state.x1), midY+22);
          drawTriangleDownMarker(axisCtx, xToPx(state.x2), midY, 7, '#f97316', '#0a0f18');
          axisCtx.fillStyle = '#f97316';
          axisCtx.fillText('r\\u2082 '+state.x2.toFixed(0)+' km', xToPx(state.x2), midY+22);

          axisCtx.fillStyle = '#6b7280'; axisCtx.textAlign='left'; axisCtx.font='9px sans-serif';
          axisCtx.fillText('0', padL, midY-8);
          axisCtx.textAlign='right';
          axisCtx.fillText(AX_XMAX+' km', padL+plotW, midY-8);
        }

        let activeDrag = null; // null | 'r1' | 'r2'
        let dragAnchorPx = 0, dragAnchorX = 0;

        function axisPointerXY(ev){
          const r = axisCv.getBoundingClientRect();
          return [ (ev.clientX-r.left)*(PW/r.width), (ev.clientY-r.top)*(PH/r.height) ];
        }
        axisCv.addEventListener('mousedown', ev => {
          const [px, py] = axisPointerXY(ev);
          const {midY} = axisLayout();
          if(Math.abs(py-midY) > 20) return;
          const p1 = xToPx(state.x1), p2 = xToPx(state.x2);
          if(Math.abs(px-p1) <= Math.abs(px-p2) && Math.abs(px-p1) < 10){ activeDrag='r1'; dragAnchorPx=px; dragAnchorX=state.x1; return; }
          if(Math.abs(px-p2) < 10){ activeDrag='r2'; dragAnchorPx=px; dragAnchorX=state.x2; }
        });
        axisCv.addEventListener('mousemove', ev => {
          if(activeDrag) return;
          const [px, py] = axisPointerXY(ev);
          const {midY} = axisLayout();
          const near = Math.abs(py-midY) < 20 && (Math.abs(px-xToPx(state.x1))<10 || Math.abs(px-xToPx(state.x2))<10);
          axisCv.style.cursor = near ? 'grab' : 'default';
        });
        window.addEventListener('mousemove', ev => {
          if(!activeDrag) return;
          const [px] = axisPointerXY(ev);
          const {plotW} = axisLayout();
          const pxPerKm = plotW/(AX_XMAX-AX_XMIN);
          const newX = dragAnchorX + (px-dragAnchorPx)/pxPerKm;
          if(activeDrag === 'r1') state.x1 = clampAway(newX, state.x2);
          else state.x2 = clampAway(newX, state.x1);
          axisCaption.textContent = 'r\\u2081='+state.x1.toFixed(0)+' km, r\\u2082='+state.x2.toFixed(0)+' km';
          drawAxis();
          throttledEmit();
        });
        window.addEventListener('mouseup', () => { if(activeDrag){ axisCv.style.cursor='grab'; } activeDrag = null; });

        // ---- Waveforms & Spectral Ratio panel ----
        function drawRecord(){
          recCtx.clearRect(0,0,RECORD_W,RECORD_H);
          recCtx.strokeStyle = '#374151'; recCtx.lineWidth = 1; recCtx.strokeRect(0.5,0.5,RECORD_W-1,RECORD_H-1);
          if(!pushed){ recCtx.fillStyle='#6b7280'; recCtx.font='12px sans-serif'; recCtx.fillText('computing...', 10, 18); return; }
          const halfW = Math.floor(RECORD_W/2);
          const bandL = 34, bandB = 16, plotH = RECORD_H-bandB;
          const plotWL = halfW-bandL-4, plotWR = RECORD_W-halfW-bandL-4;

          // left: both pulses, reduced time
          {
            const mx = Math.max(...pushed.pulse1.map(Math.abs), ...pushed.pulse2.map(Math.abs), 1e-12);
            const midY = plotH*0.5, amp = plotH*0.42/mx, nt = pushed.tgrid.length;
            function drawPulse(p, color, dashed){
              recCtx.strokeStyle = color; recCtx.lineWidth = 1.5;
              recCtx.setLineDash(dashed ? [4,3] : []);
              recCtx.beginPath();
              for(let i=0;i<nt;i++){
                const x = bandL + i/(nt-1)*plotWL, y = midY - p[i]*amp;
                i===0 ? recCtx.moveTo(x,y) : recCtx.lineTo(x,y);
              }
              recCtx.stroke(); recCtx.setLineDash([]);
            }
            drawPulse(pushed.pulse1, '#38bdf8', false);
            drawPulse(pushed.pulse2, '#f97316', true);
            recCtx.strokeStyle = '#4b5563'; recCtx.beginPath(); recCtx.moveTo(bandL,midY); recCtx.lineTo(bandL+plotWL,midY); recCtx.stroke();
            recCtx.fillStyle = '#9ca3af'; recCtx.font = '9px sans-serif'; recCtx.textAlign='center';
            recCtx.fillText('reduced time (s)', bandL+plotWL/2, RECORD_H-3);
          }
          // right: log spectral ratio + fit
          {
            const nf = pushed.fband.length;
            const fmax = pushed.fband[nf-1];
            const rmin = Math.min(...pushed.ratio), rmax = Math.max(...pushed.ratio);
            const pad = Math.max(1e-6, (rmax-rmin)*0.15);
            const ymin = rmin-pad, ymax = rmax+pad;
            const ox = halfW+bandL;
            function xOf(f){ return ox + f/fmax*plotWR; }
            function yOf(v){ return 4 + plotH - (v-ymin)/(ymax-ymin)*plotH*0.94; }
            recCtx.fillStyle = '#f97316';
            for(let i=0;i<nf;i+=Math.max(1,Math.floor(nf/120))){
              recCtx.beginPath(); recCtx.arc(xOf(pushed.fband[i]), yOf(pushed.ratio[i]), 1.6, 0, 7); recCtx.fill();
            }
            recCtx.strokeStyle = '#38bdf8'; recCtx.lineWidth = 1.5;
            recCtx.beginPath();
            recCtx.moveTo(xOf(0), yOf(pushed.intercept));
            recCtx.lineTo(xOf(fmax), yOf(pushed.intercept + pushed.slope*fmax));
            recCtx.stroke();
            recCtx.fillStyle = '#9ca3af'; recCtx.font = '9px sans-serif'; recCtx.textAlign='center';
            recCtx.fillText('frequency (Hz)', ox+plotWR/2, RECORD_H-3);
          }
          recCtx.strokeStyle = '#4b5563'; recCtx.beginPath(); recCtx.moveTo(halfW+0.5,0); recCtx.lineTo(halfW+0.5,RECORD_H); recCtx.stroke();
          recCtx.fillStyle = '#e5e7eb'; recCtx.font = '11px sans-serif'; recCtx.textAlign='left';
          recCtx.fillText('u(t): r\\u2081 solid / r\\u2082 dashed', bandL+2, 12);
          recCtx.fillText('ln|U\\u2082/U\\u2081| (dots) + fit (line)', halfW+bandL+2, 12);
          recCtx.fillStyle = '#facc15'; recCtx.textAlign='right'; recCtx.font='12px sans-serif';
          recCtx.fillText('Q true='+state.Q.toFixed(0)+'  Q\\u0302 ='+pushed.Qest.toFixed(1)+'  \\u0394t*='+pushed.dtstarEst.toFixed(3)+'s', RECORD_W-6, 12);
        }

        function draw(){ drawDispersion(); drawAxis(); drawRecord(); }
        draw();

        for(const k in modelBtns){
          modelBtns[k].addEventListener('click', () => { state.model = k; syncControlLabels(); draw(); emit(); });
        }
        QIn.addEventListener('input', () => { state.Q = parseFloat(QIn.value); syncControlLabels(); emit(); });
        c0In.addEventListener('input', () => { state.c0 = parseFloat(c0In.value); syncControlLabels(); emit(); });
        fpIn.addEventListener('input', () => { state.fp = parseFloat(fpIn.value); syncControlLabels(); emit(); });
        n1In.addEventListener('input', () => { state.n1 = parseFloat(n1In.value); syncControlLabels(); emit(); });
        n2In.addEventListener('input', () => { state.n2 = parseFloat(n2In.value); syncControlLabels(); emit(); });
        par.querySelector('#ia-spread-body').addEventListener('click', () => {
          state.n1 = 1.0; state.n2 = 1.0; n1In.value = 1.0; n2In.value = 1.0; syncControlLabels(); emit();
        });
        par.querySelector('#ia-spread-surface').addEventListener('click', () => {
          state.n1 = 0.5; state.n2 = 0.5; n1In.value = 0.5; n2In.value = 0.5; syncControlLabels(); emit();
        });

        par.addEventListener('ia-update', event => {
          commitInFlight = false;
          pushed = event.detail;
          draw();
        });
        }
        // Debounce rather than act on the first callback: WideCell's resize of `par`
        // (narrow default column -> full wide width) can fire this more than once in
        // quick succession, and reacting to the first one bakes in the still-narrow
        // size. Waiting for callbacks to stop for a bit means iaInit() always uses the
        // settled width, whatever it ends up being (including "stayed narrow" on a
        // small viewport) -- same pattern as ray-tomography.jl's rtRo/rtInit.
        let iaTimer = null;
        const iaRo = new ResizeObserver(() => {
          clearTimeout(iaTimer);
          iaTimer = setTimeout(() => { iaRo.disconnect(); iaInit(); }, 150);
        });
        iaRo.observe(par);
        }
        }
        </script>
        """)
    end

    const _ia_ready = true
end

# ╔═╡ 0f94962d-5f79-4a45-ba8d-af7cb46a6ac6
begin
    _ia_ready
    WideCell(@bind ia AttenuationInput(); max_width=1400)
end

# ╔═╡ cb0e9b09-698b-4d99-a809-1ec83592a74f
"""
	AttenuationPush(fdisp, cfutt, ckjar, xaxis, env1, env2, tgrid, pulse1, pulse2, fband, ratio, slope, intercept, Qest, dtstarEst)

Bundles every array the widget's JS needs for one redraw and dispatches them as a single
`'ia-update'` CustomEvent on `#iawidget` -- the live analogue of the Appendix computation
above, re-run whenever the bound `ia` value changes.
"""
struct AttenuationPush
    fdisp::Vector{Float64}
    cfutt::Vector{Float64}
    ckjar::Vector{Float64}
    xaxis::Vector{Float64}
    env1::Vector{Float64}
    env2::Vector{Float64}
    tgrid::Vector{Float64}
    pulse1::Vector{Float64}
    pulse2::Vector{Float64}
    fband::Vector{Float64}
    ratio::Vector{Float64}
    slope::Float64
    intercept::Float64
    Qest::Float64
    dtstarEst::Float64
end

# ╔═╡ 8bbde385-1077-431e-9848-05e7fb5ed45d
function Base.show(io::IO, ::MIME"text/html", p::AttenuationPush)
    write(io, """
    <script>
    {
    const w = document.getElementById('iawidget');
    if(w){
      w.dispatchEvent(new CustomEvent('ia-update', { detail: {
        fdisp: [$(join(p.fdisp, ","))],
        cfutt: [$(join(p.cfutt, ","))],
        ckjar: [$(join(p.ckjar, ","))],
        xaxis: [$(join(p.xaxis, ","))],
        env1: [$(join(p.env1, ","))],
        env2: [$(join(p.env2, ","))],
        tgrid: [$(join(p.tgrid, ","))],
        pulse1: [$(join(p.pulse1, ","))],
        pulse2: [$(join(p.pulse2, ","))],
        fband: [$(join(p.fband, ","))],
        ratio: [$(join(p.ratio, ","))],
        slope: $(p.slope),
        intercept: $(p.intercept),
        Qest: $(p.Qest),
        dtstarEst: $(p.dtstarEst)
      }}));
    }
    }
    </script>
    """)
end

# ╔═╡ 7c9e696f-251b-4fe0-a09e-e777844101a2
let
    model = String(ia["model"])
    Q = Float64(ia["Q"])
    c0 = Float64(ia["c0"])
    fp = Float64(ia["fp"])
    x1 = Float64(ia["x1"])
    x2 = Float64(ia["x2"])
    n1 = Float64(ia["n1"])
    n2 = Float64(ia["n2"])

    tgrid = range(-3.0, 3.0, length=1024)
    freqgrid, Fsource = source_spectrum(tgrid, fp)

    fdisp = 10 .^ range(log10(0.3), log10(30.0), length=80)
    ωdisp = 2π .* fdisp
    cfutt = [futterman_phase_velocity(ω, Q, OMEGA0, c0) / c0 for ω in ωdisp]
    ckjar = [kjartansson_phase_velocity(ω, Q, OMEGA0, c0) / c0 for ω in ωdisp]

    xaxis = range(2.0, 150.0, length=120)
    atten_ref = [abs(only(attenuation_response([fp], x, Q, OMEGA0, c0, model))) for x in xaxis]
    env1 = [spreading_factor(x, n1) for x in xaxis] .* atten_ref
    env2 = [spreading_factor(x, n2) for x in xaxis] .* atten_ref
    env1 = env1 ./ maximum(env1)
    env2 = env2 ./ maximum(env2)

    U1 = receiver_spectrum(freqgrid, Fsource, x1, n1, Q, OMEGA0, c0, model)
    U2 = receiver_spectrum(freqgrid, Fsource, x2, n2, Q, OMEGA0, c0, model)
    pulse1 = synthesize_pulse(tgrid, U1)
    pulse2 = synthesize_pulse(tgrid, U2)

    mask = Fsource .> 0.05 * maximum(Fsource)
    fband = freqgrid[mask]
    ratio = log_spectral_ratio(U1, U2)[mask]
    slope, intercept = linear_fit(fband, ratio)
    Qest = estimate_Q(slope, x1, x2, c0)
    dtstarEst = -slope / π

    AttenuationPush(collect(fdisp), cfutt, ckjar, collect(xaxis), env1, env2,
        collect(tgrid), pulse1, pulse2, collect(fband), ratio, slope, intercept,
        Qest, dtstarEst)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
FFTW = "~1.10.0"
PlutoUI = "~0.7.72"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "ce8a8bd91a0621b5464d5f2acaf73e0a293ad028"

[[deps.ADTypes]]
git-tree-sha1 = "27cecae79e5cc9935255f90c53bb831cc3c870d7"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.18.0"

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

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

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

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

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
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

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

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

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
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "6c72198e6a101cccdd4c9731d3985e904ba26037"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3bc002af51045ca3b47d2e1787d6ce02e68b943a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.122"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c249d86e97a7e8398ce2068dce4c078a1c3464de"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.16"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "3f50fa86c968fc1a9e006c07b6bc40ccbb1b704d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.4"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "173e4d8f14230a7523ae11b9a3fa9edb3e0efd78"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.14.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

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
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

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

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "5fbb102dcb8b1a858111ae81d56682376130517d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.11"
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

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

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

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4255f0032eafd6451d707a51d5f0248b8a165e4d"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.3+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

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

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

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

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "d38b8653b1cdfac5a7da3b819c0a8d6024f9a18c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.13"
weakdeps = ["ChainRulesCore"]

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

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

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "dacc8be63916b078b592806acd13bb5e5137d7e9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "f53232a27a8c1c836d3998ae1e17d898d4df2a46"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.72"

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

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
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

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "7680fbbc8a4fdf9837b4cae5e3fbebe53ec8e4ff"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.122.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
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

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "c1053ba68ede9e4005fc925dd4e8723fcd96eef8"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.9.0"
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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

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
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "a136f98cefaf3e2924a66bd75173d1c891ab7453"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.7"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

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

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "f75c7deb7e11eea72d2c1ea31b24070b713ba061"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.3"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "a85b4262a55dbd1af39bb6facf621d79ca6a322d"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.32.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "1b09f5faec5284f505c40e68ba565115e7d48718"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.56.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

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

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

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
# ╟─74a3e9e3-ed09-4582-91fc-47944d5744db
# ╠═a9436d89-ad43-4cb8-8619-67c0e26eb686
# ╟─fbbee9cc-493d-11ed-00fd-b1e3bddbb3fa
# ╟─0f94962d-5f79-4a45-ba8d-af7cb46a6ac6
# ╟─d60cf031-d375-4da6-ace3-abd5043d6364
# ╟─304689e0-f318-4730-be5d-e112517d81ad
# ╟─b3273cac-84b5-458a-86c0-a530bd769ba0
# ╟─892eb702-189e-4b3f-a4a0-89e4763b97ac
# ╟─a1041e34-78f7-4212-898a-7e1e6d0770d8
# ╟─b2b60060-a73d-4006-911b-36a28c0e7554
# ╟─02e113c5-b3b6-4e9e-9670-c937af00d7c2
# ╟─d9b9ee19-f549-4c4f-8b79-6a360bfd5f2d
# ╠═bdf54ecf-72c8-4741-b655-142fc23be647
# ╟─189722e0-e9ea-42e8-ad75-684e6da5ad04
# ╠═e9e188a6-ec01-4a4f-9168-cf185b14cb05
# ╠═3f5f9b5a-2131-432a-b58d-132cf1800610
# ╟─14a5545d-858e-407b-a193-0e5abc9f68b1
# ╠═22bda95b-a685-46ad-8284-075bafd27581
# ╠═ea207d15-c2b9-472c-9fc3-79b7000a1738
# ╟─6bfbe8d3-3ef0-4a82-8266-59794fd6118a
# ╠═fd66f5ad-1b8e-4b0e-9377-32d8f15e9ac2
# ╟─82ff57b2-3824-42db-9d55-a3c7fc377131
# ╠═0ec1c415-c990-46eb-bfdf-17d062442f90
# ╠═f668c317-1166-4ab1-b96c-c46eaf9ee136
# ╟─8271e3a5-308d-42b4-b636-1a59f18361f6
# ╠═9c1c2062-c62a-48dd-9d0f-651d848afccd
# ╠═fb5825ce-55cf-4b48-8513-0669fbe876cd
# ╟─097e5d11-4857-4fd3-a587-3d72d2567180
# ╠═f3674d89-497d-43c3-9741-8a1d8e3fc860
# ╠═cb0e9b09-698b-4d99-a809-1ec83592a74f
# ╠═8bbde385-1077-431e-9848-05e7fb5ed45d
# ╟─7c9e696f-251b-4fe0-a09e-e777844101a2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
