### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Seismic Anisotropy"
#> layout = "layout.jlhtml"
#> tags = ["planewaves"]
#> description = "This notebook helps us visualize non-circular wavefronts in a homogenoeus Earth medium due to the presence of anisotropy."

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

# ╔═╡ 5f12d3e9-4d5c-434c-bdef-bcc74b79705e
begin
    using PlutoUI
    using PlutoTeachingTools
    using Symbolics
    using SymbolicUtils
    using LinearAlgebra
    using Einsum
    using Unitful
end

# ╔═╡ b48fc46c-47b2-4090-83fb-710b1974c2a2
ChooseDisplayMode()

# ╔═╡ 6e9dc83f-5d4f-4f21-a495-d43bc64f6041
TableOfContents()

# ╔═╡ 17942206-4ea3-11ed-30fa-d38c9c47f282
md"""
# Anisotropy
It is important to consider deviations from isotropy when imaging the Earth. Although there could be up to 21 independent linear elastic constants, we attach the term *anisotropy* to a situation where we use more than two elastic constants to describe the medium. 

* __Lattice-preferred anisotropy (LPO)__: homogeneous material, but there is preferred crystal orientation, e.g., olivine.
* __Shape-preferred anisotropy (SPO)__: a stack of rock layers with different isotropic properties cause seismic velocities to differ in different directions, or preferred orientation of cracks in the medium.

It is often difficult to distinguish the effects of anisotropy and those of medium heterogeneity. The presence of anisotropy results in non-circular wavefronts, even though the elastic constants are homogeneous. We can visualize wavefronts after projecting them along different planes by interacting with the plot below.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""


# ╔═╡ e28c94fd-f04d-442c-8b19-be84d565d5c7
nothing

# ╔═╡ 99e819ab-a88c-41a0-9bbf-27d49fe19557
begin
    # Superseded by the Julia-driven canvas widget near the top of this notebook.
    nothing
end

# ╔═╡ 56efaf7d-2b8b-4540-b0e8-b722edd1d8f8
md"""
## Read the velocity surface

The widget draws the full 3-D wavefront of a homogeneous VTI medium as one
rotatable surface: qP (amber), qS₁ (blue) and qS₂ (cyan) are drawn together,
radius equal to velocity times the chosen travel time. Drag inside the canvas
to orbit; a sphere would mean speed is independent of direction, so any
stretching or dimpling away from a sphere *is* the anisotropy.

Try first making `A = C` and `L = N`: all three surfaces collapse toward
spheres, the isotropic reference. Then restore contrast between the
horizontal (`A`, `N`) and vertical (`C`, `L`) constants and watch the surfaces
elongate along the symmetry axis while staying circular in cross-section
around it — that circular cross-section (rotate to look straight down the
axis) is what makes this medium VTI rather than fully general anisotropy.

Two more things worth toggling:
* **Phase vs. group velocity.** The segmented control switches between the
  phase-velocity surface (loci of equal wave-phase, what the Christoffel
  eigenproblem below actually solves for) and the group-velocity — or *ray* —
  surface (where the wave's energy has actually travelled to). They agree
  for isotropic media and can diverge sharply for strong anisotropy; the
  group surface is the one that answers "where does the wavefront actually
  go after time `t`."
* **Crystal orientation.** The Plunge/Trend sliders tilt the medium's
  symmetry axis (its LPO fabric axis) away from vertical. Turn on the
  "crystal axis" legend item to see the axis itself as a magenta arrow, and
  watch the whole velocity surface tilt rigidly with it — this is the direct
  visual link between *aligned crystals* and *directional wave speed*.
"""

# ╔═╡ 3bfab7cf-8b35-408d-8d25-ce6a37a0f19e
md"""
| **Love Constant** | **Equivalent \(c_{ij}\)** | **Physical Meaning** |
|--------------------|---------------------------|----------------------|
| \(A\) | \(c_{11}\) | P-wave modulus in the horizontal (x–y) plane |
| \(C\) | \(c_{33}\) | P-wave modulus along the symmetry axis (z) |
| \(F\) | \(c_{13}\) | Coupling between vertical and horizontal strain |
| \(L\) | \(c_{44}\) | Shear modulus for SV motion in the x–z plane |
| \(N\) | \(c_{66}\) | Shear modulus for SH motion in the x–y plane |
"""

# ╔═╡ d2027f44-361b-4d00-9efb-86fe2c7e9b68
aside(tip(md"**Shear-wave splitting** As the two quasi-shear waves travel with different speeds in the case of an anisotropic medium, it leads to shear-wave splitting. Which means, the two orthogonal components of the shear wave reach at slightly different times."))

# ╔═╡ 7d5bcebe-3f1c-485a-a094-78eb3b976352
md"""
## Transverse Anisotropy; Voigt Notation
A transversely anisotropic material can be characterized by five independent elastic constants. The Voigt matrix is a six-dimensional symmetric matrix denoted using `C`.
We follow the Love convention, where the elastic constants are represented using the following symbols.
"""

# ╔═╡ bba0ec5e-32bd-4ab1-84d1-76e18c88dc5a
md"""
$(LocalResource("../assets/images/rock_strata.jpg", :width => 300))

*Example of SPO: transverse anisotropy because of layered rock materials.*
"""

# ╔═╡ d000b26c-e67a-4e20-a2b8-b417b596ccd8
md"""
| **Type** | **Symmetry** | **Typical Cause** | **Example** |
|-----------|---------------|-------------------|--------------|
| **Transverse isotropy (TI)** | Axis of symmetry | Layered media or aligned cracks | Sedimentary basins |
| **Vertical TI (VTI)** | Axis vertical | Bedding or layering | Marine sediments |
| **Horizontal TI (HTI)** | Axis horizontal | Vertical cracks | Reservoirs |
| **Azimuthal anisotropy** | Depends on propagation direction | Aligned fractures | Upper crust |
| **Shear-wave splitting** | Two orthogonal shear waves with different speeds | Mantle fabrics | SKS phases |
"""

# ╔═╡ 8a6982d7-9d05-4d06-bcba-33cd5b79034a
md"""
In global seismology, SKS splitting is the most common measurement:
* The incoming SKS wave (nearly vertically incident) passes through azimuthal anisotropy (often HTI-like) in the upper mantle.
* Fast polarization direction \phi → indicates mantle fabric / flow direction.
* Delay time \delta t → measures strength of anisotropy or layer thickness.
"""

# ╔═╡ 1724771c-d1ef-4c1c-93dc-e889a74d80e9
md"""
| **Observable** | **Measured from** | **Related Constant** |
|-----------------|-------------------|----------------------|
| \(V_{P0}\) | Vertical P-wave | \(C\) |
| \(V_{S0}\) | Vertical S-wave | \(L\) |
| \(V_{PH}\) | Horizontal P-wave | \(A\) |
| \(V_{SHH}\) | Horizontal SH-wave | \(N\) |
"""

# ╔═╡ 846457f5-d684-470f-bba5-144b41f8dd14
md"""
| **Type** | **Example** | **Splitting?** | **Explanation** |
|-----------|-------------|----------------|-----------------|
| **VTI** | Layered sediments | Only for oblique incidence | No splitting for vertical propagation |
| **HTI** | Aligned vertical cracks | ✅ Yes, clear splitting for horizontal propagation | Fast direction = crack alignment |
| **General anisotropy (orthorhombic, triclinic)** | Complex fabrics | ✅ Yes | Two shear modes with distinct speeds |
"""

# ╔═╡ 4b2d3084-e981-4ab2-bf40-fea4335082ad
@variables A::Real N::Real F::Real C::Real L::Real

# ╔═╡ 70ebdd20-3a02-4493-a0d6-718fa5b51675
md"In the case of transverse anisotropy, with rock layers oriented along in the $xy$-plane, the Love constants can be used from the Voigt matrix as follows. In this case, the axis of symmetry is along $z$."

# ╔═╡ 4f7b0a05-004a-417d-81d2-e4b1176909d6
Ctrans = [[A, A - 2N, F, 0, 0, 0];; [A - 2N, A, F, 0, 0, 0];; [F, F, C, 0, 0, 0];; [0, 0, 0, L, 0, 0];; [0, 0, 0, 0, L, 0];; [0, 0, 0, 0, 0, N]]

# ╔═╡ 1de3540c-6e07-4e4e-b411-7a5667b5c742
aside(tip(md"For layered rocks, the seismic waves travel faster in the direction parallel to the layers, as opposed to the perpendicular direction. This is because the waves can *choose* to travel in the fast layers in the former case."))

# ╔═╡ f5f78c1e-9062-4a5b-8bc3-4efb0f1883ca
md"As an example of LPO, we now construct the Voigt matrix for olivine, where the constants are in $(u\"GPa\")."

# ╔═╡ 5c3cea52-3154-4aa3-824c-5be50a302be5
Colivine = [[192, 66, 60, 0, 0, 0];; [66, 160, 56, 0, 0, 0];; [60, 56, 272, 0, 0, 0];; [0, 0, 0, 60, 0, 0];; [0, 0, 0, 0, 62, 0];; [0, 0, 0, 0, 0, 49]]

# ╔═╡ 88fb688e-d768-4895-8fe9-822f5db335b0
md"""
For an isotropic medium, there are only two independent elastic constants $\lambda$ and $\mu$. In this case, the Voigt matrix can be derived using the following substitutions.
"""

# ╔═╡ b64d6d57-c38c-4f54-b86f-0d254d27c3ed
@variables λ::Real μ::Real

# ╔═╡ 1fa4ecc5-8c5d-49f6-a023-bb2c4b4d5c27
Ciso = map(Ctrans) do x
    substitute(x, [A => λ + 2μ, C => λ + 2 * μ, N => μ, F => λ, L => μ])
end

# ╔═╡ 12a42fbb-bfb7-4852-a638-2dec237c032d
md"""
## Elastic tensor `cijkl`
The Voigt matrix is not a tensor and no longer preserves the mathematical properties of the elastic tensor. Given the elements of the Voigt matrix, we can construct the tensor `cijkl`, for both isotropic and transverse anisotropy cases, using the following method. 
"""

# ╔═╡ 0ea0a70c-72b2-47da-902d-00cd6279fbb5
function get_cijkl(C)
    [C[i*(isequal(i, j))+(1-isequal(i, j))*(9-i-j), k*isequal(k, l)+(1-isequal(k, l))*(9-k-l)] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
end

# ╔═╡ 2d4ab726-4257-4fc4-aac0-2cc90ce76e53
ciso = get_cijkl(Ciso); # for isotropic 

# ╔═╡ eb482986-f6b0-4900-9e4e-441351691e12
ciso[:,1,:,1]

# ╔═╡ f2c3eeca-4e03-4360-b14b-4fd98e24d530
ctrans = get_cijkl(Ctrans); # for transverse anisotropy

# ╔═╡ 4ac951fc-0054-4119-9296-fa8feacdd4c2
colivine = get_cijkl(Colivine); # for anisotropy

# ╔═╡ 0028f23a-0a92-42d8-9077-cc5fd585cce6
md"Let's take a moment to realize the power of our notation using `Symbolics` and `Einsum` packages. We shall now write down something that we were familiar in the isotropic world: the stress-strain relation."

# ╔═╡ 8aa7f639-3e76-4db4-bd44-37b78f7f3781
@variables e[1:3, 1:3] # strain tensor

# ╔═╡ c6eefa83-676f-464b-85f1-4c7442100cbf
@einsum σiso[i, j] := ciso[i, j, k, l] * e[k, l]

# ╔═╡ 183cc6c5-a4ae-47c6-8c48-2c8ff2e1cbc6
md"""
Doesn't it look easy? If you were confused about the Voigt matrix above, we can simply understand it by writing the stress-strain relation in the transversely anisotropic medium.
"""

# ╔═╡ 4af3c346-f657-4b59-9a79-ee9ff0b79bbd
@einsum σtrans[i, j] := ctrans[i, j, k, l] * e[k, l]

# ╔═╡ 6bbea2f3-da33-4abe-bf71-e6211275ce41
@einsum σolivine[i, j] := colivine[i, j, k, l] * e[k, l]

# ╔═╡ 24553b91-78cb-4cde-b0d6-8dbbd7e671d6
@variables ρ::Real ω::Real t::Real

# ╔═╡ f940ed1a-1f89-4f7b-bb41-05615e2351a0
md"""
## Planewaves in Anisotropic Media
It is important to realize that we are going to work with homogeneous media, where all the elastic constants in `cikjl` don't vary with spatial coordinates. We shall now analyze plane wave solutions.
"""

# ╔═╡ 54657fec-e8e7-4040-8108-a29ae8df1c24
@variables u[1:3] # vector displacement field x[1:3] g[1:3] p::Real

# ╔═╡ 108996f0-8980-4aff-8334-5f9f9f6f09d2
@variables x[1:3] # position vector ~ [x, y, z]

# ╔═╡ 5b1f23d8-33d1-4d25-92e3-ca590fc0118d
@variables g[1:3] # amplitude vector (we already used A for one of the elastic constants)

# ╔═╡ 7f0b2df9-6b33-44cb-a1c0-d4842cb666da
@variables s[1:3] # unit slowness vector  (use shat instead?)

# ╔═╡ ffeb806f-875f-4bc8-a15f-3d8db4252e96
@variables p # slowness magnitude

# ╔═╡ 32fdd4ec-afef-4d03-8081-71a8201c3fcf
@variables v # velocity (inverse of slowness)

# ╔═╡ c9cc3ef9-54ed-41bd-9c0b-c4c38b8c613b
@syms ı # imaginary number

# ╔═╡ 05c9f91b-06ee-4cf9-9cba-ae40aad8bd12
D = [Differential(x[1]), Differential(x[2]), Differential(x[3])]

# ╔═╡ f022c469-6ab4-4b2b-aee7-70e4a84d5a5a
Dt = Differential(t)

# ╔═╡ 7d38e878-0136-46c6-b265-85b7dc85153e
utrail1 = exp(ı * ω * (t - p * dot(s, x)))

# ╔═╡ 6c964379-2169-4567-9505-745198609d49
utrail = Symbolics.scalarize(g .* utrail1)

# ╔═╡ 0b2cf81d-98c3-42f9-a688-7069358003d6
@einsum uddot[i] := ctrans[i, j, k, l] * D[j](D[l](u[k]))

# ╔═╡ dbbc72ec-64a8-4566-b044-44fbbd1d4cbf
@einsum uddot_trail1[i] := ctrans[i, j, k, l] * expand_derivatives(D[j](D[l](utrail[k])))

# ╔═╡ 12d0f4cd-ffe1-44e2-80e5-25dc5ab0dfd5
@einsum uddot_trail2[i] := expand_derivatives(Dt(Dt(utrail[i]))) * ρ

# ╔═╡ 2a2fc196-10a9-4c97-b1c1-e2a68e21119e
md"The equation of motion is satisfied only when the following two expressions are equal to each other."

# ╔═╡ 77f96830-e81b-438d-b5ae-4ac59cb031df
ga = simplify_fractions.(uddot_trail1 ./ utrail1 ./ ω^2 ./ ı^2 / p^2)

# ╔═╡ 27f97bab-bef0-4705-bc22-8e9b087444b2
gb = simplify_fractions.(uddot_trail2 ./ utrail1 ./ ω^2 ./ ı^2 * v^2)

# ╔═╡ 366ee997-e2fc-4b48-bb6e-33b666ea38f5
md"""
## Christoffel Matrix `M`
Now we would like to work towards the construction of a transformation matrix `M` such that `ga=M*g`. Notice that `gb` is already in the direction of the amplitude vector, so once we construct `M`, all we need to solve is an eigenvalue problem
```math
M*g = ρv^2*g,
```
where the eigenvalues will give us the phase velocities. Note that `M` depends on the direction of the slowness vector `s`, therefore leading to anisotropy. As `M` is symmetric, we have three solutions for general anisotropic case: a single quasi-P wave (qP) and two quasi-S waves (qS).
"""

# ╔═╡ f6675c52-740c-4098-98d6-9c9038b4f6f4
M1 = map(ga) do x
    simplify_fractions(substitute(x, [g[2] => 0, g[3] => 0]) / g[1])
end

# ╔═╡ 97da0c08-b971-4104-88fe-aa1d1f2c68b0
M2 = map(ga) do x
    simplify_fractions(substitute(x, [g[1] => 0, g[3] => 0]) / g[2])
end

# ╔═╡ cbdb003e-09a2-419b-8811-7d0cd290968d
M3 = map(ga) do x
    simplify_fractions(substitute(x, [g[2] => 0, g[1] => 0]) / g[3])
end

# ╔═╡ edfc88be-4bef-4478-b154-ec920782c295
# M = symbolics_to_sympy.([M1;; M2;; M3])

# ╔═╡ b22a4fcd-1f0c-410f-95b8-f67671d3aa63
M = ([M1;; M2;; M3]) # concat all the columns to get the Christoffel matrix

# ╔═╡ b5997b61-2ca2-4c2f-acd6-16c121f0aa13
md"""
We shall now analyze plane waves traveling along the cartesian axis individually. Towards that end, we will make necessary substitutions for the slowness vector.
The Chistoffel matrices for plane waves traveling in `x[1]`, `x[2]`, and `x[3]` directions are given by the following expressions respectively.
"""

# ╔═╡ 2bdba346-7a64-4de6-9b71-7034b9d68986
aside(correct(md"For an isotropic medium, note that the Christoffel matrices didn't change with the direction of the slowness vector. Obviously. Its eigenvalues correspond to P, SV and SH waves."))

# ╔═╡ 78449122-6762-4191-9a8d-b1f8bee69bd8
map(M) do x
    simplify(substitute(x, [s[1] => 1, s[2] => 0, s[3] => 0]))
end

# ╔═╡ 3bbbfa0d-c8ce-44f5-8546-16c9a1742e23
map(M) do x
    simplify(substitute(x, [s[1] => 0, s[2] => 1, s[3] => 0]))
end

# ╔═╡ 747ca67e-60d7-44a4-bac1-78cfa9b20cc7
map(M) do x
    simplify(substitute(x, [s[1] => 0, s[2] => sqrt(0.5), s[3] => sqrt(0.5)]))
end

# ╔═╡ 4a66af73-6b29-475e-9c71-8e62ee69802a
md"""
## References
* Jules Thomas Browaeys, Sébastien Chevrot, Decomposition of the elastic tensor and geophysical applications, Geophysical Journal International, Volume 159, Issue 2, November 2004, Pages 667–678, https://doi.org/10.1111/j.1365-246X.2004.02415.x
	"""

# ╔═╡ ebf4953d-abc5-4530-af51-8b68096a0113
md"## Appendix"

# ╔═╡ 0989b0f0-06c6-4353-9933-ec713563f0c6
md"""
## What this model shows—and its limits

This is a homogeneous VTI medium: every difference in wavefront shape comes
from the stiffness tensor and its orientation, not from a spatial velocity
anomaly — there is no heterogeneity anywhere in this widget. Within that
scope the 3-D view is now fairly complete: it shows phase velocity *and* the
group/ray-velocity surface, *and* per-mode particle-motion (polarization)
directions, *and* lets the symmetry axis itself be tilted to stand in for a
rotated crystal fabric.

What it still does not show: wave amplitude, attenuation, and mode
conversion at interfaces (that needs a boundary, which this whole-space
model doesn't have) — see the reflection/transmission notebooks for that —
and it is restricted to transverse isotropy (5 constants), not fully general
triclinic anisotropy (21 constants), so it cannot represent, e.g., the
orthorhombic asymmetry of a single olivine crystal beyond the VTI
approximation used here.
"""

# ╔═╡ aa111111-1111-4111-8111-111111111111
"""
	vti_wavefront_payload_3d(A, C, L, N, F, density, time, plunge_deg, trend_deg;
		ntheta=17, nphi=33, arrow_step=4)

Sample the full 3-D phase-velocity surface (qP, qS1, qS2) of a homogeneous VTI
medium over a `(theta, phi)` grid on the unit sphere, plus the corresponding
group-velocity (ray) surface and per-mode particle-motion (polarization)
directions. Love constants are in GPa, `density` in g/cm³, so `sqrt(C/density)`
comes out in km/s; multiplying by `time` gives plotted distance in km.

The VTI symmetry axis is tilted away from vertical `z` by `plunge_deg` (angle
from vertical) along azimuth `trend_deg`, modelling a lattice-preferred
orientation whose fabric axis is not vertical. Phase velocity in a lab-frame
direction `n` is evaluated by rotating `n` into the (fixed) crystal frame
before contracting with the VTI stiffness tensor — a physically exact
shortcut since phase velocity is a frame-invariant scalar of direction.

Group velocity uses the standard tangential-gradient construction
``U = v\\,\\hat n + \\nabla_S v`` (``\\nabla_S`` the gradient along the unit
sphere), evaluated with finite differences on the sampled grid: it is the
group/ray velocity, i.e. where a wave packet's *energy* actually travels,
which for a genuinely anisotropic medium is not generally parallel to the
phase-velocity direction `n` itself.

Returns a named tuple with the sampled grid size, flattened Cartesian vertex
arrays for the phase and group-velocity surfaces (dims `(mode, itheta, iphi,
xyz)`), subsampled polarization arrow data, the rotated crystal-axis triad,
and stability/anisotropy summary scalars.
"""
function vti_wavefront_payload_3d(
    A::Float64, C::Float64, L::Float64, N::Float64, F::Float64,
    density::Float64, time::Float64,
    plunge_deg::Float64, trend_deg::Float64;
    ntheta::Int=17, nphi::Int=33, arrow_step::Int=4,
)
    stiffness = [A A-2N F 0.0 0.0 0.0;
                 A-2N A F 0.0 0.0 0.0;
                 F F C 0.0 0.0 0.0;
                 0.0 0.0 0.0 L 0.0 0.0;
                 0.0 0.0 0.0 0.0 L 0.0;
                 0.0 0.0 0.0 0.0 0.0 N]
    pairs = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
    index(i, j) = findfirst(p -> (i == p[1] && j == p[2]) || (i == p[2] && j == p[1]), pairs)
    c = Array{Float64}(undef, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        c[i, j, k, l] = stiffness[index(i, j), index(k, l)]
    end

    θ0, φ0 = deg2rad(plunge_deg), deg2rad(trend_deg)
    cθ0, sθ0 = cos(θ0), sin(θ0)
    cφ0, sφ0 = cos(φ0), sin(φ0)
    Ry = [cθ0 0.0 sθ0; 0.0 1.0 0.0; -sθ0 0.0 cθ0]
    Rz = [cφ0 -sφ0 0.0; sφ0 cφ0 0.0; 0.0 0.0 1.0]
    R = Rz * Ry # crystal symmetry axis in the lab frame is R*[0,0,1]
    Rt = R'

    θs = collect(range(0.0, π, length=ntheta))
    φs = [(i - 1) * 2π / nphi for i in 1:nphi]

    vphase = Array{Float64}(undef, 3, ntheta, nphi)
    polvec = Array{Float64}(undef, 3, 3, ntheta, nphi)

    stable = true
    for iθ in 1:ntheta, iφ in 1:nphi
        θ, φ = θs[iθ], φs[iφ]
        sθ, cθ = sincos(θ)
        sφ, cφ = sincos(φ)
        n_lab = (sθ * cφ, sθ * sφ, cθ)
        n_cry = (Rt[1, 1] * n_lab[1] + Rt[1, 2] * n_lab[2] + Rt[1, 3] * n_lab[3],
                 Rt[2, 1] * n_lab[1] + Rt[2, 2] * n_lab[2] + Rt[2, 3] * n_lab[3],
                 Rt[3, 1] * n_lab[1] + Rt[3, 2] * n_lab[2] + Rt[3, 3] * n_lab[3])
        Γ = zeros(3, 3)
        for i in 1:3, k in 1:3, j in 1:3, l in 1:3
            Γ[i, k] += c[i, j, k, l] * n_cry[j] * n_cry[l]
        end
        E = eigen(Symmetric(Γ))
        stable &= minimum(E.values) > 0.0
        for m in 1:3
            λ = max(E.values[m], 0.0)
            vphase[m, iθ, iφ] = sqrt(λ / density)
            polvec[m, :, iθ, iφ] .= R * E.vectors[:, m] # rotate polarization back to lab frame
        end
    end

    dθ = θs[2] - θs[1]
    dφ = φs[2] - φs[1]
    phase_xyz = Array{Float64}(undef, 3, ntheta, nphi, 3)
    group_xyz = Array{Float64}(undef, 3, ntheta, nphi, 3)
    for m in 1:3, iθ in 1:ntheta, iφ in 1:nphi
        θ = θs[iθ]
        sθ, cθ = sincos(θ)
        φ = φs[iφ]
        sφ, cφ = sincos(φ)
        n = (sθ * cφ, sθ * sφ, cθ)
        θhat = (cθ * cφ, cθ * sφ, -sθ)
        φhat = (-sφ, cφ, 0.0)

        dvdθ = if iθ == 1
            (vphase[m, 2, iφ] - vphase[m, 1, iφ]) / dθ
        elseif iθ == ntheta
            (vphase[m, ntheta, iφ] - vphase[m, ntheta-1, iφ]) / dθ
        else
            (vphase[m, iθ+1, iφ] - vphase[m, iθ-1, iφ]) / (2dθ)
        end
        iφprev = iφ == 1 ? nphi : iφ - 1
        iφnext = iφ == nphi ? 1 : iφ + 1
        dvdφ = (vphase[m, iθ, iφnext] - vphase[m, iθ, iφprev]) / (2dφ)
        invsinθ = 1.0 / max(sθ, 1e-9)

        v = vphase[m, iθ, iφ]
        U = v .* n .+ dvdθ .* θhat .+ (invsinθ * dvdφ) .* φhat

        phase_xyz[m, iθ, iφ, :] .= time .* v .* n
        group_xyz[m, iθ, iφ, :] .= time .* U
    end

    arrow_pts = NTuple{6,Float64}[]
    arrow_mode = Int[]
    for m in 1:3, iθ in 1:arrow_step:ntheta, iφ in 1:arrow_step:nphi
        p = @view phase_xyz[m, iθ, iφ, :]
        d = @view polvec[m, :, iθ, iφ]
        push!(arrow_pts, (p[1], p[2], p[3], d[1], d[2], d[3]))
        push!(arrow_mode, m)
    end

    qP = @view vphase[3, :, :]
    anisotropy_percent = 100 * (maximum(qP) - minimum(qP)) / max(sum(qP) / length(qP), eps())

    return (
        ntheta=ntheta, nphi=nphi,
        phase_xyz=phase_xyz, group_xyz=group_xyz,
        arrow_pts=arrow_pts, arrow_mode=arrow_mode,
        axis_sym=R * [0.0, 0.0, 1.0], axis_a=R * [1.0, 0.0, 0.0], axis_b=R * [0.0, 1.0, 0.0],
        stable=stable, anisotropy_percent=anisotropy_percent,
    )
end

# ╔═╡ aa333333-3333-4333-8333-333333333333
"""
	ani_push_message(payload)

Serialize a [`vti_wavefront_payload_3d`](@ref) result into the flat JSON the
widget's `<script>` expects: mode-major flattened vertex arrays for the phase
and group-velocity surfaces (so a JS index `((mode*ntheta+itheta)*nphi+iphi)*3`
recovers each `x,y,z`), the polarization arrow list, and the crystal-axis
triad.
"""
function ani_push_message(payload)
    number(x) = isfinite(x) ? string(round(Float64(x), digits=6)) : "0"
    arr(values) = "[" * join(number.(values), ",") * "]"
    # (mode,itheta,iphi,xyz) -> flatten with xyz fastest, then iphi, then itheta, then mode slowest
    flat(A) = vec(permutedims(A, (4, 3, 2, 1)))

    arrows_json = join([
        "{\"mode\":$(m),\"p\":$(arr([p[1], p[2], p[3]])),\"d\":$(arr([p[4], p[5], p[6]]))}"
        for (m, p) in zip(payload.arrow_mode, payload.arrow_pts)
    ], ",")

    return string(
        "{\"ntheta\":", payload.ntheta, ",\"nphi\":", payload.nphi,
        ",\"phaseXYZ\":", arr(flat(payload.phase_xyz)),
        ",\"groupXYZ\":", arr(flat(payload.group_xyz)),
        ",\"arrows\":[", arrows_json, "]",
        ",\"axisSym\":", arr(payload.axis_sym),
        ",\"axisA\":", arr(payload.axis_a),
        ",\"axisB\":", arr(payload.axis_b),
        ",\"stable\":", payload.stable,
        ",\"anisotropyPercent\":", number(payload.anisotropy_percent),
        "}",
    )
end

# ╔═╡ aa222222-2222-4222-8222-222222222222
begin
    struct VTIWavefrontInput
        time::Float64
        density::Float64
        A::Float64
        C::Float64
        L::Float64
        N::Float64
        F::Float64
        plunge::Float64
        trend::Float64
    end

    VTIWavefrontInput(; time=0.8, density=4.0, A=272.0, C=160.0, L=60.0, N=50.0, F=60.0, plunge=0.0, trend=0.0) =
        VTIWavefrontInput(Float64(time), Float64(density), Float64(A), Float64(C), Float64(L), Float64(N), Float64(F), Float64(plunge), Float64(trend))

    Base.get(w::VTIWavefrontInput) = Dict{String,Any}(
        "time" => w.time, "density" => w.density, "A" => w.A, "C" => w.C,
        "L" => w.L, "N" => w.N, "F" => w.F, "plunge" => w.plunge, "trend" => w.trend,
    )

    function Base.show(io::IO, ::MIME"text/html", w::VTIWavefrontInput)
        write(io, """
        <div id="ani-widget"><style>
        #ani-widget{width:100%;max-width:1400px;margin:auto;color:#e5e7eb;font:14px system-ui,sans-serif}#ani-widget *{box-sizing:border-box}
        #ani-widget .ani-title{padding:10px 14px;margin-bottom:10px;text-align:center;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px}#ani-widget .ani-title b{font-size:17px}#ani-widget .ani-hint{margin-top:3px;color:#9ca3af;font-size:13px}
        #ani-widget .ani-workspace{display:block}#ani-widget .ani-panel,#ani-widget .ani-group{min-width:0;padding:10px;background:#050505;border:1px solid #2f3744;border-radius:6px}
        #ani-widget .ani-panel-title{font-size:15px;font-weight:700;color:#f3f4f6}#ani-widget .ani-caption{min-height:20px;margin:4px 0 7px;color:#9ca3af;font-size:13px;line-height:1.3}
        #ani-widget canvas{display:block;width:100%;height:480px;background:#000;border:1px solid #374151;border-radius:4px;cursor:grab}#ani-widget canvas:active{cursor:grabbing}
        #ani-widget .ani-legend{display:flex;gap:14px;flex-wrap:wrap;align-items:center;margin-top:9px;font-size:12.5px}#ani-widget .ani-key{display:inline-flex;align-items:center;gap:5px;cursor:pointer;user-select:none}.ani-dot{width:18px;height:3px;border-radius:2px;background:currentColor}#ani-widget .ani-key input{accent-color:currentColor;margin:0}
        #ani-widget .ani-seg{display:inline-flex;border:1px solid #4b5563;border-radius:4px;overflow:hidden;margin-left:auto}#ani-widget .ani-seg button{border:none;background:#111827;color:#9ca3af;padding:5px 10px;font-size:12.5px;cursor:pointer}#ani-widget .ani-seg button.active{background:#3b5c85;color:#fff}
        #ani-widget .ani-controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:8px;margin-top:10px}#ani-widget .ani-group-title{margin-bottom:7px;font-size:16px;font-weight:700}#ani-widget .ani-row{display:grid;grid-template-columns:minmax(70px,115px) minmax(70px,1fr) minmax(42px,62px);gap:7px;align-items:center;margin:7px 0;color:#d1d5db}#ani-widget input[type=range]{width:100%;min-width:0;accent-color:#f59e0b}#ani-widget .ani-value{color:#fbbf24;text-align:right;font-variant-numeric:tabular-nums}#ani-widget button.ani-btn{border:1px solid #9ca3af;border-radius:4px;padding:6px 12px;background:#606060;color:#f3f4f6;font-size:14px;cursor:pointer}#ani-widget .ani-status{margin-top:8px;color:#9ca3af;font-size:13px}
        @media(max-width:760px){#ani-widget canvas{height:340px}}
        </style><div class="ani-title"><b>Directional stiffness turns a spherical wavefront into a velocity surface.</b><div class="ani-hint">Drag the 3-D view to orbit &middot; tilt the crystal &middot; compare phase and ray (group) velocity surfaces</div></div>
        <div class="ani-workspace"><section class="ani-panel"><div class="ani-panel-title">3-D wavefront</div><div class="ani-caption" id="ani-caption">radius = velocity &times; time, in the direction shown</div><canvas id="ani-3d"></canvas>
        <div class="ani-legend">
        <label class="ani-key" id="ani-key-qP" style="color:#f59e0b"><input type="checkbox" id="ani-show-qP" checked><i class="ani-dot"></i>qP</label>
        <label class="ani-key" id="ani-key-qS1" style="color:#60a5fa"><input type="checkbox" id="ani-show-qS1" checked><i class="ani-dot"></i>qS&#8321;</label>
        <label class="ani-key" id="ani-key-qS2" style="color:#22d3ee"><input type="checkbox" id="ani-show-qS2" checked><i class="ani-dot"></i>qS&#8322;</label>
        <label class="ani-key" style="color:#e879f9"><input type="checkbox" id="ani-show-axes"><i class="ani-dot"></i>crystal axis</label>
        <label class="ani-key" style="color:#d1d5db"><input type="checkbox" id="ani-show-pol"><i class="ani-dot"></i>polarization</label>
        <div class="ani-seg" id="ani-seg-surface"><button data-v="phase" class="active">Phase velocity</button><button data-v="group">Group velocity (ray)</button></div>
        </div>
        </section></div>
        <div class="ani-controls">
        <section class="ani-group"><div class="ani-group-title">Propagation</div><label class="ani-row"><span>Time</span><input id="ani-time" type="range" min="0.2" max="1.2" step="0.1" value="$(w.time)"><span id="ani-time-v" class="ani-value"></span></label><label class="ani-row"><span>Density</span><input id="ani-density" type="range" min="2.5" max="5" step="0.1" value="$(w.density)"><span id="ani-density-v" class="ani-value"></span></label></section>
        <section class="ani-group"><div class="ani-group-title">P-wave Love constants</div><label class="ani-row"><span>A</span><input id="ani-A" type="range" min="50" max="300" step="1" value="$(w.A)"><span id="ani-A-v" class="ani-value"></span></label><label class="ani-row"><span>C</span><input id="ani-C" type="range" min="50" max="300" step="1" value="$(w.C)"><span id="ani-C-v" class="ani-value"></span></label><label class="ani-row"><span>F</span><input id="ani-F" type="range" min="-100" max="150" step="1" value="$(w.F)"><span id="ani-F-v" class="ani-value"></span></label></section>
        <section class="ani-group"><div class="ani-group-title">Shear Love constants</div><label class="ani-row"><span>L</span><input id="ani-L" type="range" min="20" max="140" step="1" value="$(w.L)"><span id="ani-L-v" class="ani-value"></span></label><label class="ani-row"><span>N</span><input id="ani-N" type="range" min="20" max="140" step="1" value="$(w.N)"><span id="ani-N-v" class="ani-value"></span></label></section>
        <section class="ani-group"><div class="ani-group-title">Crystal orientation (LPO fabric axis)</div><label class="ani-row"><span>Plunge</span><input id="ani-plunge" type="range" min="0" max="90" step="1" value="$(w.plunge)"><span id="ani-plunge-v" class="ani-value"></span></label><label class="ani-row"><span>Trend</span><input id="ani-trend" type="range" min="0" max="360" step="1" value="$(w.trend)"><span id="ani-trend-v" class="ani-value"></span></label><button id="ani-reset" type="button" class="ani-btn" style="margin-top:4px">Reset VTI example</button></section>
        </div><div id="ani-status" class="ani-status"></div></div>
        <script>
        (function(){
        const root=document.currentScript.previousElementSibling;
        const by=id=>root.querySelector('#'+id);
        const state={time:$(w.time),density:$(w.density),A:$(w.A),C:$(w.C),L:$(w.L),N:$(w.N),F:$(w.F),plunge:$(w.plunge),trend:$(w.trend)};
        const ids={"ani-time":"time","ani-density":"density","ani-A":"A","ani-C":"C","ani-L":"L","ani-N":"N","ani-F":"F","ani-plunge":"plunge","ani-trend":"trend"};
        let data=null;
        let surfaceMode='phase';
        let showAxes=false, showPol=false;
        let visible={qP:true,qS1:true,qS2:true};
        let cam={az:0.6,el:0.35};
        let camD=10; // camera distance in DATA units; rescaled per-draw to the data's own extent
        function labels(){
          for(const [id,key] of Object.entries(ids)){
            const v=state[key];
            by(id+'-v').textContent = key==='time' ? v.toFixed(1)+' s' : key==='density' ? v.toFixed(1)+' g/cm³' : (key==='plunge'||key==='trend') ? v.toFixed(0)+'°' : v.toFixed(0)+' GPa';
          }
        }
        function emit(){ root.value={...state}; root.dispatchEvent(new CustomEvent('input')); }
        function setup(c){
          const r=c.getBoundingClientRect(), d=devicePixelRatio||1;
          c.width=Math.round(r.width*d); c.height=Math.round(r.height*d);
          const x=c.getContext('2d'); x.setTransform(d,0,0,d,0,0);
          return [x, r.width, r.height];
        }
        function rot(p, az, el){
          let [x,y,z]=p;
          let x1=x*Math.cos(az)+z*Math.sin(az);
          let z1=-x*Math.sin(az)+z*Math.cos(az);
          let y2=y*Math.cos(el)-z1*Math.sin(el);
          let z2=y*Math.sin(el)+z1*Math.cos(el);
          return [x1,y2,z2];
        }
        function project(p, W, H, scale){
          const [x,y,z]=rot(p, cam.az, cam.el);
          const f=1/(1+z/camD);
          return {sx:W/2+x*scale*f, sy:H/2-y*scale*f, depth:z, f};
        }
        const MODE_COLOR={1:'#22d3ee',2:'#60a5fa',3:'#f59e0b'};
        const MODE_KEY={1:'qS2',2:'qS1',3:'qP'};
        function idxXYZ(m0, ith, iph, nphi, ntheta){ return (((m0*ntheta)+ith)*nphi+iph)*3; }
        function autoScale(){
          if(!data) return 1;
          let mx=1e-9;
          for(const v of data.phaseXYZ) mx=Math.max(mx, Math.abs(v));
          for(const v of data.groupXYZ) mx=Math.max(mx, Math.abs(v));
          return mx;
        }
        function drawAxes(x, W, H, scale, R){
          const pts=[[R,0,0],[0,R,0],[0,0,R]], labs=['x','y','z'], cols=['#ef4444','#22c55e','#3b82f6'];
          for(let i=0;i<3;i++){
            const o=project([0,0,0],W,H,scale), t=project(pts[i],W,H,scale);
            x.strokeStyle=cols[i]; x.globalAlpha=0.85; x.lineWidth=1.4;
            x.beginPath(); x.moveTo(o.sx,o.sy); x.lineTo(t.sx,t.sy); x.stroke();
            x.fillStyle=cols[i]; x.globalAlpha=1; x.font='700 13px sans-serif';
            x.fillText(labs[i], t.sx+4, t.sy-4);
            const nt=project([-pts[i][0],-pts[i][1],-pts[i][2]],W,H,scale);
            x.globalAlpha=0.35; x.beginPath(); x.moveTo(o.sx,o.sy); x.lineTo(nt.sx,nt.sy); x.stroke();
          }
          x.globalAlpha=1;
        }
        function drawArrowHead2D(x, sx, sy, dx, dy, color){
          const len=Math.hypot(dx,dy); if(len<1e-6) return;
          const ux=dx/len, uy=dy/len, ah=6, aw=3.5;
          const bx=sx-ux*ah, by=sy-uy*ah;
          x.beginPath();
          x.moveTo(sx,sy);
          x.lineTo(bx-uy*aw, by+ux*aw);
          x.lineTo(bx+uy*aw, by-ux*aw);
          x.closePath();
          x.fillStyle=color; x.fill();
        }
        function drawCrystalAxes(x, W, H, scale, R){
          if(!data) return;
          const arrows=[[data.axisA, '#a78bfa', 0.55],[data.axisB,'#a78bfa',0.55],[data.axisSym,'#e879f9',1.0]];
          for(const [v,color,alpha] of arrows){
            const tip=[v[0]*R*0.9, v[1]*R*0.9, v[2]*R*0.9];
            const o=project([0,0,0],W,H,scale), t=project(tip,W,H,scale);
            x.globalAlpha=alpha; x.strokeStyle=color; x.lineWidth=2.2;
            x.beginPath(); x.moveTo(o.sx,o.sy); x.lineTo(t.sx,t.sy); x.stroke();
            drawArrowHead2D(x, t.sx, t.sy, t.sx-o.sx, t.sy-o.sy, color);
          }
          x.globalAlpha=1;
        }
        function drawPolarization(x, W, H, scale){
          if(!data) return;
          for(const a of data.arrows){
            if(!visible[MODE_KEY[a.mode]]) continue;
            const p=[a.p[0],a.p[1],a.p[2]];
            const r=Math.hypot(p[0],p[1],p[2]);
            const arrowLen=Math.max(r*0.16, 0.02);
            const tip=[p[0]+a.d[0]*arrowLen, p[1]+a.d[1]*arrowLen, p[2]+a.d[2]*arrowLen];
            const base=[p[0]-a.d[0]*arrowLen, p[1]-a.d[1]*arrowLen, p[2]-a.d[2]*arrowLen];
            const o=project(base,W,H,scale), t=project(tip,W,H,scale);
            x.globalAlpha=0.85; x.strokeStyle=MODE_COLOR[a.mode]; x.lineWidth=1.6;
            x.beginPath(); x.moveTo(o.sx,o.sy); x.lineTo(t.sx,t.sy); x.stroke();
            drawArrowHead2D(x, t.sx, t.sy, t.sx-o.sx, t.sy-o.sy, MODE_COLOR[a.mode]);
          }
          x.globalAlpha=1;
        }
        function buildSegments(arr, ntheta, nphi){
          const segs=[];
          for(let m0=0;m0<3;m0++){
            const key=MODE_KEY[m0+1];
            if(!visible[key]) continue;
            const pt=(ith,iph)=>{const b=idxXYZ(m0,ith,iph,nphi,ntheta); return [arr[b],arr[b+1],arr[b+2]];};
            for(let iph=0; iph<nphi; iph++){
              for(let ith=0; ith<ntheta-1; ith++){
                segs.push([pt(ith,iph), pt(ith+1,iph), m0+1]);
              }
            }
            for(let ith=0; ith<ntheta; ith++){
              for(let iph=0; iph<nphi; iph++){
                const iph2=(iph+1)%nphi;
                segs.push([pt(ith,iph), pt(ith,iph2), m0+1]);
              }
            }
          }
          return segs;
        }
        function draw(){
          const c=by('ani-3d');
          const [x,W,H]=setup(c);
          x.fillStyle='#000'; x.fillRect(0,0,W,H);
          if(!data){ return; }
          const R=autoScale();
          camD=R*3.5; // keep the camera well outside the largest sampled radius (any mode)
          const scale=Math.min(W,H)*0.34/R;
          const arr = surfaceMode==='phase' ? data.phaseXYZ : data.groupXYZ;
          const segs=buildSegments(arr, data.ntheta, data.nphi);
          const withDepth=segs.map(s=>{
            const o=project(s[0],W,H,scale), t=project(s[1],W,H,scale);
            const depth=(o.depth+t.depth)/2;
            return {o,t,depth,mode:s[2]};
          });
          withDepth.sort((a,b)=>a.depth-b.depth);
          let dmin=Infinity, dmax=-Infinity;
          for(const s of withDepth){ dmin=Math.min(dmin,s.depth); dmax=Math.max(dmax,s.depth); }
          const drange=Math.max(dmax-dmin, 1e-6);
          x.lineWidth=1.1;
          for(const s of withDepth){
            const t=(s.depth-dmin)/drange;
            const alpha=0.22+0.68*t;
            x.globalAlpha=alpha;
            x.strokeStyle=MODE_COLOR[s.mode];
            x.beginPath(); x.moveTo(s.o.sx,s.o.sy); x.lineTo(s.t.sx,s.t.sy); x.stroke();
          }
          x.globalAlpha=1;
          drawAxes(x,W,H,scale,R*1.15);
          if(showAxes) drawCrystalAxes(x,W,H,scale,R*1.15);
          if(showPol) drawPolarization(x,W,H,scale);
          by('ani-caption').textContent = (surfaceMode==='phase' ? 'Phase-velocity surface' : 'Group-velocity (ray/energy) surface') + ' — radius = velocity × '+state.time.toFixed(1)+' s. Drag to orbit.';
          by('ani-status').textContent = data.stable ? 'Vertical qP anisotropy: '+data.anisotropyPercent.toFixed(1)+'%' : 'Warning: this stiffness set is not mechanically stable in every sampled direction.';
        }
        let dragging=false, lastX=0, lastY=0;
        const canvas=by('ani-3d');
        canvas.addEventListener('mousedown', e=>{ dragging=true; lastX=e.clientX; lastY=e.clientY; });
        window.addEventListener('mouseup', ()=>{ dragging=false; });
        window.addEventListener('mousemove', e=>{
          if(!dragging) return;
          const dx=e.clientX-lastX, dy=e.clientY-lastY;
          lastX=e.clientX; lastY=e.clientY;
          cam.az += dx*0.008;
          cam.el = Math.max(-1.5, Math.min(1.5, cam.el + dy*0.008));
          draw();
        });
        canvas.addEventListener('touchstart', e=>{ if(e.touches.length===1){ dragging=true; lastX=e.touches[0].clientX; lastY=e.touches[0].clientY; } }, {passive:true});
        canvas.addEventListener('touchmove', e=>{
          if(!dragging || e.touches.length!==1) return;
          const dx=e.touches[0].clientX-lastX, dy=e.touches[0].clientY-lastY;
          lastX=e.touches[0].clientX; lastY=e.touches[0].clientY;
          cam.az += dx*0.008; cam.el=Math.max(-1.5,Math.min(1.5,cam.el+dy*0.008));
          draw();
        }, {passive:true});
        canvas.addEventListener('touchend', ()=>{ dragging=false; });
        function change(e){
          const key=ids[e.target.id]; if(!key) return;
          e.stopImmediatePropagation();
          state[key]=Number(e.target.value);
          labels(); emit();
        }
        root.addEventListener('input', change, true);
        root.addEventListener('change', change, true);
        by('ani-reset').onclick=()=>{
          Object.assign(state,{time:.8,density:4,A:272,C:160,L:60,N:50,F:60,plunge:0,trend:0});
          for(const [id,key] of Object.entries(ids)) by(id).value=state[key];
          labels(); emit();
        };
        by('ani-show-qP').addEventListener('change', e=>{ visible.qP=e.target.checked; draw(); });
        by('ani-show-qS1').addEventListener('change', e=>{ visible.qS1=e.target.checked; draw(); });
        by('ani-show-qS2').addEventListener('change', e=>{ visible.qS2=e.target.checked; draw(); });
        by('ani-show-axes').addEventListener('change', e=>{ showAxes=e.target.checked; draw(); });
        by('ani-show-pol').addEventListener('change', e=>{ showPol=e.target.checked; draw(); });
        root.querySelectorAll('#ani-seg-surface button').forEach(btn=>{
          btn.addEventListener('click', ()=>{
            surfaceMode=btn.dataset.v;
            root.querySelectorAll('#ani-seg-surface button').forEach(b=>b.classList.toggle('active', b===btn));
            draw();
          });
        });
        window.addEventListener('ani-results', e=>{ data = e.detail ? JSON.parse(e.detail) : null; draw(); });
        window.addEventListener('resize', draw);
        labels(); draw();
        })();
        </script></div>
        """)
    end
    const _ani_ready = true
end

# ╔═╡ fbe96b0a-84f8-4cae-905c-ba093b1b4340
begin
    # `VTIWavefrontInput` is defined at the end of the Appendix. The bare
    # reference makes the dependency explicit on a cold Pluto load.
    _ani_ready
    PlutoUI.WideCell(@bind _ani VTIWavefrontInput(); max_width=1400)
end

# ╔═╡ 109f4df4-aa1b-400c-b9d4-ff3604afa7d3
begin
    ani_safe = _ani isa AbstractDict ? _ani : Dict{String,Any}()
    ani_time = clamp(Float64(get(ani_safe, "time", 0.8)), 0.2, 1.2)
    ani_density = clamp(Float64(get(ani_safe, "density", 4.0)), 2.5, 5.0)
    ani_A = clamp(Float64(get(ani_safe, "A", 272.0)), 50.0, 300.0)
    ani_C = clamp(Float64(get(ani_safe, "C", 160.0)), 50.0, 300.0)
    ani_L = clamp(Float64(get(ani_safe, "L", 60.0)), 20.0, 140.0)
    ani_N = clamp(Float64(get(ani_safe, "N", 50.0)), 20.0, 140.0)
    ani_F = clamp(Float64(get(ani_safe, "F", 60.0)), -100.0, 150.0)
    ani_plunge = clamp(Float64(get(ani_safe, "plunge", 0.0)), 0.0, 90.0)
    ani_trend = clamp(Float64(get(ani_safe, "trend", 0.0)), 0.0, 360.0)
end

# ╔═╡ b99cb15b-b51b-4bd6-ac88-d5f6b99bf10b
let
    payload = vti_wavefront_payload_3d(ani_A, ani_C, ani_L, ani_N, ani_F, ani_density, ani_time, ani_plunge, ani_trend)
    message = ani_push_message(payload)
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('ani-results', {detail: $(repr(message))}));
    </script>""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Einsum = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Einsum = "~0.4.1"
PlutoTeachingTools = "~0.4.6"
PlutoUI = "~0.7.73"
SymbolicUtils = "~3.32.0"
Symbolics = "~6.57.0"
Unitful = "~1.25.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "3d0c311a717a73fee64aac659d8f1176df9b45bd"

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

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

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

[[deps.Einsum]]
deps = ["Compat"]
git-tree-sha1 = "4a6b3eee0161c89700b6c1949feae8b851da5494"
uuid = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
version = "0.4.1"

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
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

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
git-tree-sha1 = "d922b4d80d1e12c658da7785e754f4796cc1d60d"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.36"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "dacc8be63916b078b592806acd13bb5e5137d7e9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3faff84e6f97a7f18e0dd24373daa229fd358db5"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.73"

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
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLLogging", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "7614a1b881317b6800a8c66eb1180c6ea5b986f3"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.124.0"

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
git-tree-sha1 = "5a026f5549ad167cda34c67b62f8d3dc55754da3"
uuid = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
version = "1.3.1"

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
git-tree-sha1 = "8206e177903a41519145f577cb7f3793f3b7c960"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.57.0"

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

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "83360bda12f61c250835830cc40b64f487cc2230"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.25.1"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

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

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═b48fc46c-47b2-4090-83fb-710b1974c2a2
# ╠═6e9dc83f-5d4f-4f21-a495-d43bc64f6041
# ╟─17942206-4ea3-11ed-30fa-d38c9c47f282
# ╟─fbe96b0a-84f8-4cae-905c-ba093b1b4340
# ╟─109f4df4-aa1b-400c-b9d4-ff3604afa7d3
# ╟─b99cb15b-b51b-4bd6-ac88-d5f6b99bf10b
# ╟─e28c94fd-f04d-442c-8b19-be84d565d5c7
# ╟─99e819ab-a88c-41a0-9bbf-27d49fe19557
# ╟─56efaf7d-2b8b-4540-b0e8-b722edd1d8f8
# ╟─3bfab7cf-8b35-408d-8d25-ce6a37a0f19e
# ╟─d2027f44-361b-4d00-9efb-86fe2c7e9b68
# ╟─7d5bcebe-3f1c-485a-a094-78eb3b976352
# ╟─bba0ec5e-32bd-4ab1-84d1-76e18c88dc5a
# ╟─d000b26c-e67a-4e20-a2b8-b417b596ccd8
# ╟─8a6982d7-9d05-4d06-bcba-33cd5b79034a
# ╟─1724771c-d1ef-4c1c-93dc-e889a74d80e9
# ╟─846457f5-d684-470f-bba5-144b41f8dd14
# ╠═4b2d3084-e981-4ab2-bf40-fea4335082ad
# ╟─70ebdd20-3a02-4493-a0d6-718fa5b51675
# ╠═4f7b0a05-004a-417d-81d2-e4b1176909d6
# ╠═1fa4ecc5-8c5d-49f6-a023-bb2c4b4d5c27
# ╟─1de3540c-6e07-4e4e-b411-7a5667b5c742
# ╟─f5f78c1e-9062-4a5b-8bc3-4efb0f1883ca
# ╠═5c3cea52-3154-4aa3-824c-5be50a302be5
# ╟─88fb688e-d768-4895-8fe9-822f5db335b0
# ╠═b64d6d57-c38c-4f54-b86f-0d254d27c3ed
# ╟─12a42fbb-bfb7-4852-a638-2dec237c032d
# ╠═0ea0a70c-72b2-47da-902d-00cd6279fbb5
# ╠═eb482986-f6b0-4900-9e4e-441351691e12
# ╠═2d4ab726-4257-4fc4-aac0-2cc90ce76e53
# ╠═f2c3eeca-4e03-4360-b14b-4fd98e24d530
# ╠═4ac951fc-0054-4119-9296-fa8feacdd4c2
# ╟─0028f23a-0a92-42d8-9077-cc5fd585cce6
# ╠═8aa7f639-3e76-4db4-bd44-37b78f7f3781
# ╠═c6eefa83-676f-464b-85f1-4c7442100cbf
# ╟─183cc6c5-a4ae-47c6-8c48-2c8ff2e1cbc6
# ╠═4af3c346-f657-4b59-9a79-ee9ff0b79bbd
# ╠═6bbea2f3-da33-4abe-bf71-e6211275ce41
# ╠═24553b91-78cb-4cde-b0d6-8dbbd7e671d6
# ╟─f940ed1a-1f89-4f7b-bb41-05615e2351a0
# ╠═54657fec-e8e7-4040-8108-a29ae8df1c24
# ╠═108996f0-8980-4aff-8334-5f9f9f6f09d2
# ╠═5b1f23d8-33d1-4d25-92e3-ca590fc0118d
# ╠═7f0b2df9-6b33-44cb-a1c0-d4842cb666da
# ╠═ffeb806f-875f-4bc8-a15f-3d8db4252e96
# ╠═32fdd4ec-afef-4d03-8081-71a8201c3fcf
# ╠═c9cc3ef9-54ed-41bd-9c0b-c4c38b8c613b
# ╠═05c9f91b-06ee-4cf9-9cba-ae40aad8bd12
# ╠═f022c469-6ab4-4b2b-aee7-70e4a84d5a5a
# ╠═7d38e878-0136-46c6-b265-85b7dc85153e
# ╠═6c964379-2169-4567-9505-745198609d49
# ╠═0b2cf81d-98c3-42f9-a688-7069358003d6
# ╠═dbbc72ec-64a8-4566-b044-44fbbd1d4cbf
# ╠═12d0f4cd-ffe1-44e2-80e5-25dc5ab0dfd5
# ╟─2a2fc196-10a9-4c97-b1c1-e2a68e21119e
# ╠═77f96830-e81b-438d-b5ae-4ac59cb031df
# ╠═27f97bab-bef0-4705-bc22-8e9b087444b2
# ╟─366ee997-e2fc-4b48-bb6e-33b666ea38f5
# ╠═f6675c52-740c-4098-98d6-9c9038b4f6f4
# ╠═97da0c08-b971-4104-88fe-aa1d1f2c68b0
# ╠═cbdb003e-09a2-419b-8811-7d0cd290968d
# ╠═edfc88be-4bef-4478-b154-ec920782c295
# ╠═b22a4fcd-1f0c-410f-95b8-f67671d3aa63
# ╟─b5997b61-2ca2-4c2f-acd6-16c121f0aa13
# ╟─2bdba346-7a64-4de6-9b71-7034b9d68986
# ╠═78449122-6762-4191-9a8d-b1f8bee69bd8
# ╠═3bbbfa0d-c8ce-44f5-8546-16c9a1742e23
# ╠═747ca67e-60d7-44a4-bac1-78cfa9b20cc7
# ╟─4a66af73-6b29-475e-9c71-8e62ee69802a
# ╟─ebf4953d-abc5-4530-af51-8b68096a0113
# ╠═5f12d3e9-4d5c-434c-bdef-bcc74b79705e
# ╟─0989b0f0-06c6-4353-9933-ec713563f0c6
# ╠═aa111111-1111-4111-8111-111111111111
# ╠═aa333333-3333-4333-8333-333333333333
# ╠═aa222222-2222-4222-8222-222222222222
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
