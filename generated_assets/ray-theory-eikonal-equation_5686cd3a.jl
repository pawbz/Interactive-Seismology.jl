### A Pluto.jl notebook ###
# v0.20.13

#> [frontmatter]
#> title = "Ray Theory and The Eikonal Equation"
#> tags = ["raytheory"]
#> layout = "layout.jlhtml"
#> description = "This interactive notebook explores the Eikonal equation under the high-frequency approximation and its applications in tracing seismic rays in heterogeneous Earth models."

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

# ╔═╡ 22e38218-34cf-11ed-1808-97f785a5c673
begin
    using Symbolics
    using SymbolicUtils
    using LinearAlgebra
    using PlutoTeachingTools
    using PlutoPlotly
    using PlutoTest
    using Interpolations
    using Statistics
    using ProgressLogging
    using Latexify
    using PlutoUI
    using ImageFiltering
    using ForwardDiff
    using Einsum
    using Eikonal
    using DrWatson: savename
end;

# ╔═╡ 5053a4a4-312c-4a33-9c4f-79eb7bda2019
TableOfContents()

# ╔═╡ f715731e-7d18-423f-8ddf-75ae6b084e2c
ChooseDisplayMode()

# ╔═╡ d2dcd687-3623-433d-b591-cc8c2b8403eb
md"""
# Ray Theory and The Eikonal Equation
Seismic ray theory is a fundamental tool for understanding how high-frequency seismic waves propagate through a smoothly varying Earth medium. It is essential for understanding the imaging of the Earth's interior using seismic data.
This interactive notebook explores the Eikonal equation under the high-frequency approximation and its applications in tracing seismic rays in heterogeneous Earth models. The derivation of the governing equations of ray theoretical methods and their numerical solutions are presented, with a focus on two-dimensional scenarios that can be extended to 3D. An exciting feature of this notebook is that users can input any slowness function of
and leverage Julia's automatic differentiation capability to compute the gradient of this function, which is crucial for ray tracing.



[Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: Pawan Bharadwaj, Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ d2447315-f975-4447-ab92-5d5e267eaac5
md"""
Enter the name of the slowness-generating function 
$(@bind slowness_fn_name confirm(TextField(default="slowness_prograde")))
"""

# ╔═╡ dd079922-edcb-4a62-a0bc-c1e86c961dfb
md"""
Velocity increases linearly with depth. To use this function type `slowness_linear` above.
"""

# ╔═╡ b6447c75-4205-4c51-8cdc-552cdd841354
function slowness_zlinear(z, x)
    return inv(2000.0 + abs(z) * 20.0)
end

# ╔═╡ d3a75387-b9df-4fd2-b414-3ee662af813b
function slowness_gaussian(z, x)
    return inv(2000.0 + 4000.0 * exp(-(x - 250)^2 / 1e4) * exp(-(z - 50)^2 / 1e4))
end

# ╔═╡ a10d01b8-8ccd-4eff-b138-2b4eefb9d4da
function slowness_custom(z, x)
    # you can write a custom function, e.g.,
    return inv(2000.0 + abs(z) + abs(x))
end

# ╔═╡ d41b77fc-32de-4ae6-848f-8bf55772a8f0
function slowness_homogeneous(z, x)
    return inv(2000.0)
end

# ╔═╡ f332175d-a0b2-4c1e-8c80-fc48e751fdde
function slowness_LVZ(z, x)
    velocity = if z < 30.0
        2000.0 + z * 20.0
    elseif 30.0 <= z < 40.0
        2000.0
    else
        2000.0 + z * 20.0
    end
    return inv(velocity)
end

# ╔═╡ 7b2a5f24-273a-41e7-83af-cfa9ea4a6f1f
function slowness_prograde(z, x)
    velocity = if z <= 50.0
        2000 + z * 10.0
    elseif 50.0 < z <= 60.0
        -6000.0 + z * 170.0
    elseif z > 60.0
        3100.0 + z * 20.0
    end
    return inv(velocity)
end

# ╔═╡ 15d8f9d0-aa5f-4c62-868a-2486127e5800
md"""
## Ansatz 
We shall assume a solution of the form 
"""

# ╔═╡ 32a89a32-de52-481e-a017-1ae10c2c4bc4
@syms U(t)

# ╔═╡ 8d1ff0cb-78da-481c-b588-b664bf84828a
@syms ı

# ╔═╡ 7ae51402-9b74-462d-9603-2344a9766dd5
md"...with variables"

# ╔═╡ b949b16b-aad7-43ce-9ddf-769e92cd2eb3
@syms x::Real z::Real t::Real ω::Real y::Real

# ╔═╡ b6b2d21e-ebf7-4359-8b11-2ba8f4a8a4a5
@syms β(x, z) ρ(x, z) μ(x, z)

# ╔═╡ b3aa15cf-6baa-4209-9352-cebcb84dc3c5
@syms T(x, z) A(x, z) # travel-time and amplitude A(x,z)

# ╔═╡ 946772f8-3229-41cd-9a56-feae400ad11b
ϕ = A(x, z) * exp(ı * ω * (t - T(x, z))) # trail solution

# ╔═╡ ff0a4ccf-7518-4662-9203-ff639cb2bce2
md"...and the following operators."

# ╔═╡ f234fd11-93ba-4cbd-ab43-19997408be00
begin
    Dx = Differential(x)
    Dz = Differential(z)
    Dy = Differential(y)
    Dt = Differential(t)
end

# ╔═╡ 8f398848-ad04-453d-9762-d76e2ecc6c7d
D = [Dx, Dy, Dz]

# ╔═╡ 6b82596f-3d05-427d-83cc-f893f533458a
L(ϕ) = ρ(x, z) * Dt(Dt(ϕ)) - (Dx(μ(x, z) * Dx(ϕ)) + Dz(μ(x, z) * Dz(ϕ)))

# ╔═╡ 2a131f55-a6fd-489b-b4d5-205f19de5ce4
md"Here, T denotes the travel time and $A$ denotes the amplitude. Now, let's check the terms after substituting the solution $\phi$ into the scalar homogeneous wave-equation."

# ╔═╡ 0f457a47-4342-4582-b5e6-fe6fad263f2e
args_Lϕ = arguments(simplify(substitute(simplify(expand_derivatives(L(ϕ))), [ı^2 => -1])))

# ╔═╡ 28b95ced-5252-403b-8ff2-c92ab743103d
md"""
## High-frequency Approximation
This approximation ignores the terms that are sufficiently small in the far-field, i.e. when the high frequencies are considered. In this case, we shall divide the arguments obtained after substituting ϕ₂ in the 2-D wave equation with $ω^2$, and ignore $\frac{1}{ω^2}$ terms when $\omega$ is sufficiently large.
"""

# ╔═╡ 8ef1a5e5-aaad-41e0-bcc7-299822554714
map(simplify, args_Lϕ ./ ω^2)

# ╔═╡ 43d09ae6-5e46-4961-8c0a-d9f8ad0addba
md"""
This results in the 2D Eikonal equation.
"""

# ╔═╡ 78dcaf9d-84fc-49d6-b6ba-9c6f7b2bf688
begin
    eik(T) = Dx(T(x, z))^2 + Dz(T(x, z))^2 - 1 / β(x, z)^2
    eik(T) ~ 0
end

# ╔═╡ 5cb46243-8562-44ed-af21-d62a993f97c4
md"... and"

# ╔═╡ 9c30d437-507a-4300-921d-0c36a0911247
begin
    eikA(A) = 2 * dot(Symbolics.gradient(A(x, z), [x, z]), Symbolics.gradient(T(x, z), [x, z])) + A(x, z) * (Dx(Dx(T(x, z))) + Dz(Dz(T(x, z))))
    eikA(A) ~ 0
end

# ╔═╡ aba0fd7b-3cb6-49f3-9257-eb610d7a47dc
tip(md"""
Note the similarity of the above Eikonal equation to the dispersion relation we derived using a plane-wave solution of the homogeneous scalar wave equation.  The magnitude of the wavenumber vector $\vec{k}$ should be equal to $ω^2/α^2$, which also constrains the magnitude of the slowness vector $\vec{s}$ to $1/α^2$.""")

# ╔═╡ 2e1cec4b-11ff-4811-862f-1dea2c52efa6
tip(md"""... when solving the Eikonal equation, we think of a *local plane wave* at $(x,z)$ with corresponding slowness vector $\vec{s}(x,z)$. Take a moment to realize that $(Dx(T)) gives the x-component of the slowness vector.
""")

# ╔═╡ b38c287e-5fa2-4442-94cb-7ea89f34af3d
md"""
## Wavefronts
The function $T(x,y,z) = const.$ defines surfaces called wavefronts. For a medium with homogeneous elastic properties, these wavefronts are spherical. They can have arbitrary shapes in the presence of medium heterogeneities.
"""

# ╔═╡ a311274d-322e-4b87-964d-1d8db379c218
md"""
## Raypath Tracing
Lines perpendicular to the wavefronts i.e. $T(x,y,z) = const.$ surfaces are termed rays. In other words, rays are parallel to the gradient of the travel-time function 
$∇T(x,y,z)$.

If $\hat{s}$ denoted the direction along $\nabla T$, then  
```math
\nabla T(x,y,z) = \vec{s}(x,z) = |\vec{s}(x,z)|\hat{s}(x,z), \qquad\qquad (1)
```
where $|\vec{s}|=\frac{1}{α}$ is the local slowness.

We shall start this section by declaring the necessary symbols.
"""

# ╔═╡ 48472040-7d98-437c-8d80-97313f674446
begin
    @syms s(x::Real, z::Real)::Real # slowness
    @syms sx(x::Real, z::Real)::Real # slowness vector; x component
    @syms sz(x::Real, z::Real)::Real # ...
end

# ╔═╡ 663efee6-b542-4f06-b468-13eb61622dd5
svec = [sx(x, z), sz(x, z)] # slowness vector

# ╔═╡ 252804ea-2774-41ad-a832-59a997e3daab
# The Jacobian Matrix
J = Symbolics.jacobian(svec, [x, z])

# ╔═╡ 2093a743-0dd5-4766-8fee-a8607d70a675
md"""
The journey along a ray begins at the source, where we choose the direction of the outgoing plane wavefronts $\hat{s}_0$. Denoting the source position with $\vec{p}_0$, we simply move along $\hat{s}_0$ by an incremental length $dl$, and the updated position is given by:
```math
\vec{p}_1 = \vec{p}_0 + \hat{s}_0 \,dl.   \qquad\qquad(2)
```
In order to trace the ray path further, we need to estimate the change in $s$ as we moved from $\vec{p}_0$ to $\vec{p}_1$.  Note that the plane wave that we are riding gets transformed as it propagates in the medium due to changes in the slowness $s=|\vec{s}|$. 

In order to estimate this change, we shall first consider a Jacobian Matrix 
`J`=$J.

Then, we will simplify the derivative of $\vec{s}$ with respect to the length along the ray path ($l$), which is `dsdl`=`J`$\hat{s}$, using the Eikonal equation, to get `dsdl`=$\nabla\,s$.

Finally, the slowness vector at the new position $\vec{p}_1$ is given by 
```math
\vec{s}_1 = \vec{s}_0 + \nabla\,s\,dl. \qquad\qquad (3)
```
This notebook solves equations (2) and (3) numerically to trace the ray path in 2-D heterogeneous media.
"""

# ╔═╡ a2559f67-481b-4139-ac44-653b35b71f46
dsdl1 = J * (svec ./ s(x, z))

# ╔═╡ 7098ee62-8c43-4bfa-83be-472aed997975
md"""
Using equation (1), we shall now substitute `sx` and `sz` in `dsdl` and simplify.
"""

# ╔═╡ e35ea643-e180-4743-a46a-38090b397071
dsdl2 = broadcast(dsdl1) do ⋅
    simplify(substitute(⋅, [sx(x, z) => Dx(T(x, z)), sz(x, z) => Dz(T(x, z))]))
end

# ╔═╡ fd6d6a1b-b038-40be-af5b-864d693ab32c
r1 = @acrule (Dx(Dx(~T)) * Dx(~T) + ~B) / ~A => (1 / 2 * (Dx(Dx(~T) * Dx(~T))) + ~B) / ~A

# ╔═╡ 7e1d125a-cc07-4713-8e21-f5d49ce41797
r2 = @acrule (~B + Dz(Dx(~T)) * Dz(~T)) / ~A => (~B + 1 / 2 * (Dx(Dz(~T) * Dz(~T)))) / ~A

# ╔═╡ 1f9e802b-5523-4000-9486-60d77034f808
r3 = @acrule (Dz(Dz(~T)) * Dz(~T) + ~B) / ~A => (1 / 2 * (Dz(Dz(~T) * Dz(~T))) + ~B) / ~A

# ╔═╡ 80d2c246-2c18-41d2-b744-ea40d87da87f
r4 = @acrule (~B + Dx(Dz(~T)) * Dx(~T)) / ~A => (~B + 1 / 2 * (Dz(Dx(~T) * Dx(~T)))) / ~A

# ╔═╡ 3ffbf6e6-58d9-4b45-80b7-af1accad6b34
dsdl4 = [r1(r2(dsdl2[1])), r4(r3(dsdl2[2]))]

# ╔═╡ ebb5334f-8619-4e2c-b329-20a4818da3cc
dsdl3 = [r1(r2(dsdl2[1])), r4(r3(dsdl2[2]))]

# ╔═╡ 6c400f08-0262-4ab2-adea-283a43afbc7f
r5 = @acrule (~A * Dx(~B) + ~A * Dx(~C)) / ~D => (~A * Dx(~B + ~C)) / ~D

# ╔═╡ 24eaebcb-0020-4916-85af-fe1a3a80f8f6
r6 = @acrule (~A * Dz(~B) + ~A * Dz(~C)) / ~D => (~A * Dz(~B + ~C)) / ~D

# ╔═╡ 8042d4f4-89f6-41c5-9fed-38ccc648dd61
md"""
We can now finally use the Eikonal equation to derive equation (3).
"""

# ╔═╡ 123ef679-307b-4043-9318-96c91fe0ff18
dsdl = expand_derivatives.(substitute.(dsdl4, Dx(T(x, z)) * Dx(T(x, z)) + Dz(T(x, z)) * Dz(T(x, z)) => s(x, z) * s(x, z)))

# ╔═╡ bbd33fc8-e9b0-418a-a1f3-10015d8dec6f
tip(md"From this derivation, it is obvious that `dsdl` determines how the horizontal and vertical components of the slowness vector change along the raypath. For example, if $(Dx(s(x,z))) is zero in a layered Earth medium, the horizontal slowness remains constant!")

# ╔═╡ 9fa624d8-013a-4f4f-b440-a349a023dc47
@test iszero(dsdl .- Symbolics.gradient(s(x, z), [x, z]))

# ╔═╡ 7d2d4e9c-e440-472e-9900-8d3266bdeb89
md"""
## Amplitudes
"""

# ╔═╡ 1fbd1c3d-c84c-4052-ae3f-714d87a1d6e6
eikA(A) ~ 0

# ╔═╡ aec4dafd-8686-4d6d-b80a-c22490d5c429
begin
    eikA_arg(A) = simplify(expand_derivatives(substitute(eikA(A), [Dx(T(x, z)) => sx(x, z), Dz(T(x, z)) => sz(x, z)])))
    eikA_arg(A)
end

# ╔═╡ 8dfecdba-6350-485a-ae9b-27105948b3fd
md"""
Lets assume a solution of the form `A(x, z) = exp(Ã(x, z))`.
"""

# ╔═╡ fc7beea4-0a99-4624-8407-f4b00c9e61b2
@syms Ã(x::Real, z::Real)

# ╔═╡ 012f25c3-4ce2-4cb0-99e7-3ede25005427
eikÃ1 = expand_derivatives(substitute(eikA_arg(A), A(x, z) => exp(Ã(x, z))))

# ╔═╡ 0be23d53-8c21-4702-bfa1-31f420c73c1f
eikÃ2 = simplify(eikÃ1 / exp(Ã(x, z)))

# ╔═╡ c89f84e0-e5fd-4f1e-bc7b-e220e0123f72
eikÃ = sum(arguments(eikÃ2)[1:2]) / 2 / s(x, z) ~ -sum(arguments(eikÃ2)[3:4]) / s(x, z) / 2

# ╔═╡ 5cfb3c7d-c182-4431-b99e-7964b07255f7
md"""
The LHS of the above equation corresponds to the projection of the gradient of A along the ray path. We can now integrate along the ray path to obtain the amplitude.
"""

# ╔═╡ 7e94debf-3f99-4eb9-8950-0c50462edbd1
@syms ∫ₚₐₜₕ(x)

# ╔═╡ d9f0ccd2-b902-4e78-95a9-3078c01354bf
exp(∫ₚₐₜₕ(-eikÃ.lhs))

# ╔═╡ e31506b1-2fe9-44a9-9e79-88f1464fae90
md"""
Intuition: Divergence of the slowness field.
"""

# ╔═╡ 97307d52-d30c-46f9-8d55-9a0626879360
md"""
## Appendix
"""

# ╔═╡ fcef78b7-7c31-449f-b620-251249f83eb6
md"""
### UI
"""

# ╔═╡ 0c2a78b6-e859-4085-a5ad-1f742e5c70ac
function layered_medium_input(n) # n is the number of layers

    return PlutoUI.combine() do Child

        inputs = [
            md""" Layer $(string(i)): $(
Child(string("L",i), Slider(1000:10000, default=1000+i*1000, show_value=true))
            )"""
            for i in 1:n
        ]

        md"""
        #### Layered Earth
        These sliders are only active if you input `slowness_layered` below.
        Adjust the seismic velocities ∈ [1, 10] km/s.
        $(inputs)
                       		
        """
    end
end

# ╔═╡ 633f5b9a-77da-48e5-b6b3-00a5bc3e42d4
md"""
### Medium
"""

# ╔═╡ 690a6780-5169-4377-a7f1-795d89362c08
begin
    dx = 0.5 # choosing dx other than 1.0 needs adjustment for Eikonal package
    zgrid = range(0.0, stop=150.0, step=dx)
    xgrid = range(0.0, stop=500.0, step=dx)
end;

# ╔═╡ b4685924-854c-4058-af0a-bd7937f669b6
function source_input()
    return PlutoUI.combine() do Child
        dinput = [
            md""" θ∈[1, 360]$(
              Child("θ", RangeSlider(range(0, stop=360, step=5); left=0, right=60, show_value=true))
   )"""
            for i in 1:1
        ]

        linput = [
            md""" $(x): $(
            	Child(string(x, "pos"), Slider(grid, default=10, show_value=true))
            )"""
            for (x, grid) in zip(["x", "z"], [xgrid, zgrid])
        ]

        md"""
#### Source Parameters
Adjust the position from which the ray originates.
        $(linput...)
The direction of the outgoing slowness vector.
        $(dinput...)
        """
    end
end

# ╔═╡ a66ab3cd-c293-45ce-9e58-36b93712dbf2
TwoColumn(md"""$(@bind lmedium confirm(layered_medium_input(5)))""",
    md"""$(@bind source confirm(source_input()))
    """)

# ╔═╡ 2d79a52e-a11b-4841-b588-f1eddb1be8d5
begin
    layers = [getindex(lmedium, k) for k in Symbol.(filter(x -> occursin("L", x), string.(keys(lmedium))))]
    # convert the input velocity values to a slowness field 
    zlayer = collect(range(zgrid[1], stop=zgrid[end], length=length(layers)))
    xlayer = [xgrid[1], xgrid[end]]
    slayer = inv.(collect(layers))
    slayer_itp = extrapolate(interpolate((zlayer, xlayer), hcat(slayer, slayer), Gridded(Interpolations.Linear())), Flat())
end;

# ╔═╡ 412d0a5d-d4df-4c37-9fe3-90441bfcb32a
slowness_layered(z, x) = slayer_itp[z, x]

# ╔═╡ b8d6ef4d-d567-4128-89e3-b529ca6a3e3b
slowness(z, x) = eval(Symbol(slowness_fn_name))(z, x)

# ╔═╡ 307721e0-00fd-4b7e-9e17-2be8460d76b3
slowness_grid = [slowness(z, x) for z in zgrid, x in xgrid];

# ╔═╡ 10f25f52-b85c-47a4-89ea-68d94dd2912b
begin
    slowness_x(z, x) = ForwardDiff.derivative(x -> slowness(z, x), x)
    slowness_z(z, x) = ForwardDiff.derivative(z -> slowness(z, x), z)
end

# ╔═╡ ef77b591-6c37-46c9-a419-277367e48c68
md"### Trace"

# ╔═╡ 1670fc06-6e3f-4d0b-9202-f3cbac21386d
pa = (; slowness, slowness_x, slowness_z, xgrid, zgrid);

# ╔═╡ 25eda9c1-524c-46d1-a31d-90e1793bb8f1
ds = 1.0;

# ╔═╡ e0619921-389e-4351-8799-02431574a01d
function get_raypaths(N, ds, Xsource, initial_slowness_vectors, pa)
    (; slowness, slowness_x, slowness_z, xgrid, zgrid) = pa

    Xrays = map(initial_slowness_vectors) do slowness_vector
        # keep the direction of S_init, but adjust the magnitude to match the slowness at the source
        slowness_vector = (slowness_vector ./ norm(slowness_vector)) .* norm([slowness(Xsource[1], Xsource[2])])

        Xraysave = Array{Any}(missing, 2, N)
        Xray = deepcopy(Xsource)
        S = deepcopy(slowness_vector)
        for i = 1:N
            Xs = view(Xraysave, :, i)
            copyto!(Xs, Xray)
            # equation 2
            Xray[1] = Xray[1] + (S[1] / slowness(Xray[1], Xray[2])) * ds
            Xray[2] = Xray[2] + (S[2] / slowness(Xray[1], Xray[2])) * ds
            # equation 3
            S[1] = S[1] + ds * slowness_z(Xray[1], Xray[2])
            S[2] = S[2] + ds * slowness_x(Xray[1], Xray[2])

            # exit, if the ray reaches the edge of the medium
            !(Xray[2] >= xgrid[1] && Xray[2] <= xgrid[end] && Xray[1] >= zgrid[1] && Xray[1] <= zgrid[end]) && break
        end
        return Xraysave
    end
    return Xrays
end

# ╔═╡ c05a5082-0175-4a24-9aeb-de26cb22e6c6
raypaths = get_raypaths(500, ds, [source.zpos, source.xpos], [[sind(θ), cosd(θ)] for θ in source.θ], pa);

# ╔═╡ edb59274-e19a-4c64-8bab-2b3258a9a6fa
md"Now lets integrate slowness along the raypaths to get traveltimes."

# ╔═╡ 1de3052f-9543-4025-b92b-b75a73effc3d
function get_raytraveltimes(N, Xrays, pa)
    (; slowness, slowness_x, slowness_z, xgrid, zgrid) = pa
    traveltimes = map(Xrays) do Xray
        Xray = filter(x -> !any(ismissing.(x)), eachslice(Xray, dims=2))
        traveltime = Array{Any}(missing, N)
        traveltime[1] = 0.0
        for i = 2:length(Xray)
            # this is supposed to be constant, as we are moving a fixed length along the ray during ray tracing, but anyway...
            distance = sqrt(sum(abs2.(Xray[i] - Xray[i-1])))
            traveltime[i] = traveltime[i-1] + (distance * slowness(Xray[i-1][1], Xray[i-1][2]))
        end
        return traveltime
    end
end

# ╔═╡ 52a6ac1a-11ef-4993-9fb5-b0a51de58aae
raytraveltimes = get_raytraveltimes(500, raypaths, pa);

# ╔═╡ 613623ed-324d-4eae-8a93-888f205b83d4
begin
    # findout if the ray intersects the surface (z=0), if yes, return the index along the ray
    Iz0 = map(raypaths, raytraveltimes) do raypath, raytt
        Xray = filter(x -> !any(ismissing.(x)), eachslice(raypath, dims=2))
        i = findlast(x -> (x[1] < 1), Xray)
    end
    # return traveltime, when ray intersects the surface
    rayXtraveltimes = map(raytraveltimes, Iz0) do raytt, i
        if (i === nothing)
            return missing
        else
            return raytt[i]
        end
    end
    # return distance to the point where ray intersects the surface
    rayX = map(raypaths, Iz0) do raypath, i
        if (i === nothing)
            return missing
        else
            return raypath[2, i]
        end
    end
end;

# ╔═╡ 58338184-fdf5-4a03-8900-650cdbe36c1c
md"### Eikonal"

# ╔═╡ 9bd17627-8819-4d03-a02d-968a16b6dd9b
begin
    fastsweep = FastSweeping(slowness_grid)
    init!(fastsweep, (argmin(abs.(source.zpos .- zgrid)), argmin(abs.(source.xpos .- xgrid))))
    sweep!(fastsweep, verbose=false)
end;

# ╔═╡ eb401f31-e74d-46d4-b348-aaae693e8c15
TT_grid = fastsweep.t * dx; # traveltime to each (x,z)

# ╔═╡ e4d79480-0caf-4d5d-a01b-1f28e5690519
function diff2_z(u, dx)
    u2 = zero(u)
    for ix in 1:size(u, 2)
        u2[2:end-1, ix] .= (u[1:end-2, ix] .- 2u[2:end-1, ix] .+ u[3:end, ix]) ./ dx^2
    end
    return u2
end

# ╔═╡ a34e9b61-b240-45e2-98a9-6f465f756afa
function diff2_x(u, dx)
    u2 = zero(u)
    for iz in 1:size(u, 1)
        u2[iz, 2:end-1] .= (u[iz, 1:end-2] .- 2u[iz, 2:end-1] .+ u[iz, 3:end]) ./ dx^2
    end
    return u2
end

# ╔═╡ f4bdc936-f4aa-4947-82b8-1424fcc99724
div_s_grid = diff2_x(TT_grid, dx) + diff2_z(TT_grid, dx);

# ╔═╡ d073f00e-57db-4c4b-8074-02f0668b9362
md"### Amplitudes"

# ╔═╡ b1705abe-56bb-42f0-9aa3-e998424a5662
function get_amplitudes(raypaths)
    amplitudes = map(raypaths) do raypath
        amps = Array{Any}(missing, size(raypath, 2))
        amp = 0.0 # initialize
        for (i, p) in enumerate(eachslice(raypath, dims=2))
            if (any(.!ismissing.(p)))
                ipz = argmin(abs.(p[1] .- zgrid))
                ipx = argmin(abs.(p[2] .- xgrid))

                amp = amp - (2.0 * div_s_grid[ipz, ipx] / slowness(p[1], p[2]))
                amps[i] = amp
            end
        end
        return exp.(amps)
    end
    return amplitudes
end

# ╔═╡ 46cca7b7-f81a-47ad-8810-0cb5c77b98f4
amplitudes = get_amplitudes(raypaths)

# ╔═╡ 0a76470f-ffe4-4ae8-8dd6-f6886ac77454
md"""
### Plots
"""

# ╔═╡ abbe4e4e-6f0d-4f23-ad72-2930118c1ffe
begin
    source
    # create a plot object `rayplot` that will be updated later
    rayplot = Plot(Layout(showlegend=false, yaxis_autorange="reversed", height=250, width=700, title=attr(font_size=12,),
        legend=attr(
            x=-0.6,
            y=0.0,), font=attr(
            size=10), yaxis=attr(scaleanchor="x"), Subplots(horizontal_spacing=0.3, rows=1, cols=1, subplot_titles=["2-D (x, z) Ray Tracing" ""])))
    # add velocity model to `rayplot`
    add_trace!(rayplot, heatmap(
            x=xgrid,
            y=zgrid,
            z=inv.(slowness_grid), colorscale="Cividis", colorbar_title="Velocity<br>(m/s)", colorbar_x=1), row=1, col=1)
end;

# ╔═╡ df7f0572-50cd-4a84-96ba-9c91cae9605d
# update `rayplot` by adding a raypath that was just traced
function update_rayplot!(rayplot, raypaths)
    map(raypaths, amplitudes) do raypath, amplitude
        add_trace!(rayplot, scatter(
                x=raypath[2, :],
                y=raypath[1, :],
                mode="markers", marker=attr(size=1, color=log.(amplitude), colorscale="Blackbody"), line=attr(color="black", width=1.5), name=savename(source),), row=1, col=1)
    end
    return nothing
end

# ╔═╡ d3f909d1-1843-4580-ae75-de1c461dd433
begin
    update_rayplot!(rayplot, raypaths)
    PlutoPlotly.plot(rayplot)
end

# ╔═╡ b9e26366-01cc-4031-99d1-daac67b8f39d
function plot_tt_contours()
    plot(contour(x=xgrid,
            y=zgrid, z=TT_grid, showlabels=true,
            colorscale="Viridis",
            contours=attr(
                coloring="heatmap",
                showlabels=true, # show labels on contours
                labelfont=attr( # label font properties
                    size=12,
                    color="white",
                )
            )), Layout(title="First Arrival Traveltimes", yaxis_autorange="reversed", height=250, width=700, colorbar=false,
            yaxis=attr(scaleanchor="x"),))
end

# ╔═╡ f9054319-1963-4bba-92c2-c5b39753f5b5
plot_tt_contours()

# ╔═╡ 79065ef9-15fc-4cfd-966a-9cbf3d1b4f25
function plot_traveltimes()
    trPlot = Plot(Layout(showlegend=false, margin=0.5, height=300, width=680, title=attr(font_size=12,), xaxis=attr(range=(0, maximum(xgrid))),
        legend=attr(
            x=-0.6,
            y=0.0,), font=attr(
            size=10), Subplots(horizontal_spacing=0.1, rows=1, cols=2, subplot_titles=["First Arrival Traveltime" "Raytracing Traveltime"], x_title="Distance", y_title="Traveltime", shared_yaxes=true)))
    add_trace!(trPlot, scatter(x=xgrid, y=TT_grid[1, :], mode="markers", marker=attr(size=2)), col=1, row=1)
    add_trace!(trPlot, scatter(x=rayX, y=rayXtraveltimes, mode="markers+lines", marker=attr(size=5)), col=2, row=1)
    return plot(trPlot)
end

# ╔═╡ e69b3755-515b-4df7-ac7e-9c5745f5fc73
plot_traveltimes()

# ╔═╡ ee179fd5-c5c0-42f5-8bb8-b6a4acabb70c
md"## TODO"

# ╔═╡ c012fbb8-d696-403d-8752-61773c4f6d86
md"""
- Amplitudes!
- Prove Fermat Principle using Variational Calculus
"""

# ╔═╡ e4aaf1ea-f2f0-4083-bd4c-1069d98ee298
md"""
## Fermat's principle

Consider two points $A$ and $B$. We would like to show that the ray function minimizes the total travel time from $A$ to $B$.

In variational calculus, we are trying to find a function that minimizes a functional.
The travel-time is given by the integral 
```math
\mathbb{T} = \int_A^B s(x)\,dl
```
Notice that this equation of $\mathbb{T}$ is nonlinear in the slowness field
as the integration path depends on the velocity.

In order to consider all other paths, let's consider some
$\eta(x)$ that is an arbitrary path and scale it by a factor $\epsilon$. $\epsilon\eta(x)$ is the variation of $s(x)$. $\eta$ is twice differentiable. 


```math
\bar{s}(x) = s(x) + \epsilon\eta(x)
```
```math
\eta(A) = \eta(B) = 0
```

Now we are going to set the derivative of $I$ w.r.t. to $\epsilon$ be zero.

```math
\frac{d\mathbb{T}}{d\epsilon}|_{\epsilon=0}
```

```math
\int_A^B \eta(x)
```
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DrWatson = "634d3b9d-ee7a-5ddf-bec9-22491ea816e1"
Eikonal = "a6aab1ba-8f88-4217-b671-4d0788596809"
Einsum = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
ImageFiltering = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
DrWatson = "~2.12.7"
Eikonal = "~0.1.1"
Einsum = "~0.4.1"
ForwardDiff = "~0.10.36"
ImageFiltering = "~0.7.8"
Interpolations = "~0.14.7"
Latexify = "~0.16.1"
PlutoPlotly = "~0.3.9"
PlutoTeachingTools = "~0.2.13"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.52"
ProgressLogging = "~0.1.4"
SymbolicUtils = "~1.4.0"
Symbolics = "~5.5.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "f8e1f3a7be0479cb7ce7c773a522056fe006136d"

[[deps.ADTypes]]
git-tree-sha1 = "016833eb52ba2d6bea9fcb50ca295980e728ee24"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.7"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "c3c29bf6363b3ac3e421dc8b2ba8e33bdacbd245"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.32.5"

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
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c5aeb516a84459e0318a02507d2261edad97eb75"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.1"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "d8b0439d2be438a5f2cd68ec158fe08a7b2595b7"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "5a97e67919535d6841172016c9530fd69494e5ec"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.6"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

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
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

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

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d8a9c0b6ac2d9081bf76324b39c78ca3ce4f0c98"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.6"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a692f5e257d332de1e554e4566a4e5a8a72de2b2"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.4"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

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

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "0b4190661e8a4e51a842070e7dd4fae440ddb7f4"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.118"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "51b4b84d33ec5e0955b55ff4b748b99ce2c3faa9"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.7"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DrWatson]]
deps = ["Dates", "FileIO", "JLD2", "LibGit2", "MacroTools", "Pkg", "Random", "Requires", "Scratch", "UnPack"]
git-tree-sha1 = "d79f58f511d90e721496f1cdfd2ef74a313f226b"
uuid = "634d3b9d-ee7a-5ddf-bec9-22491ea816e1"
version = "2.12.7"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "30a1848c4f4fc35d1d4bbbd125650f6a11b5bc6c"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.7"

[[deps.Eikonal]]
deps = ["DataStructures", "Images", "LinearAlgebra", "PrecompileTools", "Printf"]
git-tree-sha1 = "ac89a6cf8c89a741448deb8692aaacba745ecee0"
uuid = "a6aab1ba-8f88-4217-b671-4d0788596809"
version = "0.1.1"

[[deps.Einsum]]
deps = ["Compat"]
git-tree-sha1 = "4a6b3eee0161c89700b6c1949feae8b851da5494"
uuid = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
version = "0.4.1"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "7de7c78d681078f027389e067864a8d53bd7c3c9"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+3"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43ba3d3c82c18d88471cfd2924931658838c9d8f"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+4"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1dc470db8b1131cfc7fb4c115de89fe391b9e780"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.12.0"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "ExprTools", "Logging", "MultivariatePolynomials", "Primes", "Random", "SIMD", "SnoopPrecompile"]
git-tree-sha1 = "44f595de4f6485ab5ba71fe257b5eadaa3cf161e"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.4.4"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6df9cd6ee79fc59feab33f63a1b3c9e95e2461d5"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.2"

[[deps.HistogramThresholding]]
deps = ["ImageBase", "LinearAlgebra", "MappedArrays"]
git-tree-sha1 = "7194dfbb2f8d945abdaf68fa9480a965d6661e69"
uuid = "2c695a8d-9458-5d45-9878-1b8a99cf7853"
version = "0.3.1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8e070b599339d622e9a081d17230d74a5c473293"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.17"

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

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageBinarization]]
deps = ["HistogramThresholding", "ImageCore", "LinearAlgebra", "Polynomials", "Reexport", "Statistics"]
git-tree-sha1 = "33485b4e40d1df46c806498c73ea32dc17475c59"
uuid = "cbc4b850-ae4b-5111-9e64-df94c024a13d"
version = "0.3.1"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageCorners]]
deps = ["ImageCore", "ImageFiltering", "PrecompileTools", "StaticArrays", "StatsBase"]
git-tree-sha1 = "24c52de051293745a9bad7d73497708954562b79"
uuid = "89d5987c-236e-4e32-acd0-25bd6bd87b70"
version = "0.1.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "33cb509839cc4011beb45bde2316e64344b0f92b"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.9"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "c5c5478ae8d944c63d6de961b19e6d3324812c35"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.4.0"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fa01c98985be12e5d75301c4527fff2c46fa3e0e"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "7.1.1+1"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.ImageMorphology]]
deps = ["DataStructures", "ImageCore", "LinearAlgebra", "LoopVectorization", "OffsetArrays", "Requires", "TiledIteration"]
git-tree-sha1 = "6f0a801136cb9c229aebea0df296cdcd471dbcd1"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.4.5"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "3db3bb9f7014e86f13692581fa2feb6460bdee7e"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.8.4"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "e0884bdf01bbbb111aea77c348368a86fb4b5ab6"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.10.1"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageBinarization", "ImageContrastAdjustment", "ImageCore", "ImageCorners", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "a49b96fd4a8d1a9a718dfd9cde34c154fc84fcd5"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.26.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "b842cbff3f44804a84fda409745cc8f04c029a20"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.6"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
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
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "Requires", "TranscodingStreams"]
git-tree-sha1 = "89e1e5c3d43078d42eed2306cab2a11b13e5c6ae"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.54"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "a434e811d10e7cbf4f0674285542e697dca605d0"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.42"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "d1f981fba6eb3ec393eede4821bca3f2b7592cd4"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.15.1"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd714447457c660382fe634710fb56eb255ee42e"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.6"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

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

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "8be878062e0ffa2c3f67bb58a595375eda5de80b"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.0+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "ff3b4b9d35de638936a525ecd36e86a8bb919d11"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "df37206100d39f79b3376afb6b9cee4970041c61"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.51.1+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

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

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "8084c25a250e00ae427a379a5b607e7aed96a2dd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.171"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "688d6d9e098109051ae33d126fcfc88c4ce4a021"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.1.0"

[[deps.MIMEs]]
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "e9650bea7f91c3397eb9ae6377343963a22bf5b8"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.8.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "8d39779e29f80aa6c071e7ac17101c6e31f075d7"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "491bdcdc943fcbc4c005900d7463c9f216aabf4c"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "cc0a5deefdb12ab3a096f00a6d42133af4560d71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "8a3271d8309285f4db73b4f662b1b290c715e85e"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.21"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "a414039192a155fb38c4599a60110f0018c6ec82"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.16.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "966b85253e959ea89c53a9abebbf2e964fbf593b"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.32"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

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
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Colors", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "90af5c9238c1b3b25421f1fdfffd1e8fca7a7133"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.20"

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

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Colors", "Dates", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "PackageExtensionCompat", "PlotlyBase", "PlutoUI", "Reexport"]
git-tree-sha1 = "9a77654cdb96e8c8a0f1e56a053235a739d453fe"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.3.9"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "555c272d20fc80a2658587fb9bbda60067b93b7c"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.19"

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

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "6c62ce45f268f3f958821a1e5192cf91c75ae89c"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.24"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "8b3fc30bc0390abdce15f8822c889f669baed73d"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "994cc27cdacca10e68feb291673ec3a76aa2fae9"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "b8a399e95663485820000f26b6a43c794e166a49"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.4"

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
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d7087c013e8a496ff396bae843b1e16d9a30ede8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.10"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "9bb80533cb9769933954ea4ffbecb3025a783198"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.7.2"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "456f610ca2fbd1c14f5fcf31c6bfadc55e7d66e0"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.43"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "ChainRulesCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "916b8a94c0d61fa5f7c5295649d3746afb866aff"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.98.1"

    [deps.SciMLBase.extensions]
    ZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "6149620767866d4b0f0f7028639b6e661b6a1e44"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.12"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

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

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "3e5f165e58b18204aed03158664c4982d691f454"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.5.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "f737d444cb0ad07e61b3c1bef8eb91203c321eff"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.2.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

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
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
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
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "2f3fa844bcd33e40d8c29de5ee8dded7a0a70422"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.4.0"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "ac7f8825d029b568f82dbf2cb49da9cebcadaffb"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.5.3"

    [deps.Symbolics.extensions]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

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
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

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

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "f21231b166166bebc73b99cea236071eb047525b"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.3"

[[deps.TiledIteration]]
deps = ["OffsetArrays", "StaticArrayInterface"]
git-tree-sha1 = "1176cc31e867217b06928e2f140c90bd1bc88283"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.5.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f57facfd1be61c42321765d3551b3df50f7e09f6"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.28"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

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

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "4ab62a49f1d8d9548a1c8d1a75e5f55cf196f64e"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.71"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "5f24e158cf4cee437052371455fe361f526da062"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.6"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "7d1671acbe47ac88e981868a078bd6b4e27c5191"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.42+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "9dafcee1d24c4f024e7edc92603cedba72118283"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+3"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e9216fdcd8514b7072b43653874fd688e4c6c003"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.12+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "89799ae67c17caa5b3b5a19b8469eeee474377db"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.5+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d7155fea91a4123ef59f42c4afb5ab3b4ca95058"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+3"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c57201109a9e4c0585b208bb408bc41d205ac4e9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.2+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "1a74296303b6524a0472a8cb12d3d87a78eb3612"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6dba04dbfb72ae3ebe5418ba33d087ba8aa8cb00"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.1+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "434b3de333c75fc446aa0d19fc394edafd07ab08"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.7"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "068dfe202b0a05b8332f1e8e6b4080684b9c7700"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.47+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "ccbb625a89ec6195856a50aa2b668a5c08712c94"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.4.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═5053a4a4-312c-4a33-9c4f-79eb7bda2019
# ╠═f715731e-7d18-423f-8ddf-75ae6b084e2c
# ╟─d2dcd687-3623-433d-b591-cc8c2b8403eb
# ╟─a66ab3cd-c293-45ce-9e58-36b93712dbf2
# ╟─d3f909d1-1843-4580-ae75-de1c461dd433
# ╟─f9054319-1963-4bba-92c2-c5b39753f5b5
# ╟─e69b3755-515b-4df7-ac7e-9c5745f5fc73
# ╟─d2447315-f975-4447-ab92-5d5e267eaac5
# ╟─dd079922-edcb-4a62-a0bc-c1e86c961dfb
# ╠═b6447c75-4205-4c51-8cdc-552cdd841354
# ╠═d3a75387-b9df-4fd2-b414-3ee662af813b
# ╠═a10d01b8-8ccd-4eff-b138-2b4eefb9d4da
# ╠═d41b77fc-32de-4ae6-848f-8bf55772a8f0
# ╠═f332175d-a0b2-4c1e-8c80-fc48e751fdde
# ╠═7b2a5f24-273a-41e7-83af-cfa9ea4a6f1f
# ╟─15d8f9d0-aa5f-4c62-868a-2486127e5800
# ╠═32a89a32-de52-481e-a017-1ae10c2c4bc4
# ╠═946772f8-3229-41cd-9a56-feae400ad11b
# ╠═8d1ff0cb-78da-481c-b588-b664bf84828a
# ╟─7ae51402-9b74-462d-9603-2344a9766dd5
# ╠═b949b16b-aad7-43ce-9ddf-769e92cd2eb3
# ╠═b6b2d21e-ebf7-4359-8b11-2ba8f4a8a4a5
# ╠═b3aa15cf-6baa-4209-9352-cebcb84dc3c5
# ╟─ff0a4ccf-7518-4662-9203-ff639cb2bce2
# ╠═f234fd11-93ba-4cbd-ab43-19997408be00
# ╠═8f398848-ad04-453d-9762-d76e2ecc6c7d
# ╠═6b82596f-3d05-427d-83cc-f893f533458a
# ╟─2a131f55-a6fd-489b-b4d5-205f19de5ce4
# ╠═0f457a47-4342-4582-b5e6-fe6fad263f2e
# ╟─28b95ced-5252-403b-8ff2-c92ab743103d
# ╠═8ef1a5e5-aaad-41e0-bcc7-299822554714
# ╟─43d09ae6-5e46-4961-8c0a-d9f8ad0addba
# ╠═78dcaf9d-84fc-49d6-b6ba-9c6f7b2bf688
# ╟─5cb46243-8562-44ed-af21-d62a993f97c4
# ╠═9c30d437-507a-4300-921d-0c36a0911247
# ╟─aba0fd7b-3cb6-49f3-9257-eb610d7a47dc
# ╟─2e1cec4b-11ff-4811-862f-1dea2c52efa6
# ╟─b38c287e-5fa2-4442-94cb-7ea89f34af3d
# ╟─a311274d-322e-4b87-964d-1d8db379c218
# ╠═48472040-7d98-437c-8d80-97313f674446
# ╠═663efee6-b542-4f06-b468-13eb61622dd5
# ╟─2093a743-0dd5-4766-8fee-a8607d70a675
# ╠═252804ea-2774-41ad-a832-59a997e3daab
# ╠═a2559f67-481b-4139-ac44-653b35b71f46
# ╟─7098ee62-8c43-4bfa-83be-472aed997975
# ╠═e35ea643-e180-4743-a46a-38090b397071
# ╠═fd6d6a1b-b038-40be-af5b-864d693ab32c
# ╠═7e1d125a-cc07-4713-8e21-f5d49ce41797
# ╠═1f9e802b-5523-4000-9486-60d77034f808
# ╠═80d2c246-2c18-41d2-b744-ea40d87da87f
# ╠═3ffbf6e6-58d9-4b45-80b7-af1accad6b34
# ╠═ebb5334f-8619-4e2c-b329-20a4818da3cc
# ╠═6c400f08-0262-4ab2-adea-283a43afbc7f
# ╠═24eaebcb-0020-4916-85af-fe1a3a80f8f6
# ╟─8042d4f4-89f6-41c5-9fed-38ccc648dd61
# ╠═123ef679-307b-4043-9318-96c91fe0ff18
# ╟─bbd33fc8-e9b0-418a-a1f3-10015d8dec6f
# ╠═9fa624d8-013a-4f4f-b440-a349a023dc47
# ╟─7d2d4e9c-e440-472e-9900-8d3266bdeb89
# ╠═1fbd1c3d-c84c-4052-ae3f-714d87a1d6e6
# ╠═aec4dafd-8686-4d6d-b80a-c22490d5c429
# ╟─8dfecdba-6350-485a-ae9b-27105948b3fd
# ╠═fc7beea4-0a99-4624-8407-f4b00c9e61b2
# ╠═012f25c3-4ce2-4cb0-99e7-3ede25005427
# ╠═0be23d53-8c21-4702-bfa1-31f420c73c1f
# ╠═c89f84e0-e5fd-4f1e-bc7b-e220e0123f72
# ╟─5cfb3c7d-c182-4431-b99e-7964b07255f7
# ╠═7e94debf-3f99-4eb9-8950-0c50462edbd1
# ╠═d9f0ccd2-b902-4e78-95a9-3078c01354bf
# ╟─e31506b1-2fe9-44a9-9e79-88f1464fae90
# ╟─97307d52-d30c-46f9-8d55-9a0626879360
# ╠═22e38218-34cf-11ed-1808-97f785a5c673
# ╟─fcef78b7-7c31-449f-b620-251249f83eb6
# ╠═0c2a78b6-e859-4085-a5ad-1f742e5c70ac
# ╠═b4685924-854c-4058-af0a-bd7937f669b6
# ╟─633f5b9a-77da-48e5-b6b3-00a5bc3e42d4
# ╠═690a6780-5169-4377-a7f1-795d89362c08
# ╠═2d79a52e-a11b-4841-b588-f1eddb1be8d5
# ╠═412d0a5d-d4df-4c37-9fe3-90441bfcb32a
# ╠═b8d6ef4d-d567-4128-89e3-b529ca6a3e3b
# ╠═307721e0-00fd-4b7e-9e17-2be8460d76b3
# ╠═10f25f52-b85c-47a4-89ea-68d94dd2912b
# ╟─ef77b591-6c37-46c9-a419-277367e48c68
# ╠═1670fc06-6e3f-4d0b-9202-f3cbac21386d
# ╠═25eda9c1-524c-46d1-a31d-90e1793bb8f1
# ╠═c05a5082-0175-4a24-9aeb-de26cb22e6c6
# ╠═e0619921-389e-4351-8799-02431574a01d
# ╟─edb59274-e19a-4c64-8bab-2b3258a9a6fa
# ╠═52a6ac1a-11ef-4993-9fb5-b0a51de58aae
# ╠═1de3052f-9543-4025-b92b-b75a73effc3d
# ╠═613623ed-324d-4eae-8a93-888f205b83d4
# ╟─58338184-fdf5-4a03-8900-650cdbe36c1c
# ╠═9bd17627-8819-4d03-a02d-968a16b6dd9b
# ╠═eb401f31-e74d-46d4-b348-aaae693e8c15
# ╠═e4d79480-0caf-4d5d-a01b-1f28e5690519
# ╠═a34e9b61-b240-45e2-98a9-6f465f756afa
# ╠═f4bdc936-f4aa-4947-82b8-1424fcc99724
# ╟─d073f00e-57db-4c4b-8074-02f0668b9362
# ╠═b1705abe-56bb-42f0-9aa3-e998424a5662
# ╠═46cca7b7-f81a-47ad-8810-0cb5c77b98f4
# ╟─0a76470f-ffe4-4ae8-8dd6-f6886ac77454
# ╠═abbe4e4e-6f0d-4f23-ad72-2930118c1ffe
# ╠═df7f0572-50cd-4a84-96ba-9c91cae9605d
# ╠═b9e26366-01cc-4031-99d1-daac67b8f39d
# ╠═79065ef9-15fc-4cfd-966a-9cbf3d1b4f25
# ╟─ee179fd5-c5c0-42f5-8bb8-b6a4acabb70c
# ╟─c012fbb8-d696-403d-8752-61773c4f6d86
# ╟─e4aaf1ea-f2f0-4083-bd4c-1069d98ee298
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
