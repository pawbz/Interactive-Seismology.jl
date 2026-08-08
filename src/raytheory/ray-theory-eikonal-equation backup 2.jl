### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ a66ab3cd-c293-45ce-9e58-36b93712dbf2
TwoColumn(md"""$(@bind lmedium confirm(layered_medium_input(5)))""",
    md"""$(@bind source confirm(source_input()))
    """)

# ╔═╡ d3f909d1-1843-4580-ae75-de1c461dd433
begin
    update_rayplot!(rayplot, raypaths)
    PlutoPlotly.plot(rayplot)
end

# ╔═╡ f9054319-1963-4bba-92c2-c5b39753f5b5
plot_tt_contours()

# ╔═╡ e69b3755-515b-4df7-ac7e-9c5745f5fc73
plot_traveltimes()

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

# ╔═╡ 946772f8-3229-41cd-9a56-feae400ad11b
ϕ = A(x, z) * exp(ı * ω * (t - T(x, z))) # trail solution

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

# ╔═╡ 252804ea-2774-41ad-a832-59a997e3daab
# The Jacobian Matrix
J = Symbolics.jacobian(svec, [x, z])

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

# ╔═╡ c05a5082-0175-4a24-9aeb-de26cb22e6c6
raypaths = get_raypaths(500, ds, [source.zpos, source.xpos], [[sind(θ), cosd(θ)] for θ in source.θ], pa);

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

# ╔═╡ edb59274-e19a-4c64-8bab-2b3258a9a6fa
md"Now lets integrate slowness along the raypaths to get traveltimes."

# ╔═╡ 52a6ac1a-11ef-4993-9fb5-b0a51de58aae
raytraveltimes = get_raytraveltimes(500, raypaths, pa);

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
