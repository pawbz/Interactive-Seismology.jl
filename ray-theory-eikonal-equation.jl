### A Pluto.jl notebook ###
# v0.19.27

#> [frontmatter]
#> title = "Ray Theory and The Eikonal Equation"
#> description = "This interactive notebook explores the Eikonal equation under the high-frequency approximation and its applications in tracing seismic rays in heterogeneous Earth models."

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
    using ForwardDiff
    using Einsum
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
$(@bind slowness_fn_name confirm(TextField(default="slowness_layered")))
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
    return inv(2000.0 + exp(-(x - 250)^2 / 1e5) * exp(-(z - 50)^2 / 1e5))
end

# ╔═╡ a10d01b8-8ccd-4eff-b138-2b4eefb9d4da
function slowness_custom(z, x)
	# you can write a custom function, e.g.,
    return inv(2000.0 + abs(z) + abs(x))
end

# ╔═╡ 15d8f9d0-aa5f-4c62-868a-2486127e5800
md"""
## Ansatz 
We shall assume a solution of the form 
"""

# ╔═╡ 7ae51402-9b74-462d-9603-2344a9766dd5
md"...with variables"

# ╔═╡ b949b16b-aad7-43ce-9ddf-769e92cd2eb3
@syms x::Real z::Real t::Real ω::Real α::Real y::Real

# ╔═╡ b3aa15cf-6baa-4209-9352-cebcb84dc3c5
@syms T(x, z) A(x, z) # travel-time and amplitude A(x,z)

# ╔═╡ 946772f8-3229-41cd-9a56-feae400ad11b
ϕ = A(x, z) * sin(ω * (t - T(x, z))) # trail solution

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
L(ϕ, α) = Dt(Dt(ϕ)) / α^2 - (Dx(Dx(ϕ)) + Dz(Dz(ϕ)))

# ╔═╡ 2a131f55-a6fd-489b-b4d5-205f19de5ce4
md"Here, T denotes the travel time and $A$ denotes the amplitude. Now, let's check the terms after substituting the solution $\phi$ into the scalar homogeneous wave-equation."

# ╔═╡ 0f457a47-4342-4582-b5e6-fe6fad263f2e
args_Lϕ = arguments(expand_derivatives(L(ϕ, α)))

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
    eik(T) = Dx(T(x, z))^2 + Dz(T(x, z))^2 - 1 / α^2
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
\hat{s}_1 = \hat{s}_0 + \nabla\,s\,dl. \qquad\qquad (3)
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

# ╔═╡ ebb5334f-8619-4e2c-b329-20a4818da3cc
dsdl3 = [r2(r1(dsdl2[1])), r4(r3(dsdl2[2]))]

# ╔═╡ 6c400f08-0262-4ab2-adea-283a43afbc7f
r5 = @acrule (~A * Dx(~B) + ~A * Dx(~C)) / ~D => (~A * Dx(~B + ~C)) / ~D

# ╔═╡ 24eaebcb-0020-4916-85af-fe1a3a80f8f6
r6 = @acrule (~A * Dz(~B) + ~A * Dz(~C)) / ~D => (~A * Dz(~B + ~C)) / ~D

# ╔═╡ b42d1cd0-6610-4d3f-8950-27b13f136368
dsdl4 = [r5(dsdl3[1]), r6(dsdl3[2])]

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
@syms Ã(x::Real, z::Real)::Real

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
exp(∫ₚₐₜₕ(eikÃ.rhs))

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
    zgrid = range(0, stop=100, length=512)
    xgrid = range(0, stop=500, length=512)
end

# ╔═╡ b4685924-854c-4058-af0a-bd7937f669b6
function source_input()
    return PlutoUI.combine() do Child
        dinput = [
            md""" θ∈[1, 360]$(
              Child("θ", Slider(range(1, stop=360, step=1), default=45, show_value=true))
   )"""
            for i in 1:1
        ]

        linput = [
            md""" $(x): $(
            	Child(string(x, "pos"), Slider(grid, default=50, show_value=true))
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
    slayer_itp = extrapolate(interpolate((zlayer, xlayer), hcat(slayer, slayer), Gridded(Linear())), Flat())
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

# ╔═╡ e0619921-389e-4351-8799-02431574a01d
function get_raypath(N, ds, Xinit, Sinit, pa)
    (; slowness, slowness_x, slowness_z, xgrid, zgrid) = pa
    # keep the direction of S_init, but adjust the magnitude to match the slowness at the source
    Sinit = (Sinit ./ norm(Sinit)) .* norm([slowness(Xinit[1], Xinit[2])])

    Xsave = Array{Any}(missing, 2, N)
    X = deepcopy(Xinit)
    S = deepcopy(Sinit)
    for i = 1:N
        Xs = view(Xsave, :, i)
        copyto!(Xs, X)
        # equation 2
        X[1] = X[1] + (S[1] / slowness(X[1], X[2])) * ds
        X[2] = X[2] + (S[2] / slowness(X[1], X[2])) * ds
        # equation 3
        S[1] = S[1] + ds * slowness_z(X[1], X[2])
        S[2] = S[2] + ds * slowness_x(X[1], X[2])

        # exit, if the ray reaches the edge of the medium
        (((X[2] - xgrid[1]) * (xgrid[end] - X[2]) * (X[1] - zgrid[1]) * (zgrid[end] - X[1])) < 0.0) && break
    end
    return Xsave
end

# ╔═╡ c05a5082-0175-4a24-9aeb-de26cb22e6c6
raypath = get_raypath(500, 1, [source.zpos, source.xpos], [sind(source.θ), cosd(source.θ)], pa);

# ╔═╡ 0a76470f-ffe4-4ae8-8dd6-f6886ac77454
md"""
### Plots
"""

# ╔═╡ abbe4e4e-6f0d-4f23-ad72-2930118c1ffe
begin
    # create a plot object `rayplot` that will be updated later
    rayplot = Plot(Layout(showlegend=true, yaxis_autorange="reversed", height=250, width=700, title=attr(font_size=12,), 
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
function update_rayplot!(rayplot, raypath)
    add_trace!(rayplot, scatter(
            x=raypath[2, :],
            y=raypath[1, :],
            mode="lines", line=attr(width=2), name=savename(source),), row=1, col=1)
    return nothing
end

# ╔═╡ d3f909d1-1843-4580-ae75-de1c461dd433
begin
    update_rayplot!(rayplot, raypath)
    PlutoPlotly.plot(rayplot)
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Einsum = "b7d42ee7-0b51-5a75-98ca-779d3107e4c0"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
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
Einsum = "~0.4.1"
ForwardDiff = "~0.10.36"
Interpolations = "~0.14.7"
Latexify = "~0.16.1"
PlutoPlotly = "~0.3.9"
PlutoTeachingTools = "~0.2.13"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.52"
ProgressLogging = "~0.1.4"
SymbolicUtils = "~1.2.0"
Symbolics = "~5.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "cf992dc0b9d815b5cc3f8d81c852565956e8d468"

[[deps.ADTypes]]
git-tree-sha1 = "a4c8e0f8c09d4aa708289c1a5fc23e2d1970017a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.1"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "41b45260ad00b666558c16f7cd9db0e8af408ea1"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.31.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

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

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "a1296f0fe01a4c3f9bf0dc2934efbf4416f5db31"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "d9a8f86737b665e15a9641ecbac64deef9ce6724"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.23.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

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
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fe2838a593b5f776e1597e086dcd47560d94e816"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.3"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "938fe2981db009f531b6332e31c58e9584a2f9bd"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.100"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

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

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "9b05a8bc04d7a9a9c7e9ee3b0fdf1584857b65dc"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.2"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "a20eaa3ad64254c61eeb5f230d9306e937405434"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.6.1"
weakdeps = ["SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
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

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random", "SIMD", "SnoopPrecompile"]
git-tree-sha1 = "17c81e8dd04cc8f70bccf8c50aa2d39f5ff2a3ec"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.4.2"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "8e59ea773deee525c99a8018409f64f19fb719e6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.7"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "81dc6aefcbe7421bd62cb6ca0e700779330acff8"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.25"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

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

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "f9978f23952b52b8d958b72f8b5368f84254dc02"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.1"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "5879579adbd9bc3017fe63c766371aace1a1d641"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.3.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "f9b1e033c2b1205cf30fd119f4e50881316c1923"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.1"
weakdeps = ["Requires", "TOML"]

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

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
git-tree-sha1 = "542de5acb35585afcf202a6d3361b430bc1c3fbd"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.13"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "4c9f306e5d6603ae203c2000dd460d81a5251489"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "7ed35fb5f831aaf09c2d7c8736d44667a1afdcb0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.7"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "1e597b93700fa4045d7189afa7c004e0584ea548"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.3"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0e270732477b9e551d884e6b07e23bb2ec947790"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.5"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "ChainRulesCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "a1e295a85fd6fae9a7b2072dd4577001210f512d"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.95.0"

    [deps.SciMLBase.extensions]
    ZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "65c2e6ced6f62ea796af251eb292a0e131a3613b"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.6"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "9cabadf6e7cd2349b6cf49f1915ad2028d65e881"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.2"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "9704a1d0ecbbf5d65700418dcf1bb4c680c790bb"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.2.0"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "f1d43a0dbb553890195e49fb599ea51d0e97a5ef"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.5.1"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

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

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "21c8fc7cd598ef49f11bc9e94871f5d7740e34b9"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.5"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═5053a4a4-312c-4a33-9c4f-79eb7bda2019
# ╠═f715731e-7d18-423f-8ddf-75ae6b084e2c
# ╟─d2dcd687-3623-433d-b591-cc8c2b8403eb
# ╟─a66ab3cd-c293-45ce-9e58-36b93712dbf2
# ╟─d3f909d1-1843-4580-ae75-de1c461dd433
# ╟─d2447315-f975-4447-ab92-5d5e267eaac5
# ╟─dd079922-edcb-4a62-a0bc-c1e86c961dfb
# ╠═b6447c75-4205-4c51-8cdc-552cdd841354
# ╠═d3a75387-b9df-4fd2-b414-3ee662af813b
# ╠═a10d01b8-8ccd-4eff-b138-2b4eefb9d4da
# ╟─15d8f9d0-aa5f-4c62-868a-2486127e5800
# ╠═946772f8-3229-41cd-9a56-feae400ad11b
# ╟─7ae51402-9b74-462d-9603-2344a9766dd5
# ╠═b949b16b-aad7-43ce-9ddf-769e92cd2eb3
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
# ╠═ebb5334f-8619-4e2c-b329-20a4818da3cc
# ╠═6c400f08-0262-4ab2-adea-283a43afbc7f
# ╠═24eaebcb-0020-4916-85af-fe1a3a80f8f6
# ╠═b42d1cd0-6610-4d3f-8950-27b13f136368
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
# ╠═c05a5082-0175-4a24-9aeb-de26cb22e6c6
# ╠═e0619921-389e-4351-8799-02431574a01d
# ╟─0a76470f-ffe4-4ae8-8dd6-f6886ac77454
# ╠═abbe4e4e-6f0d-4f23-ad72-2930118c1ffe
# ╠═df7f0572-50cd-4a84-96ba-9c91cae9605d
# ╟─ee179fd5-c5c0-42f5-8bb8-b6a4acabb70c
# ╟─c012fbb8-d696-403d-8752-61773c4f6d86
# ╟─e4aaf1ea-f2f0-4083-bd4c-1069d98ee298
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
