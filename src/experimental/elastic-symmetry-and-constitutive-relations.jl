### A Pluto.jl notebook ###
# v0.7.0

#> [frontmatter]
#> title = "Elastic Symmetry & Constitutive Relations"
#> layout = "layout.jlhtml"
#> tags = ["elasticity"]
#> description = "How the symmetry a material survives forces relations among its elastic constants, and reshapes how stiff it is in every direction."

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

# ╔═╡ 160890b8-404a-4c53-a4d7-06089d85c09d
begin
    using PlutoUI
    using LinearAlgebra
end

# ╔═╡ dff060e3-b550-4a33-bc48-465a10f4451c
md"""
# Elastic Symmetry & Constitutive Relations

A block of material resists being squeezed or sheared according to its **stiffness tensor** —
81 numbers relating every component of stress to every component of strain, collapsed by minor
and major symmetries into a 6×6 **Voigt matrix** with at most 21 independent entries. Every
crystal or rock fabric has *some* underlying symmetry — mirror planes, rotation axes — and that
symmetry is not a decoration on top of the elastic constants: it **forces** relations among them.

The one sentence this notebook is built around:

```math
\text{more symmetry operations survived} \;\Rightarrow\; \text{more forced relations in } C
\;\Rightarrow\; \text{fewer independent constants} \;\Rightarrow\; \text{a rounder stiffness surface.}
```

A fully generic (**triclinic**) block has all 21 constants independent and a lumpy, direction-
dependent stiffness. An **isotropic** block — the same under every possible rotation — has only 2.
In between sit **monoclinic** (13), **orthorhombic** (9), **hexagonal/transverse-isotropic** (5),
and **cubic** (3): each one degree more constrained than the last, each one a step toward the
sphere. The widget below lets you rotate a block, click a candidate symmetry operation, and watch
the Voigt matrix's cells collapse into families or zero in real time, then watch a
directional-stiffness surface round out to match.

Try the two moments the whole design is built to protect: **drag the block or the surface** — they
turn together, because the surface *is* that block's own stiffness, seen from outside. Then **click
a symmetry operation** — if the block would truly look the same afterward, watch which Voigt-matrix
cells light up to explain why.
"""

# ╔═╡ 4fcde378-5707-4ee5-9dce-9a9ec5ac78d2
PlutoUI.TableOfContents()

# ╔═╡ 729629f2-70b5-4a34-b190-0fa86e5f70e6
md"""
## Reading the Widget

**The block (left)** is not a literal crystal unit cell — its shape and any texture were chosen so
the shape's *own* symmetry group matches the class's *elastic* symmetry group exactly: a sphere
(isotropic) is unchanged by any rotation at all; a cube with the same dot pattern on every face
survives 90° and 120°-about-a-corner rotations; a circular cylinder survives spinning about its
axis; a brick with three different edge lengths survives every 180° turn and mirror but not a 90°
one (its footprint visibly isn't square); a sheared prism survives only the two operations that
respect its unique axis; the fully tilted block survives nothing but inversion. Every green/red
verdict a button gives is the same computation ([`is_symmetry_operation`](@ref)) that built the
shape in the first place — nothing about "which operations should work" is hardcoded per class.

**The matrix (middle)** colors cells by *why* they take the value they do: grey for zero, a shared
color for a family forced equal by some symmetry, a mottled cell for one *derived* from others by a
formula (hexagonal's `C66=(C11-C12)/2` is the only one that appears among the six classes here).
Clicking an operation lights up the cells *that specific operation* would force, regardless of
which class is currently selected — a property of the operation alone. Hovering a cell instead asks
the reverse question, "what would straining this way *feel* like?", and draws the answer on the
block; clicking a cell commits to that probe and actually deforms the block along it.

**The surface (right)** is the directional Young's modulus (or, toggled, the three body-wave
speeds) computed from the *same* matrix, on the *same* sampling grid, sharing the block's own
orientation. A sphere means the material resists every direction equally; the more the surface
bulges and dents, the more directions the material actually distinguishes. Toggling to wave speed
splits the two shear branches apart wherever the surface is anisotropic — the same physical effect
real shear-wave splitting (e.g. SKS phases crossing an anisotropic mantle) measures in the field.

A few honest limits: the cubic, monoclinic, and triclinic presets are illustrative synthetic
numbers, not a specific measured mineral (isotropic is parameterized directly; hexagonal and
orthorhombic reuse real values from `anisotropy.jl`'s VTI defaults and olivine matrix). The probe
deformation's size is a fixed display amount with no physical scale — it shows *which direction*
a strain acts along, not how much a real material would actually move under it.
"""

# ╔═╡ c7b1cf0e-0733-480f-9d46-8c96d3449e26
md"""
## Appendix
"""

# ╔═╡ c5e30d76-f237-421d-8174-35cad72344fc
md"""
### Voigt Notation and the Elastic Tensor

The Voigt matrix `C` is a compressed notation for the elastic tensor `c_ijkl`: pairs of tensor
indices `(i,j)` collapse to one Voigt index via `11→1, 22→2, 33→3, 23→4, 13→5, 12→6`. Going from
the tensor to the matrix is a plain substitution (stiffness needs no extra factor — that's a
compliance-only subtlety, see below); going back just reads the matrix at one representative
`(i,j)` pair per Voigt index, since a physical elastic tensor is already symmetric under
`(i,j)↔(j,i)`, `(k,l)↔(l,k)`, and `(i,j)↔(k,l)`.
"""

# ╔═╡ 3c8a53dd-00c3-4a79-b78e-be46713efc12
"""
    get_cijkl(C)

Expand a 6×6 Voigt stiffness matrix `C` into the full 3×3×3×3 elastic tensor `c_ijkl`, via the
standard index map `11→1, 22→2, 33→3, 23→4, 13→5, 12→6`. General for any symmetry class — reused
verbatim from `src/Planewave Propagation/anisotropy.jl:204-206`, where it was already validated
against Okada/Christoffel machinery independent of any particular symmetry class.
"""
function get_cijkl(C)
    idx(i, j) = i == j ? i : 9 - i - j
    [C[idx(i, j), idx(k, l)] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
end

# ╔═╡ b95dce52-53ec-4aae-84f1-a41ac1c45c0e
"""
    tensor_to_voigt(c)

The inverse of [`get_cijkl`](@ref): read a 3×3×3×3 elastic tensor back into its 6×6 Voigt matrix
by evaluating `c` at one representative `(i,j)` pair per Voigt index. Only well-defined for a
tensor that already carries the physical index symmetries (see [`check_tensor_symmetries`](@ref));
no such function existed anywhere else in this repo.
"""
function tensor_to_voigt(c)
    ij_rep = [(1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2)]
    C = zeros(Float64, 6, 6)
    for a in 1:6, b in 1:6
        (i, j) = ij_rep[a]
        (k, l) = ij_rep[b]
        C[a, b] = c[i, j, k, l]
    end
    C
end

# ╔═╡ 4479f375-688e-4272-9129-c1861f4fb0f9
"""
    check_tensor_symmetries(c; tol=1e-10)

Assert that a 3×3×3×3 array `c` has the minor symmetries (`c_ijkl = c_jikl = c_ijlk`) and major
symmetry (`c_ijkl = c_klij`) every physical elastic tensor must satisfy. Used below as a self-check
after every rotation, since a coding mistake in [`bond_rotate_tensor`](@ref) would most likely
break one of these first.
"""
function check_tensor_symmetries(c; tol=1e-10)
    ok = true
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        ok &= isapprox(c[i, j, k, l], c[j, i, k, l]; atol=tol)
        ok &= isapprox(c[i, j, k, l], c[i, j, l, k]; atol=tol)
        ok &= isapprox(c[i, j, k, l], c[k, l, i, j]; atol=tol)
    end
    ok
end

# ╔═╡ a3513ce2-b9af-42c8-a71a-9ed8afe47487
md"""
### Verifying Voigt ↔ Tensor
"""

# ╔═╡ 22c406f2-4061-43a0-8d22-0a934609c4d1
let
    n_ok = 0
    for trial in 1:20
        C = zeros(Float64, 6, 6)
        for a in 1:6, b in a:6
            C[a, b] = C[b, a] = randn()
        end
        @assert isapprox(tensor_to_voigt(get_cijkl(C)), C; atol=1e-12)
        n_ok += 1
    end
    md"""
    !!! correct "Self-check"
        `tensor_to_voigt(get_cijkl(C)) == C` exactly, for **$(n_ok)** random symmetric 6×6
        matrices ✓
    """
end

# ╔═╡ 5f4c619b-a736-4909-a487-1c32de318c05
md"""
### Rotating the Elastic Tensor: the Bond Transform

A material's response can't depend on which coordinate frame we happen to describe it in — so
rotating the *tensor* by `R` (`c'_ijkl = R_ia R_jb R_kc R_ld c_abcd`, the **Bond transform**) must
give the exact same physical description in the rotated frame. If `R` is a genuine symmetry of the
material, rotating `C` by `R` reproduces `C` exactly; that single equality test —
[`is_symmetry_operation`](@ref), below — is the only thing that ever decides whether a block "looks
the same" after an operation. No per-class table of which operations are "allowed" is hardcoded
anywhere in this notebook; every green/red judgment is this one computation.
"""

# ╔═╡ 84810524-9c74-40ce-a1f6-54a890960903
"""
    bond_rotate_tensor(c, R)

Rotate an elastic tensor by the orthogonal matrix `R`: `c'_ijkl = R_ia R_jb R_kc R_ld c_abcd`, the
Bond transform in full tensor (not Voigt-matrix) form. A brute-force quadruple sum over 3×3×3×3
outputs each needing a 3×3×3×3 inner contraction (6561 multiply-adds total) — trivially fast at
this size, so there's no need for a tensor-contraction package.
"""
function bond_rotate_tensor(c, R)
    cp = zeros(Float64, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        s = 0.0
        for a in 1:3, b in 1:3, cc in 1:3, d in 1:3
            s += R[i, a] * R[j, b] * R[k, cc] * R[l, d] * c[a, b, cc, d]
        end
        cp[i, j, k, l] = s
    end
    cp
end

# ╔═╡ f9d818c3-0022-4ae0-9e87-5548dc303905
"""
    bond_transform_voigt(C, R)

[`bond_rotate_tensor`](@ref) followed by [`tensor_to_voigt`](@ref) — the Bond transform taking a
6×6 Voigt matrix straight to its rotated 6×6 Voigt matrix, via [`get_cijkl`](@ref) to enter tensor
form. This one function is the engine behind everything else in this notebook: symmetry testing,
matrix-structure derivation, and class presets all reduce to calling it.
"""
bond_transform_voigt(C, R) = tensor_to_voigt(bond_rotate_tensor(get_cijkl(C), R))

# ╔═╡ 5ceab72b-2c3e-42ba-bc16-d7d84c633b41
"""
    symmetry_generator(op::Symbol; theta=nothing)

The fixed catalog of 10 candidate symmetry operations tested throughout this notebook: 90°/180°
rotations about the coordinate axes, mirrors across the coordinate planes, inversion, a 120°
rotation about `[1,1,1]` (cyclic permutation of x,y,z — the 3-fold axis a cube survives), and a
continuous rotation about z at angle `theta` (needed for hexagonal/isotropic symmetry, which no
*single* discrete rotation captures). Every discrete operation here is an exact signed permutation
matrix; `contz` is the only one that needs a `theta` argument.
"""
function symmetry_generator(op::Symbol; theta=nothing)
    if op == :rot90z
        [0.0 -1 0; 1 0 0; 0 0 1]
    elseif op == :rot180z
        [-1.0 0 0; 0 -1 0; 0 0 1]
    elseif op == :rot180x
        [1.0 0 0; 0 -1 0; 0 0 -1]
    elseif op == :rot180y
        [-1.0 0 0; 0 1 0; 0 0 -1]
    elseif op == :mirror_xy
        [1.0 0 0; 0 1 0; 0 0 -1]
    elseif op == :mirror_xz
        [1.0 0 0; 0 -1 0; 0 0 1]
    elseif op == :mirror_yz
        [-1.0 0 0; 0 1 0; 0 0 1]
    elseif op == :inversion
        -Matrix{Float64}(I, 3, 3)
    elseif op == :rot120_111
        [0.0 0 1; 1 0 0; 0 1 0]
    elseif op == :contz
        θ = theta === nothing ? 0.0 : theta
        [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
    else
        error("unknown symmetry operation $op")
    end
end

# ╔═╡ 3a98ea1b-100a-4ce8-a205-916ed6be814e
const ALL_OPS = [:rot90z, :rot180z, :rot180x, :rot180y, :mirror_xy, :mirror_xz, :mirror_yz, :inversion, :rot120_111, :contz]

# ╔═╡ 5c84645b-c28d-46b5-9a08-108943e6368c
md"""
### Verifying the Bond Transform
"""

# ╔═╡ bae63dd3-86e2-4a29-99fb-e6acb4162039
let
    discrete_ops = filter(!=(:contz), ALL_OPS)
    for op in discrete_ops
        R = symmetry_generator(op)
        @assert isapprox(R' * R, I; atol=1e-12) "not orthogonal: $op"
    end

    C = zeros(Float64, 6, 6)
    for a in 1:6, b in a:6
        C[a, b] = C[b, a] = randn()
    end
    c = get_cijkl(C)
    for op in discrete_ops
        R = symmetry_generator(op)
        cp = bond_rotate_tensor(c, R)
        @assert check_tensor_symmetries(cp) "rotated tensor lost index symmetries: $op"
    end

    md"""
    !!! correct "Self-check"
        Every catalog operation is an orthogonal matrix ✓ · rotating a random elastic tensor by
        each one preserves the required minor/major index symmetries ✓
    """
end

# ╔═╡ 1812ef55-3447-4d25-a352-c570e34b218e
md"""
### Deriving Matrix Structure from Symmetry

This is the part that turns "click a symmetry operation" into "these Voigt cells collapse" —
**derived**, not memorized. Trying to recall e.g. monoclinic's exact zero pattern from memory is
exactly the kind of thing that's easy to get subtly wrong; instead this notebook *computes* it.

The Bond transform is **linear** in the 21 independent entries of `C`. Build the 21×21 matrix `L`
representing that linear map for a given rotation `R` by evaluating `bond_transform_voigt` on each
of the 21 one-hot symmetric basis matrices (a single independent Voigt entry set to 1, its
symmetric partner too, everything else 0 — a *one-sided* 1 would silently break
[`tensor_to_voigt`](@ref)'s single-representative-entry assumption). Requiring invariance under `R`
means every physical `C` satisfies `(L - I)·vec(C) = 0`; row-reducing `(L-I)` and reading off each
resulting row classifies the constraint it represents:

- **one nonzero term → forced zero** (zero literally falls out as "forced equal to its own
  negative" — no special case needed).
- **two terms, coefficients `+1,-1` → a pure family** (two cells forced equal — feeds a union-find
  grouping for later color-coding).
- **three or more terms, or non-±1 coefficients → derived** (e.g. hexagonal's
  `C66 = (C11-C12)/2` below) — a distinct category, never folded into the family bucket.

Every *discrete* catalog operation is a signed permutation of x,y,z, so its `L` has exactly one
nonzero `±1` entry per row (self-checked below, not assumed). The **continuous** z-rotation
(`contz`) is different: for generic `θ` it is not a signed permutation, and it forces genuine
derived relations. Sampling it at 5 generic, mutually incommensurate angles and stacking their
`(L(θ)-I)` blocks before reducing is what correctly captures that (agreement at enough generic
angles forces the underlying trig-polynomial identity to hold for *every* angle).
"""

# ╔═╡ c2f0bd3c-1e5d-400a-819e-d202c8b99620
"""
    VOIGT_PAIRS

The 21 independent Voigt-matrix entries `(α,β)` with `β ≥ α`, in a fixed order shared by every
function below that builds or reads a length-21 vector of independent constants.
"""
const VOIGT_PAIRS = [(a, b) for a in 1:6 for b in a:6]

# ╔═╡ 6274320f-24f8-42b2-b72a-26cfe1a37876
"""
    one_hot_C(a, b)

The symmetric 6×6 matrix with `C[a,b]=C[b,a]=1` and every other entry 0 — one basis vector of the
21-dimensional space of independent Voigt entries, used to build the linear map `L` in
[`build_L`](@ref).
"""
function one_hot_C(a, b)
    C = zeros(Float64, 6, 6)
    C[a, b] = 1.0
    C[b, a] = 1.0
    C
end

# ╔═╡ d36b4b9d-811d-476f-8394-b08fd8c097ce
"""
    build_L(R)

The 21×21 matrix representing the (linear) Bond transform under rotation `R`, acting on the 21
independent Voigt entries in [`VOIGT_PAIRS`](@ref) order: column `m` is `bond_transform_voigt`
applied to the `m`-th one-hot basis matrix, read back out at every other Voigt pair.
"""
function build_L(R)
    n = length(VOIGT_PAIRS)
    L = zeros(Float64, n, n)
    for (col, (a, b)) in enumerate(VOIGT_PAIRS)
        Cp = bond_transform_voigt(one_hot_C(a, b), R)
        for (row, (p, q)) in enumerate(VOIGT_PAIRS)
            L[row, col] = Cp[p, q]
        end
    end
    L
end

# ╔═╡ 29dcbb36-8d47-442c-b3af-9676d4274c36
"""
    rref_tol(A, tol=1e-8)

Reduced row-echelon form of `A` via Gaussian elimination with partial pivoting and a zero cutoff
`tol` (no RREF routine exists in this repo's dependencies). Used to classify the constraints a
symmetry operation imposes on the 21 independent Voigt entries (see
[`classify_constraints`](@ref)).
"""
function rref_tol(A::AbstractMatrix{<:Real}, tol::Float64=1e-8)
    M = Float64.(copy(A))
    nrows, ncols = size(M)
    pivot_row = 1
    for col in 1:ncols
        pivot_row > nrows && break
        p = argmax(abs.(M[pivot_row:nrows, col])) + pivot_row - 1
        abs(M[p, col]) < tol && continue
        M[[pivot_row, p], :] = M[[p, pivot_row], :]
        M[pivot_row, :] ./= M[pivot_row, col]
        for r in 1:nrows
            r == pivot_row && continue
            M[r, :] .-= M[r, col] .* M[pivot_row, :]
        end
        pivot_row += 1
    end
    M[abs.(M) .< tol] .= 0.0
    M
end

# ╔═╡ 0d3beab8-5a9c-46cd-8c01-149d5d34c00f
md"""
### Verifying `rref_tol`
"""

# ╔═╡ 8843c9fd-85ff-4e01-a6fe-1bcf5b8088a2
let
    # hand-worked: row 3 = row 1 - row 2, so rank is 2 and the reduced rows should read off as
    # x - z = 0 and y + z = 0 by direct inspection of [1 1 0; 0 1 1; 1 0 -1]
    A = [1.0 1 0; 0 1 1; 1 0 -1]
    R = rref_tol(A)
    nz_rows = count(r -> any(abs.(R[r, :]) .> 1e-8), 1:3)
    @assert nz_rows == 2
    @assert isapprox(R[1, :], [1.0, 0.0, -1.0]; atol=1e-8)
    @assert isapprox(R[2, :], [0.0, 1.0, 1.0]; atol=1e-8)
    md"""
    !!! correct "Self-check"
        `rref_tol` on a hand-worked 3×3 example (row 3 = row 1 − row 2) correctly finds rank 2 and
        reduces to exactly the two independent relations `x=z`, `y=-z` ✓
    """
end

# ╔═╡ 37455fad-da49-4de4-87b7-ba03acd8da5d
"""
    OpConstraints

The result of classifying one symmetry operation's constraints on the 21 independent Voigt
entries: `zeros` (indices into [`VOIGT_PAIRS`](@ref) forced to zero), `families` (groups of indices
forced equal), and `derived` (an index forced to a linear combination of others, each entry
`(index, [(other_index, coefficient), ...])`).
"""
struct OpConstraints
    zeros::Vector{Int}
    families::Vector{Vector{Int}}
    derived::Vector{Tuple{Int,Vector{Tuple{Int,Float64}}}}
end

# ╔═╡ af2515da-bd17-4a9d-b0a0-3eef9661f67f
"""
    classify_constraints(stacked_LmI; tol=1e-6)

Row-reduce a stacked `(L(R)-I)` system (one operation, or several stacked together) and classify
each resulting nonzero row into [`OpConstraints`](@ref)'s three categories, per the rule explained
above: 1 term → zero, 2 terms `{+1,-1}` → family, 3+ terms (or non-±1 coefficients) → derived.
"""
function classify_constraints(stacked_LmI::AbstractMatrix{Float64}; tol=1e-6)
    R = rref_tol(stacked_LmI, tol)
    zeros_ = Int[]
    families = Vector{Vector{Int}}()
    derived = Vector{Tuple{Int,Vector{Tuple{Int,Float64}}}}()
    for row in 1:size(R, 1)
        nz = findall(x -> abs(x) > tol, R[row, :])
        isempty(nz) && continue
        if length(nz) == 1
            push!(zeros_, nz[1])
        elseif length(nz) == 2 && isapprox(R[row, nz[1]], 1.0; atol=1e-6) && isapprox(R[row, nz[2]], -1.0; atol=1e-6)
            a, b = nz[1], nz[2]
            merged = false
            for fam in families
                if a in fam || b in fam
                    push!(fam, a); push!(fam, b); unique!(fam)
                    merged = true
                    break
                end
            end
            merged || push!(families, [a, b])
        else
            pivot = nz[1]
            terms = [(nz[k], -R[row, nz[k]]) for k in 2:length(nz)]
            push!(derived, (pivot, terms))
        end
    end
    OpConstraints(sort(unique(zeros_)), families, derived)
end

# ╔═╡ a64f8607-3c0c-4c81-bdf1-bf6d4b567118
"""
    stacked_constraint_system(ops)

Build and vertically stack `(L(R)-I)` for every operation in `ops`, sampling `contz` at 5 generic
angles (see the derivation note above). Shared by [`operation_constraints`](@ref),
[`class_constraints`](@ref), and [`project_onto_class`](@ref) so they can never silently disagree
about which angles or which stacking order to use.
"""
function stacked_constraint_system(ops)
    n = length(VOIGT_PAIRS)
    In = Matrix{Float64}(I, n, n)
    blocks = Matrix{Float64}[]
    for op in ops
        if op == :contz
            for θ in (17.0, 53.0, 101.0, 149.0, 197.0) .* (pi / 180)
                push!(blocks, build_L(symmetry_generator(:contz; theta=θ)) .- In)
            end
        else
            push!(blocks, build_L(symmetry_generator(op)) .- In)
        end
    end
    isempty(blocks) ? zeros(Float64, 0, n) : vcat(blocks...)
end

# ╔═╡ 9945b72f-9595-47fc-bf2e-2dcc9eda8451
"""
    operation_constraints(op)

The [`OpConstraints`](@ref) a *single* symmetry operation `op` imposes on a fully generic (21
independent entries) starting matrix — a property of the operation alone, independent of which
symmetry class is under consideration.
"""
operation_constraints(op) = classify_constraints(stacked_constraint_system([op]))

# ╔═╡ a030eb8c-2050-4640-9bf8-39e34a95c736
"""
    class_constraints(ops)

Combine several symmetry operations' constraints (a class's full generator set) by stacking their
`(L-I)` systems before reducing. Returns `(constraints::OpConstraints, independent_count)` — the
count is this notebook's master self-check: it must exactly equal 2/3/5/9/13/21 for the six classes
below, or the proposed generator set is wrong and needs revisiting (never the target number itself).
"""
function class_constraints(ops::Vector{Symbol})
    n = length(VOIGT_PAIRS)
    stacked = stacked_constraint_system(ops)
    R = rref_tol(stacked)
    rank = count(r -> any(abs.(R[r, :]) .> 1e-6), 1:size(R, 1))
    (constraints=classify_constraints(stacked), independent_count=n - rank)
end

# ╔═╡ 8908d7fe-0b41-49e0-a033-995607e41749
"""
    is_symmetry_operation(C, R; tol=1e-6)

The single function deciding whether rotation/reflection `R` is a symmetry of the concrete matrix
`C`: does the Bond transform reproduce `C` exactly? This is what drives every "does the block look
the same?" judgment in the eventual widget, for every symmetry class alike — never a per-class
hardcoded list of which operations are "allowed".
"""
is_symmetry_operation(C, R; tol=1e-6) = isapprox(bond_transform_voigt(C, R), C; atol=tol)

# ╔═╡ f2b0023e-d0a6-4890-a2ab-85e4ec7e4112
"""
    project_onto_class(C_generic, ops)

Orthogonally project an arbitrary starting matrix `C_generic` onto the subspace of matrices that
exactly satisfy every operation in `ops` (the null space of [`stacked_constraint_system`](@ref)'s
system, found via SVD). Used to build mechanically-correct-by-construction synthetic presets for
symmetry classes with no real reference material at hand (monoclinic, triclinic below).
"""
function project_onto_class(C_generic, ops)
    isempty(ops) && return copy(C_generic)
    stacked = stacked_constraint_system(ops)
    F = svd(stacked; full=true)
    tol = 1e-6 * max(maximum(F.S), 1.0)
    rank_ = count(>(tol), F.S)
    Vnull = F.V[:, rank_+1:end]
    vgeneric = [C_generic[a, b] for (a, b) in VOIGT_PAIRS]
    vproj = Vnull * (Vnull' * vgeneric)
    Cproj = zeros(Float64, 6, 6)
    for (k, (a, b)) in enumerate(VOIGT_PAIRS)
        Cproj[a, b] = Cproj[b, a] = vproj[k]
    end
    Cproj
end

# ╔═╡ e3a630d5-a7b1-45d1-b0ed-b47646ef7d87
"""
    GENERATOR_SETS

The generator operations that define each of the 6 elastic symmetry classes, and `TARGET_COUNT`,
the textbook independent-constant count each set must reduce 21 down to. These generator sets were
found by trial against `class_constraints` (below) until the count matched exactly — not assumed.
"""
const GENERATOR_SETS = Dict(
    :triclinic => Symbol[],
    :monoclinic => [:mirror_xy],
    :orthorhombic => [:rot180x, :rot180y, :rot180z],
    :hexagonal => [:contz, :rot180x],
    :cubic => [:rot90z, :rot120_111],
    :isotropic => [:contz, :rot120_111],
)

# ╔═╡ b57f45f1-9f49-4568-8b23-46c818d4bcf8
const TARGET_COUNT = Dict(:triclinic => 21, :monoclinic => 13, :orthorhombic => 9, :hexagonal => 5, :cubic => 3, :isotropic => 2)

# ╔═╡ fc3ed6e9-19bd-4170-b90a-b3c013b71c93
md"""
### Verifying the Symmetry-Derivation Machinery
"""

# ╔═╡ c4b0936b-3db1-4a3f-a764-f9cf0ee4934b
let
    # every discrete op's L must be a signed permutation: exactly one +-1 per row
    for op in filter(!=(:contz), ALL_OPS)
        L = build_L(symmetry_generator(op))
        for row in 1:size(L, 1)
            nz = findall(x -> abs(x) > 1e-8, L[row, :])
            @assert length(nz) == 1 "op $op row $row: $(length(nz)) nonzero entries"
            @assert isapprox(abs(L[row, nz[1]]), 1.0; atol=1e-8) "op $op row $row: entry not ±1"
        end
    end

    # inversion must ALWAYS be a symmetry -- the elastic tensor can never see chirality
    for trial in 1:5
        C = zeros(Float64, 6, 6)
        for a in 1:6, b in a:6
            C[a, b] = C[b, a] = randn()
        end
        @assert is_symmetry_operation(C, symmetry_generator(:inversion))
    end

    counts = Dict(cls => class_constraints(ops).independent_count for (cls, ops) in GENERATOR_SETS)
    for (cls, target) in TARGET_COUNT
        @assert counts[cls] == target "class $cls: got $(counts[cls]), expected $target"
    end

    hexres = class_constraints(GENERATOR_SETS[:hexagonal])
    @assert length(hexres.constraints.derived) == 1 "expected exactly one derived relation for hexagonal"

    counts_str = join(["$(cls)→$(counts[cls])" for cls in (:triclinic, :monoclinic, :orthorhombic, :hexagonal, :cubic, :isotropic)], ", ")
    md"""
    !!! correct "Self-check"
        Every discrete operation's linear map is a signed permutation ✓ · inversion is always a
        symmetry of any elastic tensor (chirality is invisible to elasticity) ✓ · every one of the
        six proposed generator sets reduces 21 independent constants to exactly its textbook count
        — $(counts_str) ✓
        · hexagonal's generators force exactly one derived (not merely equal) relation, matching
        the textbook `C66=(C11-C12)/2` ✓
    """
end

# ╔═╡ 6cffaa93-19f8-483c-a780-ce8cd54b33ed
md"""
### From Stiffness to Compliance: Directional Young's Modulus

Directional stiffness needs the **compliance** tensor `S = C⁻¹`, and this is where a classic bug
lives: converting the Voigt (engineering) compliance matrix into the true tensor `S_ijkl` needs a
correction the stiffness conversion ([`get_cijkl`](@ref)) does **not** — engineering shear strain
is defined as `γ=2ε`, doubling shear terms, so recovering the tensor divides by 1 (normal-normal),
2 (one shear index), or 4 (both shear) *for compliance only*. Getting this backwards, or "fixing"
[`get_cijkl`](@ref) to match, is the trap; the two are asymmetric on purpose.
"""

# ╔═╡ 3c043dec-4383-4f2f-b77e-4203d2dd734a
"""
    get_sijkl(C)

The compliance tensor `S_ijkl`, derived directly: invert the Voigt stiffness matrix, then divide
each entry by the engineering shear factor (1, 2, or 4) appropriate to how many of its four tensor
indices are a shear pair.
"""
function get_sijkl(C)
    Svoigt = inv(C)
    factor(i, j) = i == j ? 1.0 : 2.0
    idx(i, j) = i == j ? i : 9 - i - j
    [Svoigt[idx(i, j), idx(k, l)] / (factor(i, j) * factor(k, l)) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
end

# ╔═╡ 08baab70-8994-4e90-9a91-98d2308a9418
"""
    get_sijkl_mandel(C)

An independent second derivation of the same compliance tensor, via Kelvin-Mandel scaling
(`D=diag(1,1,1,√2,√2,√2)`), under which stiffness and compliance genuinely are matrix inverses of
each other (`S^Mandel = (D·C·D)⁻¹`) with no ad hoc shear factor at the matrix-inversion step — the
factor still has to come back out afterward (`inv(D)·Sᴹ·inv(D)`) to recover the plain tensor, since
Mandel normalization doesn't remove it, only defers it. Used purely as a cross-check against
[`get_sijkl`](@ref): a different arithmetic path (a different matrix inversion) that must land on
the same tensor.
"""
function get_sijkl_mandel(C)
    d = [1.0, 1.0, 1.0, sqrt(2), sqrt(2), sqrt(2)]
    D = Diagonal(d)
    Cm = D * C * D
    Sm = inv(Cm)
    S_normalized = inv(D) * Sm * inv(D)
    get_cijkl(S_normalized)
end

# ╔═╡ ec8260cc-96c3-43c9-8b43-c5ee13a0fa93
"""
    directional_youngs_modulus(n, Sijkl)

The directional Young's modulus `E(n̂) = 1/(S_ijkl n_i n_j n_k n_l)` for unit direction `n` and
compliance tensor `Sijkl` (from [`get_sijkl`](@ref)) — the quantity Panel 3's stiffness surface
will plot as a radius in every direction.
"""
function directional_youngs_modulus(n, Sijkl)
    s = 0.0
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        s += Sijkl[i, j, k, l] * n[i] * n[j] * n[k] * n[l]
    end
    1.0 / s
end

# ╔═╡ c71483d4-6ee8-4dc3-9cd9-dcabdc27941c
"""
    fibonacci_sphere(n)

`n` roughly-evenly-spaced points on the unit sphere via the Fibonacci-spiral construction — used
below to sample directions for the isotropic Young's-modulus self-check, and later for Panel 3's
stiffness-surface sampling grid.
"""
function fibonacci_sphere(n)
    pts = Vector{Vector{Float64}}(undef, n)
    ga = pi * (3 - sqrt(5))
    for i in 0:n-1
        y = 1 - 2 * (i / (n - 1))
        r = sqrt(max(0.0, 1 - y^2))
        θ = ga * i
        pts[i+1] = [r * cos(θ), y, r * sin(θ)]
    end
    pts
end

# ╔═╡ 288d863f-fc57-48e0-89d5-bed10de6506a
"""
    isotropic_C(λ, μ)

The 6×6 Voigt matrix for an isotropic material from its two Lamé parameters: `C11=λ+2μ`, `C12=λ`,
`C44=μ` (and `C66=(C11-C12)/2=μ` automatically, consistent with the hexagonal derived relation
above at the isotropic limit).
"""
function isotropic_C(λ, μ)
    C = zeros(Float64, 6, 6)
    C[1, 1] = C[2, 2] = C[3, 3] = λ + 2μ
    C[1, 2] = C[2, 1] = C[1, 3] = C[3, 1] = C[2, 3] = C[3, 2] = λ
    C[4, 4] = C[5, 5] = C[6, 6] = μ
    C
end

# ╔═╡ 214a6443-67b8-40e1-90f7-c5d771b79694
md"""
### Verifying the Compliance Conversion
"""

# ╔═╡ 383498a0-b2cc-4deb-9b6a-86eef802049f
let
    # independent-derivation agreement
    maxdiff = 0.0
    for trial in 1:5
        C = zeros(Float64, 6, 6)
        for a in 1:6, b in a:6
            C[a, b] = C[b, a] = 50 + 40 * rand()
        end
        C = C + 200I
        S1, S2 = get_sijkl(C), get_sijkl_mandel(C)
        maxdiff = max(maxdiff, maximum(abs.(S1 .- S2)))
    end
    @assert maxdiff < 1e-9 "get_sijkl vs get_sijkl_mandel disagree by $maxdiff"

    # isotropic: E(n) constant in every direction, matching the closed form
    λ, μ = 40.0, 30.0
    Siso = get_sijkl(isotropic_C(λ, μ))
    E_formula = μ * (3λ + 2μ) / (λ + μ)
    Es = [directional_youngs_modulus(n, Siso) for n in fibonacci_sphere(50)]
    @assert all(e -> isapprox(e, E_formula; rtol=1e-9), Es) "isotropic E(n) not constant: $(extrema(Es))"

    # cubic: independently-derivable closed form at n=(1,1,0)/sqrt(2)
    Ccub = zeros(Float64, 6, 6)
    Ccub[1, 1] = Ccub[2, 2] = Ccub[3, 3] = 250.0
    Ccub[1, 2] = Ccub[2, 1] = Ccub[1, 3] = Ccub[3, 1] = Ccub[2, 3] = Ccub[3, 2] = 150.0
    Ccub[4, 4] = Ccub[5, 5] = Ccub[6, 6] = 100.0
    Svoigt = inv(Ccub)
    S11, S12, S44 = Svoigt[1, 1], Svoigt[1, 2], Svoigt[4, 4]
    n = [1, 1, 0] / sqrt(2)
    invE_formula = S11 - 2 * (S11 - S12 - S44 / 2) * (n[1]^2 * n[2]^2 + n[2]^2 * n[3]^2 + n[3]^2 * n[1]^2)
    E_direct = directional_youngs_modulus(n, get_sijkl(Ccub))
    @assert isapprox(1 / invE_formula, E_direct; rtol=1e-9) "cubic closed-form mismatch"

    md"""
    !!! correct "Self-check"
        Two independent derivations of the compliance tensor (direct shear-factor division vs.
        Kelvin-Mandel normalization) agree to within $(round(maxdiff, sigdigits=2)) ✓ · an
        isotropic material's `E(n)` is constant across 50 sampled directions, matching
        `E=μ(3λ+2μ)/(λ+μ)=$(round(μ*(3*40.0+2*30.0)/(40.0+30.0), digits=3))` GPa exactly ✓ · a
        cubic material's `E(n)` at `n=(1,1,0)/√2` matches the independently-derivable closed form
        `1/E=S11-2(S11-S12-S44/2)(n₁²n₂²+n₂²n₃²+n₃²n₁²)` ✓
    """
end

# ╔═╡ 982529a6-007c-4052-9e48-c664e7915dc0
md"""
### Wave Speeds: the Christoffel Equation

The three seismic body-wave speeds (one quasi-P, two quasi-S) in direction `n̂` are
`ρv² = ` the eigenvalues of the **Christoffel matrix**
`Γ_ik = Σⱼₗ c_ijkl n_j n_l`. This is the same eigenproblem `anisotropy.jl` already solves for VTI
(lines 457-467); nothing about the contraction or the eigensolve depends on symmetry class, so it's
reused verbatim.
"""

# ╔═╡ 18df75f7-0850-49e0-9838-f629c101fe5c
"""
    christoffel_velocities(n, c, ρ)

The three wave speeds `v = sqrt(eigenvalue/ρ)` and polarization vectors for direction `n̂`, elastic
tensor `c` (from [`get_cijkl`](@ref)), and density `ρ` — the Christoffel eigenproblem, reused from
`src/Planewave Propagation/anisotropy.jl:457-467`.
"""
function christoffel_velocities(n, c, ρ)
    Γ = zeros(Float64, 3, 3)
    for i in 1:3, k in 1:3
        s = 0.0
        for j in 1:3, l in 1:3
            s += c[i, j, k, l] * n[j] * n[l]
        end
        Γ[i, k] = s
    end
    E = eigen(Symmetric(Γ))
    (v=sqrt.(max.(E.values, 0.0) ./ ρ), pol=E.vectors)
end

# ╔═╡ a39cb03c-b7e0-40e4-a2d4-714388e580fe
md"""
### Verifying the Christoffel Equation
"""

# ╔═╡ 096270f3-5ce6-4a45-bd9b-2f8114cabc51
let
    ρ = 3.3
    c_iso = get_cijkl(isotropic_C(40.0, 30.0))
    for n in fibonacci_sphere(20)
        vs = sort(christoffel_velocities(n, c_iso, ρ).v)
        @assert isapprox(vs[1], vs[2]; rtol=1e-8) "isotropic shear speeds not degenerate: $vs"
    end

    Colivine = [192.0 66 60 0 0 0; 66 160 56 0 0 0; 60 56 272 0 0 0; 0 0 0 60 0 0; 0 0 0 0 62 0; 0 0 0 0 0 49]
    c_ol = get_cijkl(Colivine)
    n_offaxis = [1, 1, 1] / sqrt(3)
    vs_ol = sort(christoffel_velocities(n_offaxis, c_ol, ρ).v)
    splitting = (vs_ol[2] - vs_ol[1]) / vs_ol[1]
    @assert splitting > 0.01 "expected visible shear splitting for olivine, got $(splitting*100)%"

    md"""
    !!! correct "Self-check"
        An isotropic material's two shear-wave speeds are degenerate (no splitting) in every
        sampled direction ✓ · olivine's off-axis shear waves visibly split by
        **$(round(splitting*100, digits=1))%** — the same shear-wave-splitting phenomenon SKS
        splitting measures in the real Earth ✓
    """
end

# ╔═╡ 135da962-7b36-4b0c-9c68-f2313cfb65d7
md"""
### Presets for Each Symmetry Class

Isotropic is parameterized directly (two numbers, no lookup needed). Hexagonal and orthorhombic
reuse real, already-validated values from `anisotropy.jl` — its VTI defaults and its olivine
matrix, respectively. Cubic, monoclinic, and triclinic have no existing repo reference, so their
presets are clearly-labeled **illustrative synthetic** numbers, not a specific measured material —
monoclinic and triclinic are built by hand at exactly the zero pattern [`class_constraints`](@ref)
already derived above (never a separately-typed pattern that could silently disagree with it).
"""

# ╔═╡ fc1ed39a-818d-4ff5-8e82-dd4c1e028acd
"""
    vti_C(A, C33, L, N, F)

The 6×6 Voigt matrix for a transversely-isotropic (hexagonal) material in the same `(A,C,L,N,F)`
parameterization `anisotropy.jl` uses (`A=C11=C22`, `C33`, `L=C44=C55`, `N=C66`, `F=C13=C23`, with
`C12=A-2N`, matching the derived relation above).
"""
function vti_C(A, C33, L, N, F)
    C = zeros(Float64, 6, 6)
    C[1, 1] = C[2, 2] = A
    C[3, 3] = C33
    C[1, 2] = C[2, 1] = A - 2N
    C[1, 3] = C[3, 1] = C[2, 3] = C[3, 2] = F
    C[4, 4] = C[5, 5] = L
    C[6, 6] = N
    C
end

# ╔═╡ b2faf56b-88f6-4b70-b5b5-4129a28511bc
"""
    PRESET_C

The default 6×6 Voigt matrix for each of the six symmetry classes (GPa). `:hexagonal` reuses
`anisotropy.jl`'s own VTI defaults (A=272,C=160,L=60,N=50,F=60); `:orthorhombic` reuses its olivine
matrix (`anisotropy.jl:182`) verbatim. `:cubic`, `:monoclinic`, `:triclinic` are illustrative
synthetic values — see the module docs above, not a specific measured mineral.
"""
const PRESET_C = Dict(
    :isotropic => isotropic_C(40.0, 30.0),
    :hexagonal => vti_C(272.0, 160.0, 60.0, 50.0, 60.0),
    :orthorhombic => [192.0 66 60 0 0 0; 66 160 56 0 0 0; 60 56 272 0 0 0; 0 0 0 60 0 0; 0 0 0 0 62 0; 0 0 0 0 0 49],
    :cubic => let C = zeros(Float64, 6, 6)
        C[1,1]=C[2,2]=C[3,3]=250.0
        C[1,2]=C[2,1]=C[1,3]=C[3,1]=C[2,3]=C[3,2]=150.0
        C[4,4]=C[5,5]=C[6,6]=100.0
        C
    end,
    :monoclinic => let C = zeros(Float64, 6, 6)
        C[1,1]=220.0; C[2,2]=200.0; C[3,3]=260.0
        C[1,2]=C[2,1]=70.0; C[1,3]=C[3,1]=65.0; C[2,3]=C[3,2]=60.0
        C[4,4]=55.0; C[5,5]=58.0; C[6,6]=52.0
        C[1,6]=C[6,1]=30.0; C[2,6]=C[6,2]=-25.0; C[3,6]=C[6,3]=15.0
        C[4,5]=C[5,4]=-18.0
        C
    end,
    :triclinic => let C = zeros(Float64, 6, 6)
        C[1,1]=220.0; C[2,2]=200.0; C[3,3]=260.0
        C[1,2]=C[2,1]=70.0; C[1,3]=C[3,1]=65.0; C[2,3]=C[3,2]=60.0
        C[4,4]=55.0; C[5,5]=58.0; C[6,6]=52.0
        C[1,6]=C[6,1]=30.0; C[2,6]=C[6,2]=-25.0; C[3,6]=C[6,3]=15.0
        C[4,5]=C[5,4]=-18.0
        C[1,4]=C[4,1]=12.0; C[1,5]=C[5,1]=-9.0
        C[2,4]=C[4,2]=8.0; C[2,5]=C[5,2]=14.0
        C[3,4]=C[4,3]=-11.0; C[3,5]=C[5,3]=7.0
        C[4,6]=C[6,4]=10.0; C[5,6]=C[6,5]=-13.0
        C
    end,
)

# ╔═╡ ad3d84e0-b1fa-4d46-9213-cddb27af546b
md"""
### Verifying the Presets
"""

# ╔═╡ 6884b812-4802-49ba-a38e-0446112d3f00
let
    reports = String[]
    for (cls, ops) in GENERATOR_SETS
        C = PRESET_C[cls]
        eigs = eigen(Symmetric(C)).values
        @assert all(>(0), eigs) "$cls preset is not positive-definite: $eigs"
        for op in ops
            if op == :contz
                @assert all(θ -> is_symmetry_operation(C, symmetry_generator(:contz; theta=θ)), (0.3, 1.1, 2.4)) "$cls preset fails contz"
            else
                @assert is_symmetry_operation(C, symmetry_generator(op)) "$cls preset fails generator $op"
            end
        end
        n_holding = count(op -> op == :contz ? is_symmetry_operation(C, symmetry_generator(:contz; theta=0.53)) : is_symmetry_operation(C, symmetry_generator(op)), ALL_OPS)
        push!(reports, "$cls: stable, satisfies its $(length(ops)) generator(s), $(n_holding)/10 catalog ops hold")
    end

    # monoclinic/triclinic were hand-built at exactly class_constraints' own derived zero pattern --
    # confirm projecting them onto their own class is a no-op (they're already exactly on it)
    for cls in (:monoclinic, :triclinic)
        Cproj = project_onto_class(PRESET_C[cls], GENERATOR_SETS[cls])
        @assert isapprox(Cproj, PRESET_C[cls]; atol=1e-8) "$cls preset isn't exactly on its own symmetry subspace"
    end

    reports_str = join(["- " * r for r in reports], "\n")
    md"""
    !!! correct "Self-check"
        Every preset is positive-definite (physically stable) ✓ and satisfies exactly its own
        class's generators ✓ · the hand-built monoclinic and triclinic presets sit exactly on
        their class's derived symmetry subspace (projecting changes nothing) ✓

        $(reports_str)
    """
end

# ╔═╡ cfd61b4e-5a35-4c55-957f-3e2c7519968e
md"""
### Packaging Data for the Widget

Because every symmetry class here is a discrete preset (not a continuously-dragged parameter),
**everything the widget needs can be precomputed once, in Julia, at notebook-build time** — there
is no live `@bind` round trip anywhere in this notebook. The functions below assemble every class's
matrix, symmetry-test results, matrix-structure constraints, and stiffness/wave-speed sampling
grids into one hand-built JSON string (no JSON package is a dependency anywhere in this repo;
`anisotropy.jl`'s `ani_push_message` sets the precedent), embedded directly into the widget's
`<script>` tag below. JavaScript only ever reads this table — it never recomputes any of it.
"""

# ╔═╡ a59ccaaf-9d1f-40a1-8e85-dd737feca733
"""
    theta_phi_grid_directions(ntheta, nphi)

An `ntheta`×`nphi` grid of unit direction vectors, `θ` (polar angle from z) spanning `[0,π]` and
`φ` (azimuth) spanning `[0,2π)`. The sampling grid Panel 3's stiffness/wave-speed surfaces are
drawn on.
"""
function theta_phi_grid_directions(ntheta, nphi)
    thetas = range(0, pi; length=ntheta)
    phis = range(0, 2pi; length=nphi + 1)[1:end-1]
    [[sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)] for θ in thetas, φ in phis]
end

# ╔═╡ 3828fd1e-a4de-4f72-8cfd-115be3f731e5
"""
    youngs_modulus_grid(C; ntheta, nphi)

[`directional_youngs_modulus`](@ref) sampled over [`theta_phi_grid_directions`](@ref)'s grid, as an
`ntheta`×`nphi` matrix — Panel 3's "modulus surface" data for one class.
"""
function youngs_modulus_grid(C; ntheta, nphi)
    Sijkl = get_sijkl(C)
    [directional_youngs_modulus(n, Sijkl) for n in theta_phi_grid_directions(ntheta, nphi)]
end

# ╔═╡ 269f858c-27a3-447d-b4dd-5a8a839b2bb6
"""
    wave_speed_grid(C; ρ, ntheta, nphi)

The three sorted [`christoffel_velocities`](@ref) at every direction of
[`theta_phi_grid_directions`](@ref)'s grid, as an `ntheta`×`nphi`×3 array — Panel 3's "wave-speed
surface" data for one class (the two shear branches visibly split apart wherever the material is
anisotropic).
"""
function wave_speed_grid(C; ρ=3.3, ntheta, nphi)
    c = get_cijkl(C)
    dirs = theta_phi_grid_directions(ntheta, nphi)
    out = zeros(Float64, ntheta, nphi, 3)
    for i in 1:ntheta, j in 1:nphi
        out[i, j, :] = sort(christoffel_velocities(dirs[i, j], c, ρ).v)
    end
    out
end

# ╔═╡ 46fd8d6f-3fa2-45f2-9950-6232f85c48e6
const NTHETA, NPHI = 21, 41

# ╔═╡ 8a610175-37ef-4a3f-a919-8f9cfbcb6759
begin
    _jfloats(v) = "[" * join(round.(Float64.(v); digits=6), ",") * "]"
    _jints(v) = "[" * join(Int.(v), ",") * "]"
    _flatten_theta_phi(M) = vec(permutedims(M))  # row-major: theta outer, phi inner
    function _flatten_wave_grid(wgrid)
        ntheta, nphi, _ = size(wgrid)
        out = Float64[]
        for i in 1:ntheta, j in 1:nphi
            append!(out, wgrid[i, j, :])
        end
        out
    end
end

# ╔═╡ 91872170-0f34-4f3e-9aed-4d2335024f7b
"""
    json_constraints(oc::OpConstraints)

Serialize an [`OpConstraints`](@ref) to a JSON object string
`{"zeros":[...],"families":[[...]],"derived":[[index,[[other_index,coeff],...]],...]}`
(Voigt-pair indices are 1-based, matching [`VOIGT_PAIRS`](@ref); JS subtracts 1 on read).
"""
function json_constraints(oc::OpConstraints)
    zeros_str = _jints(oc.zeros)
    families_str = "[" * join(["[" * join(fam, ",") * "]" for fam in oc.families], ",") * "]"
    derived_str = "[" * join(
        ["[$(idx),[" * join(["[$(oi),$(round(c; digits=6))]" for (oi, c) in terms], ",") * "]]"
         for (idx, terms) in oc.derived], ",") * "]"
    "{\"zeros\":$zeros_str,\"families\":$families_str,\"derived\":$derived_str}"
end

# ╔═╡ 793cf760-0258-44c7-8cc9-f426d68e7634
"""
    json_class_data(cls)

One symmetry class's full data package for the widget: its preset matrix, independent-constant
count, which of the 10 catalog operations it satisfies, its own derived matrix-structure
constraints, and its modulus/wave-speed sampling grids.
"""
function json_class_data(cls)
    C = PRESET_C[cls]
    cc = class_constraints(GENERATOR_SETS[cls])
    valid = Dict(op => (op == :contz ?
                        all(θ -> is_symmetry_operation(C, symmetry_generator(:contz; theta=θ)), (0.31, 1.2, 2.5)) :
                        is_symmetry_operation(C, symmetry_generator(op)))
                 for op in ALL_OPS)
    valid_str = "{" * join(["\"$(op)\":$(valid[op] ? "true" : "false")" for op in ALL_OPS], ",") * "}"
    ygrid = youngs_modulus_grid(C; ntheta=NTHETA, nphi=NPHI)
    wgrid = wave_speed_grid(C; ntheta=NTHETA, nphi=NPHI)
    """{"C":$(_jfloats(vec(permutedims(C)))),"independentCount":$(TARGET_COUNT[cls]),\
"validOps":$valid_str,"constraints":$(json_constraints(cc.constraints)),\
"youngsGrid":$(_jfloats(_flatten_theta_phi(ygrid))),"waveGrid":$(_jfloats(_flatten_wave_grid(wgrid)))}"""
end

# ╔═╡ 188ae616-5423-484c-b22d-621e52e529a0
"""
    json_op_data(op)

One catalog operation's data package: the matrix-structure constraints it imposes on a generic
matrix (used for Panel 2's "incoming pulse" highlight, independent of the active class).
"""
json_op_data(op) = "{\"constraints\":$(json_constraints(operation_constraints(op)))}"

# ╔═╡ 2787de1b-57cd-4567-b46e-aeff7b68b493
const ESYM_DATA_JSON = let
    classes_str = "{" * join(["\"$(cls)\":$(json_class_data(cls))" for cls in keys(GENERATOR_SETS)], ",") * "}"
    ops_str = "{" * join(["\"$(op)\":$(json_op_data(op))" for op in ALL_OPS], ",") * "}"
    voigt_labels = ["$(a)$(b)" for (a, b) in VOIGT_PAIRS]
    labels_str = "[" * join(["\"$(l)\"" for l in voigt_labels], ",") * "]"
    pairs_str = "[" * join(["[$(a),$(b)]" for (a, b) in VOIGT_PAIRS], ",") * "]"
    """{"classes":$classes_str,"ops":$ops_str,"voigtLabels":$labels_str,"voigtPairs":$pairs_str,"ntheta":$NTHETA,"nphi":$NPHI}"""
end

# ╔═╡ 07d767d0-64cd-4ecb-9a1c-7126c483c67c
md"""
### Verifying the Packaged Data
"""

# ╔═╡ 9fbf8b5e-ec64-4e5f-ba44-b29ce5870efc
let
    @assert length(ESYM_DATA_JSON) > 1000 "packaged JSON suspiciously short"
    # balanced braces/brackets is a strong structural sanity check on hand-built JSON
    depth = 0
    minseen = 0
    for ch in ESYM_DATA_JSON
        if ch in ('{', '[')
            depth += 1
        elseif ch in ('}', ']')
            depth -= 1
        end
        minseen = min(minseen, depth)
    end
    @assert depth == 0 "unbalanced braces/brackets in packaged JSON (ended at depth $depth)"
    @assert minseen == 0 "packaged JSON closes a bracket before it was opened"
    for cls in keys(GENERATOR_SETS)
        @assert occursin("\"$(cls)\"", ESYM_DATA_JSON)
    end
    md"""
    !!! correct "Self-check"
        The packaged JSON string has balanced braces/brackets throughout (a strong structural
        sanity check on hand-built JSON, catching the most common class of bug) ✓ and every
        symmetry class's key is present ✓. Full parseability is confirmed live in the browser via
        `JSON.parse`, once the widget below runs.
    """
end

# ╔═╡ ddaa44f6-9462-4238-a06b-fc2a26a3e06c
md"""
### The Interactive Widget

`ElasticSymmetryInput` does no physics — every number it draws (the six presets, which of the 10
catalog operations each class satisfies, each operation's own matrix-structure constraints, the
modulus/wave-speed sampling grids) was already computed and self-checked above and is handed to it
as one JSON string, [`ESYM_DATA_JSON`](@ref). Because every symmetry class here is a discrete
preset rather than a continuously-dragged parameter, nothing needs to round-trip back through
Julia during an interaction — the widget has no `@bind` value at all, just three canvases reading
the same in-page JavaScript state object (`symmetryClass`, shared `orientation`, `hovered`) and
re-rendering from it. Camera orbiting and the symmetry-operation "does it look the same?" animation
are the only things computed in JavaScript, and both are pure display geometry with no physical
content of their own — whether an operation actually *is* a symmetry was already decided in Julia
(`validOps`, above) and only looked up here.
"""

# ╔═╡ 2e481c84-d654-4c7d-a8c6-e0ff3a0c56d3
begin
    struct ElasticSymmetryInput
        dataJSON::String
    end

    function Base.show(io::IO, ::MIME"text/html", w::ElasticSymmetryInput)
        write(io, """
        <div id="esymwidget">
        <style>
        #esymwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #esymwidget .esym-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #esymwidget .esym-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #esymwidget .esym-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #esymwidget .esym-classrow{display:flex;gap:6px;flex-wrap:wrap;justify-content:center;margin-bottom:12px}
        #esymwidget .esym-classbtn{border-radius:4px;border:1px solid #9ca3af;background:#374151;color:#f3f4f6;
          padding:7px 14px;font-size:13px;cursor:pointer}
        #esymwidget .esym-classbtn.active{background:#2563eb;border-color:#60a5fa}
        #esymwidget .esym-workspace{display:flex;gap:16px;flex-wrap:wrap;justify-content:center;align-items:flex-start}
        #esymwidget .esym-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #esymwidget .esym-panel-title{font-size:14px;font-weight:700;color:#e5e7eb;text-align:center;margin-bottom:4px}
        #esymwidget .esym-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:4px;min-height:32px}
        #esymwidget canvas{display:block}
        #esymwidget #esym-block, #esymwidget #esym-surface{cursor:grab}
        #esymwidget .esym-ops-grid{display:grid;grid-template-columns:1fr 1fr;gap:4px;margin-top:8px;width:260px}
        #esymwidget .esym-opbtn{border-radius:4px;border:1px solid #6b7280;background:#374151;color:#f3f4f6;
          padding:5px 4px;font-size:11px;cursor:pointer}
        #esymwidget .esym-opbtn.valid{background:#166534;border-color:#4ade80}
        #esymwidget .esym-opbtn.invalid{background:#7f1d1d;border-color:#f87171}
        #esymwidget .esym-toggle-row{display:flex;gap:6px;justify-content:center;margin-top:8px}
        #esymwidget .esym-togglebtn{border-radius:4px;border:1px solid #9ca3af;background:#374151;color:#f3f4f6;
          padding:5px 12px;font-size:12px;cursor:pointer}
        #esymwidget .esym-togglebtn.active{background:#2563eb;border-color:#60a5fa}
        #esymwidget .esym-count{font-size:22px;font-weight:700;color:#facc15;text-align:center}
        </style>

        <div class="esym-title">
          <div class="esym-title-desc">More symmetry survived &rarr; more forced relations in <i>C</i> &rarr; fewer independent constants &rarr; a rounder stiffness surface.</div>
          <div class="esym-title-hint">pick a symmetry class &middot; drag the block or the surface to orbit (they turn together) &middot; click an operation to test it &middot; hover a matrix cell</div>
        </div>

        <div class="esym-classrow" id="esym-classrow"></div>

        <div class="esym-workspace">
          <div>
            <div class="esym-panel-title">The Material Block</div>
            <div class="esym-panel"><canvas id="esym-block"></canvas></div>
            <div class="esym-caption" id="esym-block-caption"></div>
            <div class="esym-ops-grid" id="esym-ops-grid"></div>
          </div>
          <div>
            <div class="esym-panel-title">The Voigt Matrix</div>
            <div class="esym-panel"><canvas id="esym-matrix"></canvas></div>
            <div class="esym-caption" id="esym-matrix-caption">independent constants: <span class="esym-count" id="esym-count"></span></div>
          </div>
          <div>
            <div class="esym-panel-title">Directional Stiffness</div>
            <div class="esym-panel"><canvas id="esym-surface"></canvas></div>
            <div class="esym-caption" id="esym-surface-caption"></div>
            <div class="esym-toggle-row" id="esym-toggle-row"></div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        const DATA = JSON.parse('$(w.dataJSON)');
        const CLASS_ORDER = ['isotropic','cubic','hexagonal','orthorhombic','monoclinic','triclinic'];
        const CLASS_LABEL = {isotropic:'Isotropic',cubic:'Cubic',hexagonal:'Hexagonal',orthorhombic:'Orthorhombic',monoclinic:'Monoclinic',triclinic:'Triclinic'};
        const OP_ORDER = ['rot90z','rot180x','rot180y','rot180z','mirror_xy','mirror_xz','mirror_yz','inversion','rot120_111','contz'];
        const OP_LABEL = {rot90z:'90° / z',rot180x:'180° / x',rot180y:'180° / y',rot180z:'180° / z',
          mirror_xy:'mirror xy',mirror_xz:'mirror xz',mirror_yz:'mirror yz',inversion:'inversion',
          rot120_111:'120° / [111]',contz:'spin about z'};

        let state = { cls: 'orthorhombic', az: -0.55, el: 0.35, hoveredCell: null, surfaceMode: 'modulus',
          activeOp: null, opAnimStart: null, opAnimFinished: false, activeProbe: null, probeAnimStart: null };

        function currentClassData(){ return DATA.classes[state.cls]; }

        // ---------- shared camera (Panel 1 block + Panel 3 surface turn together) ----------
        function rot3(p, az, el){
          const ca=Math.cos(az), sa=Math.sin(az);
          const x1 = p[0]*ca - p[1]*sa, y1 = p[0]*sa + p[1]*ca;
          const ce=Math.cos(el), se=Math.sin(el);
          const y2 = y1*ce - p[2]*se, z2 = y1*se + p[2]*ce;
          return [x1, y2, z2];
        }
        function proj(p, cx, cy, scale, camD){
          const f = camD/(camD+p[2]);
          return [cx + p[0]*scale*f, cy - p[1]*scale*f, p[2]];
        }
        function axisAngleRotate(p, axis, angle){
          const [ax,ay,az2] = axis, c=Math.cos(angle), s=Math.sin(angle), t=1-c;
          const [x,y,z] = p;
          return [
            (t*ax*ax+c)*x + (t*ax*ay-s*az2)*y + (t*ax*az2+s*ay)*z,
            (t*ax*ay+s*az2)*x + (t*ay*ay+c)*y + (t*ay*az2-s*ax)*z,
            (t*ax*az2-s*ay)*x + (t*ay*az2+s*ax)*y + (t*az2*az2+c)*z,
          ];
        }
        // pure display geometry: how to animate a block through a candidate symmetry operation --
        // no physics here, whether it's actually a symmetry was already decided in Julia (validOps)
        function applyOpAnim(op, t, p){
          if(!op || t===null) return p;
          const HALF = Math.PI, THIRD2 = 2*Math.PI/3, QUART = Math.PI/2, FULL = 2*Math.PI;
          if(op==='rot90z') return axisAngleRotate(p, [0,0,1], QUART*t);
          if(op==='rot180x') return axisAngleRotate(p, [1,0,0], HALF*t);
          if(op==='rot180y') return axisAngleRotate(p, [0,1,0], HALF*t);
          if(op==='rot180z') return axisAngleRotate(p, [0,0,1], HALF*t);
          if(op==='rot120_111'){ const a=1/Math.sqrt(3); return axisAngleRotate(p, [a,a,a], THIRD2*t); }
          if(op==='contz') return axisAngleRotate(p, [0,0,1], FULL*t);
          const s = Math.cos(Math.PI*t);
          if(op==='mirror_xy') return [p[0],p[1],p[2]*s];
          if(op==='mirror_xz') return [p[0],p[1]*s,p[2]];
          if(op==='mirror_yz') return [p[0]*s,p[1],p[2]];
          if(op==='inversion') return [p[0]*s,p[1]*s,p[2]*s];
          return p;
        }
        // probe deformation: a purely illustrative "poke and release" pulse along one Voigt
        // strain direction (b=1..6), display geometry only -- its magnitude has no physical
        // scale, it exists only to show WHICH direction that strain column acts along, matching
        // the highlighted column in Panel 2.
        function applyStrainDeform(b, t, p){
          if(b===null || t===null) return p;
          const amt = 0.32*Math.sin(Math.PI*Math.min(1,t));
          const [x,y,z] = p;
          if(b===1) return [x*(1+amt), y, z];
          if(b===2) return [x, y*(1+amt), z];
          if(b===3) return [x, y, z*(1+amt)];
          if(b===4) return [x, y+amt*z, z];
          if(b===5) return [x+amt*z, y, z];
          if(b===6) return [x+amt*y, y, z];
          return p;
        }

        // ---------- Panel 1: the material block ----------
        // geometric correspondence, chosen so the shape's OWN symmetry group matches the class's
        // elastic symmetry group exactly: a centered, z-extruded prism (bottom/top face joined by
        // a purely-vertical e3) is automatically invariant under mirror_xy/rot180z regardless of
        // its cross-section shape; tilting e3 (triclinic) breaks that automatically.
        function ngon(n, r){ const pts=[]; for(let i=0;i<n;i++){ const a=2*Math.PI*i/n; pts.push([r*Math.cos(a), r*Math.sin(a)]); } return pts; }
        function rectCorners(hw, hh, shear){ shear = shear||0; return [[-hw-shear,-hh],[hw-shear,-hh],[hw+shear,hh],[-hw+shear,hh]]; }
        function blockGeom(cls){
          if(cls==='cubic') return {bottom: rectCorners(0.5,0.5), e3:[0,0,1.0], texture:'lattice', color:'#60a5fa'};
          if(cls==='hexagonal') return {bottom: ngon(16,0.58), e3:[0,0,0.75], texture:'stripes', color:'#34d399'};
          if(cls==='orthorhombic') return {bottom: rectCorners(0.5,0.65), e3:[0,0,0.4], texture:'plain', color:'#fbbf24'};
          if(cls==='monoclinic') return {bottom: rectCorners(0.5,0.62,0.32), e3:[0,0,0.42], texture:'plain', color:'#fb923c'};
          if(cls==='triclinic') return {bottom: rectCorners(0.48,0.58,0.28), e3:[0.22,0.18,0.42], texture:'speckle', color:'#f472b6'};
          return null; // isotropic handled separately as a sphere
        }
        const blockCv = par.querySelector('#esym-block');
        const BLOCK_SIZE = 260;
        const DPR = window.devicePixelRatio || 1;
        blockCv.style.width = BLOCK_SIZE+'px'; blockCv.style.height = BLOCK_SIZE+'px';
        blockCv.width = Math.round(BLOCK_SIZE*DPR); blockCv.height = Math.round(BLOCK_SIZE*DPR);
        const blockCtx = blockCv.getContext('2d');
        const BLOCK_SCALE = 90, BLOCK_CAMD = 4.2;

        function bilerp(p00,p10,p11,p01,u,v){
          return [ (1-u)*(1-v)*p00[0]+u*(1-v)*p10[0]+u*v*p11[0]+(1-u)*v*p01[0],
                   (1-u)*(1-v)*p00[1]+u*(1-v)*p10[1]+u*v*p11[1]+(1-u)*v*p01[1] ];
        }

        function drawBlock(){
          blockCtx.setTransform(DPR,0,0,DPR,0,0);
          blockCtx.clearRect(0,0,BLOCK_SIZE,BLOCK_SIZE);
          const cx=BLOCK_SIZE/2, cy=BLOCK_SIZE/2;
          const t = (state.opAnimStart===null) ? null : Math.min(1, (performance.now()-state.opAnimStart)/900);
          const op = state.activeOp;
          const pt = (state.probeAnimStart===null) ? null : (performance.now()-state.probeAnimStart)/950;
          function transform(p){ return rot3(applyStrainDeform(state.activeProbe, pt, applyOpAnim(op, t, p)), state.az, state.el); }

          const geom = blockGeom(state.cls);
          let faces = [];
          if(geom===null){
            // isotropic: a sphere, approximated as a shaded disk plus a handful of fixed surface
            // markers so it still visibly (if only faintly) rotates -- deliberately generic-looking
            // from every angle, since that IS the isotropic point.
            const R = 0.62;
            const seed = [[0.6,0.3,0.5],[-0.4,0.5,0.4],[0.2,-0.6,0.4],[-0.5,-0.3,-0.5],[0.5,-0.2,-0.5],[-0.2,0.6,-0.4]];
            const [_,__,zc] = rot3([0,0,1], state.az, state.el);
            const shade = 0.55 + 0.35*Math.max(0, zc);
            const [px,py] = proj([0,0,0], cx, cy, BLOCK_SCALE, BLOCK_CAMD);
            blockCtx.beginPath(); blockCtx.arc(px,py,R*BLOCK_SCALE,0,2*Math.PI);
            blockCtx.fillStyle = `rgba(167,139,250,\${shade})`; blockCtx.fill();
            blockCtx.strokeStyle = '#111827'; blockCtx.lineWidth = 1.5; blockCtx.stroke();
            seed.forEach(v=>{
              const n = Math.hypot(v[0],v[1],v[2]); const p = [v[0]/n*R, v[1]/n*R, v[2]/n*R];
              const rp = rot3(p, state.az, state.el);
              if(rp[2] < 0) return; // only draw markers on the near hemisphere
              const [mx,my] = proj(rp, cx, cy, BLOCK_SCALE, BLOCK_CAMD);
              blockCtx.beginPath(); blockCtx.arc(mx,my,2.5,0,2*Math.PI);
              blockCtx.fillStyle = '#4c1d95'; blockCtx.fill();
            });
            par.querySelector('#esym-block-caption').textContent = 'a sphere: identical from every angle';
            return;
          }

          const n = geom.bottom.length;
          const bottom3 = geom.bottom.map(([x,y])=>[x,y,-0.5*0]); // placeholder, replaced below
          const bot = geom.bottom.map(([x,y])=>[x,y,-geom.e3[2]/2 - Math.max(0,-geom.e3[2])]);
          // build bottom/top corners: bottom centered at z=-e3z/2 (offset by 0 in xy), top = bottom+e3
          const halfz = geom.e3[2]/2;
          const botC = geom.bottom.map(([x,y])=>[x,y,-halfz]);
          const topC = geom.bottom.map(([x,y])=>[x+geom.e3[0], y+geom.e3[1], -halfz+geom.e3[2]]);
          const botT = botC.map(p=>transform(p));
          const topT = topC.map(p=>transform(p));
          const botP = botT.map(p=>proj(p,cx,cy,BLOCK_SCALE,BLOCK_CAMD));
          const topP = topT.map(p=>proj(p,cx,cy,BLOCK_SCALE,BLOCK_CAMD));

          function faceDepth(idxs, pts3){ return idxs.reduce((s,i)=>s+pts3[i][2],0)/idxs.length; }
          function faceNormalShade(a,b,c){
            const u=[b[0]-a[0],b[1]-a[1],b[2]-a[2]], v=[c[0]-a[0],c[1]-a[1],c[2]-a[2]];
            const nx=u[1]*v[2]-u[2]*v[1], ny=u[2]*v[0]-u[0]*v[2], nz=u[0]*v[1]-u[1]*v[0];
            const len=Math.hypot(nx,ny,nz)||1;
            return Math.max(0.35, Math.min(1, 0.55+0.5*(nz/len)));
          }

          for(let i=0;i<n;i++){
            const j=(i+1)%n;
            const idx3 = [botT[i],botT[j],topT[j],topT[i]];
            const depth = (idx3[0][2]+idx3[1][2]+idx3[2][2]+idx3[3][2])/4;
            const shade = faceNormalShade(idx3[0],idx3[1],idx3[2]);
            faces.push({depth, pts:[botP[i],botP[j],topP[j],topP[i]], shade, kind:'side', i});
          }
          faces.push({depth: botT.reduce((s,p)=>s+p[2],0)/n, pts: botP.slice(), shade: faceNormalShade(botT[0],botT[2],botT[1]), kind:'bottom'});
          faces.push({depth: topT.reduce((s,p)=>s+p[2],0)/n, pts: topP.slice(), shade: faceNormalShade(topT[0],topT[1],topT[2]), kind:'top'});
          faces.sort((a,b)=>a.depth-b.depth);

          faces.forEach(f=>{
            blockCtx.beginPath();
            f.pts.forEach((p,k)=> k===0?blockCtx.moveTo(p[0],p[1]):blockCtx.lineTo(p[0],p[1]));
            blockCtx.closePath();
            const rgb = hexToRgb(geom.color);
            blockCtx.fillStyle = `rgb(\${Math.round(rgb[0]*f.shade)},\${Math.round(rgb[1]*f.shade)},\${Math.round(rgb[2]*f.shade)})`;
            blockCtx.fill();
            blockCtx.strokeStyle = '#111827'; blockCtx.lineWidth = 1; blockCtx.stroke();
            if(geom.texture==='lattice' && f.pts.length===4){
              blockCtx.save(); blockCtx.beginPath();
              f.pts.forEach((p,k)=> k===0?blockCtx.moveTo(p[0],p[1]):blockCtx.lineTo(p[0],p[1]));
              blockCtx.closePath(); blockCtx.clip();
              blockCtx.fillStyle = 'rgba(17,24,39,0.7)';
              for(let ui=1;ui<=3;ui++) for(let vi=1;vi<=3;vi++){
                const [dx,dy] = bilerp(f.pts[0],f.pts[1],f.pts[2],f.pts[3], ui/4, vi/4);
                blockCtx.beginPath(); blockCtx.arc(dx,dy,1.6,0,2*Math.PI); blockCtx.fill();
              }
              blockCtx.restore();
            } else if(geom.texture==='stripes' && f.kind==='side'){
              blockCtx.save(); blockCtx.beginPath();
              f.pts.forEach((p,k)=> k===0?blockCtx.moveTo(p[0],p[1]):blockCtx.lineTo(p[0],p[1]));
              blockCtx.closePath(); blockCtx.clip();
              blockCtx.strokeStyle = 'rgba(17,24,39,0.6)'; blockCtx.lineWidth = 1.3;
              [0.33,0.66].forEach(v=>{
                const p1 = bilerp(f.pts[0],f.pts[1],f.pts[2],f.pts[3], 0, v);
                const p2 = bilerp(f.pts[0],f.pts[1],f.pts[2],f.pts[3], 1, v);
                blockCtx.beginPath(); blockCtx.moveTo(p1[0],p1[1]); blockCtx.lineTo(p2[0],p2[1]); blockCtx.stroke();
              });
              blockCtx.restore();
            } else if(geom.texture==='speckle' && f.pts.length===4){
              blockCtx.save(); blockCtx.beginPath();
              f.pts.forEach((p,k)=> k===0?blockCtx.moveTo(p[0],p[1]):blockCtx.lineTo(p[0],p[1]));
              blockCtx.closePath(); blockCtx.clip();
              blockCtx.fillStyle = 'rgba(17,24,39,0.55)';
              const seedUV = [[0.2,0.3],[0.7,0.2],[0.4,0.7],[0.8,0.75],[0.15,0.8]];
              seedUV.forEach(([u,v])=>{
                const [dx,dy] = bilerp(f.pts[0],f.pts[1],f.pts[2],f.pts[3], u, v);
                blockCtx.beginPath(); blockCtx.arc(dx,dy,1.4,0,2*Math.PI); blockCtx.fill();
              });
              blockCtx.restore();
            }
          });

          // hovering a matrix cell always wins the caption -- it's a live, momentary probe, and
          // must not get stuck behind a symmetry-operation result that's still active from an
          // earlier click (activeOp is deliberately never cleared just by moving the mouse away).
          let capText = CLASS_LABEL[state.cls];
          if(state.hoveredCell){
            const [a,b] = state.hoveredCell;
            capText = 'apply strain ' + b + ' → feel stress ' + a;
          } else if(state.activeProbe!==null){
            capText = 'probing strain ' + state.activeProbe + ' — watch the block, and which row lights up in the matrix';
          } else if(op){
            // gated on opAnimFinished (a setTimeout, see the op-button handler below), not on the
            // animation's own requestAnimationFrame progress t -- rAF is throttled/paused on a
            // backgrounded tab, which would otherwise leave this stuck reading "testing..." forever.
            const valid = currentClassData().validOps[op];
            capText = state.opAnimFinished ? (OP_LABEL[op] + ': ' + (valid ? 'looks the same ✓' : 'looks different ✗')) : ('testing ' + OP_LABEL[op] + '…');
          }
          par.querySelector('#esym-block-caption').textContent = capText;

          // setTimeout, not requestAnimationFrame -- rAF is suspended on a backgrounded/hidden
          // tab, which would silently freeze this mid-animation; setTimeout keeps firing.
          if(t!==null && t<1 && !state.opAnimFinished){ setTimeout(drawBlock, 16); }
          else if(pt!==null && pt<1){ setTimeout(drawBlock, 16); }
        }
        function hexToRgb(h){ const n=parseInt(h.slice(1),16); return [(n>>16)&255,(n>>8)&255,n&255]; }

        // ---------- Panel 2: the Voigt matrix ----------
        const matCv = par.querySelector('#esym-matrix');
        const MAT_SIZE = 300;
        matCv.style.width = MAT_SIZE+'px'; matCv.style.height = MAT_SIZE+'px';
        matCv.width = Math.round(MAT_SIZE*DPR); matCv.height = Math.round(MAT_SIZE*DPR);
        const matCtx = matCv.getContext('2d');
        const MAT_ORIGIN = 30, CELL = 44;

        const FAMILY_COLORS = ['#f87171','#60a5fa','#34d399','#fbbf24','#a78bfa','#f472b6','#22d3ee'];

        function cellRole(cls, a, b){
          // a,b are 1-based Voigt indices; returns {kind, familyIdx}
          const cc = DATA.classes[cls].constraints;
          function pairIndex(a,b){
            for(let k=0;k<DATA.voigtPairs.length;k++){ const [p,q]=DATA.voigtPairs[k]; if((p===a&&q===b)) return k; }
            return -1;
          }
          const idx = pairIndex(Math.min(a,b), Math.max(a,b));
          if(idx<0) return {kind:'free'};
          if(cc.zeros.includes(idx)) return {kind:'zero'};
          for(let fi=0; fi<cc.families.length; fi++){ if(cc.families[fi].includes(idx)) return {kind:'family', fam:fi}; }
          for(const [pivot, terms] of cc.derived){ if(pivot===idx) return {kind:'derived'}; }
          return {kind:'free'};
        }

        function drawMatrix(){
          matCtx.setTransform(DPR,0,0,DPR,0,0);
          matCtx.clearRect(0,0,MAT_SIZE,MAT_SIZE);
          const C = currentClassData().C; // row-major 6x6
          matCtx.font = '11px sans-serif'; matCtx.textAlign='center'; matCtx.textBaseline='middle';
          for(let a=1;a<=6;a++){
            matCtx.fillStyle = '#9ca3af';
            matCtx.fillText(String(a), MAT_ORIGIN + (a-0.5)*CELL, 14);
            matCtx.fillText(String(a), 14, MAT_ORIGIN + (a-0.5)*CELL);
          }
          for(let a=1;a<=6;a++) for(let b=1;b<=6;b++){
            const v = C[(a-1)*6 + (b-1)];
            const role = cellRole(state.cls, a, b);
            let fill = '#1f2937';
            if(role.kind==='zero') fill = '#0b0b0b';
            else if(role.kind==='family') fill = FAMILY_COLORS[role.fam % FAMILY_COLORS.length] + '55';
            else if(role.kind==='derived') fill = '#78716c88';
            else fill = '#1e3a8a88';
            const x = MAT_ORIGIN + (b-1)*CELL, y = MAT_ORIGIN + (a-1)*CELL;
            matCtx.fillStyle = fill;
            matCtx.fillRect(x, y, CELL-2, CELL-2);
            if(state.activeProbe===b){
              // outgoing probe: the whole column lit up is "strain b feeds every one of these
              // stress rows" -- a stronger, non-zero cell in this column responds more.
              const pulse = 0.5 + 0.5*Math.sin(performance.now()/110);
              matCtx.fillStyle = `rgba(74,222,128,\${0.15+0.3*pulse})`;
              matCtx.fillRect(x, y, CELL-2, CELL-2);
            }
            if(state.activeOp!==null && DATA.ops[state.activeOp]){
              // incoming pulse: highlight cells THIS operation constrains, regardless of class --
              // persists as long as an op is selected (not just during the 900ms spin), so the
              // reveal is a stable, readable answer, not a flash the user has to catch in time.
              function pairIndex(p,q){ for(let k=0;k<DATA.voigtPairs.length;k++){ const [pp,qq]=DATA.voigtPairs[k]; if(pp===p&&qq===q) return k;} return -1; }
              const idx = pairIndex(Math.min(a,b), Math.max(a,b));
              const oc = DATA.ops[state.activeOp].constraints;
              const implicated = idx>=0 && (oc.zeros.includes(idx) || oc.families.some(f=>f.includes(idx)) || oc.derived.some(([p,t])=>p===idx));
              if(implicated){
                const pulse = 0.5 + 0.5*Math.sin(performance.now()/120);
                matCtx.fillStyle = `rgba(250,204,21,\${0.25+0.35*pulse})`;
                matCtx.fillRect(x, y, CELL-2, CELL-2);
              }
            }
            if(state.hoveredCell && ((state.hoveredCell[0]===a&&state.hoveredCell[1]===b)||(state.hoveredCell[0]===b&&state.hoveredCell[1]===a))){
              matCtx.strokeStyle = '#facc15'; matCtx.lineWidth = 2.5; matCtx.strokeRect(x+1, y+1, CELL-4, CELL-4);
            } else {
              matCtx.strokeStyle = '#374151'; matCtx.lineWidth = 1; matCtx.strokeRect(x, y, CELL-2, CELL-2);
            }
            matCtx.fillStyle = '#e5e7eb'; matCtx.font = '11px sans-serif';
            matCtx.fillText(Math.abs(v)<0.05 ? '0' : v.toFixed(0), x+(CELL-2)/2, y+(CELL-2)/2);
          }
          par.querySelector('#esym-count').textContent = currentClassData().independentCount;
          if(state.activeOp!==null || state.activeProbe!==null) setTimeout(drawMatrix, 16);
        }

        matCv.addEventListener('mousemove', e=>{
          const rect = matCv.getBoundingClientRect();
          const mx = e.clientX-rect.left, my = e.clientY-rect.top;
          const b = Math.floor((mx-MAT_ORIGIN)/CELL)+1, a = Math.floor((my-MAT_ORIGIN)/CELL)+1;
          if(a>=1&&a<=6&&b>=1&&b<=6){ state.hoveredCell=[a,b]; } else { state.hoveredCell=null; }
          drawMatrix(); drawBlock();
        });
        matCv.addEventListener('mouseleave', ()=>{ state.hoveredCell=null; drawMatrix(); drawBlock(); });
        // click a cell to actually PROBE it: deform the block along that strain column and light
        // up every stress row it could feed (a click is a deliberate "do it", unlike a hover)
        matCv.addEventListener('click', e=>{
          const rect = matCv.getBoundingClientRect();
          const mx = e.clientX-rect.left, my = e.clientY-rect.top;
          const b = Math.floor((mx-MAT_ORIGIN)/CELL)+1;
          if(b<1||b>6) return;
          state.activeOp = null; state.opAnimStart = null; state.opAnimFinished = false;
          opsGrid.querySelectorAll('.esym-opbtn').forEach(btn=>btn.classList.remove('valid','invalid'));
          state.activeProbe = b; state.probeAnimStart = performance.now();
          drawBlock(); drawMatrix();
          setTimeout(()=>{ state.activeProbe = null; state.probeAnimStart = null; drawBlock(); drawMatrix(); }, 970);
        });

        // ---------- Panel 3: directional stiffness surface ----------
        const surfCv = par.querySelector('#esym-surface');
        const SURF_SIZE = 260;
        surfCv.style.width = SURF_SIZE+'px'; surfCv.style.height = SURF_SIZE+'px';
        surfCv.width = Math.round(SURF_SIZE*DPR); surfCv.height = Math.round(SURF_SIZE*DPR);
        const surfCtx = surfCv.getContext('2d');
        const SURF_SCALE = 90, SURF_CAMD = 4.2, SURF_R = 0.85;

        function gridPoint(ntheta, nphi, i, j, radius){
          const theta = Math.PI*i/(ntheta-1), phi = 2*Math.PI*j/nphi;
          return [radius*Math.sin(theta)*Math.cos(phi), radius*Math.sin(theta)*Math.sin(phi), radius*Math.cos(theta)];
        }

        function drawWireGrid(values, ntheta, nphi, vmax, color){
          const cx=SURF_SIZE/2, cy=SURF_SIZE/2;
          const pts = [];
          for(let i=0;i<ntheta;i++){ pts.push([]); for(let j=0;j<nphi;j++){
            const v = values[i*nphi+j];
            const r = SURF_R * (v/vmax);
            const p3 = gridPoint(ntheta, nphi, i, j, r);
            pts[i].push(proj(rot3(p3, state.az, state.el), cx, cy, SURF_SCALE, SURF_CAMD));
          }}
          surfCtx.strokeStyle = color; surfCtx.lineWidth = 1;
          for(let i=0;i<ntheta;i++){ surfCtx.beginPath();
            for(let j=0;j<=nphi;j++){ const p=pts[i][j%nphi]; j===0?surfCtx.moveTo(p[0],p[1]):surfCtx.lineTo(p[0],p[1]); }
            surfCtx.stroke();
          }
          for(let j=0;j<nphi;j++){ surfCtx.beginPath();
            for(let i=0;i<ntheta;i++){ const p=pts[i][j]; i===0?surfCtx.moveTo(p[0],p[1]):surfCtx.lineTo(p[0],p[1]); }
            surfCtx.stroke();
          }
        }

        function drawSurface(){
          surfCtx.setTransform(DPR,0,0,DPR,0,0);
          surfCtx.clearRect(0,0,SURF_SIZE,SURF_SIZE);
          const cd = currentClassData();
          const ntheta = DATA.ntheta, nphi = DATA.nphi;
          if(state.surfaceMode==='modulus'){
            const vmax = Math.max(...cd.youngsGrid);
            drawWireGrid(cd.youngsGrid, ntheta, nphi, vmax, '#38bdf8');
            const vmin = Math.min(...cd.youngsGrid);
            par.querySelector('#esym-surface-caption').textContent =
              'Young\\'s modulus: ' + vmin.toFixed(0) + '–' + vmax.toFixed(0) + ' GPa' + (vmax/vmin>1.02 ? ' (anisotropic)' : ' (isotropic: a sphere)');
          } else {
            const branch = (k)=> cd.waveGrid.filter((_,idx)=>idx%3===k);
            const slow = branch(0), mid = branch(1), fast = branch(2);
            const vmax = Math.max(...cd.waveGrid);
            drawWireGrid(fast, ntheta, nphi, vmax, '#f87171');
            drawWireGrid(mid, ntheta, nphi, vmax, '#34d399');
            drawWireGrid(slow, ntheta, nphi, vmax, '#60a5fa');
            const split = (Math.max(...mid)-Math.min(...slow));
            par.querySelector('#esym-surface-caption').textContent =
              'qP (red) · qS1 (green) · qS2 (blue) — shear splitting is the two inner surfaces pulling apart';
          }
        }

        // ---------- shared drag-to-orbit (Panel 1 + Panel 3) ----------
        function wireOrbit(cv){
          let dragging=false, last=null;
          cv.addEventListener('mousedown', e=>{ dragging=true; last=[e.clientX,e.clientY]; });
          window.addEventListener('mousemove', e=>{
            if(!dragging) return;
            const dx=e.clientX-last[0], dy=e.clientY-last[1]; last=[e.clientX,e.clientY];
            state.az += dx*0.01; state.el = Math.max(-1.4, Math.min(1.4, state.el - dy*0.01));
            drawBlock(); drawSurface();
          });
          window.addEventListener('mouseup', ()=>{ dragging=false; });
        }
        wireOrbit(blockCv); wireOrbit(surfCv);

        // ---------- class picker ----------
        const classRow = par.querySelector('#esym-classrow');
        CLASS_ORDER.forEach(cls=>{
          const btn = document.createElement('button');
          btn.className = 'esym-classbtn' + (cls===state.cls ? ' active' : '');
          btn.textContent = CLASS_LABEL[cls] + ' (' + DATA.classes[cls].independentCount + ')';
          btn.dataset.cls = cls;
          btn.addEventListener('click', ()=>{
            state.cls = cls; state.activeOp = null; state.opAnimStart = null; state.opAnimFinished = false;
            state.activeProbe = null; state.probeAnimStart = null;
            classRow.querySelectorAll('.esym-classbtn').forEach(b=>b.classList.toggle('active', b.dataset.cls===cls));
            opsGrid.querySelectorAll('.esym-opbtn').forEach(b=>{ b.classList.remove('valid','invalid'); });
            drawBlock(); drawMatrix(); drawSurface();
          });
          classRow.appendChild(btn);
        });

        // ---------- symmetry-operation buttons ----------
        const opsGrid = par.querySelector('#esym-ops-grid');
        OP_ORDER.forEach(op=>{
          const btn = document.createElement('button');
          btn.className = 'esym-opbtn';
          btn.textContent = OP_LABEL[op];
          btn.addEventListener('click', ()=>{
            opsGrid.querySelectorAll('.esym-opbtn').forEach(b=>{ b.classList.remove('valid','invalid'); });
            state.activeProbe = null; state.probeAnimStart = null;
            state.activeOp = op; state.opAnimStart = performance.now(); state.opAnimFinished = false;
            drawBlock(); drawMatrix();
            // a setTimeout, not a continuation of drawBlock's own requestAnimationFrame loop --
            // rAF is throttled/paused on a backgrounded tab, but the final valid/invalid verdict
            // must still resolve even if the user switched away mid-animation.
            setTimeout(()=>{
              state.opAnimFinished = true;
              const valid = currentClassData().validOps[op];
              btn.classList.toggle('valid', valid); btn.classList.toggle('invalid', !valid);
              drawBlock();
            }, 950);
          });
          opsGrid.appendChild(btn);
        });

        // ---------- modulus / wave-speed toggle ----------
        const toggleRow = par.querySelector('#esym-toggle-row');
        [['modulus','Modulus'],['wave','Wave speed']].forEach(([mode,label])=>{
          const btn = document.createElement('button');
          btn.className = 'esym-togglebtn' + (mode===state.surfaceMode ? ' active' : '');
          btn.textContent = label;
          btn.addEventListener('click', ()=>{
            state.surfaceMode = mode;
            toggleRow.querySelectorAll('.esym-togglebtn').forEach(b=>b.classList.toggle('active', b.textContent===label));
            drawSurface();
          });
          toggleRow.appendChild(btn);
        });

        drawBlock(); drawMatrix(); drawSurface();
        }
        </script>
        """)
    end
end

# ╔═╡ 9f76c4ba-a940-41f9-b352-3509fc5382b6
ElasticSymmetryInput(ESYM_DATA_JSON)

# ╔═╡ 49b426a4-700d-442e-aea3-a1c1a605e82e
md"""
### References

- Aki, K. & Richards, P.G. (2002), *Quantitative Seismology*, 2nd ed. — Voigt notation, the Bond
  transform, and the Christoffel equation.
- The VTI defaults and olivine Voigt matrix are reused from this repo's own
  `src/Planewave Propagation/anisotropy.jl`.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "9bbb622e4cd9995606b539e0b3f2495d359cd8e4"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

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

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

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
version = "3.5.4+0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

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
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"
"""

# ╔═╡ Cell order:
# ╠═160890b8-404a-4c53-a4d7-06089d85c09d
# ╟─dff060e3-b550-4a33-bc48-465a10f4451c
# ╟─4fcde378-5707-4ee5-9dce-9a9ec5ac78d2
# ╠═9f76c4ba-a940-41f9-b352-3509fc5382b6
# ╟─729629f2-70b5-4a34-b190-0fa86e5f70e6
# ╟─c7b1cf0e-0733-480f-9d46-8c96d3449e26
# ╟─c5e30d76-f237-421d-8174-35cad72344fc
# ╠═3c8a53dd-00c3-4a79-b78e-be46713efc12
# ╠═b95dce52-53ec-4aae-84f1-a41ac1c45c0e
# ╠═4479f375-688e-4272-9129-c1861f4fb0f9
# ╟─a3513ce2-b9af-42c8-a71a-9ed8afe47487
# ╠═22c406f2-4061-43a0-8d22-0a934609c4d1
# ╟─5f4c619b-a736-4909-a487-1c32de318c05
# ╠═84810524-9c74-40ce-a1f6-54a890960903
# ╠═f9d818c3-0022-4ae0-9e87-5548dc303905
# ╠═5ceab72b-2c3e-42ba-bc16-d7d84c633b41
# ╠═3a98ea1b-100a-4ce8-a205-916ed6be814e
# ╟─5c84645b-c28d-46b5-9a08-108943e6368c
# ╠═bae63dd3-86e2-4a29-99fb-e6acb4162039
# ╟─1812ef55-3447-4d25-a352-c570e34b218e
# ╠═c2f0bd3c-1e5d-400a-819e-d202c8b99620
# ╠═6274320f-24f8-42b2-b72a-26cfe1a37876
# ╠═d36b4b9d-811d-476f-8394-b08fd8c097ce
# ╠═29dcbb36-8d47-442c-b3af-9676d4274c36
# ╟─0d3beab8-5a9c-46cd-8c01-149d5d34c00f
# ╠═8843c9fd-85ff-4e01-a6fe-1bcf5b8088a2
# ╠═37455fad-da49-4de4-87b7-ba03acd8da5d
# ╠═af2515da-bd17-4a9d-b0a0-3eef9661f67f
# ╠═a64f8607-3c0c-4c81-bdf1-bf6d4b567118
# ╠═9945b72f-9595-47fc-bf2e-2dcc9eda8451
# ╠═a030eb8c-2050-4640-9bf8-39e34a95c736
# ╠═8908d7fe-0b41-49e0-a033-995607e41749
# ╠═f2b0023e-d0a6-4890-a2ab-85e4ec7e4112
# ╠═e3a630d5-a7b1-45d1-b0ed-b47646ef7d87
# ╠═b57f45f1-9f49-4568-8b23-46c818d4bcf8
# ╟─fc3ed6e9-19bd-4170-b90a-b3c013b71c93
# ╠═c4b0936b-3db1-4a3f-a764-f9cf0ee4934b
# ╟─6cffaa93-19f8-483c-a780-ce8cd54b33ed
# ╠═3c043dec-4383-4f2f-b77e-4203d2dd734a
# ╠═08baab70-8994-4e90-9a91-98d2308a9418
# ╠═ec8260cc-96c3-43c9-8b43-c5ee13a0fa93
# ╠═c71483d4-6ee8-4dc3-9cd9-dcabdc27941c
# ╠═288d863f-fc57-48e0-89d5-bed10de6506a
# ╟─214a6443-67b8-40e1-90f7-c5d771b79694
# ╠═383498a0-b2cc-4deb-9b6a-86eef802049f
# ╟─982529a6-007c-4052-9e48-c664e7915dc0
# ╠═18df75f7-0850-49e0-9838-f629c101fe5c
# ╟─a39cb03c-b7e0-40e4-a2d4-714388e580fe
# ╠═096270f3-5ce6-4a45-bd9b-2f8114cabc51
# ╟─135da962-7b36-4b0c-9c68-f2313cfb65d7
# ╠═fc1ed39a-818d-4ff5-8e82-dd4c1e028acd
# ╠═b2faf56b-88f6-4b70-b5b5-4129a28511bc
# ╟─ad3d84e0-b1fa-4d46-9213-cddb27af546b
# ╠═6884b812-4802-49ba-a38e-0446112d3f00
# ╟─cfd61b4e-5a35-4c55-957f-3e2c7519968e
# ╠═a59ccaaf-9d1f-40a1-8e85-dd737feca733
# ╠═3828fd1e-a4de-4f72-8cfd-115be3f731e5
# ╠═269f858c-27a3-447d-b4dd-5a8a839b2bb6
# ╠═46fd8d6f-3fa2-45f2-9950-6232f85c48e6
# ╠═8a610175-37ef-4a3f-a919-8f9cfbcb6759
# ╠═91872170-0f34-4f3e-9aed-4d2335024f7b
# ╠═793cf760-0258-44c7-8cc9-f426d68e7634
# ╠═188ae616-5423-484c-b22d-621e52e529a0
# ╠═2787de1b-57cd-4567-b46e-aeff7b68b493
# ╟─07d767d0-64cd-4ecb-9a1c-7126c483c67c
# ╠═9fbf8b5e-ec64-4e5f-ba44-b29ce5870efc
# ╟─ddaa44f6-9462-4238-a06b-fc2a26a3e06c
# ╠═2e481c84-d654-4c7d-a8c6-e0ff3a0c56d3
# ╟─49b426a4-700d-442e-aea3-a1c1a605e82e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
