---
name: pluto-notebook-writing
description: "Content and prose conventions for this repo's Pluto notebooks (reference: src/misc/geoid-kernel.jl): how the notebook is organized around its interactive widget (intro before, intuition after, all Julia implementation pushed to a categorized Appendix at the end), how functions are documented with docstrings, and how to write LaTeX math inside md strings without hitting Julia's string-interpolation syntax. Use this whenever writing or restructuring the explanatory text, section order, or function documentation of a notebook under src/ -- especially when adding math, writing a docstring, or reorganizing a notebook to put an interactive widget in the middle of a narrative. For the widget's own HTML/CSS/JS, use the pluto-widget-style skill instead -- this one is about everything else in the notebook."
---

# Pluto notebook writing conventions

This is the companion to `pluto-widget-style` (which covers the interactive canvas widget's HTML/CSS/JS).
This skill covers everything else in the notebook: how the prose and the widget relate to each other, how
Julia functions are documented, and how to write LaTeX math that actually survives being embedded inside a
Julia string. `src/misc/geoid-kernel.jl` is the reference — go read the actual notebook when this doc
doesn't cover your specific case.

## Repository workflow: experimental by default

All new Pluto notebooks should be created in the staging area under `src/experimental/` unless the author
has already decided which topical section the notebook belongs to.

- `src/experimental/` is the default home for draft, exploratory, and candidate notebooks.
- Once a notebook is mature enough to belong to a topic, move it into the relevant folder such as
  `src/Elasticity and Deformation/`, `src/Planewave Propagation/`, `src/Surface Waves/`, or
  `src/Earth Structure and Internal Layers/`.
- The folder name is the section title that should appear in the site navigation and in the notebook
  registry, so the folder is the canonical source of the section label.
- Do not place a new notebook directly in a production section folder unless the notebook is already
  ready to be included in the site.

In other words: new work starts in experimental, and the user moves it as necessary into the section that
matches the underlying physics. The system should not auto-deploy every file in the repo; deployment is an
explicit editorial choice, not a default behavior.

## Deployment policy: explicit live vs static selection

The site should use the section folders as the organizational structure, but the actual inclusion of a
notebook in a live or static build must be decided by an explicit allowlist.

Recommended mechanism:

- Each section folder defines the section title.
- Each notebook is added to either the live registry or the static registry by explicit choice.
- The registry should be written as a section -> lists structure, not as a flat list of files.

Example pattern:

```yaml
sections:
  Elasticity and Deformation:
    live:
      - Elasticity and Deformation/strain-tensor.jl
      - Elasticity and Deformation/stress-tensor.jl
    static:
      - Elasticity and Deformation/simple-harmonic-motion.jl

  Planewave Propagation:
    live:
      - Planewave Propagation/reflection-transmission-SH.jl
    static:
      - Planewave Propagation/intrinsic-attenuation.jl
```

This gives the build clear rules:

- a notebook must be in the registry to appear in the site at all;
- it must be listed under a section folder to get that section name;
- it must be placed in the unified `live-notebooks.yml` under either the `live:` or `static:` list for its section.

This is the correct mechanism for a teaching site: the folder names provide the taxonomy, while the registry
file decides deployment status. The repository should not infer live/static status from file location alone,
because that is too rigid and would force the wrong notebooks onto the live server or onto the static site.

The repo now reads the unified sectioned registry from `live-notebooks.yml`, preserving section names and
flattening the nested `live:` / `static:` lists as it builds the final notebook list.

## The shape of a notebook: narrative around the widget, implementation in an Appendix

The read-through experience and the file's cell order are two different things in Pluto, and that's what
makes this structure possible: **Pluto resolves cells by dependency graph, not by their physical position
in the file** — see this repo's `CLAUDE.md`/`AGENTS.md` on the `# ╔═╡ Cell order:` block being what
controls *display* order, independent of where a cell's `# ╔═╡ <uuid>` definition physically sits. That
means the code that *defines* the widget (the `struct`, `Base.get`, `Base.show` triad) can live at the very
bottom of the file while the cell that *instantiates and displays* it (`@bind gk_blobs GeoidGlobeInput(...)`)
sits near the top of the reading order — Pluto doesn't care that the definition comes later in the file,
only that it exists somewhere. Use this to get the read-through order you actually want:

1. **Introduction** — why this question, what the reader is about to explore, pitched *before* they've
   touched anything. This is where you build the motivation (see geoid-kernel's opening "Why Does a Dense
   Blob Make a Geoid *Low*?" — it states the naive-but-wrong intuition first, then says why the notebook
   exists to correct it).
2. **The interactive widget**, displayed here even though its code lives at the bottom (see above).
3. **Detailed explanation and intuition**, written to build on what the reader just did with their hands —
   not a repeat of the introduction. This is where the real derivations, the governing equations, and the
   "aha" sections live (geoid-kernel's "Where the Geoid Comes From", "The Governing Equations", "Aha #3:
   Viscosity Controls the Flip", etc.) — each one can reference specific things the reader just painted or
   dragged.
4. **`## Appendix`** — every Julia function implementation, in display order, grouped into named
   subsections that tell the story of the computational pipeline rather than generic buckets. Don't call a
   section "Helper functions" — call it what it *does*: geoid-kernel uses `## Physical Constants`,
   `## Layer 1: The Analytic Building Blocks`, `## The Geoid Kernel`, `## Layer 3: From a Painted Picture
   to a Geoid`. Each category can (and in geoid-kernel does) end in its own small `### Verifying the ...`
   / `### Validating the ...` subsection with a visible numerical self-check — printed residuals, an
   assertion, a comparison table — so trust in the implementation is built the same way trust in the
   physics was built in the narrative above: by showing the check, not asserting correctness.
5. **The widget's own definition** (`struct`/`Base.get`/`Base.show`) is the **last thing in the Appendix**
   — after every other category of function. It's still just another subsection (e.g. `## The Interactive
   Widget`), it just comes last because everything else in the Appendix is a dependency of it or of the
   downstream cells that consume its bound value.

Don't force this order if a notebook doesn't have a widget, or has the widget as the *entire point* with
little supporting derivation — this is the shape for a notebook built around one hero interactive demo
with real physics/math underneath it, not a rule for every notebook in the repo.

## Documenting functions

Give every function that does real work a docstring — skip it only for genuine one-liners
(`Lof(ℓ) = ℓ * (ℓ + 1)`) where the name already says everything. The shape, straight from geoid-kernel:

```julia
"""
	geoid_response(s, ℓ, layers)

Solve the Stokes problem for a unit density sheet at radius `s`, degree `ℓ`.

Returns a named tuple with the three geoid contributions (m per kg/m³ per m of sheet
thickness), the surface and CMB dynamic topography, and the boundary-condition residuals
used for validation.
"""
function geoid_response(s, ℓ, layers)
    ...
```

- Tab-indented signature line, blank line, then prose — no `# Arguments` / `# Returns` header ceremony,
  just sentences. Reference parameters with single backticks (`` `ℓ` ``, `` `layers` ``).
- Say **what it computes physically and why**, not what the code does line-by-line — the code already
  says that. "Solve the Stokes problem for a unit density sheet" earns its place; "loops over layers and
  multiplies matrices" doesn't.
- If the return value's shape or meaning isn't obvious from the signature, spell it out in a `Returns ...`
  sentence.
- Cross-reference a related function with `` [`other_function`](@ref) `` when the reader would benefit
  from knowing it's the same computation viewed differently (geoid-kernel does this between
  `flow_profile` and `geoid_response`).
- Admonitions (see below) can live *inside* a docstring, not just in prose cells — geoid-kernel uses
  `!!! note "Sign convention"` inside `flow_profile`'s docstring to flag a convention that would otherwise
  cause a silent sign bug for the next person who touches that function.
- Inside the function body, use `#` comments for the *why* of a non-obvious step (a physical convention
  choice, a magic number's origin) — this is the same rule as the rest of this repo's Julia code, not
  specific to notebooks.

## The physics lives in Julia — always, no exceptions for "it's just the widget"

Every widget in this repo pairs a Julia physics layer with a JS display layer, and the line between
them is not negotiable: **any actual math — a formula, a derivative, a matrix product, a coordinate
transform that depends on the physics rather than just the camera — is a Julia function, full stop.**
JS only draws what Julia already computed. This holds even when the temptation to shortcut is strongest:
a slider drag that needs to feel instant, a per-direction sample evaluated hundreds of times a frame, a
"simple" trig formula that seems easier to just paste into the `<script>` block. None of that is an
exception. If a live interaction needs the physics to update as the user drags, round-trip it through
Julia with the same `commitInFlight`/`throttledCommit` pattern already used for camera drags — do not
reach for a client-side reimplementation to dodge the latency.

The one narrow carve-out (see `pluto-widget-style`'s "Marker shapes" and camera-drag sections) is pure
*display/camera geometry* that has no physical content of its own: which screen pixel a 3D point
projects to, which direction the mouse is dragging, how to rotate a view, where to place a text label.
That code has no formula a reader would want to see or verify against a textbook; it's plumbing. The
test is not "is this JS short/simple?" — it's "does this JS encode a physical relationship?" A rotation
matrix for the camera is plumbing. A dot product between a moment tensor and a direction vector is
physics, even though it's one line of JS and would be trivial to inline. When in doubt, put it in
Julia; a docstring costs nothing and an inlined formula in JS is exactly the kind of thing that goes
unverified and quietly drifts wrong (see `lame-theorem.jl`'s far-field `` X_0 ``-vs-`` \dot X_0 ``
bug, caught only because the physics lived in one documented, testable function instead of being
duplicated into a canvas script).

Push per-direction physics results as flat arrays via the same `CustomEvent` pattern used everywhere
else in this repo (`FieldPush`, `PfrPush`, `RayPush`, ...): Julia evaluates the quantity on a
*deterministically reproducible* direction set (a Fibonacci sphere, a `(θ,φ)` grid — see
[`fibonacci_sphere_directions`](@ref)/[`latlon_grid_directions`](@ref) in `lame-theorem.jl` for the
canonical pattern), JS regenerates the *identical* geometry independently (safe, since it's pure
geometry with no physics dependence) and zips it by index with the pushed array. This keeps the two
sides in lockstep without ever sending the direction vectors themselves over the wire.

## Performance: type stability from the `@bind` boundary down

Sources: [Modern Julia Workflows — Optimizing](https://modernjuliaworkflows.org/optimizing/),
[viralinstruction — How to optimise Julia code](https://viralinstruction.com/posts/optimise/), the
[official Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/). Full tip
lists there; this section keeps only what actually recurs in this repo's widget/Appendix pattern. Two
rules underlie almost everything below: **the compiler must be able to infer a concrete type for every
variable**, and **avoid heap allocation you don't need** (Julia's GC pauses the whole program when it
runs). Everything else is a specific way to violate — or respect — one of those two.

### The `@bind` boundary is the one place type instability is unavoidable — cut it off immediately

A widget's bound value always arrives as a loosely-typed `Dict{String,Any}` (see
`pluto_bond_value_type` in project memory), and a JS number that happens to be whole can deserialize
as `Int64` on one interaction and `Float64` on the next. This is exactly the Julia manual's "annotate
values from untyped locations" case (`x = a[1]::Int32`), and it's already this repo's established
pattern — `Born-approximation.jl`'s fallback/coercion cell is the model:

```julia
bs_safe = bs isa AbstractDict ? bs : Dict{String,Any}("vp0" => 2500.0, ...)
_bs_viewSrc = clamp(round(Int, bs_safe["viewSrc"]), 1, NSRC)   # Julia refuses to index with a Float64
_bs_dm = permutedims(reshape(Float64.(bs_safe["pert"]), NX, NZ)) ./ 100 ./ bs_safe["vp0"]^2
```

**Rule:** coerce every field out of the bond dict to its real type (`Float64.(...)`, `round(Int, ...)`)
in the cell that first reads it, then never touch the raw `Any`-typed dict again downstream — every
cell after that boundary should be working with plain, concretely-typed values. This is also why the
untyped "glue" cell (`bs_safe = bs isa AbstractDict ? ... : ...`) should do *only* the isa-check and
type coercions and immediately hand off to a real function with concrete argument types
(`acoustic_medium(bs_safe["vp0"], bs_safe["fpeak"], xgrid, zgrid)`) — this is Julia's **function-
barrier** technique: the compiler specializes at the function boundary, so pushing the
still-uncertain-until-just-now values through one function call gets you a fully type-stable,
optimized method body on the other side, instead of every downstream cell re-inferring through a
chain of `Any` lookups.

### Type-stable functions: one return type, no mid-function type changes

A function's return type shouldn't depend on the *value* of its input, and a variable shouldn't change
concrete type partway through a loop — both force the compiler to generate code that handles multiple
possible types at every subsequent use, which is exactly what boxing/dynamic dispatch is.

```julia
# unstable -- returns Int when x<0, else whatever typeof(x) is
pos(x) = x < 0 ? 0 : x
# stable -- always returns typeof(x)
pos(x) = x < 0 ? zero(x) : x

# unstable -- x starts Int64, becomes Float64 after the first /=
function foo()
    x = 1
    for i in 1:10; x /= rand(); end
    return x
end
# stable -- pick the real type up front
function foo()
    x = 1.0
    for i in 1:10; x /= rand(); end
    return x
end
```

### Preallocate, broadcast-fuse, and prefer views over copies in hot paths

`.`-broadcasting nests into one fused loop with no temporary arrays — `@. 3x^2 + 4x + 7x^3` is one
pass; `3x.^2 .+ 4x .+ 7x.^3` without the fuse macro can allocate a temporary per operation. Where a
function builds output arrays repeatedly (once per mouse-move, once per animation frame), prefer an
in-place `!`-convention version that writes into a caller-supplied buffer over allocating fresh each
call — the same instinct that already drives this repo's `commitInFlight`/`throttledCommit` pattern
(don't do more work than the interaction actually needs) applies to allocation too. Use `@view`/`@views`
instead of `x[a:b]` when slicing a large array just to read from it — a view costs nothing to create,
a slice copies.

### Loop over multidimensional arrays in column-major order

Julia arrays are stored column-major, so the innermost loop index should be the *first* index of the
array, matching how this repo already builds coordinate grids
(`[z for z in zgridM, x in xgridM]` in `Born-approximation.jl` — `z` varies fastest because it's the
first/outer comprehension variable, matching `vec`'s own column-major flatten). Getting this backwards
(looping `x` innermost against a `(z,x)`-shaped array) doesn't error, it just strides through memory
out of order and quietly costs cache misses.

### Benchmark with `BenchmarkTools`, never trust a bare `@time` on a global

`@time` at the REPL on a variable bound outside a function measures allocation *and* the cost of
looking up a non-`const` global — not the function's real cost. For anything under ~1 ms, use
`@btime`/`@benchmark` from BenchmarkTools.jl with `$` to interpolate external values:

```julia
using BenchmarkTools
@btime get_forward_operator($_bs_pa, $_bs_acq, 400.0, -15.0)   # $ avoids re-measuring global lookup
```

Do this as a **before/after** pair around any change made specifically for performance — "it feels
faster" isn't a measurement, and constant folding can make an un-interpolated benchmark of a literal
expression report an impossible near-0 ns.

### `@code_warntype` as a development-time check, same spirit as the numeric self-checks

For any Appendix function that's genuinely hot (called every mouse-move, every animation frame, or
inside a loop over the painting grid), run `@code_warntype the_function(args...)` once while writing
it and look for red/`Any`/`Union` in the output — this is the performance analogue of this skill's
`### Verifying ...` numeric self-checks: a habit for the person writing the function, not something
that becomes permanent notebook content. `JET.jl`'s `@report_opt` catches instability deeper in a call
chain than `@code_warntype` alone if the immediate function looks clean but something it calls isn't.

### `@inbounds`/`@simd` are narrow tools, not defaults

Reach for these only after profiling has identified a specific tight loop as the actual bottleneck, and
only when array indices are *provably* always in range (`@inbounds`) or loop iterations are truly
order-independent (`@simd`, which explicitly permits floating-point reassociation). Sprinkling either
across a whole notebook on the assumption it "can't hurt" trades a real correctness guarantee
(bounds checking) for a speedup that, on the array sizes typical of a widget redraw (hundreds to a few
thousand grid cells), is usually not the dominant cost — type stability and allocation almost always
are.

## LaTeX math — and the one thing that will silently break it

**Never write inline math as `$...$` or display math as `$$...$$`.** That's the Jupyter/generic-Markdown
convention, and it doesn't apply here for a very concrete reason: `$` inside *any* Julia string is
string-interpolation syntax, not a math delimiter. Writing `$x^2$` inside an `md"""..."""` block either
errors (if `x` isn't a variable in scope) or silently interpolates the wrong thing — it will not render as
math either way.

Julia's Markdown (what Pluto renders) uses different delimiters entirely:
- **Inline math**: double backticks, `` `\ell = 2` `` → renders as `` `\ell = 2` ``.
- **Display math**: a fenced block tagged `math`:
  ````
  ```math
  \nabla \cdot \mathbf{v} = 0
  ```
  ````

Neither of these uses `$`, so there's no interpolation conflict — but there's a second, easy-to-miss gotcha
about backslashes, and it depends on **what kind of string the math is sitting inside**:

- **Directly inside an `md"""..."""` block** (the normal case — most prose lives here): use a **single**
  backslash. `md"""..."""` is parsed by the `@md_str` macro, which does not interpret backslash escape
  sequences, so `\nabla`, `\ell`, `\theta` come through to the Markdown renderer exactly as typed.
- **Inside any *regular* Julia string that gets spliced into an md block via `$(...)` interpolation** —
  e.g. building a row of a table dynamically, or any string literal that isn't itself `md"""..."""` (this
  includes a function **docstring**, which is a plain `"""..."""` string, not an `md"""..."""` one) — use a
  **double** backslash. Regular strings *do* process escape sequences, so `\ell` there would actually try
  to escape the letter `l` (and either error or produce a mangled string); `\\ell` is what survives that
  string's own parsing to reach the Markdown renderer as `\ell`.

Rule of thumb: **if the text is inside something that itself accepts `$(...)` interpolation, double the
backslash; if it's static prose sitting directly in `md"""`, single is correct.** A docstring counts as
"accepts interpolation" even though you're not actually interpolating anything in it — it's still a plain
string, so the escaping rule is about the *string type*, not whether you happen to use interpolation.

```julia
# Directly in md""" — single backslash:
md"""
The geoid is ``N_\ell = \frac{4\pi G a}{g(2\ell+1)}\,\sigma\,(s/a)^{\ell+2}``.
"""

# Inside a docstring (a plain string) — double backslash:
"""
	mode_state(r, η, ℓ, p, cU, cP)

State vector ``y = [u_r,\\; u_\\theta,\\; \\tau_{rr},\\; \\tau_{r\\theta}]`` for a single mode.
"""
function mode_state(...)

# Building a table row dynamically, then interpolating it into md""" — double backslash
# in the inner string, because that inner string is a regular "..." literal:
"| ``\\vert ②/① \\vert`` | net geoid |\n|---|---:|:---|\n" * join(rows, "\n")
```

If a rendered equation shows up as a literal backslash-letter (`\ell` printed as text, not as ℓ) or Pluto
throws an `UndefVarError`/interpolation error on a cell that's pure prose, this escaping mismatch is the
first thing to check.

## Building prose from a dynamically-constructed string: two more traps

The pattern "compute a readout string in Julia, then hand it to Markdown" (e.g. a live summary line
that includes `bs_safe["vp0"]`, a source count, a view index — anything with several `$(...)`
interpolations packed together) has two failure modes beyond the escaping rule above, both silent
(no error, no exception — the cell just renders wrong) and both hit for real building
`Born-approximation.jl`'s summary line:

1. **Several tightly-packed `$(...)` interpolations directly inside `md"""..."""` can confuse
   Markdown.jl's own parser** — it can print one of the raw, unevaluated `$(...)` expressions
   verbatim (e.g. literally `$(nsrc == 1 ? "" : "s")`, rendered as garbled inline-math-looking text)
   instead of substituting its value. **Fix**: build the *entire* line as one plain Julia string
   first (ordinary `"..." * "..."` concatenation, or `"$(...)"` inside that string, is fine —
   ordinary string interpolation doesn't have this problem), then hand `md"""..."""` a single bare
   `$readout` interpolation with nothing else on that line.

2. **Once you've done that, don't indent the `$readout` line 4+ spaces inside `md"""..."""`.**
   Markdown treats 4-space (or one-tab) leading indentation as its own "indented code block" syntax
   — completely unrelated to Julia's own code indentation, but easy to introduce by accident since
   it's natural to indent `$readout` to match the surrounding `let`/`begin` block. A code block
   prints its content **verbatim**: no bold, no HTML entities, and it's wrapped in monospace styling
   — the whole line renders as raw source text instead of prose. Put `$readout` at column 0 (or
   inline on the same line as other static text) instead:
   ```julia
   readout = "Background velocity **$(round(Int, bs_safe["vp0"]))** m/s" * " · more text..."
   md"""
   $readout
   """
   ```
   If even that reads as risky to get right by eye, skip `md"""..."""` for this one cell entirely
   and call `Markdown.parse(readout)` directly (a plain function, not the string macro) — it takes
   the same indentation care but sidesteps re-checking the macro's own source-text parsing. **If you
   do this, you lose HTML-entity decoding**: `Markdown.parse` does not expand `&middot;`/`&#8320;`/etc.
   the way `md"""..."""`'s macro-time source parsing does, so those entities print as literal text.
   Use the actual Unicode character directly in the Julia string instead (`"·"`, `"₀"`) rather than
   its HTML entity name — it renders correctly either way and sidesteps the question of which code
   path is doing the decoding.

## Admonitions (callout boxes)

Julia's Markdown supports `!!! type "optional title"` blocks (indent the body one tab). Semantic usage,
as actually used in geoid-kernel:

| Type | Use for |
|---|---|
| `note` | A secondary aside or clarification that doesn't derail the main text — a sign convention, a scope caveat. |
| `tip` | A practical, actionable suggestion — "for research-grade numbers, use X instead." |
| `warning` | A specific place people's intuition or algebra commonly goes wrong — read this *before* you make the mistake. |
| `danger` | A common misconception being actively corrected, or something that looks fine but silently breaks the physics if missed. |
| `correct` | The payoff: a verification or resolution that confirms something raised earlier (in a `warning`/`danger`) actually checks out. |

Pick the type for what it does to the reader's confidence, not just how severe it sounds — `danger` is for
"this trips up almost everyone the first time," not for "this is dangerous."
