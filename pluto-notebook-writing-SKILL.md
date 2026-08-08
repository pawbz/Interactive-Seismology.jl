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
