---
name: pluto-widget-style
description: "Design system and build recipe for the dark-themed canvas + control-panel interactive widgets used in this repo's Pluto notebooks (reference implementations are src/misc/geoid-kernel.jl and src/misc/seismic-interferometry.jl). Use this whenever building, restyling, or resizing a bespoke interactive demo (a custom canvas/JS visualization driven by an @bind'd Julia struct) in ANY notebook under src/, or when asked to make a notebook widget match the style of the geoid/interferometry ones, look good on a 16:9 laptop/projector, use more of the screen width, or avoid vertical scrolling during a lecture. Not needed for plain PlutoUI sliders/plots that don't have a custom canvas."
---

# Pluto interactive-widget design system

This repo has a house style for "hero" interactive widgets: a `Base.show(io, ::MIME"text/html", w::YourInput)`
method that renders a dark-themed canvas (or canvases) plus a control panel, wired up with vanilla JS,
and bound back into Julia via Pluto's `@bind`. `src/misc/geoid-kernel.jl` and
`src/misc/seismic-interferometry.jl` are the two reference implementations — when in doubt, go read the
actual code there rather than trusting this doc to be exhaustive.

Read this whole file before touching a widget's HTML/CSS/JS. It is short on purpose; the reference
notebooks are the source of truth for anything not covered here.

## The shape of a widget

1. A Julia `struct YourInput` holding the default/initial values, constructed with sane defaults.
2. `Base.get(w::YourInput)` returning a `Dict{String,Any}` — this is the value Julia sees *before* the
   user ever touches the widget (the initial `@bind` value).
3. `Base.show(io::IO, ::MIME"text/html", w::YourInput)` that `write`s one big HTML string: a `<div id="...">`
   containing a scoped `<style>`, the markup (canvas/controls), and a `<script>` that does all the drawing,
   event handling, and calls a Pluto/AbstractPlutoDingetjes `@bind` publish function to send updates back.
4. A downstream Julia cell reads the bound dict (e.g. `_gk["visc_ratio"]`) and does the real computation.

Copy this shape from whichever reference notebook is closer to what you're building, don't invent a new
one. Keep the id prefix short and unique per widget (`gk-`, `si-`, pick your own two/three-letter prefix)
so multiple widgets' CSS/JS never collide when several notebooks are open in the same browser session.

## The one thing that will silently break on a fresh restart: bind-cell ordering

**`@bind x YourInput()` does not reliably create a Pluto dependency on the cell that defines
`YourInput`.** Pluto's reactive engine executes cells in *display* order on a fresh run/restart (not a
computed topological order), and its static analysis of a `@bind` expression does not treat a bare
constructor call like `YourInput()` as a reference to the cell defining that type. If the struct is
defined in the Appendix (this skill's own recommended layout — see `pluto-notebook-writing-SKILL.md`)
while its `@bind` cell is displayed near the top, **a fresh restart runs the bind cell before the struct
exists and throws `UndefVarError: YourInput not defined`.** Cell caching hides this the rest of the time
— once everything has run once in a session, incremental edits re-trigger correctly — so this only shows
up on `pluto-collab restart` or a real cold load, which is exactly the scenario a student hitting "run
all" or reopening the notebook will trigger. Confirmed by direct experiment: `geoid-kernel.jl`'s bind cell
*also* runs before its struct-definition cell in display order on restart, but doesn't error, because it
already carries the workaround below — it isn't a coincidence that geoid-kernel avoids this, it's a
deliberate fix already baked in that's easy to miss if you don't know to look for it.

**Fix**: add a trivial `const` flag as the very last line of the struct-definition cell, and reference it
as a bare statement — read, not used — at the top of the `@bind` cell. A bare variable reference *is*
reliably detected by Pluto's analysis, so this forces the dependency edge explicitly:

```julia
# in the struct/Base.get/Base.show cell (wherever it lives, e.g. the Appendix):
begin
    struct YourInput ... end
    Base.get(w::YourInput) = ...
    function Base.show(io::IO, ::MIME"text/html", w::YourInput) ... end

    const _yourprefix_ready = true   # last line — forces the ordering fix below
end

# in the bind cell (wherever it's displayed, e.g. near the top):
begin
    _yourprefix_ready
    @bind x YourInput()
end
```

Verify this the only way that actually catches it: `pluto-collab restart <nb.jl> --json`, then check
`cells_errored == 0` in the response (or grep the per-cell output for `UndefVarError`) — checking `status`
alone after incremental `run --cell` calls will not surface this, since those don't reproduce a cold
start's execution order.

**This isn't limited to the type name.** Any variable referenced only as a *nested argument* inside a
`@bind` expression — not the bound type itself — has the same problem: `@bind selected_phases
MultiCheckBox(phases, default=phases[1:1])` failed to see `phases` as a dependency even though it's a
bare variable reference, because that reference sits inside the `@bind` expression's arguments rather
than as its own top-level statement. Confirmed by direct experiment in `global-seismic-arrivals.jl`: the
fix is the exact same shape — hoist a bare `phases` statement above the `@bind` line in a `begin...end`
block. If a widget's own construction doesn't need any upstream data (the common case once you've pushed
computed values in via `CustomEvent` instead — see below), this whole class of bug disappears by
construction: nothing for the `@bind` cell to depend on but the struct itself.

## A simpler, official alternative to the hand-rolled wide-cell CSS

Everything in "Making it wide and scroll-free" below (the `pluto-cell:has(...)` override, the
`availW`/`heightBudget` sizing math) was hand-rolled to reach ~80vw before this was discovered:
**`PlutoUI.WideCell(content; max_width=...)`** is the actual, supported PlutoUI API for exactly this —
wrap whatever you'd otherwise `write(io, ...)` directly, e.g. `WideCell(@bind x YourInput(); max_width=1500)`
(matching `geoid-kernel.jl`'s `WideCell(@bind _gk GeoidGlobeInput())`). It uses a `ResizeObserver` to keep
the cell sized to the available editor width (up to `max_width`) continuously, including on window
resizes after the initial load — which the CSS-override approach below does not do, since it only sizes
once at script-load time. The two reference notebooks in this repo predate this discovery and still use
the manual CSS approach (verified working, not urgent to migrate), but **prefer `WideCell` for anything
new** — it's less code, officially supported, and more correctly responsive. If you do use it, remember
the bind-cell-ordering fix above still applies: `WideCell(...)` wraps the bind expression, it doesn't
change what Pluto's dependency analysis sees inside it.

## Rendering a lat/lon grid on a sphere: watch for Moiré, not just wrong data

If a widget paints a gridded dataset (a tomography model, any lat/lon map) onto a sphere as a mesh of
small quads — one quad per `(theta, phi)` cell, the same pattern `drawGlobe`'s outer-sphere loop uses in
both reference widgets — **pick a render-mesh resolution clearly finer than the data grid's own
resolution, not merely "close to" it.** Two periodic samplings with nearly-matching-but-different periods
(e.g. a 50-division render mesh over 180° ≈ 3.6°/band, sampling a real 3°-spaced data grid) produce a
Moiré beat pattern: clean, regular, *evenly-spaced* bands that look exactly like a bug in the data
(wrong lat/lon indexing, a transposed array, a reading-order mistake) but aren't — the underlying data can
be completely correct and still render as a barcode of stripes. This is easy to misdiagnose in the wrong
direction, since the symptom (suspiciously regular banding) really does look data-shaped.

Don't guess which one it is — check the data directly first, it's cheap: sample the actual array at a
fixed latitude across several longitude indices, and at a fixed longitude across several latitude indices
(e.g. `console.log(TOMO[depthIdx][30])` for one row, `TOMO[depthIdx].map(r => r[10])` for one column), and
confirm both show real, irregular variation — if they do, the data is fine and the bug is in the render
mesh, not the grid. Fix: make the render mesh's own angular step comfortably smaller than the data's
(roughly half or less) — e.g. going from a 50×100 to a 100×200 quad mesh against a 3° data grid resolved
this cleanly in `earth-internal-structure.jl`'s tomography globe, with no other code changes.

## Showing a whole family of curves: fade-all + highlight-on-hover beats a checkbox list

When a widget's Julia side produces a *variable-length* set of related curves (every seismic phase that
connects a given source/receiver, every candidate model fit, ...), resist the urge to add a checkbox per
item so the user can "select what to show" — once that list is more than a handful of items (tens of
seismic phases is typical), a checkbox panel is both a lot of UI chrome to build and worse pedagogically
than just showing everything at once. `global-seismic-arrivals.jl`'s ray-path widget draws every
computed path at a low `ctx.globalAlpha` (~0.28) so the whole family reads as a visible "bundle", then
hit-tests the mouse against all of them (`distToSegment` from the cursor to each polyline segment, cheap
at this data size) and redraws just the nearest one at full alpha, on top, with a floating tooltip (name +
value) anchored near the cursor. This also sidesteps the `@bind`-dependency trap above entirely: the
widget never needs to know the *available* item list up front, so there's nothing upstream for its
`@bind` cell to depend on — Julia just computes the full unfiltered set and pushes it in via the same
`CustomEvent` pattern used for any other computed result.

Two details that make this hold together:
- Draw the non-hovered items first, then draw the hovered one **last, in a second pass** — otherwise a
  faded neighbor drawn after the highlighted curve can visually clip over it.
- Recompute the hover hit-test and redraw on every `mousemove` (not just when the nearest item changes) —
  at this data scale (dozens of short polylines) it's cheap enough that simplicity beats the minor
  redundant-redraw cost, and it keeps the tooltip position glued to the cursor without extra bookkeeping.

## Visual design tokens

Dark theme, extracted from the two reference widgets after a 2026 readability pass for 16:9 laptop/projector
viewing. Don't go below these sizes/contrasts without a specific reason — they were bumped up once already
because the originals were too small to read from the back of a room.

| Token | Value | Used for |
|---|---|---|
| Background | `#000` (canvas), `#050505` (control-group panels), `#0b0b0b` (nested cards) | |
| Borders | `#374151` / `#2f3744` / `#1f2937` | canvas border, control-group border, card border |
| Primary text | `#e5e7eb` / `#f3f4f6` | headings, values |
| Secondary text | `#9ca3af` | labels, captions, axis ticks, hints — **not** `#6b7280`, it's too low-contrast for a projector |
| Body text | `#d1d5db` | |
| Accent blue / red | `#3b82f6` / `#ef4444` (or `#38bdf8`/`#dc2626` for canvas-drawn variants) | the two poles of a diverging quantity (dense/light, high/low, causal/acausal) |
| Base font | `14px sans-serif` | control panel body text |
| Control-title font | `20px`, `font-weight:700` | each control-group's heading |
| Canvas caption font | `13-14px sans-serif` minimum | any `ctx.font=...` text drawn on a canvas — never go below 12px, it's unreadable projected |
| Buttons | `border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px` | |

CSS class-naming convention (swap `gk`/`si` for your own prefix): `.{p}-workspace` (flex row holding the
canvas/canvases), `.{p}-controls` (the control grid), `.{p}-control-group` (one bordered panel inside the
grid), `.{p}-control-title`, `.{p}-control-row` (label / input / value, 3-column grid), `.{p}-value`
(the readout span), `.{p}-actions` (flex-wrap row of checkboxes/buttons).

## Marker shapes: source = star, receiver = downward triangle

When a widget draws a seismic **source** and **receiver** as draggable canvas markers (as opposed to a
generic drag handle for some other quantity), use a filled **star** for the source and a filled
**downward-pointing triangle** for the receiver, instead of a plain circle. This is the convention
`lame-theorem.jl`'s `PointForceRadiationInput` widget uses (`drawStarMarker`/`drawTriangleDownMarker`,
defined right next to `drawArrowHead2D`) — copy those two functions verbatim into a new widget rather than
reinventing the shape math:

```js
function drawStarMarker(cx, cy, r, fill, stroke){
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
function drawTriangleDownMarker(cx, cy, r, fill, stroke){
  ctx.beginPath();
  for(let i=0; i<3; i++){
    const ang = Math.PI/2 + i*2*Math.PI/3; // first vertex points down (canvas y grows downward)
    const x = cx + r*Math.cos(ang), y = cy + r*Math.sin(ang);
    i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y);
  }
  ctx.closePath();
  ctx.fillStyle = fill; ctx.fill();
  ctx.strokeStyle = stroke; ctx.lineWidth = 1.5; ctx.stroke();
}
```

This matches the marker shapes seismologists actually use on maps/cross-sections (a star for the
event, an inverted triangle for the station), so a reader who's seen a real epicenter map recognizes
the widget's markers instantly instead of having to learn an arbitrary circle-means-what convention.

## Panel titles, and not letting a "Readouts" box become a redundancy dump

When a widget has several canvas panels (a model view, a migration image, a data gather — see
`Born-approximation.jl`'s `BornScatteringInput`), give each one a short bold title directly above it
(`.{p}-panel-title`, `font-size:14px;font-weight:700`), distinct from the smaller `.{p}-caption` line
below the canvas. The title names *what the panel is* ("Model", "Migration Image", "Shot Gather");
the caption underneath carries *live detail specific to that panel's current state* (grid resolution,
counts, which source is being viewed). Splitting these two jobs across two lines reads faster than one
long caption trying to do both — a viewer scanning the row of panels only needs the titles; the detail
line is there when they look closer.

**Don't add a separate "Readouts" control-group that just repeats values already visible elsewhere.**
It's tempting to collect "current state" into one text box (`v₀ 2500 m/s · N sources · viewing #1 ·
mode: scatterer`), but check each line against what's *already on screen* before doing that: a slider's
own value already shows next to the slider (`.{p}-value` span), an active mode/toggle button already
shows via its own `.active` highlight, and "which source am I viewing" is naturally already part of
whichever panel actually shows that source's data (fold it into *that* panel's own title/caption
instead, e.g. `Shot Gather — Source #1`, rather than a fourth place reporting the same number).
`Born-approximation.jl`'s original Readouts group turned out to duplicate *every one* of its own four
lines once checked this way — removed entirely, with the one genuinely non-redundant fact left (source/
receiver counts) folded into the model panel's own caption instead. A readout box earns its place only
for a value that has nowhere else to naturally live.

## Spatial-axis tick labels drawn on top of a canvas: keep them off the true edge

For any canvas showing a physical x/z (or x/time) domain, drawing small tick-mark + label bands along
the bottom and left edges — directly on the canvas, as a semi-opaque strip, rather than reserving
separate space and re-plumbing your coordinate transform — is the cheapest way to add a real spatial
axis to an existing widget. Two panels sharing the same domain and the same tick step (e.g. a model
view and a migration image both spanning the same `XMAX`/`ZMAX`) read as one *shared* axis for free, as
long as both canvases are the same pixel size and you call the same drawing function with the same
`step` argument on both.

**The one thing that will clip labels at the canvas edge**: a naive `textAlign='left'` tick at `x=0` or
`textAlign='right'` tick at `x=W` places the glyph's edge *exactly* on the canvas boundary — for `left`
that's usually fine (the glyph extends inward), but for a *left-side vertical axis band*, right-aligning
each label's text so it ends near the band's own inner edge means the label's *start* can fall at a
negative x for anything wider than a couple characters (`"400 m"`, `"0.6 s"` at a narrow ~28px band
routinely start a few pixels before x=0 and get silently cut by the canvas boundary — no error, the
glyph just isn't there). **Fix**: measure the label with `ctx.measureText(label).width` and clamp its
drawn position to stay within `[2, bandSize-2]` (or `[2, W-2]` for the bottom axis) rather than trusting
`textAlign` extremes to keep it in bounds — a couple of pixels of margin on both ends is enough, and it
costs nothing to always do this rather than only where you've noticed clipping by eye.

## Multiple data panels of very different amplitude: normalize each one independently

When two (or more) canvas panels show physically related but *very differently scaled* quantities side
by side — the clearest case is a reference wavefield next to its scattered/perturbation counterpart
under a small-perturbation approximation (Born, linearized inversion, ...), where the perturbed field is
*by the approximation's own premise* much weaker than the reference — resist sharing one color/amplitude
scale computed from the combined data. A single shared `max(abs(...))` across both panels is dominated
by the larger one, and the smaller panel renders as visually flat/empty even when its own internal
structure is perfectly real and worth seeing. Compute each panel's own `max(abs(...))` from its own data
and normalize independently; if a viewer might mistake matching *colors* between panels for matching
*absolute amplitude*, say so explicitly in the caption (e.g. "each panel independently amplitude-scaled")
so the comparison isn't accidentally over-read.

## A bare `$` in embedded JS breaks the enclosing Julia string, not just LaTeX `$...$`

The well-known trap is writing math as `$x^2$` inside `md"""..."""` and having Julia try to interpolate
it (see `pluto-notebook-writing-SKILL.md`). The *same* root cause bites inside a widget's `<script>` block
too, in a form that's easy to miss because nothing about it looks like math: any JS **regex literal that
ends a pattern with `$`** (the end-of-string anchor), e.g. `` /^nm-c(\d)$/ ``, sits inside the Julia
`"""..."""` string that `write(io, """...""")` builds. Julia's parser sees that bare `$` and tries to
start an interpolation, expecting an identifier or `(...)` next — finds `/` instead — and throws
`ParseError: identifier or parenthesized expression expected after $ in string`. Confirmed on
`wave-mode-duality-1D.jl`: three regexes matching slider ids (`` /^nm-c(\d)$/ ``, `` /^nm-b(\d)$/ ``,
`` /^nm-c\d$/ ``) all needed this fix. **`Meta.parseall` on the whole file does not catch it** — the file
still parses as valid Julia overall; the error only surfaces per-cell, inside Pluto's own cell-level
parse (`pluto-collab restart` reports `mime: application/vnd.pluto.parseerror+object` with an *empty*
`output_text` for that cell, which is easy to misread as "the cell produced no output" rather than "this
cell failed to parse at all" — check the `mime` field, not just whether `output_text` is empty). **Fix:**
escape it as `\$` (a recognized Julia string escape producing a literal `$`), the same rule as any other
literal `$` you need inside a Julia string — `` /^nm-c(\\d)\$/ `` is what actually reaches the JS engine
as `` /^nm-c(\d)$/ ``. Any JS syntax that uses a bare trailing `$` (regex anchors are the common case;
`${...}` template-literal interpolation is another, less likely with `\\`-doubling already in place)
is worth a specific grep (`grep -n '\$/' file.jl`) after writing or editing widget JS, since this is
otherwise invisible until Pluto tries to load that exact cell.

## A capture-phase listener that unconditionally `stopImmediatePropagation()`s can silently disable a
## sibling element's own listener for the same event type

The delegated-listener pattern used throughout this repo's widgets — one `par.addEventListener('input',
e=>{...}, true)` (capture phase, so it intercepts before Pluto's own bond listener ever sees the event —
see "Pluto bond input bubbling") dispatching on `e.target.id` to handle every slider/checkbox in one place
— has a sharp edge: **if that handler calls `e.stopImmediatePropagation()` unconditionally for any
non-`par` target, instead of only for the ids it actually recognizes, it silently blocks every *other*
listener registered directly on that same element for the same event type**, including a listener that
element itself owns (e.g. a `<select>`'s own dedicated `addEventListener('change', ...)` for handling
layer-count changes). Capture-phase listeners on an ancestor always run *before* an element's own
target-phase listener, so the ancestor's `stopImmediatePropagation()` wins even though the target's own
listener was registered first in wall-clock time. Confirmed on `wave-mode-duality-1D.jl`: a `#nm-nlayers`
`<select>` had its own `change` listener to resize the `velocities`/`boundaries` arrays, but a sibling
capture-phase `par.addEventListener('change', e=>{ if(e.target===par) return; e.stopImmediatePropagation();
... }, true)` — meant only to intercept a handful of slider ids — stopped *every* change event at the
capture stage regardless of id, so picking a new layer count in the dropdown silently did nothing (the
`<select>`'s displayed value changed; the underlying JS array never did). Symptom is maddening to diagnose
purely by testing in a browser automation tool, because setting `select.value` directly (e.g. via a test
harness's form-fill helper) *also* doesn't fire a real `change` event — so a naive test can't tell "my
listener is blocked" apart from "my test harness never triggered the listener at all"; confirm with
`el.dispatchEvent(new Event('change', {bubbles:true}))` from `javascript_tool` to rule out the harness
before suspecting your own code. **Fix:** move `e.stopImmediatePropagation()` *inside* the `if` branch
that actually recognizes and handles the id, not before it — only stop propagation for the ids this
handler is actually the owner of, and let everything else (including a sibling element's own listener)
proceed normally.

## Title bar: orient the viewer before they touch anything

Every widget opens with a two-line header, centered, sitting above `.{p}-workspace`, boxed like the
control-groups below it (same rounded-corner language) but with a distinct accent border/background so it
still reads as "the important one" rather than blending into the row of control panels:

```html
<div class="{p}-title">
  <div class="{p}-title-desc">One sentence: what physical/scientific relationship this demo shows.</div>
  <div class="{p}-title-hint">the 1-2 interactions that matter &middot; separated by &middot;</div>
</div>
```
```css
#{id} .{p}-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
  background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
#{id} .{p}-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
#{id} .{p}-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
```

This exists so someone can tell what the widget is *for* and how to drive it before they start clicking —
without it, a viewer sees a globe or a canvas full of dots and has to read a paragraph of prose above the
cell to find out what "drag" or "click" even does. Keep both lines to one sentence each; if the hint needs
a third clause, you're describing every control instead of the one or two gestures that matter most. A
more detailed hint can still live inline near the specific control it applies to (see the `drag to rotate
...` line inside geoid-kernel's "Anomaly & physics" group) — the title bar is the *short* version for
first contact, not a replacement for it. `#0a0f18`/`#3b5c85` (a cool blue-tinted background/border) is a
deliberately different pair from the control-groups' neutral `#050505`/`#2f3744` — same box language,
different accent, so the header doesn't camouflage itself among the panels underneath it.

This title bar adds real height — **fold it into the `heightBudget` reserve** (next section) or you'll
reintroduce the exact overflow you fixed there. Don't estimate that height from "two lines of text plus a
bit of margin": **measure it live**, `getBoundingClientRect().height` on `.{p}-title`, after your actual
copy is in place. A one-sentence science description easily wraps to two lines once it's boxed at the
widget's real width, which roughly doubles the naive estimate — that gap is exactly wide enough to produce
a "barely overflows, but only sometimes" bug that looks fine on whatever one screen size you happened to
test first. (Related trap while chasing this kind of small overflow: if you resize the browser viewport
and re-measure *without reloading the page*, you're reading stale numbers — this widget's sizing math runs
once at script load, not on every resize, so a resize-without-reload will make a perfectly-fitting layout
look broken for reasons that have nothing to do with your CSS. Always reload after resizing before you
trust a measurement.)

## Making it wide and scroll-free on a 16:9 screen

The ask that motivated this section: "make the widget use ~80% of the screen width, and don't make me
scroll to see the whole thing during a lecture." Two independent problems, two independent fixes — do
both, they don't substitute for each other.

### 1. PlutoUI's wide-cell wrapper caps itself around 1000px, regardless of your screen

If the widget is wrapped in `PlutoUI.WideCellNotebook`'s wide-cell mechanism (look for `pui-big-wrapper`
in the rendered DOM), Pluto sets an **inline `width` style on the `pluto-cell` element itself**, capped
around 1000px no matter how wide the browser window is. Your own CSS can't beat that by scoping to your
widget's `#id` — you have to reach up and override the ancestor, from inside your own `<style>` block
(style tags apply document-wide, not scoped, so this is legal and it only needs to affect this one cell):

```css
pluto-cell:has(#yourwidget) { width: min(80vw, 1500px) !important; }
```

`:has()` scopes this to only the cell containing your widget — it won't touch any other cell. Verify this
assumption still holds before relying on it (inspect the live DOM — `getComputedStyle`/inline `style` on
`pluto-cell`/`pui-big-wrapper`/`main` — rather than assuming; Pluto/PlutoUI internals can change).

**This alone shifts the whole widget off-center to the right — you need a second override alongside it.**
Pluto centers a wide cell by *also* setting an inline `margin-left` on `pluto-cell`, computed for its own
~1000px cap (e.g. `margin-left: -150px` to center a 1000px cell under a 700px notebook column). Your
`width` override changes the cell's width but Pluto never recomputes that margin for the new width, so the
now-wider cell keeps the *old, too-small* negative margin and drifts right by roughly half the width
increase. Fix it by recomputing the margin in the same rule, using a `calc()` that resolves against the
*parent's* width (percentages in `margin-left` are relative to the containing block, not the element
itself, which is exactly the value you want here — no need to hardcode or measure the parent's width):

```css
pluto-cell:has(#yourwidget){
  width: min(80vw, 1500px) !important;
  margin-left: calc((100% - min(80vw, 1500px)) / 2) !important;
}
```

Check for this any time you override an ancestor's `width` from inside a widget — a sibling `margin`/`left`/
`transform` computed relative to the *old* size is a common way for that kind of override to leave things
visually off-kilter even though the override itself "worked."

### 2. Size the canvas(es) and control grid from the viewport, not a fixed pixel constant

At the top of the widget's `<script>`, before anything else, compute how much room you actually have and
derive your drawing-surface size from it:

```js
const par = currentScript.previousElementSibling   // your #yourwidget div
const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8)
const heightBudget = Math.max(<floor>, window.innerHeight - <reserveForEverythingElse>)
```

`availW` is capped by *both* 80% of the viewport *and* whatever the (now-widened) wide-cell wrapper
actually gives you — belt and suspenders. `heightBudget` is what makes "no scroll" real: it bounds your
square/tall canvas's height by the vertical space actually left after reserving room for whatever sits
above/below it (a comparison strip, the control panel, margins). **Measure that reserve empirically** —
don't guess. Load the widget, read `getBoundingClientRect().height` on the controls section and any other
canvas, and set the reserve to match. Guessing low causes the widget to overflow the screen; guessing high
just under-sizes the canvas. See "How this was tuned" below for the actual numbers that came out of doing
this for the two reference widgets.

**Then let the controls grid use the full `availW`, not the (possibly much smaller, height-limited) canvas
row's width.** These are two independent widths — a square canvas capped by *height* can end up much
narrower than the space you actually have *horizontally*. If you tie the control grid's width to the
canvas row instead of to `availW` directly, you silently throw away width you already have, and the
control panel stacks into more rows than it needs to — which is exactly the vertical space you were
trying to save. Set a CSS custom property from JS and reference it in the grid's width rule:

```js
par.style.setProperty('--totalw', Math.round(availW) + 'px')
```
```css
.{p}-controls{ width:min(var(--totalw,960px),100%); display:grid;
  grid-template-columns:repeat(auto-fit,minmax(<colFloor>px,1fr)); gap:8px; font:14px sans-serif }
```

Using `repeat(auto-fit,minmax(...))` instead of a fixed `repeat(2,...)` matters for the same reason: on a
short/wide screen, letting your 4 control-groups flow into 4 columns instead of 2 turns two cramped rows
into one, which is often the single biggest lever for eliminating a lingering scroll. Pick `<colFloor>`
by checking how many columns you actually want at your realistic minimum width — too high a floor and
`auto-fit` can't pack them in even when there's room (do the arithmetic: `numGroups * colFloor + gaps`
must fit in your typical `availW`, or you'll silently get fewer columns than you expect).

**This has a knock-on effect you have to fix too: once a control-group can legitimately be as narrow as
`<colFloor>`, every `.{p}-control-row` *inside* it needs to be able to shrink that far as well.** If a
row's own grid uses fixed-pixel columns for its label (`label-width minmax(...) value-width`), that label
width never shrinks — drop `<colFloor>` to fit 4 columns and a 200px fixed label plus a 120px slider
minimum can together demand more width than the now-narrower group actually has, and the slider's track
bleeds visibly past the group's own border into whatever sits next to it. Make every column in the row a
range too — `minmax(<min>,<max>)` instead of a bare `<max>px` — so labels can wrap onto a second line
and sliders can compress instead of overflowing:

```css
.{p}-control-row{ grid-template-columns: minmax(70px,130px) minmax(70px,1fr) minmax(36px,64px) }
```

and give the `<input type=range>` itself `min-width:0` — without it, a grid item's default automatic
minimum size is its content's intrinsic width (a range input's native min-content can be wider than you'd
expect), which silently overrides `width:100%` when the track is squeezed and produces the exact same
overflow. This is a general CSS Grid trap, not specific to sliders: any replaced element (`<input>`,
`<img>`, `<canvas>` used as a grid item) needs `min-width:0` before you can trust `minmax()` to actually
let it shrink. Catch it the same way as any other overflow — see "Verifying your changes" below — because
it looks completely fine at whatever width you happened to eyeball it at, and only shows up once a group
gets narrow enough (e.g. after adding a 4th `auto-fit` column) to be smaller than the row's true minimum.

**Any plain-HTML/CSS panel next to a canvas (a sidebar metrics/readout panel, say) needs its own width
floor, independent of your zoom/size factor.** Canvases redraw at whatever size you give them, but text
in an HTML panel does not — shrink its container width while its font stays fixed and the text just wraps
onto far more lines, so the panel gets *taller* exactly when you were trying to make the whole widget
*shorter*. Give it a `Math.max(<minPx>, ...)` floor, not a bare proportional shrink.

**Don't wrap the whole formula in an outer floor** (e.g. `Math.max(0.5, Math.min(...))`). If your height
budget legitimately wants a smaller size than the floor allows, an outer floor overrides the one
computation that was keeping you scroll-free, defeating the entire point. Put any necessary floor
*inside* one of the `Math.min(...)` branches (e.g. on `heightBudget` itself) so it can never make the
final result *larger* than what actually fits.

### 3. The one thing that can actually break: pixel coordinates baked into Julia

**Before you resize anything, check whether any Julia cell downstream hardcodes a pixel-space assumption
about the canvas** — e.g. `(I[1] - 320) / 1.5` to convert a bound coordinate back to physical units,
where `320` and `1.5` are literally the JS widget's `MID` and `SCALE` constants, copy-pasted into Julia.
`grep` the notebook for the widget's numeric constants (its canvas width, half-width, any per-km/per-unit
scale factor) outside the `<script>` block. If you find this coupling, DO NOT let those constants vary
with screen size — Julia has no way to know at parse time what the browser's viewport will be, so making
them dynamic silently corrupts every downstream physical calculation, with no error or exception to catch
it.

Instead, keep the **logical/design coordinate space fixed forever** (whatever the original constants were)
and only change the **on-screen CSS size**, via a canvas rendering-transform trick:

```js
const zoom = <your computed display-size> / <original design width>
function setupHiDPICanvas(canvas, ctx, w, h) {
  canvas.width  = Math.round(w * zoom * DPR)
  canvas.height = Math.round(h * zoom * DPR)
  canvas.style.width  = Math.round(w * zoom) + 'px'
  canvas.style.height = Math.round(h * zoom) + 'px'
  ctx.setTransform(zoom * DPR, 0, 0, zoom * DPR, 0, 0)
}
```

Every `ctx.` drawing call keeps using the *original* logical coordinates (0..640, or whatever) — the
transform does the visual scaling for you. The one thing this breaks: mouse/touch handlers read
`e.offsetX`/`e.offsetY` in the *enlarged* CSS-pixel space, so every place your code reads those, divide by
`zoom` first to get back into logical space before doing any hit-testing/dragging math:

```js
// wherever you previously wrote e.offsetX / e.offsetY, write:
(e.offsetX/zoom), (e.offsetY/zoom)
```

`seismic-interferometry.jl` needs this (its Julia cell does the `(I .- 320)/1.5` conversion above) —
`geoid-kernel.jl` does not (its widget only ever emits already-physical `(θ, φ, r_fraction, ...)` values
to Julia, never raw pixels), which is why that one just resizes its logical constants (`SEC`, `KW`, `R_OUT`)
directly instead of using the zoom-transform. **Check which situation you're in before copying either
approach** — grep for the widget's own numeric constants in any Julia cell that isn't inside the `<script>`.

## Blank space when the control panel outgrows the canvas column

A single-column `.{p}-controls` (one `{p}-control-group` stacked under the next, `flex-direction:column`)
looks fine when a widget starts with 2-3 control groups, but every widget in this repo accretes groups
over its life as features get added (a view toggle, a new background-model picker, a new slider group...).
Once the controls column is taller than the canvas+colorbar column next to it, the flex row's shorter side
(the canvas) just stops, and the empty page background below it reads as a big wasted gap — this is a
`ray-theory-eikonal-equation.jl` fix, hit at 6 control groups against a short, wide (500km × 150km) canvas.

**Fix**: once a widget's controls column exceeds ~4-5 groups, switch it from a single flex column to an
explicit 2-column CSS grid, and widen its allotted width to match:

```css
#{id} .{p}-controls{flex:0 0 540px;width:540px;display:grid;
  grid-template-columns:repeat(2, minmax(0,1fr));gap:8px;font:14px sans-serif;align-content:start}
```

This roughly halves the controls column's total height (6 groups → 3 rows instead of 6), bringing it back
in line with the canvas column instead of towering over it. Two follow-ups this rule requires, both easy
to forget:

- **Bump the JS `CONTROLS_W` constant to match** (`260+16` → `540+16` for the above) — it feeds `SEC_W =
  Math.max(<floor>, Math.min(availW - CONTROLS_W, <cap>))`, so a stale `CONTROLS_W` after widening the CSS
  either starves the canvas of width it should get, or (if left too small) lets the two columns overflow
  `availW` and wrap ungracefully.
- Prefer a fixed `repeat(2, ...)` over `repeat(auto-fit, minmax(...))` here specifically — auto-fit is the
  right call for a *family* of same-shaped items (see the control-grid-floor guidance under "Making it
  wide" below), but a widget's control groups are usually a small, known, differently-sized handful (a
  toggle-button group next to a 4-slider group next to a readout panel) where you want a predictable 2-up
  flow, not however many columns happen to fit at a given width.

Don't reach for this preemptively on a 2-3-group widget — the point is the controls/canvas height ratio,
not a fixed group count; check the actual rendered heights (or just eyeball a screenshot) before deciding
a widget needs it.

## How this was tuned (concrete numbers, so you have a starting point)

For the two reference widgets, on real testing across 1024×768 up to 1920×1080:
- `heightBudget` reserve: geoid-kernel needs ~400px reserved (title bar + its 2-column, 2-group control
  panel); seismic-interferometry needs ~540px (title bar + comparison strip + a richer 4-group control
  panel + a metrics sidebar). **Your reserve depends on how much fixed-height content sits above/below the
  canvas** — the title bar, any comparison strip, the control panel — measure, don't copy these numbers
  blindly, and re-measure any time you add content that wasn't there when you last tuned it.
- `pluto-cell` override: `min(80vw, 1500px)` — the 1500px ceiling stops it from becoming absurdly wide on
  a 4K display; drop it if you genuinely want it to keep growing.
- Control-grid column floor: 320px (geoid, 2 groups) / 260px (interferometry, 4 groups — needed a lower
  floor so 4 columns can actually fit in one row at typical laptop widths).
- Known residual limitation: on genuinely short screens (~1366×768 and below), the interferometry widget
  can still overflow by a small amount (a few percent) because its control panel + metrics sidebar have a
  content-driven minimum height that doesn't shrink further without cutting font size (which was
  explicitly *not* an acceptable trade-off — readability was the other half of this whole exercise). This
  is disclosed, not silently broken: verify your own widget's floor the same way (see below) rather than
  assuming zero-scroll is achievable at every resolution.

## Verifying your changes

You cannot eyeball this from static file edits — the whole point is what happens at different screen
sizes, which you have to actually measure in a live browser.

1. If a live Pluto/SpaceStation server is running (see this repo's `CLAUDE.md`/`AGENTS.md` for
   `pluto-collab`), edit the `.jl` directly, then `pluto-collab run <nb.jl> --stale` to apply.
2. Open the notebook in the Browser tool. Resize the viewport to a few representative sizes (a small
   laptop ~1366×768, a common laptop ~1440×900, a big monitor/projector ~1920×1080) and reload each time
   (the sizing math only runs once, at script load).
3. Measure, don't just screenshot: `document.getElementById('yourwidget').getBoundingClientRect()` for
   total width/height vs. `window.innerWidth`/`innerHeight`, and `document.body.scrollWidth >
   document.body.clientWidth` for horizontal overflow. Also check the widget is still *centered* —
   compare the midpoint of `pluto-cell`'s rect to the midpoint of `main`'s rect (see the centering fix
   above; an off-by-a-stale-margin bug looks fine in a screenshot at a glance and only shows up as "this
   feels shifted").
4. Check that stacked sections don't overlap: get `getBoundingClientRect()` for the title bar, any
   comparison strip, the workspace, and the controls grid, and confirm each one's `bottom` is above the
   next one's `top`. A tight `heightBudget` reserve makes this an easy way to accidentally overlap two
   sections instead of just running out of room at the very end — checking rects catches it where a
   screenshot at one specific viewport size might not.
5. If you touched any mouse-coordinate code (the zoom-transform case), sanity-check hit-testing by
   dispatching synthetic `MouseEvent`s at known coordinates via `javascript_tool` and checking a readout
   value changes by the expected amount — the Browser tool's drag gesture doesn't fire real intermediate
   `mousemove` events, so a plain click-and-drag test can pass or fail for the wrong reason.
6. Check `read_console_messages` for JS errors and `pluto-collab status` for errored/stale cells after
   interacting.
