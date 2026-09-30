### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> tags = ["basics"]
#> title = "How a Circular Wavefront Is Built From Plane Waves"
#> description = "Watch plane waves summed over every propagation direction collapse into a circular wavefront J0(kr) -- the converse of the Jacobi-Anger expansion."
#> layout = "layout.jlhtml"

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

# ╔═╡ 7e266e4d-e0fa-45f4-ad8b-b6c4f784a392
begin
    using Bessels
    using PlutoUI
end

# ╔═╡ ffdeeaca-a54c-47a2-9525-6a8875368703
TableOfContents()

# ╔═╡ 70042daf-6ff4-49c1-b50a-d9e8df536776
md"""
# How a Circular Wavefront Is Built From Plane Waves

A point source's wavefront is the simplest possible *circle*: a clean ring expanding outward,
identical in every direction. A plane wave is the opposite -- dead-straight, parallel stripes
marching in a single direction. So how can summing plane waves ever produce a circle?

This notebook is the deliberate converse of this repo's `Jacobi-Anger-expansion.jl`: there, a
sum over angular *order* `n` of cylindrical harmonics `Jₙ(kr)e^{inθ}` builds a single plane
wave. Here, a sum over propagation *direction* `φ` of plane waves -- all sharing the same
wavenumber `k` -- builds the circularly-symmetric standing wave `J₀(kr)`: literally the same
identity, read backwards.

```math
J_0(kr) = \frac{1}{2\pi}\int_0^{2\pi} e^{ikr\cos\phi}\,d\phi
```

Average enough plane waves over enough directions, and the criss-crossing stripes collapse into
concentric rings.
"""

# ╔═╡ 931da2f3-774e-42b6-9f92-37ccb250479c
md"""
## Reading the Panels

**Drag `N`** (how many evenly-spaced propagation directions are summed) and **`kr`** (the
aperture -- how many wavefront rings fit across the panel), or pick a station-array preset.

The **left** panel shows exactly which `N` plane-wave directions are being summed -- `N` ticks
spaced evenly around the compass. The **middle** panel is the field that summing them actually
produces,
`` S_N(x,y) = \mathrm{Re}\left[\frac1N\sum_{j=0}^{N-1} e^{i(x\cos\phi_j+y\sin\phi_j)}\right] ``,
plotted in the same *k*-scaled coordinates as the Jacobi–Anger notebook (`x=kX, y=kY`). The
**right** panel is the exact target, `` J_0(kr) ``, plotted the same way.

The dashed circle marks the aperture `kr`. Start from a small `N` (say `N=3`) and watch the
middle panel show a criss-crossed, `N`-fold-symmetric "flower" -- nothing like the smooth rings
on the right. Now grow `N` past `kr`: the flower collapses into clean concentric rings, and the
middle and right panels converge. This is exactly the same **"you need roughly `N≈kr` terms"**
rule of thumb the Jacobi–Anger notebook teaches from the opposite direction -- here it is the
number of *directions*, not the number of *orders*, that must clear the aperture.

Notice the flower never disappears completely for finite `N` -- it just shrinks. The Appendix's
self-check confirms this precisely: the discrete sum equals `J_0(kr)` *exactly*, plus a residual
of order `J_N(kr)`, which is why growing `N` (or shrinking `kr`) makes it vanish.

**Press `Play`** to watch it move. Every plane wave in the sum shares the same `ω` (they all sit
on the same slowness circle), so the field's entire time dependence is one shared rotation,
`` \mathrm{Re}[Z_N\,e^{-i\omega t}] ``, applied to the same two numbers already computed at every
pixel -- nothing needs to be re-integrated frame to frame. Set `N=1` and the "sum" is a single
plane wave: its stripes visibly slide sideways, a genuinely *propagating* wave. Now grow `N` up
past `kr`: the sliding stops, and instead the whole ring pattern brightens and dims together, in
place -- a **standing** wave, `J_0(kr)\cos(\omega t)`, pulsing rather than travelling. That
contrast *is* the answer to "why doesn't this look like an expanding circle": a superposition
confined to one wavenumber magnitude is exactly the case where every point oscillates in phase,
so nothing ever propagates outward.
"""

# ╔═╡ 62bdf2c0-7aed-467d-8f40-b5d7761a0573
md"""
## Why Doesn't This Expand? Standing vs. Outgoing Waves

Here is the precise reason, not just an analogy. At one frequency `ω`, *every* real, propagating
plane wave in a uniform medium is forced onto the same circle `|k|=ω/c` in wavenumber space --
that circle is all the dispersion relation allows. Averaging real plane waves over that circle,
exactly what this notebook does, can therefore only ever produce `J_0(kr)`: smooth, finite at
`r=0`, and *real* -- because it is built entirely from waves carrying no net energy flux in any
one direction. A real point source is not like that: it radiates a genuinely **outgoing** wave,
singular at the source, `` H_0^{(1)}(kr) `` (its 3D analogue is `` e^{ikr}/r ``, the spherical
wave from the Lamb's-problem notebook). Something has to be missing.

What's missing is captured exactly by one identity:

```math
J_0(kr) = \frac{H_0^{(1)}(kr) + H_0^{(2)}(kr)}{2}
```

`` H_0^{(1)} `` is the outgoing Hankel function, `` H_0^{(2)}=\overline{H_0^{(1)}} `` the
*incoming* one. Our circularly-symmetric standing wave is not a third, independent kind of wave
-- it is *exactly* half an outgoing ring plus half an incoming ring, superposed. Their outward and
inward energy fluxes cancel term by term, which is precisely why the animation pulses in place
instead of travelling: it already contains a wave travelling outward and an equal wave travelling
inward, perfectly overlapped.

So building the true, singular, one-way field of a point source needs more than real propagating
plane waves -- it needs **evanescent** ones too. Pick a reference plane (say `z=0`) and decompose
over the *horizontal* wavenumber `` (k_x,k_y) ``, now free to take any real value, not only values
with `` k_x^2+k_y^2 \le k^2 ``:

```math
\frac{e^{ikr}}{r} = \frac{i}{2\pi}\iint \frac{e^{\,i(k_xx+k_yy+k_z|z|)}}{k_z}\,dk_x\,dk_y,
\qquad k_z=\sqrt{k^2-k_x^2-k_y^2},\ \ \operatorname{Im}k_z\ge0.
```

(This repo's `Lamb_problem.jl` -- see its "Conical Wave" section and its `vertical_slowness`
function -- implements exactly this branch choice.) Once `` k_x^2+k_y^2 `` exceeds `` k^2 ``,
`` k_z `` turns imaginary and that term stops propagating in `z`, decaying instead as
`` e^{-|k_z||z|} ``. These **evanescent** components are the plane-wave representation of the
source's near field: a point is infinitely sharp, resolving it needs arbitrarily fine
(large-horizontal-wavenumber) structure, and nothing that fine can propagate at a fixed `ω` -- so
it stays pinned near the source and decays instead. Nothing in *this* notebook's own sum ever
goes evanescent, because every term already sits on the same real circle `|k|=ω/c`; that
restriction is exactly what makes it elementary enough to sum by eye, and exactly why it cannot
reproduce a singular, one-way-outgoing field.

This is also why seismologists build wavefields this way at all: fixing a horizontal slowness
`` p=k_x/\omega `` (conserved across interfaces by Snell's law -- see this course's plane-wave
lecture notes) lets each component reflect and transmit independently in a layered medium, then
be superposed afterward. The evanescent region `` p>1/c `` is not a footnote: head waves come
from a branch cut there, and surface waves (Rayleigh, Stoneley) are poles of the reflection
coefficient entirely inside it -- guided energy built from components that decay away from the
interface rather than ever leaving it.

The self-check `Standing/outgoing waves` in the Appendix verifies the boxed identity above
directly, alongside `` H_0^{(2)}=\overline{H_0^{(1)}} ``.
"""

# ╔═╡ 2e0c5129-76d4-4200-b2d3-471e404ef61b
md"""
## Appendix
"""

# ╔═╡ 3a1e3154-63b1-4661-ab44-601e47874271
md"""
## Summing Plane Waves Over Direction
"""

# ╔═╡ 57218488-7de4-487d-bb2f-6e94b4055458
"""
    plane_wave_directions(N)

`N` evenly-spaced propagation directions `φⱼ = 2πj/N` for `j = 0, …, N-1` -- the discrete
angular grid [`summed_plane_wave_field`](@ref) sums plane waves over.
"""
plane_wave_directions(N) = [2 * pi * j / N for j in 0:N-1]

# ╔═╡ 5dbfb284-a12c-49cf-ae60-6807513fedd1
"""
    summed_plane_wave_field_complex(X, Y, N)

The complex-valued average of `N` unit-amplitude plane waves, one per direction from
[`plane_wave_directions`](@ref), evaluated at a point `(X,Y)` given in `k`-scaled coordinates:

``Z_N(X,Y) = \\frac1N\\sum_{j=0}^{N-1} e^{i(X\\cos\\phi_j+Y\\sin\\phi_j)}``

This is the snapshot at `t=0` of a genuinely time-dependent field: each term is really
`` e^{i(X\\cos\\phi_j+Y\\sin\\phi_j-\\omega t)} `` (a plane wave travelling at the *same* `ω`, since
every direction shares the same `k`), so the field at any later time is just this same complex
number rotated by an overall phase, `` \\mathrm{Re}\\!\\left[Z_N(X,Y)\\,e^{-i\\omega t}\\right] =
\\mathrm{Re}(Z_N)\\cos(\\omega t) + \\mathrm{Im}(Z_N)\\sin(\\omega t) ``. Julia computes `Z_N` once;
the widget's `Play` button animates that closed-form combination directly in JS, with no
re-computation needed as `t` advances. [`summed_plane_wave_field`](@ref) is the `t=0` slice of
this, `\\mathrm{Re}(Z_N)`.
"""
function summed_plane_wave_field_complex(X, Y, N)
    s = zero(ComplexF64)
    for phi in plane_wave_directions(N)
        s += cis(X * cos(phi) + Y * sin(phi))
    end
    return s / N
end

# ╔═╡ 27f7ef58-1d43-4bda-be5e-a455aee66bab
"""
    summed_plane_wave_field(X, Y, N)

The plane-wave sum this notebook is built around: the real part (the `t=0` snapshot) of
[`summed_plane_wave_field_complex`](@ref),

``S_N(X,Y) = \\mathrm{Re}\\left[\\frac1N\\sum_{j=0}^{N-1} e^{i(X\\cos\\phi_j+Y\\sin\\phi_j)}\\right]``

By discrete Fourier orthogonality applied to the Jacobi–Anger series, only orders `n` with
`n ≡ 0 (mod N)` survive averaging over these `N` equally-spaced angles, so this equals *exactly*
`` J_0(r) `` (`` r=\\sqrt{X^2+Y^2} ``) plus a residual of order `` J_N(r) `` from the leading
`n=\\pm N` terms -- see the self-check below for the residual's measured size. This is the
deliberate converse of `planewave_partial_sum` in `Jacobi-Anger-expansion.jl`: there a sum over
angular *order* `n` builds a plane wave; here a sum over propagation *direction* `φ` builds a
circularly-symmetric wave.
"""
summed_plane_wave_field(X, Y, N) = real(summed_plane_wave_field_complex(X, Y, N))

# ╔═╡ 5fd8e68c-3302-4911-8831-8f65912c6e59
md"""
### Verifying the Summation
"""

# ╔═╡ 92c94385-480b-4bcb-9f0a-bb99a795e367
let
    # 1. the discrete N-direction average converges to the exact target J0(r), with the
    #    aliasing residual bounded by the leading ±N terms of the Jacobi-Anger series (a floor
    #    of 1e-14 absorbs floating-point roundoff once the true residual drops below it)
    r, theta = 8.0, 0.7
    X, Y = r * cos(theta), r * sin(theta)
    Ns = [5, 10, 20, 40, 80]
    errs = [abs(summed_plane_wave_field(X, Y, N) - besselj(0, r)) for N in Ns]
    bounds = [3 * abs(besselj(N, r)) + 1e-14 for N in Ns]
    println("Convergence at r=$r as N grows: errs=", errs, " bounds=", bounds)
    @assert all(errs .<= bounds) "error should stay within the predicted aliasing bound ~J_N(r)"
    @assert errs[end] < 1e-6

    # 2. exact N-fold rotational symmetry holds for ANY N, even small ones -- shifting the
    #    observation angle by one direction-spacing exactly permutes the summed terms
    N = 7
    shift = 2 * pi / N
    lhs = summed_plane_wave_field(r * cos(theta + shift), r * sin(theta + shift), N)
    rhs = summed_plane_wave_field(r * cos(theta), r * sin(theta), N)
    println("N-fold symmetry check (N=$N): |Δ| = ", abs(lhs - rhs))
    @assert abs(lhs - rhs) < 1e-10

    # 3. full (continuous) circular symmetry only emerges once N clears the aperture r -- a
    #    small N leaves visible angular structure, a large N does not
    function spread(N)
        vals = [summed_plane_wave_field(r * cos(t), r * sin(t), N) for t in range(0, 2 * pi / N; length=25)]
        return maximum(vals) - minimum(vals)
    end
    spread_small, spread_large = spread(5), spread(60)
    println("Angular spread at r=$r: N=5 -> $spread_small, N=60 -> $spread_large")
    @assert spread_small > 0.05
    @assert spread_large < 1e-3

    "Summation self-checks passed"
end

# ╔═╡ b4f123af-2e6e-46b4-a57e-ea7ff710d29c
md"""
## Standing Waves Are Half Outgoing, Half Incoming

`` J_0 `` is not a third kind of circular wave alongside the outgoing and incoming ones -- it is
*built* from them, exactly half of each.
"""

# ╔═╡ 2df09cda-c40f-4b8c-a254-c18a857c652e
let
    # J_0(kr) is exactly the symmetric combination of the outgoing Hankel function H_0^(1) and
    # its complex conjugate, the incoming H_0^(2) -- this is *why* averaging real plane waves
    # over one direction circle gives a standing wave: it already contains equal outgoing and
    # incoming radiation, so the net flux (and any sense of "expansion") cancels exactly.
    for r in [0.5, 3.7, 8.0, 15.3]
        h1 = hankelh1(0, r)
        h2 = hankelh2(0, r)
        j0 = besselj(0, r)
        @assert h2 ≈ conj(h1) "H_0^(2) should be the complex conjugate of H_0^(1)"
        @assert abs(j0 - real((h1 + h2) / 2)) < 1e-12
        @assert abs(imag((h1 + h2) / 2)) < 1e-12
    end

    "Standing/outgoing waves self-checks passed"
end

# ╔═╡ 3e7f180e-808c-4e66-b7cc-ee5d54213685
md"""
## Widget
"""

# ╔═╡ f88fd389-fc11-492b-b7c2-cb9708124baa
begin
    struct CircularWavefrontInput
        kr::Float64
        N::Int
    end
    CircularWavefrontInput(; kr=8.0, N=8) = CircularWavefrontInput(kr, N)

    Base.get(w::CircularWavefrontInput) = Dict{String,Any}("kr" => w.kr, "N" => w.N)

    function Base.show(io::IO, ::MIME"text/html", w::CircularWavefrontInput)
        write(io, """
        <div id="cwwidget">
        <style>
        pluto-cell:has(#cwwidget) { width: min(80vw, 1150px) !important;
          margin-left: calc((100% - min(80vw, 1150px)) / 2) !important; }
        #cwwidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #cwwidget .cw-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #cwwidget .cw-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #cwwidget .cw-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #cwwidget .cw-primary{display:flex;gap:14px;flex-wrap:wrap;justify-content:center;align-items:flex-start;margin-bottom:10px}
        #cwwidget .cw-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #cwwidget .cw-panel-title{font-size:14px;font-weight:700;color:#e5e7eb;text-align:center;margin-bottom:6px}
        #cwwidget .cw-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:6px}
        #cwwidget canvas{display:block}
        #cwwidget .cw-controls-row{display:flex;justify-content:center}
        #cwwidget .cw-control-group{flex:1 1 700px;max-width:1100px;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 14px}
        #cwwidget .cw-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #cwwidget .cw-control-row{display:grid;grid-template-columns:150px minmax(0,1fr) 90px;gap:8px;align-items:center;margin:6px 0}
        #cwwidget .cw-control-row label{font-size:13px;color:#9ca3af}
        #cwwidget .cw-control-row input[type=range]{width:100%;min-width:0}
        #cwwidget .cw-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #cwwidget select{background:#0b0b0b;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:4px;width:100%}
        #cwwidget .cw-controls-row{flex-wrap:wrap;gap:14px}
        #cwwidget .cw-btn{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;padding:6px 12px;font-size:14px;cursor:pointer}
        </style>

        <div class="cw-title">
          <div class="cw-title-desc">Drag N to sum more plane-wave directions, then press Play to watch it move.</div>
          <div class="cw-title-hint">dashed circle = the aperture (kr) &middot; N=1 travels, large N pulses in place &middot; pick a scenario below, or drag freely</div>
        </div>

        <div class="cw-primary">
          <div>
            <div class="cw-panel-title">Plane-Wave Directions</div>
            <div class="cw-panel"><canvas id="cw-directions"></canvas></div>
          </div>
          <div>
            <div class="cw-panel-title">Summed Field (N directions)</div>
            <div class="cw-panel"><canvas id="cw-summed"></canvas></div>
          </div>
          <div>
            <div class="cw-panel-title">Exact Target J0(kr)</div>
            <div class="cw-panel"><canvas id="cw-exact"></canvas></div>
          </div>
        </div>
        <div class="cw-caption" id="cw-caption"></div>

        <div class="cw-controls-row">
          <div class="cw-control-group">
            <div class="cw-control-title">Station-Array Scenario</div>
            <div class="cw-control-row"><label>preset</label>
              <select id="cw-preset">
                <option value="8,3">3 stations (triangle, N &lt; kr)</option>
                <option value="8,8" selected>8 stations (small array, N &asymp; kr)</option>
                <option value="8,20">20 stations (converged, N &gt; kr)</option>
                <option value="15,60">60 stations (dense ring array)</option>
              </select>
            </div>
            <div class="cw-control-row"><label>N (directions)</label><input type="range" id="cw-N" min="1" max="150" step="1" value="$(w.N)"><span class="cw-value" id="cw-N-v"></span></div>
            <div class="cw-control-row"><label>kr (aperture)</label><input type="range" id="cw-kr" min="0" max="1.5" step="0.01" value="$(log10(w.kr))"><span class="cw-value" id="cw-kr-v"></span></div>
          </div>
          <div class="cw-control-group" style="flex:0 1 220px;">
            <div class="cw-control-title">Time Animation</div>
            <div style="display:flex;gap:8px;margin:6px 0">
              <button id="cw-play" class="cw-btn" type="button">&#9654; Play</button>
              <button id="cw-reset-t" class="cw-btn" type="button">Reset</button>
            </div>
            <div class="cw-caption" id="cw-phase-caption" style="margin-top:2px">phase &omega;t = 0.00 rad</div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        let state = { kr: $(w.kr), N: $(w.N) };
        let pushed = null; // {halfwidth, summedRe, summedIm, exact, maxerr} from Julia
        let commitInFlight = false;
        const OMEGA = 1.3; // rad/s of animated phase -- a full 2*pi cycle takes ~4.8s
        let animT = 0, playing = false, lastTs = null, rafId = null;

        const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 1150);
        const GAP = 14, PANEL_OVERHEAD = 18;
        const PW = Math.max(220, Math.floor((availW - 2*GAP - 3*PANEL_OVERHEAD) / 3));
        const PH = PW;
        const DPR = window.devicePixelRatio || 1;

        function hidpi(canvas, ctx, w, h){
          canvas.width = Math.round(w*DPR); canvas.height = Math.round(h*DPR);
          canvas.style.width = w+'px'; canvas.style.height = h+'px';
          ctx.setTransform(DPR,0,0,DPR,0,0);
        }
        const dirCv = par.querySelector('#cw-directions'), dirCtx = dirCv.getContext('2d');
        hidpi(dirCv, dirCtx, PW, PH);
        const summedCv = par.querySelector('#cw-summed'), summedCtx = summedCv.getContext('2d');
        hidpi(summedCv, summedCtx, PW, PH);
        const exactCv = par.querySelector('#cw-exact'), exactCtx = exactCv.getContext('2d');
        hidpi(exactCv, exactCtx, PW, PH);

        function emit(){
          commitInFlight = true;
          par.value = { kr: state.kr, N: state.N };
          par.dispatchEvent(new CustomEvent('input'));
        }
        function throttledEmit(){ if(!commitInFlight) emit(); }

        // diverging colormap, blue=+1, red=-1 -- this repo's usual convention
        function velColor(v){
          const t = Math.max(-1, Math.min(1, v));
          if(t >= 0) return [Math.round(255*(1-t)), Math.round(255*(1-t)), 255];
          const s = -t;
          return [255, Math.round(255*(1-s)), Math.round(255*(1-s))];
        }

        function drawAperture(ctx, halfwidth){
          const scale = PW/(2*halfwidth);
          const cx = PW/2, cy = PH/2, r = state.kr*scale;
          ctx.strokeStyle = '#f3f4f6'; ctx.setLineDash([5,4]); ctx.lineWidth = 1.5;
          ctx.beginPath(); ctx.arc(cx, cy, r, 0, 2*Math.PI); ctx.stroke();
          ctx.setLineDash([]);
        }

        // Re[Z e^{-i*omega*t}] = Re(Z)cos(omega t) + Im(Z)sin(omega t) -- every plane wave shares
        // the same omega, so this closed-form combination is the *exact* time evolution, applied
        // fresh to the same Julia-computed grids every frame (no re-computation needed).
        function timeField(re, im, T){
          const c = Math.cos(T), s = Math.sin(T);
          const out = new Float64Array(re.length);
          for(let i=0;i<re.length;i++) out[i] = re[i]*c + (im ? im[i]*s : 0);
          return out;
        }

        function drawHeatmap(ctx, grid, halfwidth){
          ctx.clearRect(0,0,PW,PH);
          if(!grid){ return; }
          const n = Math.round(Math.sqrt(grid.length));
          const cvBuf = document.createElement('canvas'); cvBuf.width = n; cvBuf.height = n;
          const bctx = cvBuf.getContext('2d');
          const img = bctx.createImageData(n, n);
          for(let idx=0; idx<grid.length; idx++){
            const c = velColor(grid[idx]);
            const k = 4*idx;
            img.data[k]=c[0]; img.data[k+1]=c[1]; img.data[k+2]=c[2]; img.data[k+3]=255;
          }
          bctx.putImageData(img, 0, 0);
          ctx.imageSmoothingEnabled = false;
          ctx.drawImage(cvBuf, 0, 0, n, n, 0, 0, PW, PH);
          drawAperture(ctx, halfwidth);
        }

        function drawDirections(){
          const ctx = dirCtx;
          ctx.clearRect(0,0,PW,PH);
          const cx = PW/2, cy = PH/2, R = Math.min(PW,PH)*0.42;
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          ctx.beginPath(); ctx.arc(cx,cy,R,0,2*Math.PI); ctx.stroke();
          ctx.strokeStyle = '#38bdf8'; ctx.fillStyle = '#38bdf8'; ctx.lineWidth = 1.5;
          for(let j=0;j<state.N;j++){
            const phi = 2*Math.PI*j/state.N;
            const x1 = cx + 0.68*R*Math.cos(phi), y1 = cy + 0.68*R*Math.sin(phi);
            const x2 = cx + R*Math.cos(phi), y2 = cy + R*Math.sin(phi);
            ctx.beginPath(); ctx.moveTo(x1,y1); ctx.lineTo(x2,y2); ctx.stroke();
            ctx.beginPath(); ctx.arc(x2,y2,2.4,0,2*Math.PI); ctx.fill();
          }
          ctx.fillStyle = '#9ca3af'; ctx.font = '12px sans-serif'; ctx.textAlign = 'center';
          ctx.fillText(state.N + (state.N===1 ? ' direction' : ' directions'), cx, PH-8);
        }

        function draw(){
          const halfwidth = pushed ? pushed.halfwidth : Math.max(1.3*state.kr, 10.0);
          drawDirections();
          if(pushed){
            drawHeatmap(summedCtx, timeField(pushed.summedRe, pushed.summedIm, animT), halfwidth);
            drawHeatmap(exactCtx, timeField(pushed.exact, null, animT), halfwidth);
          } else {
            drawHeatmap(summedCtx, null, halfwidth);
            drawHeatmap(exactCtx, null, halfwidth);
          }
          const cap = par.querySelector('#cw-caption');
          cap.textContent = pushed ? ('max |summed - exact| at t=0, inside the aperture: ' + pushed.maxerr.toExponential(2)) : '';
          par.querySelector('#cw-phase-caption').textContent = 'phase ωt = ' + animT.toFixed(2) + ' rad';
        }

        function tick(ts){
          if(!playing) return;
          if(lastTs === null) lastTs = ts;
          const dt = Math.max(0, Math.min(0.1, (ts - lastTs) / 1000)); // clamp huge/backwards gaps
          lastTs = ts;
          animT += dt * OMEGA;
          draw();
          rafId = requestAnimationFrame(tick);
        }

        function syncControls(){
          par.querySelector('#cw-N-v').textContent = String(state.N);
          par.querySelector('#cw-kr-v').textContent = state.kr.toFixed(state.kr<10?2:1);
        }

        function onSlider(event){
          const id = event.target.id;
          if(id === 'cw-N') state.N = Math.round(Number(event.target.value));
          else if(id === 'cw-kr') state.kr = Math.pow(10, Number(event.target.value));
          else return;
          syncControls();
          draw();
          throttledEmit();
        }
        par.querySelectorAll('input[type=range]').forEach(el => el.addEventListener('input', onSlider));

        par.querySelector('#cw-preset').addEventListener('change', event => {
          const [krStr, nStr] = event.target.value.split(',');
          state.kr = Number(krStr); state.N = Number(nStr);
          par.querySelector('#cw-kr').value = Math.log10(state.kr);
          par.querySelector('#cw-N').value = state.N;
          syncControls();
          draw();
          emit();
        });

        par.addEventListener('cw-results', event => {
          pushed = event.detail || null;
          commitInFlight = false;
          draw();
        });

        par.querySelector('#cw-play').addEventListener('click', () => {
          playing = !playing;
          par.querySelector('#cw-play').innerHTML = playing ? '&#10074;&#10074; Pause' : '&#9654; Play';
          if(playing){ lastTs = null; rafId = requestAnimationFrame(tick); }
        });
        par.querySelector('#cw-reset-t').addEventListener('click', () => {
          playing = false; animT = 0; lastTs = null;
          par.querySelector('#cw-play').innerHTML = '&#9654; Play';
          draw();
        });

        syncControls();
        draw();
        }
        </script>
        """)
    end

    const _cw_ready = true
end

# ╔═╡ fa29d377-b796-4489-95f1-76181f9d9159
begin
    _cw_ready
    WideCell(@bind cw CircularWavefrontInput(); max_width=1150)
end

# ╔═╡ 8347c858-d52c-46d5-8d65-0f3e9190bbc9
begin
    struct CwPush
        halfwidth::Float64
        summedRe::String
        summedIm::String
        exact::String
        maxerr::Float64
    end
    function Base.show(io::IO, ::MIME"text/html", p::CwPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('cwwidget');
        if(w){
          w.dispatchEvent(new CustomEvent('cw-results', { detail: {
            halfwidth: $(p.halfwidth),
            summedRe: [$(p.summedRe)],
            summedIm: [$(p.summedIm)],
            exact: [$(p.exact)],
            maxerr: $(p.maxerr),
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ 5a574765-529e-49a6-a4ae-4a911e4bbc4a
begin
    local halfwidth = max(1.3 * cw["kr"], 10.0)
    local ngrid = 140
    local xs = range(-halfwidth, halfwidth; length=ngrid)
    local N = round(Int, cw["N"])
    local summed_c = [summed_plane_wave_field_complex(x, y, N) for y in xs, x in xs]
    local summed_re = real.(summed_c)
    local summed_im = imag.(summed_c)
    local exact = [besselj(0, hypot(x, y)) for y in xs, x in xs]
    local maxerr = 0.0
    for (iy, y) in enumerate(xs), (ix, x) in enumerate(xs)
        if x^2 + y^2 <= cw["kr"]^2
            maxerr = max(maxerr, abs(summed_re[iy, ix] - exact[iy, ix]))
        end
    end
    CwPush(halfwidth, join(vec(permutedims(summed_re)), ","), join(vec(permutedims(summed_im)), ","),
        join(vec(permutedims(exact)), ","), maxerr)
end

# ╔═╡ 012115ec-82c1-4dfd-909c-e593f9d042c7
md"""
## References
- Companion notebook in this repo: `Jacobi-Anger-expansion.jl` (Basics) builds a single plane
  wave from cylindrical harmonics -- the reverse of the sum performed here.
- `Lamb_problem.jl` (Seismic Sources and Faulting) uses the more general Sommerfeld/Weyl integral
  (integrating over wavenumber magnitude, not direction) to build a point source's radiation
  into a layered medium -- see its "Conical Wave" section.
- Abramowitz & Stegun, *Handbook of Mathematical Functions*, §9.1 (integral representations of
  Bessel functions) and §9.1.3-9.1.6 (Hankel functions and their relation to `J_0`, `Y_0`).
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bessels = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Bessels = "~0.2.8"
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.0"
manifest_format = "2.1"
project_hash = "da0f5d56ff09b918d2e3ac3de0d1f5c24047003e"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
registries = "General"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

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
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
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
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.18.0+1"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.103+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
registries = "General"
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
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

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
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "908fec9df6c5de98548ead82a468c95ccf6cd263"
registries = "General"
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

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.67.1+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╟─7e266e4d-e0fa-45f4-ad8b-b6c4f784a392
# ╠═ffdeeaca-a54c-47a2-9525-6a8875368703
# ╟─70042daf-6ff4-49c1-b50a-d9e8df536776
# ╟─fa29d377-b796-4489-95f1-76181f9d9159
# ╟─931da2f3-774e-42b6-9f92-37ccb250479c
# ╟─5a574765-529e-49a6-a4ae-4a911e4bbc4a
# ╟─62bdf2c0-7aed-467d-8f40-b5d7761a0573
# ╟─2e0c5129-76d4-4200-b2d3-471e404ef61b
# ╟─3a1e3154-63b1-4661-ab44-601e47874271
# ╠═57218488-7de4-487d-bb2f-6e94b4055458
# ╠═5dbfb284-a12c-49cf-ae60-6807513fedd1
# ╠═27f7ef58-1d43-4bda-be5e-a455aee66bab
# ╟─5fd8e68c-3302-4911-8831-8f65912c6e59
# ╠═92c94385-480b-4bcb-9f0a-bb99a795e367
# ╠═b4f123af-2e6e-46b4-a57e-ea7ff710d29c
# ╠═2df09cda-c40f-4b8c-a254-c18a857c652e
# ╟─3e7f180e-808c-4e66-b7cc-ee5d54213685
# ╠═f88fd389-fc11-492b-b7c2-cb9708124baa
# ╠═8347c858-d52c-46d5-8d65-0f3e9190bbc9
# ╟─012115ec-82c1-4dfd-909c-e593f9d042c7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
