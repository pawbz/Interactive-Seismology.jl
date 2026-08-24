### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Building a Plane Wave From Cylindrical Waves"
#> tags = ["basics"]
#> layout = "layout.jlhtml"
#> description = "Watch the Jacobi-Anger expansion build a plane wave out of pinwheel-shaped cylindrical harmonics, term by term."

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

# ╔═╡ 22fcbc39-b1b8-4a22-9310-eb97923ab4b0
begin
    using Bessels
    using PlutoUI
end

# ╔═╡ 4643b924-2b62-4970-8cc0-bff87395c309
TableOfContents()

# ╔═╡ a918feb0-4e91-4300-9fc3-8c18c0cfd21f
md"""
# How a Plane Wave Is Built From Cylindrical Waves

A plane wave is the simplest possible wavefront: dead-straight, parallel stripes moving in
one direction. Yet the natural building blocks centered on a point -- the cylindrical
harmonics `Jₙ(kr)e^{inθ}` used throughout this repo's `multipole-scattering.jl` notebook --
look nothing like that. Each one is a pinwheel: `2n` angular lobes riding a
radially-oscillating Bessel envelope. The **Jacobi–Anger expansion**

```math
e^{ikx} = \sum_{n=-\infty}^{\infty} i^n J_n(kr)\, e^{in\theta}
```

says a plane wave is an *exact* infinite sum of these pinwheels. Watching that sum actually
happen -- pinwheels in, straight stripes out -- is the point of this notebook.
"""

# ╔═╡ 81081809-9c45-420f-aba3-af1c07ed3e2a
md"""
## Reading the Three Panels

**Drag across the spectrum panel** (left) to choose which orders `n` get summed -- a click
selects one order, a drag selects a band `nLo:nHi`. Each stem's height is `` |J_n(kr)| ``,
the natural amplitude of order `n` at the current aperture: this *is* the plane wave's
"multipole spectrum," and its shape (rising, then collapsing once `n` exceeds `kr`) is
exactly why only finitely many terms end up mattering.

The **middle** panel shows exactly what's currently selected: a single pinwheel
`` \mathrm{Re}\big[i^nJ_n(kr)e^{in\theta}\big] `` when you've selected one order, or the sum
`` \mathrm{Re}\sum_{n=n_{\text{Lo}}}^{n_{\text{Hi}}}\big[i^nJ_n(\rho)e^{in\theta} +
i^{-n}J_{-n}(\rho)e^{-in\theta}\big] `` over the selected band when you've dragged a range
(each order `n>0` is paired with its mirror `-n`, since the two always appear together in
the real-valued sum). The **right** panel is the true plane wave, simply `cos(kx)`, plotted
in the same *k-scaled* coordinates -- horizontal axis `kx`, vertical axis `ky`, so distance
is measured in radians of phase, not physical length.

The dashed circle marks the **aperture**: a radius `kr` chosen from the preset below. Select
the band `0:N` for a growing `N` (the default behaviour) and watch a familiar picture
emerge: inside a growing central disk the middle and right panels already agree almost
perfectly, while outside it the sum still looks like a pinwheel, not a plane wave. That
disk's radius grows with `N` -- the self-check in the Appendix confirms it tracks `N`
almost exactly, matching the rule of thumb this repo's `multipole-scattering.jl` notebook
uses for its own multipole sums (**"you need roughly `N≈kr` terms"**). But you're not
limited to starting at `n=0` -- drag a band that *doesn't* include the low orders (say
`20:40`) and see that an isolated slice of the spectrum, on its own, doesn't look like a
plane wave at all, converged or not; the low-order terms are doing real, specific work.
"""

# ╔═╡ 5ab5ebd6-3223-4a8e-98be-bde38d7bd4bc
md"""
## Appendix
"""

# ╔═╡ ed0ee423-3a7d-4431-bf9a-d5f25cb94a77
md"""
## Partial-Sum Physics
"""

# ╔═╡ cf380bc3-e9e2-4a54-ba9c-c89d7e2f26e8
"""
	planewave_partial_sum(X, Y, nLo, nHi)

The Jacobi–Anger sum restricted to orders `n=nLo:nHi` and their mirror orders `-n`
(every `n>0` term always appears paired with `-n` in the real-valued sum; `n=0` is its
own mirror and is not double-counted),
`` \\mathrm{Re}\\sum_{n=n_{\\text{Lo}}}^{n_{\\text{Hi}}}\\big[i^nJ_n(\\rho)e^{in\\theta} +
i^{-n}J_{-n}(\\rho)e^{-in\\theta}\\big] ``, evaluated at a point `(X,Y)` given in
*k-scaled* coordinates (`X=kx`, `Y=ky`, so `` \\rho=\\sqrt{X^2+Y^2} `` is literally `kr`
for this point and `` \\theta=\\mathrm{atan}(Y,X) ``).

With `nLo=0`, this is the full truncated Jacobi–Anger sum through order `nHi`: as
`` n_{\\text{Hi}}\\to\\infty `` it converges to the exact plane wave `cos(X)` everywhere,
and for finite `nHi` it only converges within roughly `` \\rho\\lesssim n_{\\text{Hi}} ``
(checked below) -- exactly the "`N≈kr` terms needed" rule of thumb used throughout
`multipole-scattering.jl`, here isolated from any scattering physics so it can be seen on
its own. With `nLo>0`, this instead isolates a single *band* of the spectrum, which need
not resemble a plane wave at all -- see the widget's spectrum panel.
"""
function planewave_partial_sum(X, Y, nLo, nHi)
    rho = hypot(X, Y)
    theta = atan(Y, X)
    s = zero(ComplexF64)
    for n in nLo:nHi
        s += cis(n * pi / 2) * besselj(n, rho) * cis(n * theta)
        if n != 0
            s += cis(-n * pi / 2) * besselj(-n, rho) * cis(-n * theta)
        end
    end
    return real(s)
end

# ╔═╡ dfe6f001-3d1f-4c93-9f1a-7c1a2c9e0b3a
"""
	spectrum_magnitudes(kr, nmax)

`` |J_n(kr)| `` for `n=0:nmax` -- the plane wave's "multipole spectrum" at aperture `kr`:
how strongly each order actually contributes. Its shape (rising through small `n`, then
collapsing once `n` exceeds `kr`) is exactly why only finitely many terms end up mattering
in [`planewave_partial_sum`](@ref), and is what the widget's spectrum panel plots as a
stem chart to drag-select a band of orders from.
"""
spectrum_magnitudes(kr, nmax) = [abs(besselj(n, kr)) for n in 0:nmax]

# ╔═╡ be9945a3-d6dc-406c-b7e6-b13bfca30ec6
md"""
### Verifying the Partial Sum
"""

# ╔═╡ b51453e8-379b-4d45-a493-8836ac0db46f
let
    # 1. as N grows, a fixed point converges to the exact plane wave (band 0:N)
    rho, theta = 11.8, 0.7
    X, Y = rho * cos(theta), rho * sin(theta)
    errs = [abs(planewave_partial_sum(X, Y, 0, N) - cos(X)) for N in [5, 15, 25, 40]]
    println("Convergence at rho=$rho as N grows: ", errs)
    @assert issorted(errs, rev=true) "error should shrink monotonically as N grows past rho"
    @assert errs[end] < 1e-6

    # 2. for a fixed N, points inside rho~N are accurate and points well outside are not --
    #    the "growing disk of agreement" the notebook's prose describes
    N = 15
    err_inside = abs(planewave_partial_sum(2.0 * cos(1.1), 2.0 * sin(1.1), 0, N) - cos(2.0 * cos(1.1)))
    err_outside = abs(planewave_partial_sum(25.0 * cos(1.1), 25.0 * sin(1.1), 0, N) - cos(25.0 * cos(1.1)))
    println("N=$N: error well inside disk = $err_inside, error well outside = $err_outside")
    @assert err_inside < 1e-6
    @assert err_outside > 0.1

    # 3. an isolated band nLo>0 does not reduce to the nLo=0 case (it's a genuinely
    #    different, non-plane-wave quantity) -- sanity check the two disagree
    band_only = planewave_partial_sum(X, Y, 20, 40)
    from_zero = planewave_partial_sum(X, Y, 0, 40)
    println("band 20:40 vs 0:40 at rho=$rho differ by ", abs(band_only - from_zero))
    @assert abs(band_only - from_zero) > 0.1

    "Partial-sum self-checks passed"
end

# ╔═╡ 4eaa243a-36a2-4955-898c-5d328589008f
md"""
## Widget
"""

# ╔═╡ 92f754a7-27a2-45b2-9f32-45cfeee534f5
begin
    struct JacobiAngerInput
        kr::Float64
        selLo::Int
        selHi::Int
    end
    JacobiAngerInput(; kr=5.4978, selLo=0, selHi=10) = JacobiAngerInput(kr, selLo, selHi)

    Base.get(w::JacobiAngerInput) = Dict{String,Any}("kr" => w.kr, "selLo" => w.selLo, "selHi" => w.selHi)

    function Base.show(io::IO, ::MIME"text/html", w::JacobiAngerInput)
        write(io, """
        <div id="jawidget">
        <style>
        pluto-cell:has(#jawidget) { width: min(80vw, 1150px) !important;
          margin-left: calc((100% - min(80vw, 1150px)) / 2) !important; }
        #jawidget{font-family:sans-serif;color:#e5e7eb;width:100%;box-sizing:border-box}
        #jawidget .ja-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
          background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
        #jawidget .ja-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
        #jawidget .ja-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
        #jawidget .ja-primary{display:flex;gap:14px;flex-wrap:wrap;justify-content:center;align-items:flex-start;margin-bottom:10px}
        #jawidget .ja-panel{background:#000;border:1px solid #374151;border-radius:6px;padding:8px}
        #jawidget .ja-panel-title{font-size:14px;font-weight:700;color:#e5e7eb;text-align:center;margin-bottom:6px}
        #jawidget .ja-caption{font-size:12px;color:#9ca3af;text-align:center;margin-top:6px}
        #jawidget canvas{display:block}
        #jawidget #ja-spectrum{cursor:crosshair}
        #jawidget .ja-controls-row{display:flex;justify-content:center}
        #jawidget .ja-control-group{flex:1 1 700px;max-width:1100px;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 14px}
        #jawidget .ja-control-title{font-size:15px;font-weight:700;color:#e5e7eb;margin-bottom:6px}
        #jawidget .ja-control-row{display:grid;grid-template-columns:150px minmax(0,1fr) 90px;gap:8px;align-items:center;margin:6px 0}
        #jawidget .ja-control-row label{font-size:13px;color:#9ca3af}
        #jawidget .ja-control-row input[type=range]{width:100%;min-width:0}
        #jawidget .ja-value{font-size:13px;color:#e5e7eb;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
        #jawidget select{background:#0b0b0b;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:4px;width:100%}
        </style>

        <div class="ja-title">
          <div class="ja-title-desc">Drag across the spectrum to sum a band of pinwheel terms, and watch straight wavefronts emerge.</div>
          <div class="ja-title-hint">dashed circle = the aperture (kr) &middot; pick a real array scenario below, or drag kr freely</div>
        </div>

        <div class="ja-primary">
          <div>
            <div class="ja-panel-title">Multipole Spectrum</div>
            <div class="ja-panel"><canvas id="ja-spectrum"></canvas></div>
            <div class="ja-caption" id="ja-spectrum-caption"></div>
          </div>
          <div>
            <div class="ja-panel-title" id="ja-terms-title">Selected Terms</div>
            <div class="ja-panel"><canvas id="ja-partial"></canvas></div>
          </div>
          <div>
            <div class="ja-panel-title">Exact Plane Wave</div>
            <div class="ja-panel"><canvas id="ja-exact"></canvas></div>
          </div>
        </div>
        <div class="ja-caption" id="ja-caption"></div>

        <div class="ja-controls-row">
          <div class="ja-control-group">
            <div class="ja-control-title">Array Scenario</div>
            <div class="ja-control-row"><label>preset</label>
              <select id="ja-preset">
                <option value="1.2566,0,5">Tight local array (2 Hz S wave, 0.3 km aperture)</option>
                <option value="5.4978,0,10" selected>Broadband network (20 s Rayleigh wave, USArray TA spacing)</option>
                <option value="11.7810,0,17">Regional network (1 Hz Pn, 15 km aperture)</option>
                <option value="83.7758,0,100">Dense nodal array (5 Hz local, 8 km aperture)</option>
              </select>
            </div>
            <div class="ja-control-row"><label>kr (aperture)</label><input type="range" id="ja-kr" min="0" max="2.1" step="0.01" value="$(log10(w.kr))"><span class="ja-value" id="ja-kr-v"></span></div>
          </div>
        </div>
        </div>

        <script>
        {
        const par = currentScript.previousElementSibling;
        const NMAX_SPECTRUM = 120;
        let state = { kr: $(w.kr), selLo: $(w.selLo), selHi: $(w.selHi) };
        let pushed = null; // {halfwidth, grid, maxerr, spectrum} from Julia
        let commitInFlight = false;
        let specAnchor = null; // order index where a spectrum click/drag started, null when idle

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
        const spectrumCv = par.querySelector('#ja-spectrum'), spectrumCtx = spectrumCv.getContext('2d');
        hidpi(spectrumCv, spectrumCtx, PW, PH);
        const partialCv = par.querySelector('#ja-partial'), partialCtx = partialCv.getContext('2d');
        hidpi(partialCv, partialCtx, PW, PH);
        const exactCv = par.querySelector('#ja-exact'), exactCtx = exactCv.getContext('2d');
        hidpi(exactCv, exactCtx, PW, PH);

        function emit(){
          commitInFlight = true;
          par.value = { kr: state.kr, selLo: state.selLo, selHi: state.selHi };
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

        function drawExact(halfwidth){
          // putImageData bypasses the DPR transform set up in hidpi() and always writes in
          // physical device pixels, so the image buffer and loop must use canvas.width/height
          // (physical) here, not PW/PH (CSS pixels) -- using PW/PH left this only filling the
          // top-left quarter of the canvas at DPR=2.
          const ctx = exactCtx;
          const pw = exactCv.width, ph = exactCv.height;
          const img = ctx.createImageData(pw, ph);
          const scale = 2*halfwidth/PW;
          for(let j=0;j<ph;j++){
            for(let i=0;i<pw;i++){
              const X = (i/DPR - PW/2)*scale;
              const c = velColor(Math.cos(X));
              const k = 4*(j*pw+i);
              img.data[k]=c[0]; img.data[k+1]=c[1]; img.data[k+2]=c[2]; img.data[k+3]=255;
            }
          }
          ctx.putImageData(img, 0, 0);
          drawAperture(ctx, halfwidth);
        }

        function drawPartial(){
          const ctx = partialCtx;
          ctx.clearRect(0,0,PW,PH);
          if(!pushed){ return; }
          const n = Math.round(Math.sqrt(pushed.grid.length));
          const cvBuf = document.createElement('canvas'); cvBuf.width = n; cvBuf.height = n;
          const bctx = cvBuf.getContext('2d');
          const img = bctx.createImageData(n, n);
          for(let idx=0; idx<pushed.grid.length; idx++){
            const c = velColor(pushed.grid[idx]);
            const k = 4*idx;
            img.data[k]=c[0]; img.data[k+1]=c[1]; img.data[k+2]=c[2]; img.data[k+3]=255;
          }
          bctx.putImageData(img, 0, 0);
          ctx.imageSmoothingEnabled = false;
          ctx.drawImage(cvBuf, 0, 0, n, n, 0, 0, PW, PH);
          drawAperture(ctx, pushed.halfwidth);
        }

        function pxToN(px){
          const bandL = 30, bandR = 6;
          const plotW = PW - bandL - bandR;
          const n = Math.floor((px-bandL)/plotW*(NMAX_SPECTRUM+1));
          return Math.max(0, Math.min(NMAX_SPECTRUM, n));
        }

        function drawSpectrum(){
          const ctx = spectrumCtx, W = PW, H = PH;
          ctx.clearRect(0,0,W,H);
          ctx.strokeStyle = '#374151'; ctx.lineWidth = 1;
          ctx.strokeRect(0.5,0.5,W-1,H-1);
          if(!pushed){ return; }
          const spec = pushed.spectrum;
          const bandL = 30, bandR = 6, bandT = 6, bandB = 18;
          const plotW = W - bandL - bandR, plotH = H - bandT - bandB;
          const n = spec.length;
          const mx = Math.max(...spec, 1e-12);
          const stemW = plotW/n;
          ctx.fillStyle = '#facc15'; ctx.globalAlpha = 0.16;
          ctx.fillRect(bandL + state.selLo*stemW, bandT, (state.selHi-state.selLo+1)*stemW, plotH);
          ctx.globalAlpha = 1;
          for(let i=0;i<n;i++){
            const active = i>=state.selLo && i<=state.selHi;
            const px = bandL + (i+0.5)*stemW;
            const py = bandT + plotH - (spec[i]/mx)*plotH;
            ctx.strokeStyle = active ? '#facc15' : '#4b5563';
            ctx.lineWidth = active ? 2 : 1;
            ctx.beginPath(); ctx.moveTo(px, bandT+plotH); ctx.lineTo(px, py); ctx.stroke();
          }
          ctx.fillStyle = '#9ca3af'; ctx.font = '10px sans-serif'; ctx.textAlign = 'center';
          ctx.fillText('order n', bandL+plotW/2, H-4);
          ctx.save(); ctx.translate(10, bandT+plotH/2); ctx.rotate(-Math.PI/2);
          ctx.fillText('|J\\u2099(kr)|', 0, 0); ctx.restore();
        }

        function updateTermsTitle(){
          const el = par.querySelector('#ja-terms-title');
          el.textContent = (state.selLo === state.selHi)
            ? ('Term n=' + state.selLo)
            : ('Terms n=' + state.selLo + '\\u2013' + state.selHi);
        }

        function draw(){
          const halfwidth = pushed ? pushed.halfwidth : Math.max(1.3*state.kr, 10.0);
          drawSpectrum();
          drawExact(halfwidth);
          drawPartial();
          updateTermsTitle();
          const cap = par.querySelector('#ja-caption');
          cap.textContent = pushed ? ('difference from the exact plane wave, inside the aperture: ' + pushed.maxerr.toExponential(2)) : '';
          par.querySelector('#ja-spectrum-caption').textContent = 'click or drag to select which orders n to sum';
        }

        function syncControls(){
          par.querySelector('#ja-kr-v').textContent = state.kr.toFixed(state.kr<10?2:1);
        }

        function onSlider(event){
          const id = event.target.id;
          if(id === 'ja-kr') state.kr = Math.pow(10, Number(event.target.value));
          else return;
          syncControls();
          throttledEmit();
        }
        par.querySelectorAll('input[type=range]').forEach(el => el.addEventListener('input', onSlider));

        par.querySelector('#ja-preset').addEventListener('change', event => {
          const [krStr, loStr, hiStr] = event.target.value.split(',');
          state.kr = Number(krStr); state.selLo = Number(loStr); state.selHi = Number(hiStr);
          par.querySelector('#ja-kr').value = Math.log10(state.kr);
          syncControls();
          draw();
          emit();
        });

        function spectrumPointerX(ev){
          const rect = spectrumCv.getBoundingClientRect();
          return ev.clientX - rect.left;
        }
        spectrumCv.addEventListener('mousedown', ev => {
          const n = pxToN(spectrumPointerX(ev));
          specAnchor = n; state.selLo = n; state.selHi = n;
          draw(); throttledEmit();
        });
        window.addEventListener('mousemove', ev => {
          if(specAnchor === null) return;
          const n = pxToN(Math.max(0, Math.min(PW, spectrumPointerX(ev))));
          state.selLo = Math.min(specAnchor, n); state.selHi = Math.max(specAnchor, n);
          draw(); throttledEmit();
        });
        window.addEventListener('mouseup', () => {
          if(specAnchor === null) return;
          specAnchor = null;
          emit();
        });

        par.addEventListener('ja-results', event => {
          pushed = event.detail || null;
          commitInFlight = false;
          draw();
        });

        syncControls();
        draw();
        }
        </script>
        """)
    end

    const _ja_ready = true
end

# ╔═╡ 2e30cdbe-790b-457b-923a-e52348ac0140
begin
    _ja_ready
    WideCell(@bind ja JacobiAngerInput(); max_width=1150)
end

# ╔═╡ 6df5e7eb-7146-436d-9448-c2cabb6596da
begin
    struct JaPush
        halfwidth::Float64
        grid::String
        maxerr::Float64
        spectrum::String
    end
    function Base.show(io::IO, ::MIME"text/html", p::JaPush)
        write(io, """
        <script>
        {
        const w = document.getElementById('jawidget');
        if(w){
          w.dispatchEvent(new CustomEvent('ja-results', { detail: {
            halfwidth: $(p.halfwidth),
            grid: [$(p.grid)],
            maxerr: $(p.maxerr),
            spectrum: [$(p.spectrum)],
          }}));
        }
        }
        </script>
        """)
    end
end

# ╔═╡ b4bc0a91-dc47-49df-a5da-b2d0287e398f
begin
    local halfwidth = max(1.3 * ja["kr"], 10.0)
    local ngrid = 140
    local xs = range(-halfwidth, halfwidth; length=ngrid)
    local nLo, nHi = ja["selLo"], ja["selHi"]
    local grid = [planewave_partial_sum(x, y, nLo, nHi) for y in xs, x in xs]
    local maxerr = 0.0
    for (iy, y) in enumerate(xs), (ix, x) in enumerate(xs)
        if x^2 + y^2 <= ja["kr"]^2
            maxerr = max(maxerr, abs(grid[iy, ix] - cos(x)))
        end
    end
    local nmax_spectrum = 120
    local spectrum = spectrum_magnitudes(ja["kr"], nmax_spectrum)
    JaPush(halfwidth, join(vec(permutedims(grid)), ","), maxerr, join(spectrum, ","))
end

# ╔═╡ cc0f4687-22e2-4c5b-a688-93f540e50a07
md"""
## References
- Companion notebooks in this repo: `multipole-scattering.jl` (Scattering) uses this same
  expansion as the starting point for scattering off a cylinder; `array-beamforming.jl`
  (Imaging) is the practical flip side -- once an array's aperture needs far more multipole
  terms than is useful (the "Dense nodal array" preset above), plane-wave/ray beamforming
  is used directly instead of a multipole sum.
- Abramowitz & Stegun, *Handbook of Mathematical Functions*, §9.1.44 (Jacobi–Anger expansion).
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

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "da0f5d56ff09b918d2e3ac3de0d1f5c24047003e"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
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
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

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
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
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
# ╟─22fcbc39-b1b8-4a22-9310-eb97923ab4b0
# ╠═4643b924-2b62-4970-8cc0-bff87395c309
# ╟─a918feb0-4e91-4300-9fc3-8c18c0cfd21f
# ╟─2e30cdbe-790b-457b-923a-e52348ac0140
# ╟─81081809-9c45-420f-aba3-af1c07ed3e2a
# ╟─b4bc0a91-dc47-49df-a5da-b2d0287e398f
# ╟─5ab5ebd6-3223-4a8e-98be-bde38d7bd4bc
# ╟─ed0ee423-3a7d-4431-bf9a-d5f25cb94a77
# ╠═cf380bc3-e9e2-4a54-ba9c-c89d7e2f26e8
# ╠═dfe6f001-3d1f-4c93-9f1a-7c1a2c9e0b3a
# ╟─be9945a3-d6dc-406c-b7e6-b13bfca30ec6
# ╠═b51453e8-379b-4d45-a493-8836ac0db46f
# ╟─4eaa243a-36a2-4955-898c-5d328589008f
# ╠═92f754a7-27a2-45b2-9f32-45cfeee534f5
# ╠═6df5e7eb-7146-436d-9448-c2cabb6596da
# ╟─cc0f4687-22e2-4c5b-a688-93f540e50a07
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
