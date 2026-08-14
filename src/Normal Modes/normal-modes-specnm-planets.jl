### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "PREM Earth Free Oscillations"
#> tags = ["normalmodes"]
#> layout = "layout.jlhtml"
#> description = "Browse precomputed free-oscillation spectra for Earth, Mars, the Moon, and Europa"

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

# ╔═╡ 9f419394-528a-4bde-98a5-d62787c17fa8
using PlutoUI

# ╔═╡ cd3c5ff7-7992-41f7-9a5e-29934a1355ca
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ bdae265d-3a96-4ecc-a6dd-c166357e801c
md"""
# 🌍 Exploring Free Oscillations

**Free oscillations** describe how a whole planet -- or moon -- rings after a large
disturbance (e.g. a giant earthquake), the same way a bell rings after being struck.
This notebook browses those modes for four bodies, computed ahead of time by
**`specnm`**, a spectral-element solver for the radial gravito-elastic ODEs.

Picking a body below only *loads* its result -- nothing solves live in this notebook.
Every mode was solved once, offline, up to `fmax = 0.01` Hz, by
`scripts/precompute_specnm.py`; re-run that script (not this notebook) if a model
file changes or a different `fmax` is needed.

`specnm`:
A spectral element approach to computing normal modes,
J Kemper, M van Driel, F Munch, A Khan, D Giardini,
Geophysical Journal International, Volume 229, Issue 2, May 2022, Pages 915–932.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ f2a1c3d4-7e6b-4a5c-9d8e-1b2c3d4e5f6a
md"""
##### What are spheroidal and toroidal modes?

A planet's free oscillations split into two independent families, both indexed by
angular degree `l` (spherical-harmonic degree, more nodal lines as `l` grows) and
overtone number `n` (radial order -- `n=0` is the fundamental, higher `n` adds radial
nodes):

- **Spheroidal (`ₙSₗ`)**: motion has both radial and horizontal components, coupled
  to gravity and compression -- shown as **U** (radial displacement), **V**
  (horizontal displacement), and, for a self-gravitating solve, **P** (perturbation
  to the gravitational potential).
- **Toroidal (`ₙTₗ`)**: purely horizontal, rotational shearing motion with no radial
  component and no coupling to gravity -- shown as **W** (horizontal displacement).
  Toroidal modes don't exist in a fluid (they need shear rigidity), so they vanish
  wherever `vs = 0` in the planet's model.

The eigenfunction panel plots each component versus depth for whichever mode you
click in the spectrum -- their zero-crossings are exactly the mode's radial nodes,
and the amplitude near the surface versus near the core tells you how deeply that
mode samples the planet.
"""

# ╔═╡ 9b0a9a80-e65f-4385-a377-372e408b19ad
md"## Appendix"

# ╔═╡ dc1c3b73-27d9-40cd-a33d-604ac6990daf
"""
The four bodies precomputed by `scripts/precompute_specnm.py` and available to
browse here -- must match that script's `DEFAULT_MODELS` and the JSON files under
`src/assets/data/specnm_precomputed/`.
"""
const SPECNM_PRECOMPUTED_MODELS = ["prem_ani", "mars", "vpremoon", "europa"]

# ╔═╡ 297386a9-d8e1-4da9-9bed-aca3ef653a11
md"### Loading Precomputed Results"

# ╔═╡ 383d55aa-01cb-4068-8140-39f52a5d7884
md"### The Interactive Widgets"

# ╔═╡ 12b77c6d-5a9e-43b7-a34c-8de94eb6003d
begin
    """
        SpecnmModelPicker(; model="prem_ani")

    Picks which precomputed body to browse. There's no Compute gate here (unlike
    the live-solving widget this replaces) -- loading a body is just a file read,
    so every change to the dropdown commits immediately.
    """
    struct SpecnmModelPicker
        model::String
    end

    SpecnmModelPicker(; model="prem_ani") = SpecnmModelPicker(model)

    Base.get(w::SpecnmModelPicker) = Dict{String,Any}("model" => w.model)

    const SPECNM_MODEL_LABELS = Dict(
        "prem_ani" => "Earth (PREM, anisotropic)",
        "mars" => "Mars",
        "vpremoon" => "Moon (VPREMOON)",
        "europa" => "Europa",
    )

    function Base.show(io::IO, ::MIME"text/html", w::SpecnmModelPicker)
        options = join(["""<option value="$m"$(m == w.model ? " selected" : "")>$(get(SPECNM_MODEL_LABELS, m, m))</option>""" for m in SPECNM_PRECOMPUTED_MODELS], "\n          ")
        write(io, """
<div id="spmwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#spmwidget){width:min(70vw,700px)!important;margin-left:calc((100% - min(70vw,700px))/2)!important}
    #spmwidget{width:100%;box-sizing:border-box;color:#d1d5db;font:14px sans-serif}
    #spmwidget .spm-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #spmwidget .spm-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #spmwidget .spm-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #spmwidget .spm-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px;width:100%}
    #spmwidget .spm-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:16px}
    #spmwidget select{width:100%;background:#111827;color:#e5e7eb;border:1px solid #374151;border-radius:4px;padding:8px 6px;font-size:15px}
  </style>
  <div class="spm-title">
    <div class="spm-title-desc">Pick a body to browse its precomputed normal modes.</div>
    <div class="spm-title-hint">precomputed offline at fmax = 0.01 Hz &middot; switching loads instantly, nothing solves here</div>
  </div>
  <div class="spm-control-group">
    <div class="spm-control-title">Body</div>
    <select id="spm-model">
          $(options)
    </select>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  par.querySelector('#spm-model').addEventListener('change', e=>{
    par.value = {model: e.target.value}
    par.dispatchEvent(new CustomEvent('input'))
  })
</script>
""")
    end

    # Forces the correct execution order on a fresh restart (Pluto's static
    # dependency analysis doesn't reliably detect that the bind cell above depends
    # on this cell defining SpecnmModelPicker) -- same pattern used throughout this
    # repo's other widgets, see pluto-widget-SKILL.md.
    const _spm_ready = true
end

# ╔═╡ 3508d355-720a-410a-bbd1-50e3779e77ee
begin
	_spm_ready
	@bind spm_pick SpecnmModelPicker()
end

# ╔═╡ 473e3e64-609c-41d2-ac77-b17943011eb2
model_name = spm_pick isa AbstractDict ? spm_pick["model"] : "prem_ani"

# ╔═╡ 1b958497-3e50-46a1-9237-022ddbe8aabb
"""
Read the precomputed JSON straight off disk as text -- it's already in exactly the
shape [`SpecnmBrowseView`](@ref) expects (`spheroidal`/`toroidal`/`earth_model`/
`model_name` keys), built once by `scripts/precompute_specnm.py`, so there's nothing
left to parse or reshape here.
"""
precomputed_json = let
	path = joinpath(@__DIR__, "..", "assets", "data", "specnm_precomputed", model_name * ".json")
	isfile(path) ? read(path, String) : nothing
end

# ╔═╡ f018bd73-8f46-448c-9afb-ec47c7ea705b
# Push the precomputed payload to SpecnmBrowseView unchanged -- the widget only ever
# JSON.parses this string client-side, it never round-trips back through Julia.
if precomputed_json === nothing
	Markdown.MD(Markdown.Admonition("danger", "No precomputed data", [Markdown.Paragraph(
		"No precomputed JSON found for `$model_name`. Run `/opt/miniconda3/bin/python3 scripts/precompute_specnm.py $model_name` first."
	)]))
else
	HTML("""<script>window.dispatchEvent(new CustomEvent('specnm-browse-data', {detail: $(repr(precomputed_json))}));</script>""")
end

# ╔═╡ a9bbb105-3138-481d-9a73-0ec2b8303c36
begin
    """
        SpecnmBrowseView()

    The cheap, live-browsing widget. Entirely fed by a `CustomEvent` pushed once
    after a body is picked (see the push cell above); every interaction below
    (spheroidal/toroidal toggle, hover, click-to-select a mode) is then pure
    client-side array indexing with zero further Julia calls. Three canvases: the
    mode spectrum (angular order `l` vs frequency, colored by overtone `n`), the
    selected mode's eigenfunction(s) vs depth, and the solved body's density/
    velocity profile vs depth.
    """
    struct SpecnmBrowseView end

    function Base.show(io::IO, ::MIME"text/html", w::SpecnmBrowseView)
        write(io, """
<div id="sbvwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#sbvwidget){width:min(90vw,1400px)!important;margin-left:calc((100% - min(90vw,1400px))/2)!important}
    #sbvwidget{width:100%;box-sizing:border-box;color:#d1d5db;font:14px sans-serif}
    #sbvwidget .sbv-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #sbvwidget .sbv-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #sbvwidget .sbv-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #sbvwidget .sbv-workspace{display:flex;gap:12px;align-items:flex-start;justify-content:center;width:100%;flex-wrap:wrap}
    #sbvwidget canvas{background:#000;border:1px solid #374151;border-radius:6px;display:block;cursor:crosshair}
    #sbvwidget .sbv-controls{width:100%;margin-top:8px;display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:8px;font:14px sans-serif}
    #sbvwidget .sbv-control-group{box-sizing:border-box;background:#050505;border:1px solid #2f3744;border-radius:6px;padding:10px 12px}
    #sbvwidget .sbv-control-title{font-weight:700;color:#e5e7eb;margin-bottom:8px;font-size:18px}
    #sbvwidget .sbv-actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    #sbvwidget label{color:#d1d5db}
  </style>
  <div class="sbv-title">
    <div class="sbv-title-desc">Every mode from the last Compute, browsed for free -- no re-solve needed.</div>
    <div class="sbv-title-hint">click a point in the spectrum to see its eigenfunction &middot; hover to identify a mode &middot; toggle spheroidal/toroidal</div>
  </div>
  <div class="sbv-workspace">
    <canvas id="sbvspectrum" width="440" height="440"></canvas>
    <canvas id="sbveigen" width="260" height="440"></canvas>
    <canvas id="sbvmodel" width="260" height="440"></canvas>
  </div>
  <div class="sbv-controls">
    <div class="sbv-control-group">
      <div class="sbv-control-title">Mode Family</div>
      <div class="sbv-actions">
        <label><input type="radio" name="sbv-type" id="sbv-type-s" checked> Spheroidal</label>
        <label><input type="radio" name="sbv-type" id="sbv-type-t"> Toroidal</label>
      </div>
    </div>
    <div class="sbv-control-group">
      <div class="sbv-control-title">Status</div>
      <div id="sbv-status" style="font-size:13px;color:#9ca3af">waiting for a body to load</div>
    </div>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.9, par.clientWidth || window.innerWidth*0.9, 1400)
  const totalW = Math.max(760, availW)
  const SPW = Math.round(Math.min(totalW*0.5, 560)), SIDEW = Math.round(Math.min(totalW*0.28, 340))
  // Math.min alone has no floor -- on a short viewport (embedded preview panes,
  // small laptop windows) window.innerHeight-380 can collapse toward zero,
  // squeezing all three canvases down to a sliver where curves are effectively
  // invisible even though they're being drawn correctly. Math.max(360, ...) keeps
  // them usable; the notebook scrolls instead of the widget disappearing.
  const H = Math.max(360, Math.round(Math.min(window.innerHeight - 380, 480, 440)))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)

  let data = null       // {spheroidal:{...}, toroidal:{...}, earth_model:{...}, model_name}
  let curType = 'spheroidal'
  let selIdx = -1, hoverIdx = -1
  let eigenHidden = new Set()   // component indices (0=U,1=V,2=P) toggled off via legend click
  let modelHidden = new Set()   // series keys ('rho'/'vp'/'vs') toggled off via legend click
  let eigenLegendRects = [], modelLegendRects = []  // recomputed each draw, for click hit-testing

  const spCvs = par.querySelector('#sbvspectrum'), sctx = spCvs.getContext('2d')
  const eiCvs = par.querySelector('#sbveigen'), ectx = eiCvs.getContext('2d')
  const moCvs = par.querySelector('#sbvmodel'), mctx = moCvs.getContext('2d')
  function hidpi(cv, cx, w, h){ cv.width=Math.round(w*DPR); cv.height=Math.round(h*DPR); cv.style.width=w+'px'; cv.style.height=h+'px'; cx.setTransform(DPR,0,0,DPR,0,0) }
  hidpi(spCvs, sctx, SPW, H); hidpi(eiCvs, ectx, SIDEW, H); hidpi(moCvs, mctx, SIDEW, H)

  const PAD_L = 42, PAD_R = 12, PAD_T = 10, PAD_B = 26

  // Categorical palette cycling by overtone number, readable on black.
  const PALETTE = ['#38bdf8','#f97316','#4ade80','#facc15','#f472b6','#a78bfa','#fb7185','#2dd4bf','#fbbf24','#818cf8']
  function branchColor(n){ return PALETTE[((n%PALETTE.length)+PALETTE.length)%PALETTE.length] }

  function cur(){ return data ? data[curType] : null }

  function drawSpectrum(){
    sctx.clearRect(0,0,SPW,H)
    const d = cur()
    if(!d){
      sctx.fillStyle = '#6b7280'; sctx.font = '13px sans-serif'
      sctx.fillText('pick a body above to load its modes', 12, 20)
      return
    }
    const x0=PAD_L, x1=SPW-PAD_R, y0=PAD_T, y1=H-PAD_B
    const lmax = Math.max(...d.l), fmax = Math.max(...d.f)*1.05
    const X = l => x0 + (l/lmax)*(x1-x0)
    const Y = f => y1 - (f/fmax)*(y1-y0)

    sctx.strokeStyle = '#1f2937'; sctx.fillStyle = '#9ca3af'; sctx.font = '10px sans-serif'; sctx.textAlign = 'right'
    for(let i=0;i<=5;i++){ const f = fmax*i/5, y = Y(f); sctx.beginPath(); sctx.moveTo(x0,y); sctx.lineTo(x1,y); sctx.stroke(); sctx.fillText(f.toFixed(1), x0-6, y+3) }
    sctx.textAlign = 'center'
    for(let i=0;i<=5;i++){ const l = lmax*i/5, x = X(l); sctx.fillText(Math.round(l)+'', x, y1+14) }
    sctx.fillStyle = '#e5e7eb'; sctx.font = '12px sans-serif'
    sctx.fillText((curType==='spheroidal'?'Spheroidal':'Toroidal') + ' spectrum -- ' + (data.model_name||''), (x0+x1)/2, 16)
    sctx.font = '11px sans-serif'; sctx.fillStyle = '#9ca3af'
    sctx.fillText('angular degree l', (x0+x1)/2, H-4)
    sctx.save(); sctx.translate(12, (y0+y1)/2); sctx.rotate(-Math.PI/2); sctx.fillText('frequency (mHz)', 0, 0); sctx.restore()

    for(let i=0;i<d.l.length;i++){
      const isSel = i===selIdx, isHov = i===hoverIdx
      sctx.beginPath()
      sctx.arc(X(d.l[i]), Y(d.f[i]), isSel?5:(isHov?4:2.5), 0, 2*Math.PI)
      sctx.fillStyle = branchColor(d.n[i])
      sctx.globalAlpha = (selIdx===-1 || isSel || isHov) ? 0.9 : 0.35
      sctx.fill()
      sctx.globalAlpha = 1
      if(isSel){ sctx.strokeStyle = '#f5f3ef'; sctx.lineWidth = 1.4; sctx.stroke() }
    }

    if(hoverIdx>=0){
      const label = (curType==='spheroidal'?'S':'T') + ': l='+d.l[hoverIdx]+' n='+d.n[hoverIdx]+' f='+d.f[hoverIdx].toFixed(3)+'mHz'
      sctx.font='12px sans-serif'
      const tw = sctx.measureText(label).width
      const tx = Math.min(X(d.l[hoverIdx])+8, SPW-tw-14), ty = 34
      sctx.fillStyle = 'rgba(11,18,32,0.9)'; sctx.fillRect(tx-5, ty-13, tw+10, 18)
      sctx.strokeStyle = '#374151'; sctx.strokeRect(tx-5, ty-13, tw+10, 18)
      sctx.fillStyle = '#e5e7eb'; sctx.textAlign='left'; sctx.fillText(label, tx, ty)
    }
  }

  function nearestIndex(mx, my){
    const d = cur(); if(!d) return -1
    const x0=PAD_L, x1=SPW-PAD_R, y0=PAD_T, y1=H-PAD_B
    const lmax = Math.max(...d.l), fmax = Math.max(...d.f)*1.05
    const X = l => x0 + (l/lmax)*(x1-x0), Y = f => y1 - (f/fmax)*(y1-y0)
    let best=-1, bestD=10
    for(let i=0;i<d.l.length;i++){
      const dd = Math.hypot(mx-X(d.l[i]), my-Y(d.f[i]))
      if(dd<bestD){ bestD=dd; best=i }
    }
    return best
  }

  function drawEigen(){
    ectx.clearRect(0,0,SIDEW,H)
    const d = cur()
    const x0=PAD_L, x1=SIDEW-PAD_R, y0=PAD_T+14
    if(!d || selIdx<0){
      ectx.fillStyle = '#6b7280'; ectx.font = '12px sans-serif'
      ectx.fillText('click a mode', 10, 18)
      return
    }
    const nr = d.nradial, depth = d.depth
    // depth[] runs center-to-surface (decreasing, depth[last]===0) -- the true max
    // is depth[0], but Math.max(...) doesn't assume an ordering that might change.
    // Using depth[depth.length-1] (===0) here made every Y(dep>0) a division-by-
    // zero (NaN/Infinity), which is why the curve never actually rendered even
    // though the underlying eigenfunction data itself was correct.
    const depthMax = Math.max(...depth)
    // -12 reserves a strip below the depth axis for the amplitude x-axis ticks
    // added below -- eigenfunction scale varies by many orders of magnitude
    // between modes (and between U/V/P within one mode), so it has to be drawn
    // fresh each time rather than assumed from a fixed axis.
    const y1 = H - PAD_B - 12
    // Each component (U/V/P) gets its OWN amplitude scale, always its own true
    // peak -- these can differ by many orders of magnitude within a single mode
    // (surface-wave-like overtones are evanescent toward the interior), so a
    // single shared scale collapses the smaller component(s) to an invisible
    // flat line even though their own internal structure is perfectly real. The
    // "always" (not just over currently-visible components) means toggling one
    // component off via the legend never changes another's own scale either.
    const vmaxOf = comp => {
      let v = 1e-30
      for(let i=0;i<nr;i++) v = Math.max(v, Math.abs(comp[selIdx*nr+i]))
      return v
    }
    const vmaxes = d.ef.map(vmaxOf)
    const Xc = c => v => (x0+x1)/2 + (v/vmaxes[c])*(x1-x0)*0.48
    const Y = dep => y0 + (dep/depthMax)*(y1-y0)

    ectx.strokeStyle = '#1f2937'; ectx.fillStyle = '#9ca3af'; ectx.font='10px sans-serif'; ectx.textAlign='right'
    for(let i=0;i<=5;i++){ const dep=depthMax*i/5, y=Y(dep); ectx.beginPath(); ectx.moveTo(x0,y); ectx.lineTo(x1,y); ectx.stroke(); ectx.fillText(Math.round(dep)+'', x0-6, y+3) }
    ectx.strokeStyle = '#374151'; ectx.beginPath(); ectx.moveTo((x0+x1)/2,y0); ectx.lineTo((x0+x1)/2,y1); ectx.stroke()

    const colors = ['#38bdf8','#f97316','#4ade80']
    d.ef.forEach((comp, c) => {
      if(eigenHidden.has(c)) return
      const X = Xc(c)
      ectx.strokeStyle = colors[c%colors.length]; ectx.lineWidth = 1.8
      ectx.beginPath()
      for(let i=0;i<nr;i++){ const px=X(comp[selIdx*nr+i]), py=Y(depth[i]); i===0?ectx.moveTo(px,py):ectx.lineTo(px,py) }
      ectx.stroke()
    })
    // Each visible component gets its own pair of amplitude ticks (+-peak, drawn
    // in that component's own color) at its own scale -- a single shared numeric
    // axis would be meaningless once the curves are independently normalized.
    ectx.font = '9px sans-serif'; ectx.textAlign = 'center'; ectx.textBaseline = 'top'
    d.ef.forEach((comp, c) => {
      if(eigenHidden.has(c)) return
      const X = Xc(c)
      ectx.strokeStyle = colors[c%colors.length]; ectx.fillStyle = colors[c%colors.length]
      for(const s of [-1,1]){
        const px = X(vmaxes[c]*s)
        ectx.beginPath(); ectx.moveTo(px,y1); ectx.lineTo(px,y1+4); ectx.stroke()
        ectx.fillText(vmaxes[c].toExponential(1), px, y1+5)
      }
    })
    ectx.textBaseline = 'alphabetic'

    ectx.textAlign='center'; ectx.fillStyle='#e5e7eb'; ectx.font='12px sans-serif'
    const nS = (curType==='spheroidal') ? 'S' : 'T'
    ectx.fillText(d.n[selIdx]+nS+d.l[selIdx]+'  f='+d.f[selIdx].toFixed(3)+'mHz', (x0+x1)/2, 16)
    ectx.font='11px sans-serif'
    eigenLegendRects = []
    let lx = x0+4
    d.labels.forEach((lab,c)=>{
      const hidden = eigenHidden.has(c)
      const text = lab + ' (±' + vmaxes[c].toExponential(1) + ')'
      ectx.fillStyle = hidden ? '#4b5563' : colors[c%colors.length]
      ectx.textAlign = 'left'
      ectx.fillText(text, lx, H-6)
      const tw = ectx.measureText(text).width
      eigenLegendRects.push({x:lx-2, y:H-16, w:tw+4, h:14, idx:c})
      if(hidden){ ectx.strokeStyle = '#4b5563'; ectx.beginPath(); ectx.moveTo(lx,H-10); ectx.lineTo(lx+tw,H-10); ectx.stroke() }
      lx += tw + 14
    })
    ectx.save(); ectx.translate(12, (y0+y1)/2); ectx.rotate(-Math.PI/2); ectx.fillStyle='#9ca3af'; ectx.fillText('depth (km)', 0, 0); ectx.restore()
  }

  function drawModel(){
    mctx.clearRect(0,0,SIDEW,H)
    const em = data ? data.earth_model : null
    const x0=PAD_L, x1=SIDEW-PAD_R, y0=PAD_T+14, y1=H-PAD_B
    if(!em){
      mctx.fillStyle = '#6b7280'; mctx.font = '12px sans-serif'
      mctx.fillText('no body loaded yet', 10, 18)
      return
    }
    const depthMax = Math.max(...em.depth)  // see the matching comment in drawEigen
    // -12 reserves a strip below the depth axis for the amplitude x-axis ticks below.
    const yAxisBottom = y1
    const y1m = yAxisBottom - 12
    const Y = dep => y0 + (dep/depthMax)*(y1m-y0)
    const series = [['rho','#f97316','density (kg/m3)'], ['vp','#38bdf8','vp (m/s)'], ['vs','#4ade80','vs (m/s)']]

    mctx.strokeStyle = '#1f2937'; mctx.fillStyle = '#9ca3af'; mctx.font='10px sans-serif'; mctx.textAlign='right'
    for(let i=0;i<=5;i++){ const dep=depthMax*i/5, y=Y(dep); mctx.beginPath(); mctx.moveTo(x0,y); mctx.lineTo(x1,y); mctx.stroke(); mctx.fillText(Math.round(dep)+'', x0-6, y+3) }

    // rho (kg/m3) and vp/vs (m/s) sit in the same numeric range for real planetary
    // bodies, so a single shared 0..vmax scale (over all series, not just the
    // currently-visible ones -- toggling a curve off shouldn't rescale the others)
    // keeps the three curves comparable and gives one coherent x-axis.
    let vmax = 1
    series.forEach(([key])=>{ vmax = Math.max(vmax, ...em[key]) })
    vmax *= 1.05
    const X = v => x0 + (v/vmax)*(x1-x0)

    series.forEach(([key,color])=>{
      if(modelHidden.has(key)) return
      const vals = em[key]
      mctx.strokeStyle = color; mctx.lineWidth = 1.8
      mctx.beginPath()
      for(let i=0;i<vals.length;i++){ const px=X(vals[i]), py=Y(em.depth[i]); i===0?mctx.moveTo(px,py):mctx.lineTo(px,py) }
      mctx.stroke()
    })

    // Amplitude x-axis -- rebuilt every draw from the shared vmax above.
    mctx.strokeStyle = '#4b5563'; mctx.fillStyle = '#9ca3af'; mctx.font = '9px sans-serif'
    mctx.textAlign = 'center'; mctx.textBaseline = 'top'
    for(let i=0;i<=5;i++){
      const v = vmax*i/5, px = X(v)
      mctx.beginPath(); mctx.moveTo(px,y1m); mctx.lineTo(px,y1m+4); mctx.stroke()
      mctx.fillText(Math.round(v)+'', px, y1m+5)
    }
    mctx.textBaseline = 'alphabetic'

    mctx.textAlign='center'; mctx.fillStyle='#e5e7eb'; mctx.font='12px sans-serif'
    mctx.fillText('Body model', (x0+x1)/2, y0-2)

    // Legend lives at the bottom (not top-left) so it never collides with the
    // centered title above -- same placement as the eigenfunction panel's legend.
    mctx.font='11px sans-serif'; mctx.textAlign='left'
    modelLegendRects = []
    let lx = x0+4
    series.forEach(([key,color,label])=>{
      const hidden = modelHidden.has(key)
      mctx.fillStyle = hidden ? '#4b5563' : color
      mctx.fillText(label, lx, H-6)
      const tw = mctx.measureText(label).width
      modelLegendRects.push({x:lx-2, y:H-16, w:tw+4, h:14, key})
      if(hidden){ mctx.strokeStyle = '#4b5563'; mctx.beginPath(); mctx.moveTo(lx,H-10); mctx.lineTo(lx+tw,H-10); mctx.stroke() }
      lx += tw + 14
    })

    mctx.save(); mctx.translate(12, (y0+y1m)/2); mctx.rotate(-Math.PI/2); mctx.fillStyle='#9ca3af'; mctx.font='11px sans-serif'; mctx.fillText('depth (km)', 0, 0); mctx.restore()
  }

  function redraw(){ drawSpectrum(); drawEigen(); drawModel() }
  redraw()

  spCvs.addEventListener('mousemove', e=>{
    const i = nearestIndex(e.offsetX, e.offsetY)
    if(i !== hoverIdx){ hoverIdx = i; drawSpectrum() }
  })
  spCvs.addEventListener('mouseleave', ()=>{ if(hoverIdx!==-1){ hoverIdx=-1; drawSpectrum() } })
  spCvs.addEventListener('click', e=>{
    const i = nearestIndex(e.offsetX, e.offsetY)
    if(i>=0){ selIdx = i; drawSpectrum(); drawEigen() }
  })

  function hitLegend(rects, mx, my){
    for(const r of rects){ if(mx>=r.x && mx<=r.x+r.w && my>=r.y && my<=r.y+r.h) return r }
    return null
  }
  eiCvs.addEventListener('click', e=>{
    const hit = hitLegend(eigenLegendRects, e.offsetX, e.offsetY)
    if(!hit) return
    eigenHidden.has(hit.idx) ? eigenHidden.delete(hit.idx) : eigenHidden.add(hit.idx)
    drawEigen()
  })
  moCvs.addEventListener('click', e=>{
    const hit = hitLegend(modelLegendRects, e.offsetX, e.offsetY)
    if(!hit) return
    modelHidden.has(hit.key) ? modelHidden.delete(hit.key) : modelHidden.add(hit.key)
    drawModel()
  })
  eiCvs.addEventListener('mousemove', e=>{
    eiCvs.style.cursor = hitLegend(eigenLegendRects, e.offsetX, e.offsetY) ? 'pointer' : 'default'
  })
  moCvs.addEventListener('mousemove', e=>{
    moCvs.style.cursor = hitLegend(modelLegendRects, e.offsetX, e.offsetY) ? 'pointer' : 'default'
  })

  // eigenHidden is index-based (0/1/2 = U/V/P for spheroidal, 0 = W for toroidal) --
  // those indices mean different physical quantities across families, so clear it
  // on a family switch rather than carrying, say, "U hidden" over onto W.
  par.querySelector('#sbv-type-s').addEventListener('change', ()=>{ curType='spheroidal'; selIdx=-1; hoverIdx=-1; eigenHidden.clear(); redraw() })
  par.querySelector('#sbv-type-t').addEventListener('change', ()=>{ curType='toroidal'; selIdx=-1; hoverIdx=-1; eigenHidden.clear(); redraw() })

  window.addEventListener('specnm-browse-data', e=>{
    data = JSON.parse(e.detail)
    par.querySelector('#sbv-status').textContent = 'body: ' + (data.model_name||'') + '  (' + data.spheroidal.l.length + ' spheroidal, ' + data.toroidal.l.length + ' toroidal modes)'
    selIdx = -1; hoverIdx = -1
    redraw()
  })
</script>
""")
    end
end

# ╔═╡ bf3349da-ad56-4c4b-9b5c-c2e138abd1c0
SpecnmBrowseView()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.83"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "40c9f1cac973d64f8ca3ef3a09f769ff947e80f3"

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
git-tree-sha1 = "5253f44481f18cd938d4559d5e44fa82198408a6"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.3"

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
# ╠═cd3c5ff7-7992-41f7-9a5e-29934a1355ca
# ╟─bdae265d-3a96-4ecc-a6dd-c166357e801c
# ╟─f2a1c3d4-7e6b-4a5c-9d8e-1b2c3d4e5f6a
# ╟─3508d355-720a-410a-bbd1-50e3779e77ee
# ╟─bf3349da-ad56-4c4b-9b5c-c2e138abd1c0
# ╟─9b0a9a80-e65f-4385-a377-372e408b19ad
# ╠═9f419394-528a-4bde-98a5-d62787c17fa8
# ╠═dc1c3b73-27d9-40cd-a33d-604ac6990daf
# ╟─297386a9-d8e1-4da9-9bed-aca3ef653a11
# ╠═473e3e64-609c-41d2-ac77-b17943011eb2
# ╠═1b958497-3e50-46a1-9237-022ddbe8aabb
# ╟─f018bd73-8f46-448c-9afb-ec47c7ea705b
# ╟─383d55aa-01cb-4068-8140-39f52a5d7884
# ╠═12b77c6d-5a9e-43b7-a34c-8de94eb6003d
# ╟─a9bbb105-3138-481d-9a73-0ec2b8303c36
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
