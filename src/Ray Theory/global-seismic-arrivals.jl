### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Global Bodywave Arrivals"
#> date = "2025-08-05"
#> tags = ["raytheory"]
#> description = "What is the arrival you are looking for?"
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

# ╔═╡ 47b2c09a-2ae8-49f0-ba73-ddb6868417b1
begin
    # TauP is installed in PythonCall's managed CPython environment.
    using CondaPkg
    CondaPkg.add("python"; version=">=3.11,<3.14")
    CondaPkg.add_pip("obspy")
    CondaPkg.resolve()
    using PythonCall
    using PlutoUI
end

# ╔═╡ 025b2827-ed43-45f5-a981-56dd599c72cb
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ 6b3bba88-b693-4e39-8866-8166dfc55c30
md"""
# Global Bodywave Arrivals
This notebook interactively visualizes global seismic bodywave arrivals using the TauP toolkit via Python in Julia. Users can select source depth and receiver distance to explore ray paths, phases, and travel times for various seismic waves.


##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)


Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ 8967b290-ec9f-4f8d-bca7-91d2c8c8ff18
begin
    """A draggable Earth cross-section: the source (depth) and receiver (angular
    distance) are set by dragging directly on the circle, not by sliders. The bound
    value only updates on release -- dragging alone never triggers a TauP recompute."""
    struct RayGeometryInput
        distance_deg::Float64
        depth_km::Float64
    end

    RayGeometryInput(; receiver_distance=120.0, source_depth=20.0) =
        RayGeometryInput(Float64(receiver_distance), Float64(source_depth))

    Base.get(w::RayGeometryInput) = Dict{String,Any}(
        "receiver_distance" => w.distance_deg,
        "source_depth" => w.depth_km,
    )

    function Base.show(io::IO, ::MIME"text/html", w::RayGeometryInput)
        write(io, """
<div id="rgwidget" style="display:flex;flex-direction:column;align-items:center;width:100%;color:#9ca3af">
  <style>
    pluto-cell:has(#rgwidget) {
      width: min(80vw, 900px) !important;
      margin-left: calc((100% - min(80vw, 900px)) / 2) !important;
    }
    #rgwidget { width: 100%; box-sizing: border-box; color: #d1d5db; font: 14px sans-serif; }
    #rgwidget .rg-title { width: 100%; box-sizing: border-box; text-align: center; margin-bottom: 10px;
      background: #0a0f18; border: 1px solid #3b5c85; border-radius: 6px; padding: 10px 14px; }
    #rgwidget .rg-title-desc { font-size: 17px; font-weight: 700; color: #e5e7eb; }
    #rgwidget .rg-title-hint { font-size: 13px; color: #9ca3af; margin-top: 3px; }
    #rgwidget .rg-actions { margin-top: 10px; display: flex; justify-content: center;
      align-items: center; gap: 12px; flex-wrap: wrap; }
    #rgwidget button { border-radius: 4px; border: 1px solid #9ca3af; background: #606060; color: #f3f4f6;
      padding: 6px 12px; font-size: 14px; cursor: pointer; }
  </style>
  <div class="rg-title">
    <div class="rg-title-desc">Where the earthquake and receiver sit determines which seismic phases connect them, and how fast each one travels.</div>
    <div class="rg-title-hint">drag the red source or the blue receiver &middot; release to compute the rays &middot; hover a ray to identify it</div>
  </div>
  <canvas id="rgcvs" style="background:#000;border:1px solid #374151;border-radius:6px;display:block"></canvas>
  <div class="rg-actions">
    <button id="rgreset" type="button">Reset defaults</button>
    <span style="font-size:13px;color:#9ca3af">dashed rings are real PREM discontinuity depths</span>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const availW = Math.min(window.innerWidth*0.8, par.clientWidth || window.innerWidth*0.8, 900)
  const heightBudget = Math.max(360, window.innerHeight - 300)
  const SEC = Math.round(Math.min(availW, heightBudget, 640))
  const DPR = Math.min(window.devicePixelRatio || 1, 2)
  const REARTH = 6371
  // Real PREM discontinuity depths (km), from src/assets/data/specnm_models/prem_ani.
  const DISCS = [[24.4,'Moho'],[400,'400'],[670,'670'],[2891,'CMB'],[5149.5,'ICB']]
  let distanceDeg = $(w.distance_deg), depthKm = $(w.depth_km)
  let rayPaths = []      // filled in by the 'raypath-results' push from Julia, below
  let hoverIdx = -1, hoverPos = null   // which rayPaths[] entry the cursor is over

  const cvs = par.querySelector('#rgcvs'), ctx = cvs.getContext('2d')
  function hidpi(cv, cx, w, h){
    cv.width = Math.round(w*DPR); cv.height = Math.round(h*DPR)
    cv.style.width = w+'px'; cv.style.height = h+'px'
    cx.setTransform(DPR,0,0,DPR,0,0)
  }
  hidpi(cvs, ctx, SEC, SEC)
  const CX = SEC/2, CY = SEC/2
  // Smaller than the canvas half-width so the epicentral-distance arc and the
  // source/receiver/distance labels all have room to sit outside the Earth circle
  // without being clipped by the canvas edge.
  const R = SEC*0.5*0.68

  function polarXY(thetaDeg, canvasR){
    const th = thetaDeg*Math.PI/180
    return [CX + canvasR*Math.sin(th), CY - canvasR*Math.cos(th)]
  }
  function toXY(thetaDeg, rKm){
    return polarXY(thetaDeg, (rKm/REARTH)*R)
  }
  function sourcePt(){ return toXY(0, REARTH-depthKm) }
  function receiverPt(){ return toXY(distanceDeg, REARTH) }

  // Source/receiver markers matching plot_rays' own convention: a yellow star at the
  // source, a rust-colored marker at the receiver pointing inward along the local
  // radial direction (obspy draws an arrow into the surface; a triangle is the
  // canvas-friendly equivalent, oriented per receiver angle rather than fixed "down"
  // since the receiver can sit anywhere around the circle).
  function drawStar(cx, cy, r){
    ctx.beginPath()
    for(let i=0;i<10;i++){
      const ang = -Math.PI/2 + i*Math.PI/5
      const rr = i%2===0 ? r : r*0.42
      const x = cx+rr*Math.cos(ang), y = cy+rr*Math.sin(ang)
      i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y)
    }
    ctx.closePath()
    ctx.fillStyle = '#FEF215'; ctx.fill()
    ctx.strokeStyle = '#4b5563'; ctx.lineWidth = 1.4; ctx.stroke()
  }
  function drawReceiverMarker(px, py, angDeg, r){
    const th = angDeg*Math.PI/180
    const outX = Math.sin(th), outY = -Math.cos(th)
    const tanX = Math.cos(th), tanY = Math.sin(th)
    const baseX = px + outX*r*1.6, baseY = py + outY*r*1.6
    ctx.beginPath()
    ctx.moveTo(px, py)
    ctx.lineTo(baseX + tanX*r*0.9, baseY + tanY*r*0.9)
    ctx.lineTo(baseX - tanX*r*0.9, baseY - tanY*r*0.9)
    ctx.closePath()
    ctx.fillStyle = '#C95241'; ctx.fill()
    ctx.strokeStyle = '#4b5563'; ctx.lineWidth = 1.4; ctx.stroke()
  }

  function drawArrowHead(px, py, dirX, dirY, size){
    const ang = Math.atan2(dirY, dirX)
    ctx.beginPath()
    ctx.moveTo(px, py)
    ctx.lineTo(px - size*Math.cos(ang-0.4), py - size*Math.sin(ang-0.4))
    ctx.lineTo(px - size*Math.cos(ang+0.4), py - size*Math.sin(ang+0.4))
    ctx.closePath()
    ctx.fillStyle = '#9ca3af'; ctx.fill()
  }

  // Epicentral distance measured the way it's actually defined: along the surface,
  // from the source's surface projection (theta=0) to the receiver -- not the
  // straight chord between them. Drawn just outside the Earth circle as a
  // double-headed dimension arc, like a technical drawing's angle callout.
  function drawEpicentralArc(){
    const arcR = R + 16
    const steps = Math.max(8, Math.round(Math.abs(distanceDeg)/3))
    ctx.beginPath()
    for(let i=0;i<=steps;i++){
      const th = distanceDeg*i/steps
      const p = polarXY(th, arcR)
      i===0 ? ctx.moveTo(p[0],p[1]) : ctx.lineTo(p[0],p[1])
    }
    ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 1.4; ctx.stroke()

    const p0 = polarXY(0, arcR), p1 = polarXY(distanceDeg, arcR)
    const th1 = distanceDeg*Math.PI/180
    drawArrowHead(p0[0], p0[1], -1, 0, 7)
    drawArrowHead(p1[0], p1[1], Math.cos(th1), Math.sin(th1), 7)
  }

  // A single legend box (bottom-left) explains all three markers and carries their
  // current values -- keeps the circle itself free of scattered text labels that
  // would otherwise compete with the ray paths for space.
  function drawLegend(){
    const rows = [
      {icon:'source', text: 'Source · depth ' + Math.round(depthKm) + ' km'},
      {icon:'receiver', text: 'Receiver'},
      {icon:'arc', text: 'Epicentral distance ' + Math.round(distanceDeg) + '°'},
      {icon:'p', text: 'P-type leg (solid)'},
      {icon:'s', text: 'S-type leg (wiggly)'},
    ]
    ctx.font = '12px sans-serif'
    const pad = 10, rowH = 20, iconW = 24
    let textW = 0
    for(const r of rows) textW = Math.max(textW, ctx.measureText(r.text).width)
    const boxW = pad*2 + iconW + textW, boxH = pad*2 + rows.length*rowH - (rowH-14)
    const bx = 10, by = SEC - boxH - 10
    ctx.fillStyle = 'rgba(11,18,32,0.9)'; ctx.fillRect(bx, by, boxW, boxH)
    ctx.strokeStyle = '#374151'; ctx.lineWidth = 1; ctx.strokeRect(bx, by, boxW, boxH)

    rows.forEach((r, i) => {
      const cy = by + pad + i*rowH + 7
      const cx = bx + pad + iconW/2
      if(r.icon === 'source'){
        drawStar(cx, cy, 6)
      } else if(r.icon === 'receiver'){
        ctx.beginPath()
        ctx.moveTo(cx, cy+6); ctx.lineTo(cx-6, cy-5); ctx.lineTo(cx+6, cy-5); ctx.closePath()
        ctx.fillStyle = '#C95241'; ctx.fill(); ctx.strokeStyle = '#4b5563'; ctx.lineWidth = 1; ctx.stroke()
      } else if(r.icon === 'arc'){
        ctx.beginPath(); ctx.arc(cx, cy+11, 12, Math.PI*1.15, Math.PI*1.85)
        ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 1.4; ctx.stroke()
      } else if(r.icon === 'p'){
        ctx.beginPath(); ctx.moveTo(cx-8, cy+3); ctx.lineTo(cx+8, cy+3)
        ctx.strokeStyle = colorFor('p'); ctx.lineWidth = 2; ctx.stroke()
      } else if(r.icon === 's'){
        ctx.beginPath()
        for(let k=-8;k<=8;k++) k===-8 ? ctx.moveTo(cx+k, cy+3+2.5*Math.sin(k*1.1)) : ctx.lineTo(cx+k, cy+3+2.5*Math.sin(k*1.1))
        ctx.strokeStyle = colorFor('s'); ctx.lineWidth = 1.5; ctx.stroke()
      }
      ctx.fillStyle = '#e5e7eb'
      ctx.fillText(r.text, bx + pad + iconW, cy+4)
    })
  }

  // Wave-family colors: warm for P-type legs, cool for S-type -- matches the
  // wiggly/solid line-style distinction below so the two encode the same thing twice.
  function colorFor(wave){
    if(wave === 's') return '#38bdf8'
    if(wave === 'diff') return '#fbbf24'
    return '#f97316'
  }

  function projectSegment(seg){
    const pts = []
    for(let i=0;i<seg.dist.length;i++) pts.push(toXY(seg.dist[i]*180/Math.PI, REARTH-seg.depth[i]))
    return pts
  }

  // Upsample a polyline by linear interpolation so the wiggle (below) has enough
  // points to look smooth even on a short segment -- same reason obspy's own
  // plot_rays interpolates before applying its sketch effect on s-wave legs.
  function densify(pts, minPoints){
    if(pts.length >= minPoints || pts.length < 2) return pts
    const steps = Math.ceil((minPoints-1)/(pts.length-1))
    const out = []
    for(let i=0;i<pts.length-1;i++){
      for(let k=0;k<steps;k++){
        const t = k/steps
        out.push([pts[i][0]+(pts[i+1][0]-pts[i][0])*t, pts[i][1]+(pts[i+1][1]-pts[i][1])*t])
      }
    }
    out.push(pts[pts.length-1])
    return out
  }

  // A clean sinusoidal perpendicular-offset wiggle -- the canvas analogue of
  // matplotlib's path.sketch effect obspy uses to mark s-wave legs.
  function drawWiggly(pts, amp, wavelength){
    if(pts.length < 2) return
    let acc = 0
    ctx.beginPath()
    for(let i=0;i<pts.length;i++){
      if(i>0) acc += Math.hypot(pts[i][0]-pts[i-1][0], pts[i][1]-pts[i-1][1])
      const p0 = pts[Math.max(0,i-1)], p1 = pts[Math.min(pts.length-1,i+1)]
      let tx = p1[0]-p0[0], ty = p1[1]-p0[1]
      const tl = Math.hypot(tx,ty) || 1
      tx /= tl; ty /= tl
      const nx = -ty, ny = tx
      const s = amp*Math.sin(2*Math.PI*acc/wavelength)
      const x = pts[i][0]+nx*s, y = pts[i][1]+ny*s
      i===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y)
    }
    ctx.stroke()
  }

  function drawSegment(seg, alpha, isHover){
    ctx.globalAlpha = alpha
    ctx.strokeStyle = colorFor(seg.wave)
    if(seg.wave === 's'){
      ctx.lineWidth = isHover ? 2.4 : 1.5
      drawWiggly(densify(projectSegment(seg), 60), 3, 14)
    } else {
      ctx.lineWidth = isHover ? 3 : 2
      const pts = projectSegment(seg)
      ctx.beginPath()
      pts.forEach((xy,i)=> i===0?ctx.moveTo(xy[0],xy[1]):ctx.lineTo(xy[0],xy[1]))
      ctx.stroke()
    }
    ctx.globalAlpha = 1
  }

  function distToSegment(px,py, x1,y1, x2,y2){
    const dx=x2-x1, dy=y2-y1
    const len2 = dx*dx+dy*dy
    let t = len2 ? ((px-x1)*dx+(py-y1)*dy)/len2 : 0
    t = Math.max(0, Math.min(1, t))
    return Math.hypot(px-(x1+t*dx), py-(y1+t*dy))
  }

  function nearestPathIndex(mx,my){
    let best = -1, bestD = 8
    rayPaths.forEach((p, idx) => {
      for(const seg of p.segments){
        const pts = projectSegment(seg)
        for(let i=0;i<pts.length-1;i++){
          const d = distToSegment(mx,my, pts[i][0],pts[i][1], pts[i+1][0],pts[i+1][1])
          if(d < bestD){ bestD = d; best = idx }
        }
      }
    })
    return best
  }

  function redraw(){
    ctx.clearRect(0,0,SEC,SEC)
    ctx.beginPath(); ctx.arc(CX,CY,R,0,2*Math.PI)
    ctx.fillStyle = '#0b1220'; ctx.fill()
    ctx.strokeStyle = '#374151'; ctx.lineWidth = 1.4; ctx.stroke()

    // Angled off to the left of straight-up so labels don't sit inside the dense
    // near-vertical bundle of rays leaving the source.
    ctx.font = '12px sans-serif'
    ctx.textAlign = 'right'
    for(const [d,label] of DISCS){
      const rf = ((REARTH-d)/REARTH)*R
      ctx.beginPath(); ctx.setLineDash([3,4])
      ctx.arc(CX,CY,rf,0,2*Math.PI); ctx.strokeStyle = '#2f3744'; ctx.lineWidth = 1; ctx.stroke()
      ctx.setLineDash([])
      const lp = toXY(-24, REARTH-d)
      ctx.fillStyle = '#6b7280'; ctx.fillText(label, lp[0]-4, lp[1]+3)
    }
    ctx.textAlign = 'left'

    // Every ray drawn faded at once (so the whole family of arrivals is visible),
    // then the hovered one redrawn last, full-strength, on top of the others.
    // Fade the rest out further while something is hovered, for more contrast.
    const fadeAlpha = hoverIdx >= 0 ? 0.12 : 0.28
    rayPaths.forEach((p, idx) => {
      if(idx === hoverIdx) return
      for(const seg of p.segments) drawSegment(seg, fadeAlpha, false)
    })
    if(hoverIdx >= 0 && rayPaths[hoverIdx]){
      for(const seg of rayPaths[hoverIdx].segments) drawSegment(seg, 1.0, true)
    }

    drawEpicentralArc()

    const sp = sourcePt(), rp = receiverPt()
    drawStar(sp[0], sp[1], 9)
    drawReceiverMarker(rp[0], rp[1], distanceDeg, 8)

    drawLegend()

    // A stable top-right readout mirrors the floating cursor tooltip below, so the
    // hovered phase's name/time can be read without following the cursor.
    ctx.font = '12px sans-serif'
    ctx.textAlign = 'right'
    if(hoverIdx >= 0 && rayPaths[hoverIdx]){
      ctx.fillStyle = '#e5e7eb'
      ctx.fillText(rayPaths[hoverIdx].name + '   ' + rayPaths[hoverIdx].time.toFixed(1) + ' s', SEC-10, 16)
    } else if(rayPaths.length){
      ctx.fillStyle = '#6b7280'
      ctx.fillText(rayPaths.length + ' phases — hover a ray to identify it', SEC-10, 16)
    }
    ctx.textAlign = 'left'

    if(hoverIdx >= 0 && rayPaths[hoverIdx] && hoverPos){
      const p = rayPaths[hoverIdx]
      const label = p.name + '   ' + p.time.toFixed(1) + ' s'
      ctx.font = '13px sans-serif'
      const tw = ctx.measureText(label).width
      const tx = Math.min(hoverPos[0]+12, SEC-tw-16), ty = Math.max(hoverPos[1]-12, 16)
      ctx.fillStyle = 'rgba(11,18,32,0.9)'; ctx.fillRect(tx-6, ty-14, tw+12, 20)
      ctx.strokeStyle = '#374151'; ctx.lineWidth = 1; ctx.strokeRect(tx-6, ty-14, tw+12, 20)
      ctx.fillStyle = '#e5e7eb'; ctx.fillText(label, tx, ty)
    }
  }

  function emit(){
    par.value = {receiver_distance: distanceDeg, source_depth: depthKm}
    par.dispatchEvent(new CustomEvent('input'))
  }

  function hitTest(mx, my){
    const sp = sourcePt(), rp = receiverPt()
    const ds = Math.hypot(mx-sp[0], my-sp[1])
    const dr = Math.hypot(mx-rp[0], my-rp[1])
    if(ds < 12 && ds <= dr) return 'source'
    if(dr < 12) return 'receiver'
    return null
  }

  let dragging = null
  cvs.addEventListener('mousedown', e=>{ dragging = hitTest(e.offsetX, e.offsetY) })
  cvs.addEventListener('mousemove', e=>{
    if(dragging === 'source'){
      let rf = CY - e.offsetY
      rf = Math.max(((REARTH-700)/REARTH)*R, Math.min(R, rf))
      depthKm = REARTH - (rf/R)*REARTH
      redraw()
    } else if(dragging === 'receiver'){
      let ang = Math.atan2(e.offsetX-CX, -(e.offsetY-CY)) * 180/Math.PI
      distanceDeg = Math.max(0, Math.min(180, ang))
      redraw()
    } else {
      const h = hitTest(e.offsetX, e.offsetY)
      cvs.style.cursor = h ? 'grab' : 'default'
      hoverPos = [e.offsetX, e.offsetY]
      hoverIdx = h ? -1 : nearestPathIndex(e.offsetX, e.offsetY)
      redraw()
    }
  })
  cvs.addEventListener('mouseleave', ()=>{
    if(hoverIdx !== -1){ hoverIdx = -1; redraw() }
  })
  // Only publish the bound value on release -- dragging is purely local/visual, so
  // TauP (a Python round-trip) recomputes once per gesture, not once per pixel.
  window.addEventListener('mouseup', ()=>{
    if(dragging) emit()
    dragging = null
  })

  par.querySelector('#rgreset').addEventListener('click', ()=>{
    distanceDeg = 120; depthKm = 20; redraw(); emit()
  })

  window.addEventListener('raypath-results', e=>{
    const d = e.detail ? JSON.parse(e.detail) : null
    if(!d) return
    rayPaths = d.paths || []
    hoverIdx = -1
    redraw()
  })

  redraw(); emit()
</script>
""")
    end

    const _rg_ready = true
end

# ╔═╡ 0426c6fd-4bb8-413f-b552-0112434d907c
begin
    _rg_ready
    @bind arrival_controls RayGeometryInput()
end

# ╔═╡ e0c1ab0d-32c5-47b4-9d89-e2e51c2fe0dd
begin
    geometry = (
        receiver_distance=Float64(arrival_controls["receiver_distance"]),
        source_depth=Float64(arrival_controls["source_depth"]),
    )
end

# ╔═╡ 9b0b4e7e-fd8f-4573-af80-ed76ff2848f5
md"## Appendix"

# ╔═╡ dd4cb9d8-8d6e-4ea8-b6bf-545631fecff8
taup = pyimport("obspy.taup")

# ╔═╡ 1e7a3c9a-6c2f-4b3a-9c5f-2a6f7e8b9d10
# `split_ray_path` classifies each leg of a ray as p/s/diff, exactly what
# plot_rays(indicate_wave_type=true) uses to decide which legs to draw wiggly.
taup_utils = pyimport("obspy.taup.utils")

# ╔═╡ 23f2f44c-144a-4fd1-a429-656bf0af4cca
model = taup.TauPyModel(model="iasp91")

# ╔═╡ 7568eb42-fe1b-44bc-86a1-dda9200bb49b
arrivals = model.get_ray_paths(geometry.source_depth, geometry.receiver_distance)

# ╔═╡ 460863de-314c-4196-b72b-eb6382cce6f7
# Push every computed ray path (not just a filtered subset) straight into the
# RayGeometryInput widget above -- it stays mounted across reruns of this cell, same
# CustomEvent pattern geoid-kernel uses to push its geoid/topography maps back to the
# already-rendered globe. The widget draws them all faded and highlights on hover, so
# there is no separate phase-selection step: showing the whole family of arrivals is
# the point.
let
    num(x) = isfinite(x) ? string(round(x, sigdigits=6)) : "0"
    jsonarr(v) = "[" * join(num.(v), ",") * "]"
    entries = String[]
    for arrival in arrivals
        name = pyconvert(String, arrival.name)
        t = pyconvert(Float64, arrival.time)
        paths_py, waves_py = taup_utils.split_ray_path(arrival.path, model.model)
        n_seg = pyconvert(Int, paths_py.__len__())
        seg_entries = String[]
        for si in 0:(n_seg-1)
            seg = paths_py[si]
            wave = pyconvert(String, waves_py[si])
            dist = pyconvert(Vector{Float64}, seg["dist"])
            depth = pyconvert(Vector{Float64}, seg["depth"])
            push!(seg_entries, string(
                "{\"wave\":\"", wave, "\",\"dist\":", jsonarr(dist), ",\"depth\":", jsonarr(depth), "}",
            ))
        end
        push!(entries, string(
            "{\"name\":\"", name, "\",\"time\":", num(t), ",\"segments\":[", join(seg_entries, ","), "]}",
        ))
    end
    payload = "{\"paths\":[" * join(entries, ",") * "]}"
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('raypath-results', {detail: $(repr(payload))}));
    </script>""")
end

# ╔═╡ 5d94e67e-5334-4e1c-9838-749b2318c66d
md"""
## Credits
- [https://www.seis.sc.edu/taup/](https://www.seis.sc.edu/taup/)
- [https://docs.obspy.org/packages/obspy.taup.html](https://docs.obspy.org/packages/obspy.taup.html)
"""




# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"

[compat]
CondaPkg = "~0.2.36"
PlutoUI = "~0.7.83"
PythonCall = "~0.9.35"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "7131236ddec56d88ba0b0a5a94dd48f90fe89dfc"

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

[[deps.CondaPkg]]
deps = ["JSON", "Markdown", "MicroMamba", "Pidfile", "Pkg", "Preferences", "Scratch", "TOML", "pixi_jll"]
git-tree-sha1 = "2b1afb8ae65a0758795b00adafb37f97e67ef0e9"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.36"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c89d196f5ffb64bfbf80985b699ea913b0d2c211"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.6.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

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

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MicroMamba]]
deps = ["Pkg", "Scratch", "micromamba_jll"]
git-tree-sha1 = "535656ce55266bfed0575cd051acc4f36dc869a0"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.15"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "05f45c2e0de6259db764adbfd2f1dc6d3f8de13c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "32a4e09c5f29402573d673901778a0e03b0807b9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.6"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Preferences", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "2b67e030054dd9438a00e3d7f59927e839b00569"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.35"

    [deps.PythonCall.extensions]
    CategoricalArraysExt = "CategoricalArrays"
    PyCallExt = "PyCall"

    [deps.PythonCall.weakdeps]
    CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"

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

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

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

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "82bee338d650aa515f31866c460cb7e3bcef90b8"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

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
git-tree-sha1 = "0f38a06c83f0007bbab3cf911262841c9a0f07e0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.13.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "3b0738bd7c5645641845da25cbd99800b8718689"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnsafePointers]]
git-tree-sha1 = "c81331b3b2e60a982be57c046ec91f599ede674a"
uuid = "e17b2a0c-0bdf-430a-bd0c-3a23cae4ff39"
version = "1.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.micromamba_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "717df6f6892af4ee13279a73aa58474e58a88667"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "2.3.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.pixi_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "3667b0931a7fe50f0a5554c61af00e5640019e21"
uuid = "4d7b5844-a134-5dcd-ac86-c8f19cd51bed"
version = "0.63.2+0"
"""

# ╔═╡ Cell order:
# ╠═025b2827-ed43-45f5-a981-56dd599c72cb
# ╟─6b3bba88-b693-4e39-8866-8166dfc55c30
# ╟─0426c6fd-4bb8-413f-b552-0112434d907c
# ╠═e0c1ab0d-32c5-47b4-9d89-e2e51c2fe0dd
# ╠═8967b290-ec9f-4f8d-bca7-91d2c8c8ff18
# ╟─9b0b4e7e-fd8f-4573-af80-ed76ff2848f5
# ╠═47b2c09a-2ae8-49f0-ba73-ddb6868417b1
# ╠═dd4cb9d8-8d6e-4ea8-b6bf-545631fecff8
# ╠═1e7a3c9a-6c2f-4b3a-9c5f-2a6f7e8b9d10
# ╠═23f2f44c-144a-4fd1-a429-656bf0af4cca
# ╠═7568eb42-fe1b-44bc-86a1-dda9200bb49b
# ╠═460863de-314c-4196-b72b-eb6382cce6f7
# ╟─5d94e67e-5334-4e1c-9838-749b2318c66d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
