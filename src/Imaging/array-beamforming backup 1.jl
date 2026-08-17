### A Pluto.jl notebook ###
# v0.20.19

#> [frontmatter]
#> title = "Array Beamforming"
#> tags = ["imaging", "arrays", "plane waves"]
#> layout = "layout.jlhtml"
#> description = "Build a seismic array and see how its geometry focuses a monochromatic plane wave in slowness space."

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

# ╔═╡ c1000001-0001-4001-8001-000000000001
using PlutoUI

# ╔═╡ c1000002-0002-4002-8002-000000000002
PlutoUI.TableOfContents()

# ╔═╡ c1000003-0003-4003-8003-000000000003
md"""
# Array Beamforming

A seismic array does not merely record the same wave many times. Its geometry lets us test many candidate plane-wave slownesses and identify the one that aligns the recordings most coherently. Build an array below by clicking to add receivers; click an existing receiver to remove it.

The synthetic signal is one monochromatic plane wave. The beam image is calculated in Julia from the receiver locations, selected period, velocity, and propagation direction.

##### [Interactive Seismology Notebooks](https://pawbz.github.io/Interactive-Seismology.jl/)

Instructor: *Pawan Bharadwaj*,
Indian Institute of Science, Bengaluru, India
"""

# ╔═╡ c1000004-0004-4004-8004-000000000004
begin
    # The widget definition lives at the end of the Appendix. This explicit
    # dependency keeps cold Pluto loads in the correct order.
    _bf_ready
    PlutoUI.WideCell(@bind _bf ArrayBeamformingInput(); max_width=1400)
end

# ╔═╡ c1000005-0005-4005-8005-000000000005
begin
    bf_safe = _bf isa AbstractDict ? _bf : Dict{String,Any}()
    bf_period = clamp(Float64(get(bf_safe, "period", 1.0)), 0.25, 4.0)
    bf_velocity = clamp(Float64(get(bf_safe, "velocity", 3.5)), 1.5, 6.0)
    bf_azimuth = mod(Float64(get(bf_safe, "azimuth", 40.0)), 360.0)
    bf_receivers = receiver_points(get(bf_safe, "receivers", default_receivers()))
end

# ╔═╡ c1000006-0006-4006-8006-000000000006
let
    payload = beamforming_payload(bf_receivers, bf_period, bf_velocity, bf_azimuth)
    number(x) = isfinite(x) ? string(round(Float64(x), digits = 7)) : "0"
    array(values) = "[" * join(number.(values), ",") * "]"
    message = string(
        "{\"power\":", array(payload.power), ",\"receiver_x\":", array(payload.receiver_x),
        ",\"receiver_y\":", array(payload.receiver_y), ",\"pmax\":", number(payload.pmax),
        ",\"nslowness\":", payload.nslowness, ",\"true_px\":", number(payload.true_px),
        ",\"true_py\":", number(payload.true_py), ",\"direction_x\":", number(payload.direction_x),
        ",\"direction_y\":", number(payload.direction_y), ",\"peak_px\":", number(payload.peak_px),
        ",\"peak_py\":", number(payload.peak_py), ",\"aperture\":", number(payload.aperture),
        ",\"wavelength\":", number(payload.wavelength), ",\"receiver_count\":", payload.receiver_count, "}",
    )
    HTML("""<script>
      window.dispatchEvent(new CustomEvent('beamforming-results', {detail: $(repr(message))}));
    </script>""")
end

# ╔═╡ c1000007-0007-4007-8007-000000000007
md"""
## What the image means

The right panel is a beam power map in horizontal-slowness space. Its coordinates are `p_x` and `p_y` in seconds per kilometre. The cyan marker is the slowness used to create the synthetic wave; the warm peak is what the array recovers.

Try a compact cluster first, then spread receivers around the map. A larger aperture narrows the peak, while an irregular geometry often breaks symmetric sidelobes. Increase the period and observe that the wavelength grows relative to the array, reducing directional resolution.

!!! warning "A beam is not automatically a unique source"
	A monochromatic array response can have sidelobes and ambiguity. Real processing uses bandwidth, time windows, and multiple frequencies to distinguish a physical arrival from an array artifact.
"""

# ╔═╡ c1000008-0008-4008-8008-000000000008
md"""
## The delay-and-sum test

For a plane wave travelling with horizontal slowness vector `\mathbf{p}_0`, the signal at receiver `i` is delayed by `\mathbf{p}_0 \cdot \mathbf{r}_i`. For a selected angular frequency `\omega`, delay-and-sum beamforming tests a candidate `\mathbf{p}` using

```math
B(\mathbf{p}) = \left|\frac{1}{N}\sum_{i=1}^{N}
\exp\left[i\omega(\mathbf{p}-\mathbf{p}_0)\cdot\mathbf{r}_i\right]\right|^2.
```

When the candidate slowness equals the true slowness, every phase is aligned and `B = 1`. Elsewhere, the receiver phases interfere destructively—except where the finite array geometry permits a sidelobe.
"""

# ╔═╡ c1000009-0009-4009-8009-000000000009
md"## Appendix"

# ╔═╡ c1000010-0010-4010-8010-000000000010
md"### Receiver geometry"

# ╔═╡ c1000011-0011-4011-8011-000000000011
"""
	default_receivers()

Return a small, asymmetric starting array in kilometres. The geometry is deliberately
not a perfect circle so students can see how array shape affects sidelobes.
"""
default_receivers() = [
    (-8.0, -3.5), (-5.0, 6.0), (-1.0, -7.0), (0.0, 0.0), (3.5, 7.0),
    (5.5, -4.0), (8.0, 2.0), (1.5, 3.0),
]

# ╔═╡ c1000012-0012-4012-8012-000000000012
"""
	receiver_points(value)

Convert a browser-provided nested array into a concrete vector of `(x, y)` receiver
positions in kilometres. Invalid points are discarded and the map bounds are enforced
at the widget boundary.
"""
function receiver_points(value)
    value isa AbstractVector || return default_receivers()
    points = NTuple{2,Float64}[]
    for point in value
        point isa AbstractVector && length(point) >= 2 || continue
        x = try Float64(point[1]) catch; continue end
        y = try Float64(point[2]) catch; continue end
        isfinite(x) && isfinite(y) || continue
        candidate = (clamp(x, -20.0, 20.0), clamp(y, -20.0, 20.0))
        any(hypot(candidate[1] - q[1], candidate[2] - q[2]) < 0.25 for q in points) && continue
        push!(points, candidate)
        length(points) == 16 && break
    end
    points
end

# ╔═╡ c1000013-0013-4013-8013-000000000013
md"### Monochromatic array response"

# ╔═╡ c1000014-0014-4014-8014-000000000014
"""
	beamforming_payload(receivers, period, velocity, azimuth; nslowness=91)

Evaluate the normalized delay-and-sum power of one plane wave over a horizontal
slowness grid. `azimuth` is clockwise from north and describes the propagation
direction. Returns display-ready arrays; the browser only renders them.
"""
function beamforming_payload(
    receivers::Vector{NTuple{2,Float64}},
    period::Float64,
    velocity::Float64,
    azimuth::Float64;
    nslowness::Int = 91,
)
    θ = deg2rad(azimuth)
    direction_x, direction_y = sin(θ), cos(θ)
    true_px, true_py = direction_x / velocity, direction_y / velocity
    pmax = max(0.55, 1.2 / velocity)
    pgrid = collect(range(-pmax, pmax; length = nslowness))
    power = zeros(Float64, nslowness, nslowness)
    ω = 2π / period
    nreceivers = length(receivers)

    if nreceivers > 0
        for iy in eachindex(pgrid), ix in eachindex(pgrid)
            phase_sum = 0.0 + 0.0im
            px, py = pgrid[ix], pgrid[iy]
            @inbounds for (x, y) in receivers
                phase_sum += cis(ω * ((px - true_px) * x + (py - true_py) * y))
            end
            power[iy, ix] = abs2(phase_sum) / nreceivers^2
        end
    end

    peak = argmax(power)
    peak_py, peak_px = pgrid[peak[1]], pgrid[peak[2]]
    aperture = nreceivers < 2 ? 0.0 : maximum(
        hypot(receivers[i][1] - receivers[j][1], receivers[i][2] - receivers[j][2])
        for i in eachindex(receivers) for j in i + 1:length(receivers)
    )
    return (
        power = [power[iy, ix] for iy in eachindex(pgrid) for ix in eachindex(pgrid)],
        receiver_x = first.(receivers), receiver_y = last.(receivers), pmax = pmax,
        nslowness = nslowness, true_px = true_px, true_py = true_py,
        direction_x = direction_x, direction_y = direction_y,
        peak_px = peak_px, peak_py = peak_py, aperture = aperture,
        wavelength = velocity * period, receiver_count = nreceivers,
    )
end

# ╔═╡ c1000015-0015-4015-8015-000000000015
md"### Verifying the focused beam"

# ╔═╡ c1000016-0016-4016-8016-000000000016
let
    check = beamforming_payload(default_receivers(), 1.0, 3.5, 40.0)
    grid_step = 2check.pmax / (check.nslowness - 1)
    peak_error = hypot(check.peak_px - check.true_px, check.peak_py - check.true_py)
    # The displayed grid need not contain the exact continuous slowness, but
    # its nearest sample must still retain a strongly coherent beam.
    @assert maximum(check.power) > 0.98
    @assert peak_error <= sqrt(2) * grid_step
    (normalised_peak_power = round(maximum(check.power), digits = 6), peak_slowness_error_s_per_km = round(peak_error, digits = 5))
end

# ╔═╡ c1000017-0017-4017-8017-000000000017
md"### The interactive widget"

# ╔═╡ c1000018-0018-4018-8018-000000000018
begin
    struct ArrayBeamformingInput
        period::Float64
        velocity::Float64
        azimuth::Float64
        receivers::Vector{NTuple{2,Float64}}
    end

    ArrayBeamformingInput(; period=1.0, velocity=3.5, azimuth=40.0, receivers=default_receivers()) =
        ArrayBeamformingInput(Float64(period), Float64(velocity), Float64(azimuth), receiver_points(receivers))

    Base.get(w::ArrayBeamformingInput) = Dict{String,Any}(
        "period" => w.period, "velocity" => w.velocity, "azimuth" => w.azimuth,
        "receivers" => [[x, y] for (x, y) in w.receivers],
    )

    function Base.show(io::IO, ::MIME"text/html", w::ArrayBeamformingInput)
        points = "[" * join(("[" * string(x) * "," * string(y) * "]" for (x, y) in w.receivers), ",") * "]"
        write(io, """
        <div id="bf-widget"><style>
        #bf-widget{width:100%;max-width:1400px;margin:auto;color:#e5e7eb;font:14px system-ui,sans-serif}#bf-widget *{box-sizing:border-box}
        #bf-widget .bf-title{padding:10px 14px;margin-bottom:10px;text-align:center;background:#0a0f18;border:1px solid #3b5c85;border-radius:6px}#bf-widget .bf-title-desc{font-size:17px;font-weight:700}#bf-widget .bf-title-hint{margin-top:3px;color:#9ca3af;font-size:13px}
        #bf-widget .bf-workspace{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}#bf-widget .bf-panel,#bf-widget .bf-group{min-width:0;padding:10px;background:#050505;border:1px solid #2f3744;border-radius:6px}#bf-widget .bf-panel-title{font-size:15px;font-weight:700;color:#f3f4f6}#bf-widget .bf-caption{min-height:34px;margin:4px 0 7px;color:#9ca3af;font-size:13px;line-height:1.3}
        #bf-widget canvas{display:block;width:100%;height:310px;background:#000;border:1px solid #374151;border-radius:4px;touch-action:manipulation}#bf-widget .bf-legend{display:flex;gap:12px;flex-wrap:wrap;margin-top:7px;font-size:12px}.bf-key{display:inline-flex;align-items:center;gap:5px}.bf-swatch{width:18px;height:3px;border-radius:2px;background:currentColor}.bf-dot{width:10px;height:10px;border-radius:50%;background:currentColor}
        #bf-widget .bf-controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:8px;margin-top:10px}.bf-group-title{margin-bottom:7px;font-size:16px;font-weight:700;color:#f3f4f6}.bf-row{display:grid;grid-template-columns:minmax(78px,125px) minmax(70px,1fr) minmax(46px,70px);gap:7px;align-items:center;margin:7px 0;color:#d1d5db}#bf-widget input{width:100%;min-width:0;accent-color:#f59e0b}.bf-value{color:#fbbf24;text-align:right;font-variant-numeric:tabular-nums}#bf-widget button{border:1px solid #9ca3af;border-radius:4px;padding:6px 12px;background:#606060;color:#f3f4f6;font-size:14px;cursor:pointer;margin:3px 5px 0 0}#bf-widget button:hover{background:#737373}.bf-status{margin-top:8px;color:#d1d5db;font-size:13px;line-height:1.45}.bf-status b{color:#fbbf24}@media(max-width:760px){#bf-widget .bf-workspace{grid-template-columns:1fr}#bf-widget canvas{height:260px}}
        </style><div class="bf-title"><div class="bf-title-desc">An array focuses a plane wave where trial delays match the receiver delays.</div><div class="bf-title-hint">Click the array map to add receivers &middot; click a receiver to remove it &middot; inspect the focused slowness peak</div></div>
        <div class="bf-workspace"><section class="bf-panel"><div class="bf-panel-title">Receiver array</div><div class="bf-caption">Map extent: ±20 km. The arrow is the synthetic wave-propagation direction.</div><canvas id="bf-array"></canvas><div class="bf-legend"><span class="bf-key" style="color:#fbbf24"><i class="bf-dot"></i>receiver</span><span class="bf-key" style="color:#38bdf8"><i class="bf-swatch"></i>plane-wave direction</span></div></section><section class="bf-panel"><div class="bf-panel-title">Delay-and-sum beam power</div><div class="bf-caption">Warm colours mark coherent trial slownesses. Cyan: true synthetic slowness.</div><canvas id="bf-beam"></canvas><div class="bf-legend"><span class="bf-key" style="color:#22d3ee"><i class="bf-dot"></i>true slowness</span><span class="bf-key" style="color:#fbbf24"><i class="bf-dot"></i>recovered peak</span></div></section></div>
        <div class="bf-controls"><section class="bf-group"><div class="bf-group-title">Incident plane wave</div><label class="bf-row"><span>Period</span><input id="bf-period" type="range" min="0.25" max="4" step="0.05" value="$(w.period)"><span id="bf-period-v" class="bf-value"></span></label><label class="bf-row"><span>Velocity</span><input id="bf-velocity" type="range" min="1.5" max="6" step="0.1" value="$(w.velocity)"><span id="bf-velocity-v" class="bf-value"></span></label><label class="bf-row"><span>Azimuth</span><input id="bf-azimuth" type="range" min="0" max="359" step="1" value="$(w.azimuth)"><span id="bf-azimuth-v" class="bf-value"></span></label></section><section class="bf-group"><div class="bf-group-title">Array editing</div><div class="bf-status">Click open map space to add a receiver. Click a nearby gold receiver to remove it. Up to 16 receivers are allowed.</div><button id="bf-reset" type="button">Reset array</button><button id="bf-clear" type="button">Clear array</button><div id="bf-summary" class="bf-status"></div></section></div></div>
        <script>const currentScript=document.currentScript,par=currentScript.previousElementSibling,by=id=>par.querySelector('#'+id);const state={period:$(w.period),velocity:$(w.velocity),azimuth:$(w.azimuth),receivers:$points};const ids={'bf-period':'period','bf-velocity':'velocity','bf-azimuth':'azimuth'};const LIMIT=20;let result=null;function labels(){by('bf-period-v').textContent=state.period.toFixed(2)+' s';by('bf-velocity-v').textContent=state.velocity.toFixed(1)+' km/s';by('bf-azimuth-v').textContent=state.azimuth.toFixed(0)+'°'}function emit(){par.value={...state,receivers:state.receivers.map(q=>[q[0],q[1]])};par.dispatchEvent(new CustomEvent('input'))}function setup(canvas){const r=canvas.getBoundingClientRect(),d=devicePixelRatio||1,W=Math.max(1,r.width),H=Math.max(1,r.height);canvas.width=Math.round(W*d);canvas.height=Math.round(H*d);const x=canvas.getContext('2d');x.setTransform(d,0,0,d,0,0);x.fillStyle='#000';x.fillRect(0,0,W,H);return[x,W,H]}function mapAxes(x,W,H,title,xlabel,ylabel){const m={l:38,r:14,t:28,b:30},px=v=>m.l+(v+LIMIT)/(2*LIMIT)*(W-m.l-m.r),py=v=>H-m.b-(v+LIMIT)/(2*LIMIT)*(H-m.t-m.b);x.strokeStyle='#1f2937';x.lineWidth=1;for(let i=0;i<5;i++){const q=m.l+i*(W-m.l-m.r)/4,r=m.t+i*(H-m.t-m.b)/4;x.beginPath();x.moveTo(q,m.t);x.lineTo(q,H-m.b);x.moveTo(m.l,r);x.lineTo(W-m.r,r);x.stroke()}x.strokeStyle='#4b5563';x.beginPath();x.moveTo(px(0),m.t);x.lineTo(px(0),H-m.b);x.moveTo(m.l,py(0));x.lineTo(W-m.r,py(0));x.stroke();x.fillStyle='#e5e7eb';x.font='600 13px sans-serif';x.fillText(title,m.l,18);x.fillStyle='#9ca3af';x.font='12px sans-serif';x.fillText(xlabel,W-75,H-8);x.save();x.translate(12,m.t+45);x.rotate(-Math.PI/2);x.fillText(ylabel,0,0);x.restore();return{m,px,py}}function triangle(x,cx,cy,r,fill){x.beginPath();for(let i=0;i<3;i++){const a=Math.PI/2+i*2*Math.PI/3,qx=cx+r*Math.cos(a),qy=cy+r*Math.sin(a);i?x.lineTo(qx,qy):x.moveTo(qx,qy)}x.closePath();x.fillStyle=fill;x.fill();x.strokeStyle='#111827';x.lineWidth=1;x.stroke()}function drawArray(){const [x,W,H]=setup(by('bf-array')),a=mapAxes(x,W,H,'east–west (km)','E','N');if(result){const cx=W/2,cy=H/2,dx=result.direction_x,dy=-result.direction_y;x.strokeStyle='#38bdf8';x.lineWidth=2;x.beginPath();x.moveTo(cx-42*dx,cy-42*dy);x.lineTo(cx+42*dx,cy+42*dy);x.stroke();x.fillStyle='#38bdf8';x.beginPath();x.arc(cx+42*dx,cy+42*dy,4,0,2*Math.PI);x.fill();for(let i=0;i<result.receiver_x.length;i++)triangle(x,a.px(result.receiver_x[i]),a.py(result.receiver_y[i]),7,'#fbbf24')}}function colour(v){const t=Math.max(0,Math.min(1,v));const r=Math.round(20+235*t),g=Math.round(28+160*Math.pow(t,.65)),b=Math.round(58*(1-t));return[r,g,b]}function drawBeam(){if(!result)return;const [x,W,H]=setup(by('bf-beam')),m={l:40,r:14,t:28,b:30},n=result.nslowness,iw=W-m.l-m.r,ih=H-m.t-m.b,off=document.createElement('canvas');off.width=n;off.height=n;const ox=off.getContext('2d'),im=ox.createImageData(n,n);for(let j=0;j<n;j++)for(let i=0;i<n;i++){const c=colour(result.power[j*n+i]),k=4*(j*n+i);im.data[k]=c[0];im.data[k+1]=c[1];im.data[k+2]=c[2];im.data[k+3]=255}ox.putImageData(im,0,0);x.imageSmoothingEnabled=true;x.drawImage(off,m.l,m.t,iw,ih);const pmax=result.pmax,px=v=>m.l+(v+pmax)/(2*pmax)*iw,py=v=>H-m.b-(v+pmax)/(2*pmax)*ih;x.strokeStyle='rgba(229,231,235,.5)';x.lineWidth=1;x.beginPath();x.moveTo(px(0),m.t);x.lineTo(px(0),H-m.b);x.moveTo(m.l,py(0));x.lineTo(W-m.r,py(0));x.stroke();x.fillStyle='#e5e7eb';x.font='600 13px sans-serif';x.fillText('trial slowness p (s/km)',m.l,18);x.fillStyle='#9ca3af';x.font='12px sans-serif';x.fillText('pₓ',W-32,H-8);x.save();x.translate(13,m.t+34);x.rotate(-Math.PI/2);x.fillText('pᵧ',0,0);x.restore();for(const q of [[result.true_px,result.true_py,'#22d3ee'],[result.peak_px,result.peak_py,'#fbbf24']]){x.fillStyle=q[2];x.beginPath();x.arc(px(q[0]),py(q[1]),5,0,2*Math.PI);x.fill();x.strokeStyle='#111827';x.stroke()}}function draw(){drawArray();drawBeam();if(result){const resolution=result.aperture>0?state.period/result.aperture:Infinity;by('bf-summary').innerHTML='<b>'+result.receiver_count+'</b> receivers · aperture <b>'+result.aperture.toFixed(1)+' km</b><br>wavelength <b>'+result.wavelength.toFixed(1)+' km</b> · nominal slowness resolution '+(isFinite(resolution)?'<b>'+resolution.toFixed(3)+' s/km</b>':'needs two receivers');}}function onControl(event){const key=ids[event.target.id];if(!key)return;event.stopImmediatePropagation();state[key]=Number(event.target.value);labels();emit()}function arrayClick(event){const c=by('bf-array'),r=c.getBoundingClientRect(),mx=event.clientX-r.left,my=event.clientY-r.top,W=r.width,H=r.height,m={l:38,r:14,t:28,b:30},x=(mx-m.l)/(W-m.l-m.r)*2*LIMIT-LIMIT,y=(H-m.b-my)/(H-m.t-m.b)*2*LIMIT-LIMIT;let nearest=-1,best=12;for(let i=0;i<state.receivers.length;i++){const q=state.receivers[i],dx=(q[0]-x)/(2*LIMIT)*(W-m.l-m.r),dy=(q[1]-y)/(2*LIMIT)*(H-m.t-m.b);const d=Math.hypot(dx,dy);if(d<best){best=d;nearest=i}}if(nearest>=0)state.receivers.splice(nearest,1);else if(state.receivers.length<16)state.receivers.push([Math.max(-LIMIT,Math.min(LIMIT,x)),Math.max(-LIMIT,Math.min(LIMIT,y))]);emit()}par.addEventListener('input',onControl,true);par.addEventListener('change',onControl,true);by('bf-array').addEventListener('click',arrayClick);by('bf-reset').onclick=()=>{state.receivers=[[-8,-3.5],[-5,6],[-1,-7],[0,0],[3.5,7],[5.5,-4],[8,2],[1.5,3]];emit()};by('bf-clear').onclick=()=>{state.receivers=[];emit()};window.addEventListener('beamforming-results',event=>{result=event.detail?JSON.parse(event.detail):null;draw()});new ResizeObserver(draw).observe(par);labels();draw();</script></div>
        """)
    end

    const _bf_ready = true
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.60"
"""

# ╔═╡ Cell order:
# ╠═c1000001-0001-4001-8001-000000000001
# ╠═c1000002-0002-4002-8002-000000000002
# ╟─c1000003-0003-4003-8003-000000000003
# ╟─c1000004-0004-4004-8004-000000000004
# ╟─c1000005-0005-4005-8005-000000000005
# ╟─c1000006-0006-4006-8006-000000000006
# ╟─c1000007-0007-4007-8007-000000000007
# ╟─c1000008-0008-4008-8008-000000000008
# ╟─c1000009-0009-4009-8009-000000000009
# ╟─c1000010-0010-4010-8010-000000000010
# ╠═c1000011-0011-4011-8011-000000000011
# ╠═c1000012-0012-4012-8012-000000000012
# ╟─c1000013-0013-4013-8013-000000000013
# ╠═c1000014-0014-4014-8014-000000000014
# ╟─c1000015-0015-4015-8015-000000000015
# ╠═c1000016-0016-4016-8016-000000000016
# ╟─c1000017-0017-4017-8017-000000000017
# ╠═c1000018-0018-4018-8018-000000000018
# ╟─00000000-0000-0000-0000-000000000001
