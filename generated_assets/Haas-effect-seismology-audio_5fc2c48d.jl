### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Why Seismology? — Hearing Echoes"
#> tags = ["introduction", "waves", "seismology"]
#> layout = "layout.jlhtml"
#> description = "Hear how delayed copies turn a simple sound into a seismic wavefield."

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

# ╔═╡ 9e9c44d5-7940-4cc0-9f68-f79a315c5a27
begin
    using PlutoUI
end

# ╔═╡ 5c5c34cc-b71d-4bdc-a857-8763ae803e8c
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ a63f99b1-a70d-4a17-bcdb-91f0ef0b3963
md"""
# There is no Haas effect in seismology. Our data are complex.

A seismogram does not report one tidy pulse. It records the direct wave **and** every
reflection and multiple that reaches the instrument. Your ears receive exactly the same
kind of superposition in a room — but your brain normally makes it seem simple.
"""

# ╔═╡ 6a4cd45d-4e9e-4fd7-8856-22d77bcac1b7
md"""
## Make one sound become many

Use the sentence:

> “The Earth has many seismic reflections.”

You can record it directly in the notebook, or load any short audio file. Then play it
with a delayed copy. Start with the built-in voiced test signal if you want to explore
the timing before recording.

!!! note "Microphone permission"
    The notebook asks for microphone access only after you choose **Record sentence**.
    Nothing is uploaded; the clip stays in the browser tab.
"""

# ╔═╡ 7f0c9d2e-3b4a-4a4b-8e7f-1c2d3e4f5a6b
HTML("""
<figure id="re-room-earth-sketch" style="margin:18px 0;color:#d1d5db">
  <style>
    #re-room-earth-sketch{background:#050505;border:1px solid #2f3744;border-radius:7px;padding:12px;box-sizing:border-box}
    #re-room-earth-sketch .re-sketch-title{color:#f3f4f6;text-align:center;font:700 18px sans-serif;margin:0 0 8px}
    #re-room-earth-sketch svg{display:block;width:100%;height:auto;background:#000;border:1px solid #374151;border-radius:5px}
    #re-room-earth-sketch figcaption{color:#9ca3af;font:13px sans-serif;text-align:center;margin-top:8px}
  </style>
  <div class="re-sketch-title">The same wavefield idea at two scales</div>
  <svg viewBox="0 0 1200 410" role="img" aria-label="Room and Earth conceptual sketch showing direct and reflected waves">
    <defs>
      <marker id="re-blue-arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#38bdf8"/></marker>
      <marker id="re-red-arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#ef4444"/></marker>
    </defs>
    <rect x="18" y="20" width="554" height="350" rx="9" fill="#050505" stroke="#2f3744" stroke-width="2"/>
    <rect x="628" y="20" width="554" height="350" rx="9" fill="#050505" stroke="#2f3744" stroke-width="2"/>
    <text x="295" y="58" text-anchor="middle" fill="#f3f4f6" font-family="sans-serif" font-size="25" font-weight="700">Room</text>
    <text x="905" y="58" text-anchor="middle" fill="#f3f4f6" font-family="sans-serif" font-size="25" font-weight="700">Earth</text>

    <line x1="72" y1="112" x2="518" y2="112" stroke="#9ca3af" stroke-width="5"/>
    <line x1="72" y1="305" x2="518" y2="305" stroke="#9ca3af" stroke-width="5"/>
    <text x="75" y="96" fill="#9ca3af" font-family="sans-serif" font-size="16">wall</text>
    <text x="75" y="334" fill="#9ca3af" font-family="sans-serif" font-size="16">wall</text>
    <circle cx="125" cy="208" r="19" fill="#e5e7eb"/><path d="M145 190 q26 18 0 36" fill="none" stroke="#e5e7eb" stroke-width="4"/>
    <text x="125" y="251" text-anchor="middle" fill="#e5e7eb" font-family="sans-serif" font-size="16">speaker</text>
    <rect x="475" y="183" width="17" height="50" rx="8" fill="#e5e7eb"/><circle cx="483" cy="174" r="14" fill="#e5e7eb"/>
    <text x="483" y="263" text-anchor="middle" fill="#e5e7eb" font-family="sans-serif" font-size="16">microphone</text>
    <path d="M155 208 L464 208" fill="none" stroke="#38bdf8" stroke-width="5" marker-end="url(#re-blue-arrow)"/>
    <path d="M151 199 L285 122 L465 194" fill="none" stroke="#ef4444" stroke-width="4" marker-end="url(#re-red-arrow)"/>
    <path d="M151 217 L300 295 L465 222" fill="none" stroke="#ef4444" stroke-width="4" marker-end="url(#re-red-arrow)"/>
    <text x="310" y="194" text-anchor="middle" fill="#38bdf8" font-family="sans-serif" font-size="17">direct sound</text>
    <text x="325" y="142" text-anchor="middle" fill="#ef4444" font-family="sans-serif" font-size="17">wall reflection</text>

    <circle cx="875" cy="220" r="130" fill="#111827" stroke="#e5e7eb" stroke-width="5"/>
    <circle cx="875" cy="220" r="121" fill="#172554" stroke="#38bdf8" stroke-width="3"/>
    <circle cx="875" cy="220" r="80" fill="#3f3f46" stroke="#9ca3af" stroke-width="3"/>
    <circle cx="875" cy="220" r="46" fill="#7f1d1d" stroke="#ef4444" stroke-width="3"/>
    <text x="775" y="275" text-anchor="middle" fill="#f3f4f6" font-family="sans-serif" font-size="34">★</text>
    <text x="742" y="334" text-anchor="middle" fill="#e5e7eb" font-family="sans-serif" font-size="16">earthquake</text>
    <path d="M958 91 L992 91 L975 122 z" fill="#f3f4f6"/>
    <text x="1072" y="104" text-anchor="middle" fill="#e5e7eb" font-family="sans-serif" font-size="16">seismometer</text>
    <path d="M789 265 L960 122" fill="none" stroke="#38bdf8" stroke-width="5" marker-end="url(#re-blue-arrow)"/>
    <path d="M789 267 Q850 314 908 263 Q944 205 960 122" fill="none" stroke="#ef4444" stroke-width="4" marker-end="url(#re-red-arrow)"/>
    <path d="M789 259 Q870 102 958 118" fill="none" stroke="#ef4444" stroke-width="4" marker-end="url(#re-red-arrow)"/>
    <text x="1040" y="174" fill="#e5e7eb" font-family="sans-serif" font-size="16">crust + surface</text>
    <line x1="1027" y1="170" x2="990" y2="170" stroke="#38bdf8" stroke-width="2"/>
    <text x="1040" y="212" fill="#e5e7eb" font-family="sans-serif" font-size="16">mantle</text>
    <line x1="1027" y1="208" x2="947" y2="208" stroke="#9ca3af" stroke-width="2"/>
    <text x="1040" y="250" fill="#e5e7eb" font-family="sans-serif" font-size="16">core</text>
    <line x1="1027" y1="246" x2="918" y2="246" stroke="#ef4444" stroke-width="2"/>
    <text x="1040" y="288" fill="#fca5a5" font-family="sans-serif" font-size="16">CMB</text>
    <line x1="1027" y1="284" x2="946" y2="284" stroke="#ef4444" stroke-width="2" stroke-dasharray="5 4"/>
  </svg>
  <figcaption><span style="color:#38bdf8">Blue</span> = direct arrival &nbsp;•&nbsp; <span style="color:#ef4444">red</span> = a reflected path. Treat the surface and CMB as the room's two giant walls.</figcaption>
</figure>
""")

# ╔═╡ 28b6cf3a-61c9-40d4-9d0f-4f290d21874b
md"""
## One idea, two laboratories

| A room: hearing echoes | Earth: interpreting seismograms |
|:---|:---|
| A **speaker** launches a sound | An **earthquake** launches seismic energy |
| **Acoustic** waves travel through air | **Elastic** waves travel through rock |
| Two walls make paths and echoes | The **surface** and **CMB** act as two great reflecting boundaries |
| A microphone records air-pressure variations | A seismometer records ground motion |
| Direct sound, wall reflections, multiple echoes | Direct P, PP/reflections, and multiples |
| The auditory system fuses or suppresses reflections | A seismologist identifies and interprets phases |

!!! danger "The incoming wavefield is not the perception"
    A microphone and a seismometer report what arrives. Ears and seismologists must infer
    what produced those arrivals.
"""

# ╔═╡ 47ee7135-61df-4c58-bcfc-07d23a65d6e6
md"""
## The question to ask

> **“If our ears behaved like a seismometer, every room would sound like a cave.”**

Why do rooms usually sound simple instead? The brain must infer a source and suppress or
reinterpret reflections — the same broad task a seismologist faces when separating direct
arrivals, reflected phases, and multiples in a seismogram.
"""

# ╔═╡ 1de07b55-f4c8-4b6b-b090-f20eced1bea0
md"""
## Appendix

### The interactive room-echo widget
"""

# ╔═╡ 28ef5c14-eadc-4ff3-996a-006e34c4a33d
begin
    """
        RoomEchoInput(; delay_milliseconds=0)

    Control state for the room-echo demonstration. The bound dictionary reports the
    delay between the direct sound and one reflected copy in milliseconds.
    """
    struct RoomEchoInput
        delay_milliseconds::Int
    end

    RoomEchoInput(; delay_milliseconds=0) =
        RoomEchoInput(clamp(Int(delay_milliseconds), 0, 300))

    Base.get(w::RoomEchoInput) = Dict{String,Any}("delay_ms" => w.delay_milliseconds)

    function Base.show(io::IO, ::MIME"text/html", w::RoomEchoInput)
        write(io, """
<div id="rewidget">
  <style>
    #rewidget{width:100%;box-sizing:border-box;color:#d1d5db;font:14px sans-serif}
    #rewidget .re-title{width:100%;box-sizing:border-box;text-align:center;margin-bottom:10px;
      background:#0a0f18;border:1px solid #3b5c85;border-radius:6px;padding:10px 14px}
    #rewidget .re-title-desc{font-size:17px;font-weight:700;color:#e5e7eb}
    #rewidget .re-title-hint{font-size:13px;color:#9ca3af;margin-top:3px}
    #rewidget .re-workspace{background:#000;border:1px solid #374151;border-radius:6px;padding:10px}
    #rewidget canvas{display:block;width:100%;height:auto}
    #rewidget .re-controls{margin-top:8px;display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:8px}
    #rewidget .re-control-group{background:#050505;border:1px solid #2f3744;border-radius:6px;padding:12px}
    #rewidget .re-control-title{color:#f3f4f6;font-size:20px;font-weight:700;margin-bottom:8px}
    #rewidget .re-control-row{display:grid;grid-template-columns:minmax(80px,130px) minmax(70px,1fr) minmax(48px,72px);
      align-items:center;gap:8px;margin:8px 0;color:#9ca3af}
    #rewidget input[type=range]{width:100%;min-width:0;accent-color:#3b82f6}
    #rewidget .re-value{color:#e5e7eb;text-align:right;font-variant-numeric:tabular-nums}
    #rewidget .re-actions{display:flex;flex-wrap:wrap;gap:7px}
    #rewidget button,#rewidget .re-file-label{border-radius:4px;border:1px solid #9ca3af;background:#606060;color:#f3f4f6;
      padding:6px 12px;font-size:14px;cursor:pointer}
    #rewidget #re-direct{background:#075985;border-color:#38bdf8;color:#f3f4f6}
    #rewidget #re-direct:hover{background:#0c4a6e}
    #rewidget #re-echo{background:linear-gradient(110deg,#075985 0%,#075985 47%,#991b1b 53%,#991b1b 100%);
      border-color:#f3f4f6;color:#f3f4f6}
    #rewidget #re-echo:hover{background:linear-gradient(110deg,#0c4a6e 0%,#0c4a6e 47%,#b91c1c 53%,#b91c1c 100%)}
    #rewidget button[disabled]{opacity:.55;cursor:not-allowed}
    #rewidget .re-file-label{display:inline-block}
    #rewidget input[type=file]{display:none}
    #rewidget .re-presets{display:flex;flex-wrap:wrap;gap:6px;margin-top:10px}
    #rewidget .re-presets button{background:#1f2937}
    #rewidget .re-presets button.active{border-color:#38bdf8;background:#164e63}
    #rewidget .re-status{min-height:20px;color:#9ca3af;margin-top:10px}
  </style>
  <div class="re-title">
    <div class="re-title-desc">A reflection is a delayed copy of the direct wave.</div>
    <div class="re-title-hint">Record or load a sentence · choose a delay · hear the reflected arrival</div>
  </div>
  <div class="re-workspace">
    <canvas id="re-canvas" width="960" height="300" aria-label="Timeline of direct, reflected, and received sound waveforms"></canvas>
  </div>
  <div class="re-controls">
    <section class="re-control-group">
      <div class="re-control-title">Reflection delay</div>
      <label class="re-control-row" for="re-delay"><span>Delay</span>
        <input id="re-delay" type="range" min="0" max="300" step="1" value="$(w.delay_milliseconds)">
        <output id="re-delay-value" class="re-value"></output>
      </label>
      <div class="re-presets" aria-label="Delay presets">
        <button type="button" data-delay="0">0 ms</button>
        <button type="button" data-delay="5">5 ms</button>
        <button type="button" data-delay="10">10 ms</button>
        <button type="button" data-delay="20">20 ms</button>
        <button type="button" data-delay="40">40 ms</button>
        <button type="button" data-delay="80">80 ms</button>
        <button type="button" data-delay="300">300 ms</button>
      </div>
    </section>
    <section class="re-control-group">
      <div class="re-control-title">Source and playback</div>
      <div class="re-actions">
        <button id="re-direct" type="button">Play direct only</button>
        <button id="re-echo" type="button">Play direct + echo</button>
        <button id="re-record" type="button">Record sentence</button>
        <button id="re-stop" type="button" disabled>Stop recording</button>
        <label class="re-file-label" for="re-file">Load audio file</label>
        <input id="re-file" type="file" accept="audio/*">
      </div>
      <div id="re-status" class="re-status">Source: built-in voiced test signal</div>
    </section>
  </div>
</div>
<script>
  const par = currentScript.previousElementSibling
  const canvas = par.querySelector('#re-canvas')
  const ctx2d = canvas.getContext('2d')
  const plotWidth = 960
  const plotHeight = 300
  const delayInput = par.querySelector('#re-delay')
  const delayValue = par.querySelector('#re-delay-value')
  const status = par.querySelector('#re-status')
  const directButton = par.querySelector('#re-direct')
  const echoButton = par.querySelector('#re-echo')
  const recordButton = par.querySelector('#re-record')
  const stopButton = par.querySelector('#re-stop')
  const fileInput = par.querySelector('#re-file')
  const presetButtons = [...par.querySelectorAll('[data-delay]')]
  let audioContext = null
  let sourceBuffer = null
  let fallbackBuffer = null
  let recorder = null
  let recordingStream = null
  let chunks = []
  let playheadSeconds = null
  let animationFrame = null

  function resizeCanvas() {
    const cssWidth = canvas.clientWidth || plotWidth
    const dpr = Math.min(window.devicePixelRatio || 1, 2)
    canvas.width = Math.round(cssWidth * dpr)
    canvas.height = Math.round(cssWidth * plotHeight / plotWidth * dpr)
    canvas.style.height = Math.round(cssWidth * plotHeight / plotWidth) + 'px'
    ctx2d.setTransform(canvas.width / plotWidth, 0, 0, canvas.height / plotHeight, 0, 0)
  }

  function audio() {
    if (!audioContext) audioContext = new (window.AudioContext || window.webkitAudioContext)()
    return audioContext
  }

  function voicedTestSignal(context) {
    const rate = 24000
    const duration = 2.1
    const result = context.createBuffer(1, Math.floor(rate * duration), rate)
    const data = result.getChannelData(0)
    const centers = [0.20, 0.48, 0.76, 1.08, 1.40, 1.72]
    for (let i = 0; i < data.length; i++) {
      const t = i / rate
      let envelope = 0
      for (const center of centers) {
        const x = (t - center) / 0.11
        envelope += Math.exp(-x * x * 7)
      }
      const carrier = Math.sin(2 * Math.PI * 145 * t) +
        0.38 * Math.sin(2 * Math.PI * 290 * t) +
        0.16 * Math.sin(2 * Math.PI * 435 * t)
      data[i] = 0.16 * envelope * carrier
    }
    return result
  }

  function activeBuffer() {
    if (sourceBuffer) return sourceBuffer
    if (!fallbackBuffer) fallbackBuffer = voicedTestSignal(audio())
    return fallbackBuffer
  }

  function drawTrace(data, sampleRate, startX, baseline, timeScale, color) {
    ctx2d.strokeStyle = color
    ctx2d.lineWidth = 1.4
    ctx2d.beginPath()
    const points = 1200
    for (let i = 0; i < points; i++) {
      const sample = data[Math.min(data.length - 1, Math.floor(i * data.length / points))]
      const x = startX + i * data.length * timeScale / (points * sampleRate)
      const y = baseline - sample * 95
      if (i === 0) ctx2d.moveTo(x, y); else ctx2d.lineTo(x, y)
    }
    ctx2d.stroke()
  }

  function sampleAt(data, sampleRate, time) {
    const index = Math.floor(time * sampleRate)
    return index >= 0 && index < data.length ? data[index] : 0
  }

  function drawReceivedTrace(data, sampleRate, startX, baseline, timeScale, delay) {
    const totalDuration = data.length / sampleRate + delay / 1000
    const points = 1600
    ctx2d.strokeStyle = '#f3f4f6'
    ctx2d.lineWidth = 1.5
    ctx2d.beginPath()
    for (let i = 0; i < points; i++) {
      const time = i * totalDuration / (points - 1)
      const direct = sampleAt(data, sampleRate, time)
      const reflected = delay > 0 ? sampleAt(data, sampleRate, time - delay / 1000) : 0
      const x = startX + time * timeScale
      const y = baseline - (direct + reflected) * 82
      if (i === 0) ctx2d.moveTo(x, y); else ctx2d.lineTo(x, y)
    }
    ctx2d.stroke()
  }

  function draw() {
    const delay = Number(delayInput.value)
    const buffer = activeBuffer()
    const data = buffer.getChannelData(0)
    // Reserve a dedicated label column so text never collides with a waveform.
    const left = 250
    const right = 32
    const duration = buffer.duration + delay / 1000
    const scale = (plotWidth - left - right) / Math.max(duration, 0.4)
    ctx2d.fillStyle = '#000'
    ctx2d.fillRect(0, 0, plotWidth, plotHeight)
    ctx2d.strokeStyle = '#374151'
    ctx2d.lineWidth = 1
    ctx2d.beginPath()
    ctx2d.moveTo(left, 70); ctx2d.lineTo(plotWidth - right, 70)
    ctx2d.moveTo(left, 155); ctx2d.lineTo(plotWidth - right, 155)
    ctx2d.moveTo(left, 245); ctx2d.lineTo(plotWidth - right, 245)
    ctx2d.stroke()
    ctx2d.fillStyle = '#9ca3af'
    ctx2d.font = '16px sans-serif'
    ctx2d.fillText('direct arrival', 12, 74)
    ctx2d.fillText(delay === 0 ? 'no reflection' : 'reflected copy', 12, 159)
    ctx2d.fillText(delay > 0 ? 'receiver: direct + reflection' : 'receiver: direct only', 12, 249)
    if (delay > 0) {
      ctx2d.fillStyle = '#ef4444'
      ctx2d.textAlign = 'right'
      ctx2d.fillText(delay + ' ms reflection', plotWidth - right, 28)
      ctx2d.textAlign = 'left'
    }
    drawTrace(data, buffer.sampleRate, left, 70, scale, '#38bdf8')
    if (delay > 0) {
      drawTrace(data, buffer.sampleRate, left + delay * scale / 1000, 155, scale, '#ef4444')
    }
    drawReceivedTrace(data, buffer.sampleRate, left, 245, scale, delay)
    if (playheadSeconds !== null) {
      const x = left + Math.min(playheadSeconds, duration) * scale
      ctx2d.strokeStyle = '#fbbf24'
      ctx2d.lineWidth = 2
      ctx2d.setLineDash([6, 4])
      ctx2d.beginPath()
      ctx2d.moveTo(x, 18); ctx2d.lineTo(x, 280)
      ctx2d.stroke()
      ctx2d.setLineDash([])
    }
  }

  function sync() {
    const delay = Number(delayInput.value)
    delayValue.textContent = delay + ' ms'
    presetButtons.forEach(button => button.classList.toggle('active', Number(button.dataset.delay) === delay))
    draw()
  }

  function publish() {
    par.value = {delay_ms: Number(delayInput.value)}
    par.dispatchEvent(new CustomEvent('input'))
  }

  function schedule(buffer, when, gain) {
    const node = audio().createBufferSource()
    const volume = audio().createGain()
    node.buffer = buffer
    volume.gain.value = gain
    node.connect(volume).connect(audio().destination)
    node.start(when)
  }

  function animatePlay(start, totalDuration) {
    if (animationFrame) cancelAnimationFrame(animationFrame)
    function frame() {
      playheadSeconds = Math.max(0, audio().currentTime - start)
      draw()
      if (playheadSeconds < totalDuration) {
        animationFrame = requestAnimationFrame(frame)
      } else {
        playheadSeconds = null
        animationFrame = null
        draw()
      }
    }
    animationFrame = requestAnimationFrame(frame)
  }

  async function play(withEcho) {
    const context = audio()
    await context.resume()
    const buffer = activeBuffer()
    const delay = Number(delayInput.value)
    const start = context.currentTime + 0.06
    schedule(buffer, start, withEcho && delay > 0 ? 0.52 : 0.82)
    if (withEcho && delay > 0) schedule(buffer, start + delay / 1000, 0.52)
    animatePlay(start, buffer.duration + (withEcho ? delay / 1000 : 0))
    status.textContent = withEcho && delay > 0 ? 'Playing direct sound + reflected copy at ' + delay + ' ms.' : 'Playing the direct sound only.'
  }

  async function decode(arrayBuffer, label) {
    try {
      sourceBuffer = await audio().decodeAudioData(arrayBuffer)
      status.textContent = 'Source: ' + label
      draw()
    } catch (_) {
      status.textContent = 'This browser could not decode that audio file. Try a WAV, MP3, M4A, or WebM recording.'
    }
  }

  async function startRecording() {
    if (!navigator.mediaDevices || !window.MediaRecorder) {
      status.textContent = 'Recording is not available here. Load a short audio file instead.'
      return
    }
    try {
      recordingStream = await navigator.mediaDevices.getUserMedia({audio:true})
      recorder = new MediaRecorder(recordingStream)
      chunks = []
      recorder.addEventListener('dataavailable', event => { if (event.data.size) chunks.push(event.data) })
      recorder.addEventListener('stop', async () => {
        const clip = new Blob(chunks, {type: recorder.mimeType})
        recordingStream.getTracks().forEach(track => track.stop())
        await decode(await clip.arrayBuffer(), 'recorded sentence')
        recordButton.disabled = false
        stopButton.disabled = true
      }, {once:true})
      recorder.start()
      recordButton.disabled = true
      stopButton.disabled = false
      status.textContent = 'Recording… say: “The Earth has many seismic reflections.”'
    } catch (_) {
      status.textContent = 'Microphone access was not granted. You can still load an audio file or use the built-in test signal.'
    }
  }

  delayInput.addEventListener('input', () => { sync(); publish() })
  presetButtons.forEach(button => button.addEventListener('click', () => {
    delayInput.value = button.dataset.delay
    sync()
    publish()
  }))
  directButton.addEventListener('click', () => { play(false) })
  echoButton.addEventListener('click', () => { play(true) })
  recordButton.addEventListener('click', startRecording)
  stopButton.addEventListener('click', () => { if (recorder && recorder.state === 'recording') recorder.stop() })
  fileInput.addEventListener('change', async () => {
    const file = fileInput.files && fileInput.files[0]
    if (file) await decode(await file.arrayBuffer(), file.name)
  })

  resizeCanvas()
  const resizeObserver = new ResizeObserver(() => { resizeCanvas(); draw() })
  resizeObserver.observe(canvas)
  if (typeof invalidation !== 'undefined') invalidation.then(() => resizeObserver.disconnect())
  sync()
  par.value = {delay_ms: Number(delayInput.value)}
</script>
""")
    end

    const _re_ready = true
end

# ╔═╡ 5ef3b597-6e38-4f4d-8428-d4813844e55b
begin
    _re_ready
    PlutoUI.WideCell(@bind room_echo RoomEchoInput(); max_width=1100)
end

# ╔═╡ 06fafab8-a76c-4c99-9016-10dc73d8c99c
begin
    echo_delay = Int(room_echo["delay_ms"])
    hearing_description =
        echo_delay == 0 ? "one sound" :
        echo_delay <= 20 ? "one sound — the copy fuses with the direct arrival" :
        echo_delay <= 40 ? "the beginning of an echo" :
        "a clear echo"
end

# ╔═╡ b505b02d-12ce-4a89-b609-7c7ffd5e8f7c
md"""
### What students should hear at **$(echo_delay) ms**

**$(hearing_description)**

| Delay | Typical perception |
|---:|:---|
| 0 ms | one sound |
| 5 ms | one sound |
| 10 ms | one sound |
| 20 ms | one sound |
| 40 ms | beginning of echo |
| 80 ms | clear echo |

The source waveform has not changed. Only the arrival time of its copy has changed.
That is the essential move in reflection seismology.

!!! warning "Who is fooling us: the eyes or the ears?"
    Neither. The blue direct trace and red reflected trace really are different incoming
    waveforms, and the receiver really contains both. But for a short enough delay, the
    auditory system infers **one source** and fuses the arrivals into one sound. A
    seismometer does not make that inference: it leaves the direct arrival and reflection
    together for us to separate.
"""

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
# ╠═5c5c34cc-b71d-4bdc-a857-8763ae803e8c
# ╟─a63f99b1-a70d-4a17-bcdb-91f0ef0b3963
# ╟─6a4cd45d-4e9e-4fd7-8856-22d77bcac1b7
# ╟─5ef3b597-6e38-4f4d-8428-d4813844e55b
# ╟─06fafab8-a76c-4c99-9016-10dc73d8c99c
# ╟─b505b02d-12ce-4a89-b609-7c7ffd5e8f7c
# ╟─7f0c9d2e-3b4a-4a4b-8e7f-1c2d3e4f5a6b
# ╟─28b6cf3a-61c9-40d4-9d0f-4f290d21874b
# ╟─47ee7135-61df-4c58-bcfc-07d23a65d6e6
# ╟─1de07b55-f4c8-4b6b-b090-f20eced1bea0
# ╟─28ef5c14-eadc-4ff3-996a-006e34c4a33d
# ╠═9e9c44d5-7940-4cc0-9f68-f79a315c5a27
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
