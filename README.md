# Interactive Seismology Notebooks

A [collection](https://pawbz.github.io/Interactive-Seismology.jl/) of *reactive notebooks* using the [Julia](https://julialang.org/) notebook system, Pluto. These notebooks are designed to illustrate fundamental concepts in seismology in an interactive manner. They are primarily intended to support our courses on Seismology (ES218) and Inverse Problems (ES219), and are accessible to students with minimal programming experience. The notebooks contain all the essential mathematical expressions, which are easily manipulated using [Symbolics.jl](https://symbolics.juliasymbolics.org/dev/), a modern Computer Algebra System (CAS). By using these standalone notebooks, students can intuitively explore and gain a better understanding of wave phenomena in seismology.

## Deployment Options

### GitHub Pages (Default)
- View notebooks: [pawbz.github.io/Interactive-Seismology.jl](https://pawbz.github.io/Interactive-Seismology.jl/)
- Static HTML + Binder fallback
- Automatically updated on push to main

### VM Deployment (Self-Hosted)
- Run live interactive notebooks on your own infrastructure
- Configure which notebooks are "live" via `live-notebooks.yml`
- Auto-deploy via GitHub Actions with SSH

**Setup guide:** See [VM_DEPLOYMENT_SETUP.md](VM_DEPLOYMENT_SETUP.md) for complete instructions

### Local deployment preview

`live-notebooks.yml` is the single source of truth for which Pluto notebooks appear
on the site. To build and test that exact selection locally, run:

```bash
./local_serve_test
```

The static site is available at <http://127.0.0.1:8000/> and the live notebook
service at <http://127.0.0.1:1234/>. Stop the command with `Ctrl+C`. If either
port is occupied, choose replacements, for example
`STATIC_PORT=8100 LIVE_PORT=17321 ./local_serve_test`.

GitHub Pages always publishes the static variant of this site. The remote VM
scripts set `PLUTO_SLIDER_SERVER_URL` explicitly when selected notebooks should
connect to a live SliderServer endpoint.

**Quick start:**
```bash
# Edit which notebooks are live
nano live-notebooks.yml

# Commit and push - auto-deploys to your VM
git add live-notebooks.yml
git commit -m "Update live notebooks"
git push origin main
```
