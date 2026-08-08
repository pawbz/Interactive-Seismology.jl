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

**Quick start:**
```bash
# Edit which notebooks are live
nano live-notebooks.yml

# Commit and push - auto-deploys to your VM
git add live-notebooks.yml
git commit -m "Update live notebooks"
git push origin main
```