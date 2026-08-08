# VM Deployment Implementation Summary

## What's Been Created ✅

### Configuration Files

1. **`live-notebooks.yml`** (NEW)
   - Registry of notebooks to run live on your VM
   - Currently configured with: `misc/geoid-kernel.jl` (test)
   - Edit to add/remove notebooks → auto-deploys to VM

2. **`static-notebooks.yml`** (NEW)
   - Registry of notebooks for GitHub Pages (future use)
   - Currently configured with: `misc/seismic-interferometry.jl` (test)
   - Allows organizing static content separately from live

### Deployment Scripts

3. **`deploy-to-vm.sh`** (NEW)
   - Automation script that:
     - Pulls latest code from git
     - Reads `live-notebooks.yml`
     - Exports selected notebooks via PlutoSliderServer
     - Restarts systemd service
   - Run manually or via GitHub Actions

### Docker & CI/CD

4. **`Dockerfile.vm`** (NEW)
   - Simplified Docker image for VM deployment
   - Runs PlutoSliderServer on port 8080
   - Can be used instead of systemd service if preferred

5. **`.github/workflows/DeployToVM.yml`** (NEW)
   - GitHub Actions workflow for auto-deployment
   - Triggers on: push to main, manual dispatch
   - Connects to VM via SSH and runs deploy script
   - Requires 4 GitHub secrets to be configured

### Documentation

6. **`VM_DEPLOYMENT_SETUP.md`** (NEW)
   - Complete step-by-step setup guide
   - VM preparation instructions
   - Systemd service setup
   - GitHub Actions secret configuration
   - Troubleshooting guide

7. **`DEPLOYMENT_SUMMARY.md`** (THIS FILE)
   - Quick reference

8. **`README.md`** (UPDATED)
   - Added VM deployment option documentation
   - Links to setup guide

---

## Next Steps

### Phase 1: VM Preparation (Manual)

You need to do these steps on your VM:

1. **Install Julia 1.11.7+** on your VM
2. **Clone the repo** to your VM:
   ```bash
   git clone https://github.com/pawbz/Interactive-Seismology.jl.git
   cd Interactive-Seismology.jl
   ```
3. **Create systemd service** (follow `VM_DEPLOYMENT_SETUP.md`, Step 2)
   - Create `/etc/systemd/system/pluto-slider-server.service`
   - Enable and start it
4. **Test manual deployment**:
   ```bash
   chmod +x deploy-to-vm.sh
   ./deploy-to-vm.sh
   ```

### Phase 2: GitHub Actions Setup (Local)

Configure auto-deployment:

1. **Generate SSH key pair** for GitHub Actions:
   ```bash
   ssh-keygen -t ed25519 -f ~/.ssh/vm-deploy
   # Add public key to VM: ssh-copy-id -i ~/.ssh/vm-deploy.pub user@vm
   ```

2. **Add GitHub Secrets** (Repo → Settings → Secrets and variables → Actions):
   - `VM_HOST` - Your VM's IP or domain
   - `VM_USER` - SSH user (e.g., `pluto`, `ubuntu`)
   - `VM_SSH_KEY` - Private SSH key content (from `~/.ssh/vm-deploy`)
   - `VM_DEPLOY_PATH` - Path to cloned repo on VM (e.g., `/home/pluto/Interactive-Seismology.jl`)

3. **Test the workflow**:
   ```bash
   git add -A
   git commit -m "Initial VM deployment setup"
   git push origin main
   # Watch: GitHub → Actions → "Deploy Live Notebooks to VM"
   ```

### Phase 3: Gradual Rollout (Ongoing)

Add more notebooks incrementally:

```bash
# Edit to add more notebooks
nano live-notebooks.yml

# Example:
# sections:
#   Test:
#     - misc/geoid-kernel.jl
#   Fundamentals:
#     - introduction/waves-on-string.jl

git add live-notebooks.yml
git commit -m "Add more notebooks to live deployment"
git push origin main
# GitHub Actions automatically deploys
```

---

## Current Test Configuration

**Live notebooks** (running on VM):
- `misc/geoid-kernel.jl` (interactive globe visualization)

**Static notebooks** (on GitHub Pages):
- `misc/seismic-interferometry.jl` (static, Binder fallback)

Once you verify the pipeline works with these test notebooks, expand the registries to include your full course content.

---

## Architecture at a Glance

```
You push to main
    ↓
GitHub Actions triggered
    ├─→ ExportNotebooks.yml (existing - GitHub Pages)
    │   └─→ Generate static site from all notebooks
    │       └─→ Deploy to gh-pages (all static)
    │
    └─→ DeployToVM.yml (new - auto-deploy)
        ├─→ SSH into your VM
        ├─→ Run deploy-to-vm.sh
        │   ├─→ Pull latest code
        │   ├─→ Generate complete static site (like GitHub Pages)
        │   └─→ Copy to serving directory
        ├─→ Restart PlutoSliderServer
        │   ├─→ Serves all notebooks as base HTML
        │   └─→ Makes live-notebooks.yml entries interactive
        └─→ Restart complete
        
Result:
├─ GitHub Pages (all static): https://pawbz.github.io/Interactive-Seismology.jl/
│
└─ Your VM (hybrid): http://your-vm-ip:8080/
    ├─ Same site layout & navigation as GitHub Pages
    ├─ Static notebooks: served as pre-rendered HTML
    └─ Live notebooks (from live-notebooks.yml): interactive via PlutoSliderServer
```

---

## Key Files Reference

| File | Purpose | Edit? | Notes |
|------|---------|-------|-------|
| `live-notebooks.yml` | Which notebooks go live | ✏️ Yes | Edit to deploy/remove |
| `static-notebooks.yml` | Which notebooks on GitHub Pages | Later | For future site organization |
| `deploy-to-vm.sh` | Deployment automation | Usually no | Only if changing deploy logic |
| `Dockerfile.vm` | Docker image for VM | Optional | Alternative to systemd |
| `.github/workflows/DeployToVM.yml` | GitHub Actions workflow | Rarely | Only if changing CI/CD |
| `VM_DEPLOYMENT_SETUP.md` | Setup instructions | Reference only | Comprehensive guide |

---

## Verification Checklist

Before considering this complete:

- [ ] VM has Julia 1.11.7+ installed
- [ ] Repository cloned on VM
- [ ] systemd service created and running (`pluto-slider-server`)
- [ ] Manual deployment script works: `./deploy-to-vm.sh`
- [ ] Static site generated in `_site/` directory
- [ ] PlutoSliderServer running on port 8080
- [ ] GitHub secrets configured (4 secrets)
- [ ] First push triggers DeployToVM workflow successfully
- [ ] VM site accessible at: `http://vm-ip:8080/`
- [ ] Static notebook appears as HTML: `http://vm-ip:8080/misc/seismic-interferometry/`
- [ ] Live notebook is interactive: `http://vm-ip:8080/misc/geoid-kernel/` (widgets work)
- [ ] Navigation and layout identical to GitHub Pages
- [ ] GitHub Pages still works: https://pawbz.github.io/Interactive-Seismology.jl/

---

## Troubleshooting

**Workflow fails?** Check:
- GitHub secrets are set correctly and non-empty
- VM is online and reachable
- SSH key has access to VM
- `deploy-to-vm.sh` is executable on VM

**Notebooks not appearing?** Check:
- `live-notebooks.yml` syntax is correct
- Notebook paths exist: `ls src/misc/geoid-kernel.jl`
- PlutoSliderServer service is running: `sudo systemctl status pluto-slider-server`
- Port 8080 is accessible

**See detailed troubleshooting in `VM_DEPLOYMENT_SETUP.md`**

---

## Next Phase (After Validation)

Once the test notebooks work:

1. **Expand registries**: Add all your notebooks to both YAML files
2. **Organize by sections**: Group notebooks into teaching topics
3. **Set up domain**: Use custom domain + SSL instead of IP:port
4. **Add monitoring**: Uptime checks, error logging
5. **Optimize**: Cache notebook outputs, consider performance tuning

---

Questions? Refer to:
- VM setup: [VM_DEPLOYMENT_SETUP.md](VM_DEPLOYMENT_SETUP.md)
- Project plan: [.claude/plans/crystalline-cuddling-bengio.md](.claude/plans/crystalline-cuddling-bengio.md)
- PlutoSliderServer: https://github.com/JuliaPluto/PlutoSliderServer.jl
