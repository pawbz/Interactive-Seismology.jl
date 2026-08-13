# VM Deployment Setup Guide

This guide walks you through setting up live notebook deployment to your VM using PlutoSliderServer.

## Overview

**The VM serves the exact same site as GitHub Pages, but with selective interactivity:**

- **Complete site generated**: All notebooks exported as static HTML (like GitHub Pages)
- **Selective interactivity**: Notebooks in `live-notebooks.yml` become interactive via PlutoSliderServer
- **Same navigation**: Identical layout, sections, and organization as pawbz.github.io
- **Auto-deployment**: Updates via GitHub Actions when you push to main

**Result**: Your VM is a complete mirror of GitHub Pages, but with live-editable notebooks

## Step 1: Prepare Your VM

### 1.1 Install Dependencies

```bash
# SSH into your VM
ssh user@your-vm-ip

# Install Julia 1.12.0 (if not already installed)
# See https://julialang.org/downloads/platform/ for your OS
# For Ubuntu/Debian:
# curl -fsSL https://install.julialang.org | sh
# Or install Julia manually from https://julialang.org/downloads/

# Install git
sudo apt-get update
sudo apt-get install git

# Verify installations
julia --version
git --version
```

### 1.2 Clone Repository

```bash
# Create deployment directory
mkdir -p ~/notebooks
cd ~/notebooks

# Clone the repository
git clone https://github.com/pawbz/Interactive-Seismology.jl.git
cd Interactive-Seismology.jl
```

### 1.3 Test Manual Deployment

```bash
# Make the deploy script executable
chmod +x deploy-to-vm.sh

# Test deployment (this will export live notebooks)
./deploy-to-vm.sh

# Check if it worked
echo "Status: $?"
```

## Step 2: Set Up Systemd Service

Create a systemd service to run PlutoSliderServer automatically.

### 2.1 Create service file

```bash
sudo nano /etc/systemd/system/pluto-slider-server.service
```

Paste this content (adjust paths as needed):

```ini
[Unit]
Description=Pluto Slider Server for Interactive Seismology
After=network.target

[Service]
Type=simple
User=pluto
WorkingDirectory=/home/pluto/Interactive-Seismology.jl
Environment="PATH=/home/pluto/.julia/bin:/usr/local/bin:/usr/bin"

# Run only the notebooks in live-notebooks.yml.
# Keep this port private; nginx serves the generated _site directory separately.
ExecStart=/usr/bin/julia --project=pluto-deployment-environment \
  serve_live_notebooks.jl --host 127.0.0.1 --port 1234

Restart=always
RestartSec=10s
StandardOutput=journal
StandardError=journal

[Install]
WantedBy=multi-user.target
```

### 2.2 Enable and start service

```bash
sudo systemctl daemon-reload
sudo systemctl enable pluto-slider-server
sudo systemctl start pluto-slider-server

# Check status
sudo systemctl status pluto-slider-server

# View logs
sudo journalctl -u pluto-slider-server -f
```

### 2.3 Serve the generated site with nginx

Install nginx and point its document root to the directory used by
`deploy-to-vm.sh` (`/home/pluto/_site` by default). Proxy only the two
PlutoSliderServer action endpoints; the rest of the site remains static.

```nginx
server {
    listen 80;
    server_name notebooks.example.edu;
    root /home/pluto/_site;

    location /staterequest/ {
        proxy_pass http://127.0.0.1:1234;
    }

    location /bondconnections/ {
        proxy_pass http://127.0.0.1:1234;
    }

    location / {
        try_files $uri $uri/ =404;
    }
}
```

With this configuration, `deploy-to-vm.sh` exports live pages with `/` as their
slider-server URL, so they use the same public domain while PlutoSliderServer
remains private on port 1234. Add TLS before making the site public.

## Step 3: Configure GitHub Secrets

For auto-deployment to work, add these secrets to your GitHub repository:

1. Go to: Settings → Secrets and variables → Actions
2. Add these secrets:

| Secret | Value | Example |
|--------|-------|---------|
| `VM_HOST` | Your VM's IP or domain | `192.168.1.100` or `seismology.myuniv.edu` |
| `VM_USER` | SSH user on VM | `pluto` or `ubuntu` |
| `VM_SSH_KEY` | Private SSH key for authentication | (copy your private key) |
| `VM_DEPLOY_PATH` | Path to cloned repo | `/home/pluto/Interactive-Seismology.jl` |

### Generating SSH Key (if needed)

```bash
# On your local machine
ssh-keygen -t ed25519 -f ~/.ssh/vm-deploy -C "github-actions-vm-deploy"

# Copy public key to VM
ssh-copy-id -i ~/.ssh/vm-deploy.pub user@your-vm-ip

# Add to GitHub as VM_SSH_KEY secret (copy the PRIVATE key)
cat ~/.ssh/vm-deploy
```

## Step 4: Test the Pipeline

### 4.1 Local Test

```bash
# Edit live-notebooks.yml to add/remove notebooks
# Example:
cat > live-notebooks.yml << 'EOF'
sections:
  Test:
    - misc/geoid-kernel.jl
EOF

# Commit and push
git add live-notebooks.yml
git commit -m "Test: deploy single notebook to VM"
git push origin main
```

### 4.2 Check GitHub Actions

1. Go to your GitHub repo → Actions tab
2. Find the "Deploy Live Notebooks to VM" workflow
3. Wait for it to complete (should be green ✅)

### 4.3 Verify on VM

```bash
# SSH into VM and check service is running
ssh user@your-vm-ip
sudo systemctl status pluto-slider-server

# Check if notebooks were exported
ls -la /var/lib/pluto-slider-server/_site/misc/

# Try accessing in browser
# http://your-vm-ip:8080/misc/geoid-kernel/
```

## Step 5: Managing Notebooks

### Add a Live Notebook

```bash
# Edit live-notebooks.yml
# Add the notebook path under appropriate section
git add live-notebooks.yml
git commit -m "Add notebook to live deployment"
git push origin main

# GitHub Actions will automatically deploy it
```

### Remove a Live Notebook

```bash
# Edit live-notebooks.yml and delete the line
# GitHub Actions will automatically remove it on next deployment
git add live-notebooks.yml
git commit -m "Remove notebook from live deployment"
git push origin main
```

## Troubleshooting

### Workflow Fails: "SSH deployment failed"

**Check:**
1. GitHub secrets are set correctly (`VM_HOST`, `VM_USER`, `VM_SSH_KEY`, `VM_DEPLOY_PATH`)
2. VM is reachable: `ssh -i /path/to/key user@vm-host`
3. Repository is cloned on VM at `VM_DEPLOY_PATH`
4. `deploy-to-vm.sh` exists in repo

### Service Won't Start

```bash
# Check service logs
sudo journalctl -u pluto-slider-server -n 50

# Check if port 8080 is already in use
sudo lsof -i :1234

# Try starting manually
cd /home/pluto/Interactive-Seismology.jl
julia --project=pluto-deployment-environment \
  serve_live_notebooks.jl --host 127.0.0.1 --port 1234
```

### Notebooks Not Appearing

```bash
# Check if export directory exists
ls -la /var/lib/pluto-slider-server/_site/

# Check live-notebooks.yml is being read correctly
cat live-notebooks.yml

# Re-run deploy script manually
cd /home/pluto/Interactive-Seismology.jl
./deploy-to-vm.sh
```

## Next Steps

- **Gradual rollout**: Start with 1-2 test notebooks, add more as you gain confidence
- **Domain setup**: Set up a custom domain and SSL certificate (Let's Encrypt)
- **Reverse proxy**: Use Nginx to proxy requests if behind a custom domain
- **Monitoring**: Add monitoring/alerting for uptime and performance
- **Single registry**: Use the unified `live-notebooks.yml` with per-section `live:` and `static:` entries

## Security Notes

- Keep `VM_SSH_KEY` secret secure
- Consider restricting VM access by IP (VPN)
- Use firewall rules to limit port 8080 access
- Monitor SSH logs for unauthorized access attempts

---

For questions or issues, refer to:
- [PlutoSliderServer.jl docs](https://github.com/JuliaPluto/PlutoSliderServer.jl)
- [Julia documentation](https://docs.julialang.org)
