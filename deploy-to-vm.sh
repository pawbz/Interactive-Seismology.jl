#!/bin/bash
# Deploy complete site to VM with selective notebook interactivity
# VM serves identical site to GitHub Pages, but selected notebooks are LIVE

set -e  # Exit on error

REPO_DIR="${1:-.}"
SERVICE_NAME="pluto-slider-server"
EXPORT_DIR="${HOME}/_site"

echo "🚀 Deploying site to VM (with live notebooks from live-notebooks.yml)..."

# Navigate to repo
cd "$REPO_DIR"

# Pull latest from git
echo "📦 Pulling latest from git..."
git pull origin main

# Generate the complete static site using PlutoPages
# This exports ALL notebooks from both live-notebooks.yml and static-notebooks.yml
echo "🔧 Generating complete static site via PlutoPages.jl..."
julia --project=pluto-deployment-environment -e "
    import Pkg
    Pkg.instantiate()
    include(\"generate.jl\")
" || {
    echo "❌ Site generation failed"
    exit 1
}

# Copy generated site to serving directory
echo "📂 Deploying static site..."
mkdir -p "$EXPORT_DIR"
if [ -d "_site" ]; then
    cp -r _site/* "$EXPORT_DIR/" 2>/dev/null || true
    echo "✅ Static site deployed to $EXPORT_DIR"
else
    echo "⚠️  _site directory not found"
fi

# Restart PlutoSliderServer (for live notebook interactivity)
echo "🔄 Restarting PlutoSliderServer (for live interactivity)..."
sudo systemctl restart $SERVICE_NAME || echo "⚠️  Could not restart service (may require manual restart)"

echo "✅ Deployment complete!"
echo "📍 Site available at: http://localhost:8080/"
echo "ℹ️  Live notebooks (from live-notebooks.yml) are interactive"
echo "ℹ️  Static notebooks appear as regular HTML pages"
