#!/bin/bash
# Deploy complete site to VM with selective notebook interactivity
# VM serves identical site to GitHub Pages, but selected notebooks are LIVE

set -e  # Exit on error

REPO_DIR="${1:-.}"
SERVICE_NAME="pluto-slider-server"
EXPORT_DIR="${SITE_DIR:-${HOME}/_site}"

if [ -z "$EXPORT_DIR" ] || [ "$EXPORT_DIR" = "/" ]; then
    echo "❌ SITE_DIR must name a specific directory, not the filesystem root" >&2
    exit 2
fi

echo "🚀 Deploying site to VM (with live notebooks from live-notebooks.yml)..."

# Navigate to repo
cd "$REPO_DIR"

# Pull latest from git
echo "📦 Pulling latest from git..."
git pull origin main

# Generate the complete static site using PlutoPages
# This exports all notebooks listed in the unified live-notebooks.yml registry
echo "🔧 Generating complete static site via PlutoPages.jl..."
PLUTO_SLIDER_SERVER_URL="${PLUTO_SLIDER_SERVER_URL:-/}" \
julia --project=pluto-deployment-environment -e "
    import Pkg
    Pkg.instantiate()
    include(\"generate.jl\")
" || {
    echo "❌ Site generation failed"
    exit 1
}

# Replace the served site atomically, so removing a notebook from either YAML
# registry removes its old page too.
echo "📂 Deploying static site..."
if [ -d "_site" ]; then
    STAGING_DIR="${EXPORT_DIR}.next"
    PREVIOUS_DIR="${EXPORT_DIR}.previous"
    rm -rf "$STAGING_DIR" "$PREVIOUS_DIR"
    mkdir -p "$STAGING_DIR"
    cp -R _site/. "$STAGING_DIR/"
    if [ -e "$EXPORT_DIR" ]; then
        mv "$EXPORT_DIR" "$PREVIOUS_DIR"
    fi
    mv "$STAGING_DIR" "$EXPORT_DIR"
    rm -rf "$PREVIOUS_DIR"
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
