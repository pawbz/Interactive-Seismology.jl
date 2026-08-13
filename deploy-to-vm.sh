#!/bin/bash
# Deploy complete site to VM with selective notebook interactivity
# VM serves identical site to GitHub Pages, but selected notebooks are LIVE

set -euo pipefail

REPO_DIR="${1:-.}"
SERVICE_NAME="pluto-slider-server"
EXPORT_DIR="${SITE_DIR:-/var/lib/interactive-seismology/_site-serving}"
BASE_PATH="${BASE_PATH:-/~interactive-seismology/}"
DEPOT_PATH="${JULIA_DEPOT_PATH:-/var/lib/interactive-seismology/.julia-depot}"

if [ -z "$EXPORT_DIR" ] || [ "$EXPORT_DIR" = "/" ]; then
    echo "❌ SITE_DIR must name a specific directory, not the filesystem root" >&2
    exit 2
fi

echo "🚀 Deploying site to VM (with live notebooks from live-notebooks.yml)..."

# Navigate to repo
cd "$REPO_DIR"

# Pull latest from git
echo "📦 Pulling latest from git..."
git pull --ff-only origin main

# Generate the complete static site using PlutoPages
# This exports all notebooks listed in the unified live-notebooks.yml registry
echo "🔧 Generating complete static site via PlutoPages.jl..."
JULIA_DEPOT_PATH="$DEPOT_PATH" \
PLUTO_SLIDER_SERVER_URL="${PLUTO_SLIDER_SERVER_URL:-$BASE_PATH}" \
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
sudo systemctl restart "$SERVICE_NAME"
sudo systemctl is-active --quiet "$SERVICE_NAME"

echo "✅ Deployment complete!"
echo "📍 Site available at: $BASE_PATH"
echo "ℹ️  Live notebooks (from live-notebooks.yml) are interactive"
echo "ℹ️  Static notebooks appear as regular HTML pages"
