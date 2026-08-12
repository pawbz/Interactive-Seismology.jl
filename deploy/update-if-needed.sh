#!/usr/bin/env bash

set -euo pipefail

REPO_DIR="/var/lib/interactive-seismology"
DEPLOYED_REV_FILE="$REPO_DIR/.deployed-main-revision"

cd "$REPO_DIR"
git fetch --quiet origin main

remote_revision="$(git rev-parse origin/main)"
deployed_revision="$(cat "$DEPLOYED_REV_FILE" 2>/dev/null || true)"

if [[ "$remote_revision" == "$deployed_revision" ]]; then
    echo "No update needed; origin/main is already deployed at $remote_revision"
    exit 0
fi

echo "Deploying origin/main at $remote_revision"
bash "$REPO_DIR/deploy-to-vm.sh"
deployed_revision="$(git rev-parse HEAD)"
printf '%s\n' "$deployed_revision" > "$DEPLOYED_REV_FILE"
echo "Deployment recorded at $deployed_revision"
