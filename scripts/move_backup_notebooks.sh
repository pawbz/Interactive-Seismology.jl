#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TARGET_DIR="$REPO_ROOT/backup_notebooks"
DRY_RUN=false

if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=true
elif [[ -n "${1:-}" ]]; then
  echo "Usage: $0 [--dry-run]" >&2
  exit 1
fi

mkdir -p "$TARGET_DIR"

count=0
while IFS= read -r -d '' file; do
  rel_path="${file#$REPO_ROOT/}"
  relative_in_src="${rel_path#src/}"
  dest="$TARGET_DIR/$relative_in_src"

  if [[ -e "$dest" ]]; then
    stem="${dest%.*}"
    ext="${dest##*.}"
    n=1
    while [[ -e "${stem}_${n}.${ext}" ]]; do
      ((n++))
    done
    dest="${stem}_${n}.${ext}"
  fi

  if [[ "$DRY_RUN" == true ]]; then
    echo "Would move: $rel_path -> ${dest#$REPO_ROOT/}"
  else
    mkdir -p "$(dirname "$dest")"
    mv "$file" "$dest"
    echo "Moved: $rel_path -> ${dest#$REPO_ROOT/}"
  fi

  ((count++))

done < <(find "$REPO_ROOT/src" -type f -iname '* backup *.jl' -print0)

if [[ "$count" -eq 0 ]]; then
  echo "No backup notebooks found under src/."
else
  echo "Processed $count backup notebook(s)."
fi
