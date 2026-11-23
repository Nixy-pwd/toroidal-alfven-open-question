#!/usr/bin/env bash
set -euo pipefail

SOURCE="$(cd "$(dirname "$0")" && pwd)"
REMOTE="https://github.com/Nixy-pwd/toroidal-alfven-open-question.git"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

if ! command -v git >/dev/null 2>&1; then
  echo "Git is required." >&2
  exit 1
fi
if ! command -v rsync >/dev/null 2>&1; then
  echo "rsync is required (normally available on macOS)." >&2
  exit 1
fi

git clone "$REMOTE" "$TMP/repo"
rsync -a --delete --exclude='.git' "$SOURCE/" "$TMP/repo/"
git -C "$TMP/repo" add -A
if git -C "$TMP/repo" diff --cached --quiet; then
  echo "No changes to publish."
else
  git -C "$TMP/repo" commit -m "Add complete reproducibility data and research materials"
  git -C "$TMP/repo" push origin main
fi