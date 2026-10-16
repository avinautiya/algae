#!/usr/bin/env bash
# Render all five scenes at 1080p30 and stitch them into one MP4.
# Usage: ./render.sh [media_dir]   (needs manim>=0.18 and ffmpeg; no LaTeX)
set -euo pipefail
cd "$(dirname "$0")"
MEDIA="${1:-media}"
SCENES=(S1_Smartphone S2_Extraction S3_Processing S4_Recycling S5_ClosedLoop)
for s in "${SCENES[@]}"; do
  manim -r 1920,1080 --fps 30 --media_dir "$MEDIA" battery_second_life.py "$s"
done
DIR="$MEDIA/videos/battery_second_life/1080p30"
: > "$DIR/list.txt"
for s in "${SCENES[@]}"; do echo "file '$s.mp4'" >> "$DIR/list.txt"; done
mkdir -p renders
ffmpeg -y -v error -f concat -safe 0 -i "$DIR/list.txt" -c:v libx264 -pix_fmt yuv420p -crf 20 \
  -movflags +faststart renders/from_mine_to_market.mp4
echo "Wrote renders/from_mine_to_market.mp4"
