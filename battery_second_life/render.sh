#!/usr/bin/env bash
# Render the video in two versions and stitch them.
#   renders/preview.mp4           silhouettes mark where the presenter goes (for review)
#   renders/clean_scenes_1-4.mp4  empty presenter windows ready for your footage
#   renders/S5_overlay.mov        transparent lower third to lay over full-screen presenter footage
#   renders/S5_overlay.webm       the same overlay in a smaller file
# Needs manim>=0.18 and ffmpeg. No LaTeX needed.
set -euo pipefail
cd "$(dirname "$0")"
MEDIA="${1:-media}"
mkdir -p ~/.fonts renders && cp fonts/*.ttf ~/.fonts/ && fc-cache -f >/dev/null

SCENES=(S1_Hook S2_Extraction S3_Journey S4_Recycling S5_Closing)
join() {  # join <dir> <out> <scene...>
  local dir=$1 out=$2; shift 2
  : > "$dir/list.txt"
  for s in "$@"; do echo "file '$s.mp4'" >> "$dir/list.txt"; done
  ffmpeg -y -v error -f concat -safe 0 -i "$dir/list.txt" -c:v libx264 -pix_fmt yuv420p -crf 20 -movflags +faststart "$out"
}

for s in "${SCENES[@]}"; do
  PRESENTER_PLACEHOLDER=1 manim -r 1920,1080 --fps 30 --media_dir "$MEDIA/preview" battery_second_life.py "$s"
done
join "$MEDIA/preview/videos/battery_second_life/1080p30" renders/preview.mp4 "${SCENES[@]}"

for s in "${SCENES[@]:0:4}"; do
  manim -r 1920,1080 --fps 30 --media_dir "$MEDIA/clean" battery_second_life.py "$s"
done
join "$MEDIA/clean/videos/battery_second_life/1080p30" renders/clean_scenes_1-4.mp4 "${SCENES[@]:0:4}"

manim -r 1920,1080 --fps 30 --transparent --format mov --media_dir "$MEDIA/overlay" battery_second_life.py S5_Closing
cp "$MEDIA/overlay/videos/battery_second_life/1080p30/S5_Closing.mov" renders/S5_overlay.mov
ffmpeg -y -v error -i renders/S5_overlay.mov -c:v libvpx-vp9 -pix_fmt yuva420p -b:v 2M renders/S5_overlay.webm
echo "Done. See renders/"
