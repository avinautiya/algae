# From Mine to Market: A Battery's Second Life

An 85-second explainer for a geoscience competition. It mixes Manim animation with a presenter on camera.

## Files

| File | Use |
|---|---|
| `renders/preview.mp4` | Full 1:25 cut. Grey silhouettes show where your footage goes. |
| `renders/clean_scenes_1-4.mp4` | Scenes 1 to 4 (1:10) with empty presenter windows. Drop your footage into them. |
| `renders/S5_overlay.webm` | Scene 5 lower third with a transparent background. Lay it over your full-screen closing shot. |
| `renders/S5_overlay.mov` | The same overlay as QuickTime with alpha, for editors that don't read WebM alpha. `./render.sh` creates it; it isn't committed because of its size. |
| `battery_second_life.py` | Source. One class per scene. |
| `fonts/` | Fraunces (headlines) and Nunito Sans (labels). Both are SIL Open Font License. |

To re-render everything, run `./render.sh`. You need Manim and ffmpeg but not LaTeX.

## Look

The page is light warm paper (`#F1EEE4`) under a hand-drawn landscape of soft hills, wobbly-edged trees, pines, clouds and birds. The colours are muted earth tones: ochre for raw minerals, slate teal for lithium and water, sage for recycling, and clay red for the bottleneck. There are no scene numbers and no corner titles.

## Where your footage goes (1920×1080 pixels)

| Scene | Presenter window | Left | Top | Width | Height |
|---|---|---|---|---|---|
| 1 (0:00–0:14) | Left side | 82 | 169 | 662 | 837 |
| 3 (0:30–0:48) | Left side | 76 | 203 | 486 | 756 |
| 5 (1:10–1:25) | Full screen | 0 | 0 | 1920 | 1080 |

Scenes 2 and 4 are full-screen animation with your voice over them. Scale and crop your clip so it fills each window, and put it on the track above the animation.

## Presenter script

Timed to the rendered animation, at about 160 words a minute. The full cue table with delivery notes is in `SCRIPT.md`.

## Scene storyboard

| Time | What happens on screen |
|---|---|
| 0:00 | The presenter window fades in. A phone appears, tilts and separates into five labelled layers. |
| 0:07 | Zoom into the battery. A counter climbs to **30+** with *minerals dug from the Earth* and lights up thirty element tiles. *So what happens when it dies?* |
| 0:14 | A rock cross-section with an open pit and a spodumene vein sits beside brine ponds that fill and turn from blue to yellow green under the sun. |
| 0:20 | Particles stream out of both sites and settle into **Li₂CO₃**. *Every charge begins with tonnes of rock and months of sunlight.* |
| 0:30 | A conveyor runs from Mining to Refining to Cathode making to Cell assembly. Refining is flagged as the bottleneck. |
| 0:38 | Ships cross a hand-drawn world map from Chile and Australia to China, then Korea, then Europe, and a truck carries the battery to *You*. |
| 0:43 | Battery health falls to 72%. The battery drops into a landfill. *Too many end up buried in landfills but they don't have to.* |
| 0:48 | A spent battery drops into the recycling reactor and is shredded. Lithium, cobalt and nickel come back out. |
| 0:56 | A green arc carries about 90% of the metal to a new battery factory and new batteries head back out. |
| 1:04 | *Up to 80% less energy than mining new metal.* *We can mine our past instead of the planet.* |
| 1:10 | Over the presenter: a straight Mine → Market → Recycle line bends into an infinity loop with light flowing around it. **From Mine to Market** *A sustainable and resilient future for our planet.* |

## Figures used

These are rounded, commonly cited figures. Check them against your sources before you submit.

- More than 30 critical minerals in a smartphone.
- Brine evaporation takes many months, typically about 12 to 18.
- About 65% of lithium chemical refining happens in China (IEA *Global Critical Minerals Outlook*).
- Batteries are usually retired at about 70 to 80% of their original capacity.
- Hydrometallurgical recovery of Li, Co and Ni can reach 90% or more. Some recyclers report up to about 80% lower energy use than mining new material.
