# From Mine to Market: A Battery's Second Life

An 85-second explainer animation in five scenes, written in [Manim Community](https://www.manim.community/).

- **Code:** `battery_second_life.py` (one class per scene)
- **Rendered video:** `renders/from_mine_to_market.mp4` (1920×1080, 30 fps, 1:25)

## Render

```bash
pip install manim            # needs Cairo + Pango dev libs (apt: libcairo2-dev libpango1.0-dev) and ffmpeg
./render.sh                  # renders all 5 scenes at 1080p30 and joins them into renders/from_mine_to_market.mp4
# quick preview of one scene (480p15):
manim -pql battery_second_life.py S3_Processing
```

You don't need LaTeX. All text is set in Inter through Pango. If Inter isn't installed, change `FONT` at the top of the file.

## Design system

| Token | Hex | Meaning |
|---|---|---|
| Background | `#111827` | — |
| Slate | `#64748B` / `#CBD5E1` | structure, neutral hardware |
| Electric blue | `#38BDF8` | lithium, energy, finished cells |
| Emerald green | `#10B981` | recycling, sustainability |
| Warm orange | `#F97316` | raw minerals, mining |
| Amber / red | `#FBBF24` / `#EF4444` | caution (end of life) / bottleneck |

Every scene shares the same layout. A numbered header sits top-left, a series tag sits top-right, and each scene ends with a 0.6 s fade to the background colour. That makes the cuts between scenes seamless.

## Keyframe storyboard

### Scene 1 · 0:00–0:14 · Smartphone Mineral Breakdown (`S1_Smartphone`)
| t (s) | Keyframe |
|---|---|
| 0.0 | Header fades in. A flat phone front view (glass, OLED app grid, "9:41") rises in. |
| 1.8 | The phone shears into an isometric view (rotate 50°, y × 0.55), and darker edges give each layer some thickness. |
| 3.1 | Exploded view. Housing, battery, logic board, display and glass separate vertically by 0.72 units. |
| 4.5 | Leader lines draw to labels, each showing a part and its minerals (e.g. *Li-ion battery — Li · Co · Ni · Mn · C*). |
| 6.8 | Zoom in on the battery. The other layers drop away, and the battery morphs into a large cell cutaway with cathode, separator and anode stripes. |
| 8.1 | A counter ticks from **1 → 30**. Thirty element tiles light up in sync, with the battery minerals in blue. |
| 11.5 | **"+"** pops with a flash, giving **30+ Critical Minerals**. The battery-critical tiles brighten. |

### Scene 2 · 0:14–0:30 · Extraction & Geology (`S2_Extraction`)
| t (s) | Keyframe |
|---|---|
| 0.0 | Split screen with a dashed divider. Left: hard-rock mining (spodumene in pegmatite). Right: brine and salar ponds. |
| 0.9 | Left: a strata cross-section with an orange pegmatite vein, spodumene crystals and a stepped open pit. Right: a salt crust, a brine aquifer and three ponds. |
| 2.2 | A well pipe draws in. The ponds fill one after another, blue → teal → green (rising concentration), with arrows between ponds. The sun rotates. |
| 4.4 | Evaporation wisps rise and the levels drop. A note reads *Li concentration ↑ · ~12–18 months*. |
| 5.9 | Process boxes appear: *Crush → Roast → Leach* and *Concentrate → Purify*. |
| 6.4 | 2 × 28 particles fly out of the pit and pond 3, through the process boxes, and converge at the centre. Their colour shifts from orange or green to blue. |
| 10.6 | **Li₂CO₃** resolves from the particles, followed by *Lithium carbonate · battery-grade*. |

### Scene 3 · 0:30–0:48 · Processing & Bottlenecks (`S3_Processing`)
| t (s) | Keyframe |
|---|---|
| 0.8 | A conveyor belt with moving chevrons. Four station badges: Extraction → Refining → Cathode synthesis → Battery assembly. |
| 2.1 | An ore chunk travels the belt and morphs at each station: rock → white Li₂CO₃ powder → blue cathode hexagon → battery cell. |
| 3.6 | Refining turns red, with a pulsing ring, a **BOTTLENECK** tag and *~65% of global Li refining in one country*. |
| 7.5 | A dot-matrix world map with pins: Chile and Australia (mining), China (refining), Korea/Japan (cathodes), Europe (assembly), Market. |
| 9.1 | Ships trace the routes (Chile/Australia → China → Korea → Europe), then a truck makes the final leg to Market. |
| 13.5 | The map dims. A state-of-health panel counts down 100% → 72% as the cell turns blue → amber. |
| 15.1 | A ⚠ caution icon pops: **End of life · below ~80% capacity**. |

### Scene 4 · 0:48–1:10 · Circular Recycling Loop (`S4_Recycling`)
| t (s) | Keyframe |
|---|---|
| 0.8 | A hydrometallurgical reactor builds in, with a hopper, viewing window, teal leach liquor, rising bubbles and a turning stirrer. |
| 1.8 | The spent amber battery (with ⚠) slides in and drops into the hopper. |
| 3.7 | Shredded fragments swirl inside the window. The caption reads *Shred → Leach → Separate*. |
| 5.7 | Element nodes **Li**, **Co** and **Ni** emerge from the outlet pipe onto the left of the loop. The Battery Manufacturing node appears on the right. |
| 7.5 | A glowing green arc sweeps from the nodes over the top to Manufacturing, led by a bright dot. The pill reads **≈90% of material recovered**. |
| 9.9 | 27 material units (9 per element) stream along the arc into Manufacturing. 3 units (1 per element) drop into a *~10% residue* bin. |
| 13.3 | Manufacturing flashes. A slate return arc draws: *New cells → market → end-of-life*. Pulses then circulate the whole loop continuously. |
| 15.4 | The banner pops: **⚡ Up to 80% Energy Reduction** *vs. producing from virgin ore*. |

### Scene 5 · 1:10–1:25 · Closed-Loop Supply Chain (`S5_ClosedLoop`)
| t (s) | Keyframe |
|---|---|
| 0.8 | The linear pipeline Mine → Market → Recycle appears as one line graded orange → blue → green. The caption reads *Linear: extract → use → discard*. |
| 3.4 | The line morphs into a lemniscate (∞). Mine sits at the left lobe, Market at the crossing and Recycle at the right lobe. The caption becomes *Closed loop: materials stay in circulation*. |
| 5.4 | Direction arrowheads and a soft glow fade in. 16 light pulses run endlessly around the ∞, coloured by position. |
| 7.6 | The header fades and the title rises: **From Mine to Market:** / **Sustainable & Resilient.** The final frame holds with the pulses still moving. |

## Figures used on screen

These are rounded figures, suitable for a general-audience explainer. Check them against your preferred sources before you publish.

- "30+ critical minerals in a smartphone" (prompt brief).
- Brine evaporation takes about 12–18 months.
- About 65% of lithium chemical refining happens in China (IEA, *Global Critical Minerals Outlook*).
- EV and phone cells are typically retired at about 70–80% of their original capacity.
- Hydrometallurgical recovery of Li, Co and Ni can reach ≥90%. Some recyclers report up to about 80% lower energy use than virgin supply.
