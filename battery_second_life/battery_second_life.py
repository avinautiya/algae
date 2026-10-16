"""From Mine to Market: A Battery's Second Life — 5-scene Manim explainer.

Render all scenes at 1080p30 and stitch them with ./render.sh, or one scene:
    manim -r 1920,1080 --fps 30 battery_second_life.py S1_Smartphone
Uses only Text (Pango) — no LaTeX install required.
"""
from manim import *
import numpy as np

# ---------------------------------------------------------------- palette
BG = "#111827"
INK = "#F1F5F9"
MUTED = "#94A3B8"
SLATE = "#64748B"
SLATE_M = "#334155"
SLATE_D = "#1F2937"
SLATE_L = "#CBD5E1"
BLUE = "#38BDF8"      # lithium / energy
BLUE_D = "#0C4A6E"
GREEN = "#10B981"     # recycling / sustainability
ORANGE = "#F97316"    # raw minerals / mining
ORANGE_D = "#7C2D12"
AMBER = "#FBBF24"     # caution
RED = "#EF4444"       # bottleneck
INDIGO = "#818CF8"    # cobalt
FONT = "Inter"

config.background_color = BG


# ---------------------------------------------------------------- helpers
def T(s, size=24, color=INK, weight=NORMAL):
    return Text(s, font=FONT, font_size=size, color=color, weight=weight)


def header(num, title, color):
    n = T(f"{num:02d}", 24, color, BOLD)
    bar = Line(ORIGIN, UP * 0.34, stroke_width=3, color=color)
    t = T(title, 24, INK, MEDIUM)
    return VGroup(n, bar, t).arrange(RIGHT, buff=0.22).to_corner(UL, buff=0.45)


def chem(parts, size=72, colors=None):
    """Chemical formula with real subscripts, e.g. [("Li","2"),("CO","3")]."""
    g, x = VGroup(), 0.0
    for i, (base, sub) in enumerate(parts):
        c = colors[i] if colors else INK
        b = T(base, size, c, BOLD).move_to([x, 0, 0], aligned_edge=LEFT)
        g.add(b)
        x = b.get_right()[0] + 0.03
        if sub:
            s = T(sub, size * 0.58, c, BOLD)
            s.move_to([x, b.get_bottom()[1] - 0.08, 0], aligned_edge=DL)
            g.add(s)
            x = s.get_right()[0] + 0.05
    return g


def glow(path, color, width=4, layers=5):
    g = VGroup()
    for i in range(layers, 0, -1):
        g.add(path.copy().set_fill(opacity=0).set_stroke(color, width + i * 4, opacity=0.07))
    g.add(path.copy().set_fill(opacity=0).set_stroke(color, width, opacity=1))
    return g


def glow_dot(color, r=0.07):
    g = VGroup(*[Circle(radius=r * k, stroke_width=0).set_fill(color, o)
                 for k, o in [(3.2, 0.08), (2.2, 0.16), (1.5, 0.32)]])
    g.add(Circle(radius=r, stroke_width=0).set_fill(color, 1))
    return g


def tip_at(path, color, size=0.24):
    """Filled arrowhead at the end of an arbitrary path."""
    p1 = path.point_from_proportion(1)
    p0 = path.point_from_proportion(0.98)
    d = normalize(p1 - p0)
    tri = Triangle().set_fill(color, 1).set_stroke(width=0).scale_to_fit_width(size)
    tri.rotate(angle_of_vector(d) - PI / 2)
    return tri.move_to(p1 - d * size * 0.3)


def flow(dots, paths, c0, c1, run_time, lag):
    """Move each dot along its path while blending its color c0 -> c1."""
    anims = []
    for d, p in zip(dots, paths):
        def upd(m, a, p=p):
            m.move_to(p.point_from_proportion(smooth(a)))
            m.set_color(interpolate_color(ManimColor(c0), ManimColor(c1), a))
        anims.append(UpdateFromAlphaFunc(d, upd))
    return LaggedStart(*anims, lag_ratio=lag, run_time=run_time)


# ---------------------------------------------------------------- icons
def battery_icon(w=1.0, h=0.5, color=BLUE, level=1.0, stroke=3):
    body = RoundedRectangle(width=w, height=h, corner_radius=h * 0.18,
                            stroke_color=color, stroke_width=stroke)
    nub = RoundedRectangle(width=w * 0.07, height=h * 0.4, corner_radius=0.02,
                           stroke_width=0).set_fill(color, 1).next_to(body, RIGHT, buff=0.02)
    pad = h * 0.15
    fw = max((w - 2 * pad) * level, 0.001)
    fill = Rectangle(width=fw, height=h - 2 * pad, stroke_width=0).set_fill(color, 0.9)
    fill.move_to(body.get_left() + RIGHT * (pad + fw / 2))
    return VGroup(body, nub, fill)


def caution(size=0.6):
    tri = Triangle().set_stroke(width=0).set_fill(AMBER, 1).scale_to_fit_width(size)
    tri.round_corners(size * 0.08)
    mark = T("!", 40, BG, BOLD).scale_to_fit_height(size * 0.42)
    mark.move_to(tri.get_center() + DOWN * size * 0.08)
    return VGroup(tri, mark)


def pickaxe(color=ORANGE):
    head = Arc(radius=0.45, start_angle=PI / 2 - 1.0, angle=2.0, stroke_width=7, color=color)
    handle = Line(UP * 0.42, DOWN * 0.5, stroke_width=7, color=SLATE_L)
    return VGroup(handle, head).rotate(-PI / 5).move_to(ORIGIN)


def flask(color=SLATE_L, liquid=BLUE):
    body = Polygon([-0.12, 0.45, 0], [0.12, 0.45, 0], [0.12, 0.1, 0], [0.4, -0.4, 0],
                   [-0.4, -0.4, 0], [-0.12, 0.1, 0], stroke_color=color, stroke_width=4)
    liq = Polygon([-0.22, -0.08, 0], [0.22, -0.08, 0], [0.37, -0.36, 0], [-0.37, -0.36, 0],
                  stroke_width=0).set_fill(liquid, 0.85)
    rim = Line([-0.2, 0.45, 0], [0.2, 0.45, 0], stroke_width=4, color=color)
    return VGroup(liq, body, rim).move_to(ORIGIN)


def lattice(color=BLUE):
    hexa = RegularPolygon(6, radius=0.42, stroke_color=color, stroke_width=4)
    pts = hexa.get_vertices()
    spokes = VGroup(*[Line(ORIGIN, p, stroke_width=2, color=color, stroke_opacity=0.6) for p in pts])
    dots = VGroup(*[Dot(p, radius=0.06, color=color) for p in pts], Dot(ORIGIN, radius=0.08, color=INK))
    return VGroup(spokes, hexa, dots)


def truck(color=SLATE_L):
    cargo = RoundedRectangle(width=0.7, height=0.38, corner_radius=0.04, stroke_width=0).set_fill(color, 1)
    cab = RoundedRectangle(width=0.26, height=0.28, corner_radius=0.05, stroke_width=0).set_fill(color, 1)
    cab.next_to(cargo, RIGHT, buff=0.03).align_to(cargo, DOWN)
    win = Rectangle(width=0.1, height=0.09, stroke_width=0).set_fill(BG, 1).move_to(cab.get_center() + UR * 0.05)
    wheels = VGroup(*[Circle(radius=0.07).set_fill(BG, 1).set_stroke(color, 2.5).move_to([x, cargo.get_bottom()[1], 0])
                      for x in (cargo.get_left()[0] + 0.15, cargo.get_right()[0] - 0.12, cab.get_center()[0])])
    return VGroup(cargo, cab, win, wheels).move_to(ORIGIN)


def ship(color=SLATE_L):
    hull = Polygon([-0.5, 0, 0], [0.5, 0, 0], [0.38, -0.2, 0], [-0.42, -0.2, 0],
                   stroke_width=0).set_fill(color, 1)
    boxes = VGroup(*[Rectangle(width=0.16, height=0.12).set_fill(c, 1).set_stroke(BG, 1.5)
                     for c in (ORANGE, BLUE, GREEN, SLATE)]).arrange(RIGHT, buff=0.02)
    boxes.next_to(hull, UP, buff=0).shift(RIGHT * 0.08)
    bridge = Rectangle(width=0.12, height=0.24, stroke_width=0).set_fill(color, 1)
    bridge.next_to(hull, UP, buff=0).align_to(hull, LEFT).shift(RIGHT * 0.06)
    return VGroup(hull, boxes, bridge).move_to(ORIGIN)


def sun(color=AMBER):
    c = Circle(radius=0.28, stroke_width=0).set_fill(color, 1)
    rays = VGroup(*[Line(RIGHT * 0.4, RIGHT * 0.56, stroke_width=4, color=color).rotate(a, about_point=ORIGIN)
                    for a in np.arange(0, TAU, TAU / 8)])
    return VGroup(c, rays)


def factory(color=SLATE_L):
    pts = [(-0.5, -0.3), (0.5, -0.3), (0.5, 0.12), (0.2, -0.04), (0.2, 0.12), (-0.1, -0.04),
           (-0.1, 0.12), (-0.32, 0.0), (-0.32, 0.42), (-0.46, 0.42), (-0.46, -0.0), (-0.5, 0.0)]
    body = Polygon(*[[x, y, 0] for x, y in pts], stroke_width=0).set_fill(color, 1)
    wins = VGroup(*[Square(0.09, stroke_width=0).set_fill(BG, 1).move_to([x, -0.15, 0]) for x in (-0.25, 0.0, 0.25)])
    return VGroup(body, wins).move_to(ORIGIN)


def recycle_icon(color=GREEN, r=0.3):
    g = VGroup()
    for k in range(3):
        a = Arc(radius=r, start_angle=k * TAU / 3 + 0.3, angle=TAU / 3 - 0.6, stroke_width=5, color=color)
        a.add_tip(tip_length=0.14, tip_width=0.14)
        g.add(a)
    return g


def bolt(color=GREEN, h=0.4):
    p = Polygon([0.05, 0.5, 0], [-0.25, -0.05, 0], [-0.02, -0.05, 0], [-0.08, -0.5, 0],
                [0.25, 0.08, 0], [0.02, 0.08, 0], stroke_width=0).set_fill(color, 1)
    return p.scale_to_fit_height(h)


def pin(color):
    ring = Circle(radius=0.12).set_fill(BG, 1).set_stroke(color, 3)
    return VGroup(ring, Dot(radius=0.05, color=color))


def pill(text, color, size=14, fill_op=0.14):
    t = T(text, size, color, BOLD)
    box = RoundedRectangle(width=t.width + 0.36, height=t.height + 0.2, corner_radius=0.12)
    box.set_fill(color, fill_op).set_stroke(color, 1.5).move_to(t)
    return VGroup(box, t)


# ---------------------------------------------------------------- base scene
class Base(Scene):
    DUR = 10
    FADE = 0.6

    def add_chrome(self, num, title, color):
        self.hdr = header(num, title, color)
        self.tag = T("FROM MINE TO MARKET  ·  A BATTERY'S SECOND LIFE", 13, MUTED, SEMIBOLD)
        self.tag.to_edge(RIGHT, buff=0.45).match_y(self.hdr)
        self.play(FadeIn(self.hdr, shift=RIGHT * 0.2), FadeIn(self.tag), run_time=0.8)

    def finish(self, fade=True):
        rest = self.DUR - self.renderer.time - (self.FADE if fade else 0)
        if rest > 0:
            self.wait(rest)
        if fade:
            self.play(*[FadeOut(m) for m in self.mobjects], run_time=self.FADE)


# ================================================================ SCENE 1
class S1_Smartphone(Base):
    """0:00-0:14  Exploded smartphone -> battery -> 30+ critical minerals."""
    DUR = 14

    def construct(self):
        self.add_chrome(1, "Smartphone Mineral Breakdown", BLUE)

        W, H = 2.0, 4.0

        def slab(fill, op, stroke, sw=2, w=W, h=H, cr=0.3):
            return RoundedRectangle(width=w, height=h, corner_radius=cr).set_fill(fill, op).set_stroke(stroke, sw)

        housing = VGroup(slab(SLATE_M, 1, SLATE), slab(SLATE_D, 1, SLATE, 1, W - 0.18, H - 0.18, 0.22))
        battery = VGroup(
            RoundedRectangle(width=1.55, height=2.3, corner_radius=0.12).set_fill(BLUE_D, 1).set_stroke(BLUE, 2.5),
            T("Li-ion", 22, BLUE, BOLD),
        ).shift(DOWN * 0.6)
        board_main = RoundedRectangle(width=1.75, height=1.05, corner_radius=0.1).set_fill("#172033", 1) \
            .set_stroke(SLATE_L, 1.5).move_to(UP * 1.3)
        soc = Square(0.38).set_fill(SLATE_M, 1).set_stroke(SLATE_L, 1).move_to([0.35, 1.2, 0])
        soc_t = T("SoC", 11, SLATE_L, BOLD).move_to(soc)
        mem = Rectangle(width=0.32, height=0.18).set_fill(SLATE_M, 1).set_stroke(SLATE_L, 1).move_to([0.35, 1.6, 0])
        cams = VGroup(*[Circle(radius=0.13).set_fill(BG, 1).set_stroke(SLATE_L, 1.5).move_to([-0.5, y, 0])
                        for y in (1.52, 1.1)])
        sub = RoundedRectangle(width=1.1, height=0.25, corner_radius=0.06).set_fill("#172033", 1) \
            .set_stroke(SLATE_L, 1.5).move_to(DOWN * 1.82)
        flex = Line([0.84, 0.78, 0], [0.84, -1.7, 0], stroke_width=3, color=SLATE)
        board = VGroup(board_main, soc, soc_t, mem, cams, sub, flex)

        app_cols = [BLUE, GREEN, ORANGE, SLATE, INDIGO]
        icons = VGroup(*[RoundedRectangle(width=0.28, height=0.28, corner_radius=0.07, stroke_width=0)
                         .set_fill(app_cols[(i * 3) % 5], 0.85) for i in range(20)]) \
            .arrange_in_grid(5, 4, buff=0.14).move_to(UP * 0.15)
        dock = VGroup(*[RoundedRectangle(width=0.28, height=0.28, corner_radius=0.07, stroke_width=0)
                        .set_fill(c, 0.85) for c in (GREEN, BLUE, SLATE, ORANGE)]) \
            .arrange(RIGHT, buff=0.14).move_to(DOWN * 1.55)
        clock = T("9:41", 15, INK, SEMIBOLD).move_to(UP * 1.5)
        display = VGroup(slab("#0B1324", 1, SLATE, 1.5), icons, dock, clock)
        glass = VGroup(slab(WHITE, 0.04, SLATE_L, 2),
                       RoundedRectangle(width=0.55, height=0.14, corner_radius=0.07, stroke_width=0)
                       .set_fill(BG, 1).move_to(UP * 1.78))

        layers = [housing, battery, board, display, glass]
        phone = VGroup(*layers)
        self.play(FadeIn(phone, shift=UP * 0.3), run_time=1.0)

        # -- flat -> isometric (rotate 50 deg, squash vertically)
        th = 50 * DEGREES
        M = np.diag([1, 0.55]) @ np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        base_c = np.array([-1.7, -0.3, 0])
        targets = [L.copy().apply_matrix(M, about_point=ORIGIN).shift(base_c + UP * 0.05 * i)
                   for i, L in enumerate(layers)]
        self.play(*[Transform(L, t) for L, t in zip(layers, targets)], run_time=1.3)

        edge_cols = ["#0F172A", "#082F49", "#0B1220", "#060B16", SLATE_M]
        edges = [L[0].copy().set_fill(c, 1).set_stroke(width=0).shift(DOWN * 0.07)
                 for L, c in zip(layers, edge_cols)]
        self.remove(*layers)
        stacks = [VGroup(e, L) for e, L in zip(edges, layers)]
        self.add(*stacks)

        # -- explode
        gap = 0.72
        shifts = [UP * ((i - 2) * gap - 0.05 * i) for i in range(5)]
        self.play(LaggedStart(*[s.animate.shift(v) for s, v in zip(stacks, shifts)], lag_ratio=0.08),
                  run_time=1.4)

        names = [("Housing", "Al · Mg"), ("Li-ion battery", "Li · Co · Ni · Mn · C"),
                 ("Logic board", "Au · Ag · Cu · Ta · Pd"), ("OLED display", "In · Sn · Ga"),
                 ("Cover glass", "Si · Al · K")]
        labels = VGroup()
        for i, (s, (name, sub_t)) in enumerate(zip(stacks, names)):
            pts = s[1][0].get_points()
            anchor = pts[np.argmax(pts[:, 0])] + RIGHT * 0.1
            label_y = base_c[1] + (i - 2) * 1.05 - 0.25   # evenly spaced, independent of layer shape
            col = BLUE if i == 1 else SLATE_L
            dot = Dot(anchor, radius=0.04, color=col)
            ln = VMobject().set_points_as_corners([anchor, [1.05, anchor[1], 0], [1.35, label_y, 0], [1.5, label_y, 0]]) \
                .set_stroke(col, 1.5)
            txt = VGroup(T(name, 20, INK if i != 1 else BLUE, SEMIBOLD), T(sub_t, 15, MUTED)) \
                .arrange(DOWN, aligned_edge=LEFT, buff=0.06)
            txt.next_to(ln.get_end(), RIGHT, buff=0.15)
            labels.add(VGroup(dot, ln, txt))
        self.play(LaggedStart(*[Create(l[1]) for l in labels], *[FadeIn(l[0]) for l in labels],
                              *[FadeIn(l[2], shift=RIGHT * 0.15) for l in labels], lag_ratio=0.04),
                  run_time=1.5)
        self.wait(0.8)

        # -- zoom into the battery layer
        body = RoundedRectangle(width=2.5, height=3.5, corner_radius=0.2).set_fill(BLUE_D, 1).set_stroke(BLUE, 3)
        tab = RoundedRectangle(width=0.6, height=0.18, corner_radius=0.05, stroke_width=0).set_fill(BLUE, 1) \
            .next_to(body, UP, buff=0)
        stripes = VGroup(*[RoundedRectangle(width=2.0, height=0.13, corner_radius=0.06, stroke_width=0)
                           .set_fill(c, 0.9) for c in [BLUE, SLATE_L, SLATE] * 4]).arrange(DOWN, buff=0.07)
        stripes.move_to(body).shift(DOWN * 0.25)
        cell_name = T("Li-ion cell", 22, INK, SEMIBOLD).next_to(stripes, UP, buff=0.22)
        legend = VGroup(*[VGroup(Square(0.16, stroke_width=0).set_fill(c, 1), T(s, 14, MUTED)).arrange(RIGHT, buff=0.1)
                          for c, s in [(BLUE, "Cathode"), (SLATE_L, "Separator"), (SLATE, "Anode")]]) \
            .arrange(RIGHT, buff=0.3).next_to(body, DOWN, buff=0.3)
        big = VGroup(body, tab, stripes, cell_name, legend).move_to([-3.7, -0.35, 0])

        others = [s for i, s in enumerate(stacks) if i != 1]
        self.play(*[FadeOut(o, shift=DOWN * 0.2) for o in others], FadeOut(labels),
                  ReplacementTransform(stacks[1], big), run_time=1.3)

        # -- counter 1 -> 30+
        v = ValueTracker(1)
        num_anchor = np.array([3.05, 1.4, 0])
        num = always_redraw(lambda: T(str(int(round(v.get_value()))), 110, BLUE, BOLD)
                            .move_to(num_anchor, aligned_edge=RIGHT))
        plus = T("+", 110, BLUE, BOLD).next_to(num_anchor, RIGHT, buff=0.06)
        caption = T("Critical Minerals", 36, INK, SEMIBOLD).move_to([2.6, 0.37, 0])

        syms = ["Li", "Co", "Ni", "Mn", "C", "Al", "Cu", "Au", "Ag", "Pd", "Sn", "In", "Ga", "Ta", "W",
                "Nd", "Pr", "Dy", "Tb", "Gd", "La", "Ce", "Eu", "Y", "Sb", "Bi", "Mg", "Si", "Ge", "Zn"]
        tiles = VGroup()
        for i, s in enumerate(syms):
            c = BLUE if i < 5 else SLATE_L
            sq = RoundedRectangle(width=0.46, height=0.46, corner_radius=0.06).set_fill(c, 0.14).set_stroke(c, 1.5)
            tiles.add(VGroup(sq, T(s, 16, c, BOLD).move_to(sq)))
        tiles.arrange_in_grid(3, 10, buff=0.08).move_to([2.6, -1.1, 0])
        note = T("Battery-critical:  Li · Co · Ni · Mn · Graphite", 18, BLUE).next_to(tiles, DOWN, buff=0.3)

        self.add(num)
        self.play(FadeIn(caption, shift=UP * 0.1), run_time=0.4)
        self.play(v.animate.set_value(30),
                  LaggedStart(*[FadeIn(t, scale=0.6) for t in tiles], lag_ratio=1.0),
                  run_time=3.0, rate_func=linear)
        self.play(FadeIn(plus, scale=1.6), Flash(plus, color=BLUE, line_length=0.25, flash_radius=0.55),
                  run_time=0.45)
        self.play(FadeIn(note, shift=UP * 0.1),
                  *[t[0].animate.set_fill(BLUE, 0.35) for t in tiles[:5]], run_time=0.6)
        self.finish()


# ================================================================ SCENE 2
class S2_Extraction(Base):
    """0:14-0:30  Hard-rock spodumene vs. salar brine -> Li2CO3."""
    DUR = 16

    def construct(self):
        self.add_chrome(2, "Extraction & Geology", ORANGE)
        divider = DashedLine([0, 2.85, 0], [0, -1.45, 0], dash_length=0.08, stroke_width=1.5, color=SLATE_M)
        lt = VGroup(T("Hard-rock mining", 24, INK, SEMIBOLD), T("Spodumene in pegmatite", 17, ORANGE)) \
            .arrange(DOWN, buff=0.1).move_to([-3.6, 2.45, 0])
        rt = VGroup(T("Brine extraction", 24, INK, SEMIBOLD), T("Salar evaporation ponds", 17, BLUE)) \
            .arrange(DOWN, buff=0.1).move_to([3.4, 2.45, 0])

        # ---- left: geological cross-section
        x0, x1 = -6.6, -0.6

        def band(yl0, yr0, yl1, yr1, col):
            return Polygon([x0, yl0, 0], [x1, yr0, 0], [x1, yr1, 0], [x0, yl1, 0], stroke_width=0).set_fill(col, 1)

        strata = VGroup(band(1.2, 1.2, 0.55, 0.7, "#2B3445"), band(0.55, 0.7, -0.35, -0.1, "#252D3C"),
                        band(-0.35, -0.1, -1.4, -1.4, "#1E2533"))
        vtop = lambda x: -0.9 + (x - x0) / 6 * 1.8
        vein = Polygon([x0, vtop(x0), 0], [x1, vtop(x1), 0], [x1, vtop(x1) - 0.5, 0], [x0, vtop(x0) - 0.5, 0],
                       stroke_width=0).set_fill(ORANGE_D, 1)
        rng = np.random.default_rng(3)
        crystals = VGroup()
        for _ in range(26):
            x = rng.uniform(x0 + 0.2, x1 - 0.2)
            y = vtop(x) - rng.uniform(0.1, 0.4)
            crystals.add(Rectangle(width=rng.uniform(0.12, 0.26), height=0.07, stroke_width=0)
                         .set_fill(ORANGE, 0.95).rotate(rng.uniform(0, PI)).move_to([x, y, 0]))
        pit = Polygon(*[[x, y, 0] for x, y in [(-5.4, 1.2), (-5.05, 0.75), (-4.75, 0.75), (-4.4, 0.3), (-4.1, 0.3),
                                               (-3.85, -0.1), (-2.85, -0.1), (-2.6, 0.3), (-2.3, 0.3), (-1.95, 0.75),
                                               (-1.65, 0.75), (-1.3, 1.2)]]).set_fill(BG, 1).set_stroke(SLATE_L, 1.5)
        surface = Line([x0, 1.2, 0], [x1, 1.2, 0], stroke_width=2, color=SLATE_L)
        pit_lbl = T("Open pit", 14, MUTED).move_to([-3.35, 0.95, 0])
        vein_lbl = T("Li-bearing pegmatite vein", 14, ORANGE).move_to([-4.9, -1.2, 0])
        left = VGroup(strata, vein, crystals, pit, surface, pit_lbl, vein_lbl)

        # ---- right: salar with evaporation ponds
        gx0, gx1 = 0.6, 6.6
        crust = Polygon([gx0, 0.6, 0], [gx1, 0.6, 0], [gx1, -0.3, 0], [gx0, -0.3, 0], stroke_width=0).set_fill("#2B3445", 1)
        aquifer = Polygon([gx0, -0.3, 0], [gx1, -0.3, 0], [gx1, -1.4, 0], [gx0, -1.4, 0], stroke_width=0).set_fill("#0E2A47", 1)
        aq_lbl = T("Lithium-rich brine aquifer", 14, BLUE).move_to([3.8, -0.95, 0])
        surface2 = Line([gx0, 0.6, 0], [gx1, 0.6, 0], stroke_width=2, color=SLATE_L)
        PX = [1.95, 3.6, 5.25]
        ponds = VGroup(*[Polygon([cx - 0.7, 0.6, 0], [cx - 0.55, 0.1, 0], [cx + 0.55, 0.1, 0], [cx + 0.7, 0.6, 0])
                         .set_fill(BG, 1).set_stroke(SLATE_L, 1.5) for cx in PX])
        pipe = VMobject().set_points_as_corners([[0.85, -0.85, 0], [0.85, 0.95, 0], [1.45, 0.95, 0], [1.45, 0.55, 0]]) \
            .set_stroke(SLATE_L, 4)
        right = VGroup(crust, aquifer, aq_lbl, surface2, ponds)

        levels = [ValueTracker(0.001) for _ in PX]
        pond_cols = ["#3B82F6", "#14B8A6", GREEN]

        def liquid(cx, lv, col):
            def make():
                y = 0.1 + 0.42 * lv.get_value()
                hw = 0.55 + (y - 0.1) / 0.5 * 0.15
                return Polygon([cx - 0.55, 0.1, 0], [cx + 0.55, 0.1, 0], [cx + hw, y, 0], [cx - hw, y, 0],
                               stroke_width=0).set_fill(col, 0.85)
            return always_redraw(make)

        liquids = VGroup(*[liquid(cx, lv, c) for cx, lv, c in zip(PX, levels, pond_cols)])
        transfer = VGroup(*[Arrow([a, 0.82, 0], [b, 0.82, 0], buff=0, stroke_width=2.5, color=SLATE_L,
                                  max_tip_length_to_length_ratio=0.35) for a, b in [(2.5, 3.05), (4.15, 4.7)]])
        sunny = sun().move_to([6.15, 1.95, 0])
        conc_arrow = Arrow([1.3, 1.3, 0], [5.9, 1.3, 0], buff=0, stroke_width=2, color=MUTED,
                           max_tip_length_to_length_ratio=0.04)
        conc_lbl = T("Li concentration ↑   ·   ~12–18 months", 15, MUTED).move_to([3.6, 1.58, 0])

        self.play(Create(divider), FadeIn(lt, shift=DOWN * 0.1), FadeIn(rt, shift=DOWN * 0.1), run_time=0.9)
        self.play(FadeIn(left, lag_ratio=0.05), FadeIn(right, lag_ratio=0.05), run_time=1.3)
        self.add(liquids)
        self.play(Create(pipe), FadeIn(sunny, scale=0.6), run_time=0.7)
        self.play(LaggedStart(*[lv.animate.set_value(1) for lv in levels], lag_ratio=0.45),
                  LaggedStart(*[GrowArrow(a) for a in transfer], lag_ratio=0.5),
                  Rotate(sunny[1], PI / 4), run_time=2.2)

        # evaporation wisps
        def wisp(x):
            return ParametricFunction(lambda t: np.array([0.05 * np.sin(t * 18), t, 0]), t_range=[0, 0.35],
                                      stroke_width=2, color=SLATE_L).move_to([x, 0.75, 0]).set_stroke(opacity=0.7)
        wisps = VGroup(*[wisp(cx + d) for cx in PX for d in (-0.25, 0.25)])
        self.play(FadeIn(wisps, shift=UP * 0.2), GrowArrow(conc_arrow), FadeIn(conc_lbl), run_time=0.6)
        self.play(wisps.animate.shift(UP * 0.45).set_opacity(0),
                  *[lv.animate.set_value(0.82) for lv in levels], run_time=0.9)

        # ---- processing boxes + particle flows
        def box(text, x):
            b = RoundedRectangle(width=2.5, height=0.55, corner_radius=0.12).set_fill(SLATE_D, 1).set_stroke(SLATE_L, 1.5)
            return VGroup(b, T(text, 15, SLATE_L, MEDIUM)).move_to([x, -2.05, 0])
        box_l = box("Crush → Roast → Leach", -2.3)
        box_r = box("Concentrate → Purify", 2.3)
        formula = chem([("Li", "2"), ("CO", "3")], 66, [BLUE, INK]).move_to([0, -2.75, 0])
        f_sub = T("Lithium carbonate  ·  battery-grade", 17, MUTED).next_to(formula, DOWN, buff=0.18)
        self.play(FadeIn(box_l, shift=UP * 0.1), FadeIn(box_r, shift=UP * 0.1), run_time=0.5)

        N = 28
        fl, fr = formula.get_left()[0], formula.get_right()[0]
        fy0, fy1 = formula.get_bottom()[1], formula.get_top()[1]

        def target():
            return np.array([rng.uniform(fl + 0.1, fr - 0.1), rng.uniform(fy0 + 0.1, fy1 - 0.1), 0])

        dl, pl, dr, pr = VGroup(), [], VGroup(), []
        for _ in range(N):
            s = np.array([rng.uniform(-3.7, -3.0), rng.uniform(-0.05, 0.2), 0])
            p = VMobject().set_points_smoothly([s, [s[0] + 0.3, 1.65, 0], [-1.0, 0.9, 0], [-1.4, -1.6, 0],
                                                [-1.7, -2.05, 0], target()])
            dl.add(glow_dot(ORANGE, 0.045).move_to(s))
            pl.append(p)
            s = np.array([rng.uniform(4.85, 5.65), rng.uniform(0.15, 0.4), 0])
            p = VMobject().set_points_smoothly([s, [4.9, 1.1, 0], [1.0, 0.95, 0], [1.4, -1.6, 0],
                                                [1.7, -2.05, 0], target()])
            dr.add(glow_dot(GREEN, 0.045).move_to(s))
            pr.append(p)
        self.add(dl, dr)
        self.play(flow(dl, pl, ORANGE, BLUE, 4.2, 0.03), flow(dr, pr, GREEN, BLUE, 4.2, 0.03))
        self.play(FadeIn(formula, scale=0.85), FadeOut(dl), FadeOut(dr), run_time=0.8)
        self.play(FadeIn(f_sub, shift=UP * 0.1),
                  Circumscribe(formula, color=BLUE, buff=0.15, time_width=0.6), run_time=1.1)
        self.finish()


# ================================================================ SCENE 3
class S3_Processing(Base):
    """0:30-0:48  Conveyor of value-chain stages, global routes, end-of-life."""
    DUR = 18

    def construct(self):
        self.add_chrome(3, "Processing & Bottlenecks", SLATE_L)
        BY = 0.95
        top = BY + 0.16
        belt = RoundedRectangle(width=12.4, height=0.32, corner_radius=0.16).set_fill(SLATE_D, 1) \
            .set_stroke(SLATE_L, 2).move_to([0, BY, 0])

        def chevron(x):
            return VMobject().set_points_as_corners([[x - 0.04, BY + 0.06, 0], [x + 0.04, BY, 0], [x - 0.04, BY - 0.06, 0]]) \
                .set_stroke(SLATE, 2)
        chevrons = VGroup(*[chevron(x) for x in np.arange(-5.9, 5.9, 0.45)])

        def chev_upd(m, dt):
            for c in m:
                c.shift(RIGHT * dt * 0.9)
                if c.get_x() > 5.95:
                    c.shift(LEFT * 0.45 * len(m))
        chevrons.add_updater(chev_upd)

        SX = [-4.5, -1.5, 1.5, 4.5]
        names = ["Extraction", "Refining", "Cathode synthesis", "Battery assembly"]
        cols = [ORANGE, SLATE_L, BLUE, BLUE]
        icons = [pickaxe(ORANGE).scale(0.62), flask(SLATE_L, BLUE).scale(0.62), lattice(BLUE).scale(0.62),
                 battery_icon(0.58, 0.3, BLUE)]
        badges, slabels = VGroup(), VGroup()
        for x, c, ic, n in zip(SX, cols, icons, names):
            circ = Circle(radius=0.42).set_fill(c, 0.1).set_stroke(c, 3).move_to([x, 2.15, 0])
            badges.add(VGroup(circ, ic.move_to(circ)))
            slabels.add(T(n, 18, INK, SEMIBOLD).move_to([x, 0.5, 0]))
        stage_arrows = VGroup(*[Arrow([a + 0.55, 2.15, 0], [b - 0.55, 2.15, 0], buff=0, stroke_width=2,
                                      color=SLATE, max_tip_length_to_length_ratio=0.1)
                                for a, b in zip(SX[:-1], SX[1:])])

        self.play(FadeIn(belt), FadeIn(chevrons),
                  LaggedStart(*[FadeIn(b, scale=0.7) for b in badges], lag_ratio=0.15),
                  LaggedStart(*[GrowArrow(a) for a in stage_arrows], lag_ratio=0.2),
                  LaggedStart(*[FadeIn(l, shift=UP * 0.1) for l in slabels], lag_ratio=0.15), run_time=1.3)

        rock = Polygon(*[[x, y, 0] for x, y in [(-0.2, -0.12), (0.05, -0.15), (0.22, -0.05), (0.15, 0.12),
                                                 (-0.05, 0.17), (-0.2, 0.05)]]).set_fill(ORANGE, 1).set_stroke(ORANGE_D, 2)
        powder = Polygon(*[[x, 0.2 * (1 - (x / 0.26) ** 2), 0] for x in np.linspace(-0.26, 0.26, 24)],
                         stroke_width=0).set_fill(INK, 0.95)
        cathode = RegularPolygon(6, radius=0.19).set_fill(BLUE, 0.9).set_stroke(INK, 1.5)
        cell = battery_icon(0.5, 0.28, BLUE)

        def at(m, x):
            return m.copy().move_to([x, top + m.height / 2 + 0.02, 0])

        item = at(rock, SX[0])
        self.play(FadeIn(item, shift=DOWN * 0.3), Indicate(badges[0], color=ORANGE, scale_factor=1.12), run_time=0.6)
        shapes = [powder, cathode, cell]
        for k in range(1, 4):
            self.play(item.animate.move_to([SX[k], item.get_y(), 0]), run_time=0.75)
            self.play(Transform(item, at(shapes[k - 1], SX[k])),
                      Indicate(badges[k], color=cols[k] if k != 1 else INK, scale_factor=1.12), run_time=0.5)
            if k == 1:
                ring = Circle(radius=0.55).set_stroke(RED, 3).move_to(badges[1])
                tag = pill("BOTTLENECK", RED, 13).move_to([-1.5, 2.98, 0])
                bn_note = T("~65% of global Li refining in one country", 14, RED).move_to([-1.5, 0.17, 0])
                self.play(Create(ring), FadeIn(tag, shift=DOWN * 0.1), FadeIn(bn_note),
                          badges[1][0].animate.set_stroke(RED), run_time=0.7)
                ring.add_updater(lambda m, dt: m.set_stroke(opacity=0.55 + 0.45 * np.sin(self.renderer.time * 6)))
        self.play(item.animate.move_to([6.4, item.get_y(), 0]).set_opacity(0), run_time=0.6)
        self.remove(item)

        # ---- global supply map
        frame = RoundedRectangle(width=12.4, height=2.95, corner_radius=0.2).set_stroke(SLATE_M, 1.2) \
            .move_to([0, -1.93, 0])
        def proj(lon, lat):
            return np.array([lon / 180 * 5.8, -1.95 + lat / 90 * 1.4, 0])

        # coarse continent outlines (lon, lat) -> dot-matrix world map
        continents = [
            [(-165, 65), (-140, 70), (-95, 72), (-75, 62), (-55, 52), (-80, 25), (-97, 18), (-105, 22), (-125, 40), (-130, 55)],
            [(-80, 10), (-60, 8), (-35, -7), (-40, -22), (-58, -38), (-68, -55), (-75, -45), (-71, -18), (-81, -5)],
            [(-10, 36), (-10, 44), (0, 50), (5, 58), (25, 71), (40, 68), (40, 45), (28, 41), (15, 38)],
            [(-17, 15), (-10, 35), (10, 37), (32, 31), (43, 12), (51, 11), (40, -15), (32, -28), (20, -35), (12, -18), (9, 4), (-8, 5)],
            [(40, 45), (40, 68), (70, 73), (110, 76), (140, 72), (180, 68), (160, 60), (140, 50), (130, 40), (122, 30),
             (120, 22), (108, 10), (100, 1), (95, 15), (80, 8), (72, 20), (60, 25), (48, 30), (35, 33)],
            [(114, -22), (130, -12), (142, -11), (153, -28), (147, -38), (135, -34), (115, -34)],
        ]

        def on_land(lon, lat):
            for poly in continents:
                inside, n = False, len(poly)
                for i in range(n):
                    (xa, ya), (xb, yb) = poly[i], poly[(i + 1) % n]
                    if (ya > lat) != (yb > lat) and lon < (xb - xa) * (lat - ya) / (yb - ya) + xa:
                        inside = not inside
                if inside:
                    return True
            return False

        grid = VGroup()
        for x in np.arange(-5.9, 5.95, 0.2):
            for y in np.arange(-3.25, -0.6, 0.2):
                lon, lat = x / 5.8 * 180, (y + 1.95) / 1.4 * 90
                grid.add(Dot([x, y, 0], radius=0.022 if on_land(lon, lat) else 0.012, color=SLATE)
                         .set_opacity(0.75 if on_land(lon, lat) else 0.18))

        P = {"chile": proj(-69, -24), "aus": proj(121, -28), "china": proj(110, 33),
             "korea": proj(129, 37), "eu": proj(9, 50), "mkt": proj(24, 44)}
        pin_cols = {"chile": ORANGE, "aus": ORANGE, "china": RED, "korea": BLUE, "eu": BLUE, "mkt": GREEN}
        pin_txt = {"chile": ("Chile · mining", DOWN), "aus": ("Australia · mining", DOWN),
                   "china": ("China · refining", UL), "korea": ("Korea / Japan · cathodes", DR),
                   "eu": ("Europe · assembly", UP), "mkt": ("Market", DR)}
        pins, plabels = VGroup(), VGroup()
        for k, p in P.items():
            pins.add(pin(pin_cols[k]).move_to(p))
            txt, d = pin_txt[k]
            plabels.add(T(txt, 13, pin_cols[k] if k != "mkt" else INK).next_to(p, d, buff=0.18))
        map_title = T("Global supply chain", 15, MUTED, SEMIBOLD).move_to([-4.9, -0.7, 0])
        self.play(Create(frame), FadeIn(grid, lag_ratio=0.002), FadeIn(map_title), run_time=0.8)
        self.play(LaggedStart(*[FadeIn(p, scale=0.5) for p in pins], lag_ratio=0.1),
                  LaggedStart(*[FadeIn(l) for l in plabels], lag_ratio=0.1), run_time=0.8)

        def leg(a, b, angle, color, vehicle, rt):
            path = ArcBetweenPoints(P[a], P[b], angle=angle) if angle else Line(P[a], P[b])
            path.set_stroke(color, 2.5)
            v = vehicle.move_to(P[a])
            return path, v, [Create(path), MoveAlongPath(v, path), ]

        legs = [("chile", "china", 0.6, ORANGE, ship().scale(0.42)),
                ("aus", "china", 0.5, ORANGE, ship().scale(0.42)),
                ("china", "korea", -0.9, SLATE_L, ship().scale(0.42)),
                ("korea", "eu", -0.5, BLUE, ship().scale(0.42).flip(UP)),
                ("eu", "mkt", 0, BLUE, truck().scale(0.4))]
        built = [leg(*l, rt=0) for l in legs]
        for idx, rt in [((0, 1), 1.3), ((2,), 0.6), ((3,), 1.0), ((4,), 0.6)]:
            for i in idx:
                self.add(built[i][1])
            self.play(*[a for i in idx for a in built[i][2]], run_time=rt, rate_func=smooth)
            self.play(*[FadeOut(built[i][1]) for i in idx], run_time=0.15)

        # ---- end of life
        map_all = VGroup(frame, grid, pins, plabels, map_title, *[b[0] for b in built])
        panel = RoundedRectangle(width=7.6, height=2.2, corner_radius=0.2).set_fill(BG, 0.96) \
            .set_stroke(SLATE_M, 1.5).move_to([0.3, -1.95, 0])
        soh = ValueTracker(1.0)

        def soh_col():
            return interpolate_color(ManimColor(BLUE), ManimColor(AMBER), np.clip((1 - soh.get_value()) / 0.28, 0, 1))
        bat = always_redraw(lambda: battery_icon(1.8, 0.85, soh_col(), soh.get_value(), 4).move_to([-1.9, -1.95, 0]))
        pct = always_redraw(lambda: T(f"{int(round(soh.get_value() * 100))}%", 44, soh_col(), BOLD)
                            .move_to([0.35, -1.78, 0]))
        pct_lbl = T("State of health", 15, MUTED).move_to([0.35, -2.35, 0])
        warn = caution(0.75).move_to([2.25, -1.62, 0])
        eol = VGroup(T("End of life", 20, AMBER, BOLD), T("below ~80% capacity", 14, MUTED)) \
            .arrange(DOWN, buff=0.08).next_to(warn, DOWN, buff=0.12)
        self.play(map_all.animate.fade(0.75), FadeIn(panel), FadeIn(pct_lbl), run_time=0.5)
        self.add(bat, pct)
        self.play(soh.animate.set_value(0.72), run_time=1.6, rate_func=linear)
        self.play(FadeIn(warn, scale=1.8), FadeIn(eol, shift=UP * 0.1), run_time=0.5)
        self.play(Wiggle(warn, scale_value=1.15, rotation_angle=0.05 * TAU), run_time=0.5)
        self.finish()


# ================================================================ SCENE 4
class S4_Recycling(Base):
    """0:48-1:10  Reactor -> Li/Co/Ni -> 90% back to manufacturing."""
    DUR = 22

    def construct(self):
        self.add_chrome(4, "Circular Recycling Loop", GREEN)
        clock = ValueTracker(0).add_updater(lambda m, dt: m.increment_value(dt))
        self.add(clock)

        # ---- reactor
        RX = -4.6
        body = RoundedRectangle(width=1.9, height=2.6, corner_radius=0.35).set_fill(SLATE_D, 1) \
            .set_stroke(SLATE_L, 2.5).move_to([RX, -0.4, 0])
        hopper = Polygon([RX - 0.6, 1.6, 0], [RX + 0.6, 1.6, 0], [RX + 0.25, 1.0, 0], [RX - 0.25, 1.0, 0]) \
            .set_fill(SLATE_M, 1).set_stroke(SLATE_L, 2)
        neck = Rectangle(width=0.5, height=0.12).set_fill(SLATE_M, 1).set_stroke(SLATE_L, 2).move_to([RX, 0.94, 0])
        WC = np.array([RX, -0.45, 0])
        window = RoundedRectangle(width=1.4, height=1.6, corner_radius=0.2).set_fill(BG, 1).set_stroke(SLATE, 1.5).move_to(WC)
        liq = Intersection(window, Rectangle(width=2, height=1.2).move_to(WC + DOWN * 0.5)) \
            .set_fill("#0F766E", 0.6).set_stroke(width=0)
        shaft = Line([RX, 0.9, 0], [RX, -1.0, 0], stroke_width=3, color=SLATE_L)
        blades = always_redraw(lambda: Line(LEFT, RIGHT, stroke_width=4, color=SLATE_L)
                               .scale(0.4 * max(abs(np.cos(clock.get_value() * 4)), 0.05)).move_to([RX, -1.0, 0]))
        legs = VGroup(*[Line([RX + d, -1.7, 0], [RX + d * 1.15, -2.05, 0], stroke_width=4, color=SLATE)
                        for d in (-0.65, 0.65)])
        rng = np.random.default_rng(7)
        bubbles = VGroup(*[Circle(radius=rng.uniform(0.025, 0.05)).set_stroke(INK, 1, opacity=0.6)
                           .move_to([RX + rng.uniform(-0.55, 0.55), rng.uniform(-1.15, 0.0), 0]) for _ in range(9)])

        def bub_upd(m, dt):
            for b in m:
                b.shift(UP * dt * 0.45)
                if b.get_y() > 0.05:
                    b.set_y(-1.15)
        bubbles.add_updater(bub_upd)
        pipe = Line([RX + 0.95, -0.2, 0], [-2.8, -0.2, 0], stroke_width=6, color=SLATE)
        r_lbl = T("Hydrometallurgical reactor", 16, MUTED).move_to([RX, -2.35, 0])
        reactor = VGroup(legs, pipe, body, hopper, neck, window, liq, shaft)

        self.play(FadeIn(reactor, lag_ratio=0.05), FadeIn(r_lbl), run_time=1.0)
        self.add(blades, bubbles)

        # ---- spent battery enters
        spent = VGroup(battery_icon(0.9, 0.45, AMBER, 0.25))
        spent.add(caution(0.34).next_to(spent[0], UR, buff=-0.08))
        spent.move_to([-6.6, 2.4, 0])
        self.play(FadeIn(spent, shift=RIGHT * 0.3), run_time=0.4)
        self.play(spent.animate.move_to([RX, 2.15, 0]), run_time=0.9)
        self.play(spent.animate.move_to([RX, 1.15, 0]).scale(0.3).set_opacity(0), run_time=0.6)
        self.remove(spent)

        shred_cols = [BLUE, INDIGO, SLATE_L, ORANGE, "#B87333", AMBER]
        shreds = VGroup(*[Square(rng.uniform(0.05, 0.1), stroke_width=0).set_fill(rng.choice(shred_cols), 0.95)
                          .rotate(rng.uniform(0, PI)).move_to(WC + rng.uniform(0.1, 0.55) *
                                                               np.array([np.cos(a := rng.uniform(0, TAU)), np.sin(a), 0]))
                          for _ in range(34)])
        shreds.add_updater(lambda m, dt: m.rotate(dt * 3.2, about_point=WC))
        steps = T("Shred  →  Leach  →  Separate", 15, SLATE_L).move_to([RX, -2.72, 0])
        self.play(FadeIn(shreds, scale=0.5), FadeIn(steps, shift=UP * 0.1), run_time=0.5)
        self.wait(1.5)

        # ---- element nodes on the loop
        E = lambda th: np.array([1.0 + 3.4 * np.cos(th), -0.2 + 2.2 * np.sin(th), 0])
        elems = [("Li", "Lithium", BLUE, 150), ("Co", "Cobalt", INDIGO, 180), ("Ni", "Nickel", SLATE_L, 210)]
        nodes = VGroup()
        for sym, name, col, deg in elems:
            c = Circle(radius=0.4).set_fill(BG, 1).set_stroke(col, 3)
            c2 = Circle(radius=0.4, stroke_width=0).set_fill(col, 0.18)
            n = VGroup(c, c2, T(sym, 24, col, BOLD), T(name, 13, MUTED).next_to(c, RIGHT, buff=0.12))
            n[2].move_to(c)
            nodes.add(n.shift(E(deg * DEGREES) - c.get_center()))
        pipe_end = np.array([-2.8, -0.2, 0])
        for n in nodes:
            n.save_state()
            n.scale(0.1).move_to(pipe_end).set_opacity(0)
        self.play(LaggedStart(*[Restore(n) for n in nodes], lag_ratio=0.25),
                  shreds.animate.scale(0.2).set_opacity(0), run_time=1.2)
        self.remove(shreds)

        MX = E(0)
        mfg = VGroup(Circle(radius=0.72).set_fill(SLATE_D, 1).set_stroke(GREEN, 3), factory(SLATE_L).scale(0.9)).move_to(MX)
        mfg_lbl = VGroup(T("Battery", 18, INK, SEMIBOLD), T("manufacturing", 18, INK, SEMIBOLD)) \
            .arrange(DOWN, aligned_edge=LEFT, buff=0.06).next_to(mfg, RIGHT, buff=0.2)
        self.play(FadeIn(mfg, scale=0.7), FadeIn(mfg_lbl, shift=LEFT * 0.1), run_time=0.6)

        # material units: 10 per element, ring around node
        units, unit_paths, losses, loss_paths = VGroup(), [], VGroup(), []
        bin_c = np.array([-1.7, -2.75, 0])
        for (sym, name, col, deg), n in zip(elems, nodes):
            cen = n[0].get_center()
            for j in range(10):
                p0 = cen
                d = glow_dot(col, 0.045).move_to(p0)
                if j < 9:
                    ths = np.linspace(deg, 26, 50) * DEGREES
                    path = VMobject().set_points_smoothly([p0] + [E(t) for t in ths])
                    units.add(d)
                    unit_paths.append(path)
                else:
                    path = VMobject().set_points_smoothly([p0, p0 + DOWN * 0.55 + RIGHT * 0.3,
                                                           bin_c + UP * 0.25 + RIGHT * rng.uniform(-0.2, 0.2)])
                    losses.add(d)
                    loss_paths.append(path)
        self.add(units, losses)
        self.bring_to_back(units, losses)
        self.play(*[Indicate(n[0], color=n[0].get_stroke_color(), scale_factor=1.12) for n in nodes], run_time=0.5)

        # ---- glowing loop arrow sweeps
        top = ParametricFunction(E, t_range=[218 * DEGREES, 22 * DEGREES, -0.01])
        top_glow = glow(top, GREEN, 5)
        top_tip = tip_at(top, GREEN, 0.3)
        lead = glow_dot(GREEN, 0.08).move_to(top.get_start())
        self.add(lead)
        self.bring_to_back(top_glow)
        self.play(Create(top_glow), MoveAlongPath(lead, top), run_time=1.5)
        rec = pill("≈90% of material recovered", GREEN, 15).move_to(E(PI / 2) + UP * 0.42)
        self.play(FadeIn(top_tip, scale=0.5), FadeOut(lead), FadeIn(rec, shift=DOWN * 0.1), run_time=0.5)

        bin_box = VGroup(RoundedRectangle(width=1.0, height=0.34, corner_radius=0.08).set_fill(SLATE_D, 1)
                         .set_stroke(SLATE, 1.2), T("~10% residue", 11, MUTED)).move_to(bin_c)
        bin_box[1].move_to(bin_box[0])
        self.play(FadeIn(bin_box), run_time=0.3)
        self.play(LaggedStart(*[MoveAlongPath(u, p, rate_func=smooth) for u, p in zip(units, unit_paths)],
                              lag_ratio=0.025, run_time=3.4),
                  LaggedStart(*[MoveAlongPath(u, p) for u, p in zip(losses, loss_paths)], lag_ratio=0.3, run_time=1.6))
        self.play(FadeOut(units, scale=0.3), FadeOut(losses), Indicate(mfg[0], color=GREEN, scale_factor=1.15),
                  Flash(mfg, color=GREEN, flash_radius=0.95, line_length=0.25), run_time=0.7)

        # ---- return arc: new cells -> market -> end-of-life
        bot = ParametricFunction(E, t_range=[338 * DEGREES, 236 * DEGREES, -0.01])
        bot_glow = glow(bot, SLATE_L, 3, 3).set_stroke(opacity=0.5)
        bot_glow[-1].set_stroke(SLATE_L, 3, opacity=0.9)
        bot_tip = tip_at(bot, SLATE_L, 0.26)
        newcell = battery_icon(0.42, 0.22, BLUE).move_to(bot.get_start())
        bot_lbl = T("New cells  →  market  →  end-of-life", 15, MUTED).move_to([1.2, -2.78, 0])
        self.bring_to_back(bot_glow)
        self.play(Create(bot_glow), MoveAlongPath(newcell, bot), FadeIn(bot_lbl), run_time=1.4)
        self.play(FadeIn(bot_tip, scale=0.5), FadeOut(newcell), run_time=0.3)

        # continuous pulses around the loop
        ph = ValueTracker(0).add_updater(lambda m, dt: m.increment_value(dt * 0.6))
        self.add(ph)
        pulses = VGroup(*[glow_dot(GREEN, 0.05) for _ in range(8)])

        def pul_upd(m):
            for i, d in enumerate(m):
                th = (-ph.get_value() + i * TAU / len(m)) % TAU
                d.move_to(E(th))
                d.set_color(GREEN if np.sin(th) > 0 else SLATE_L)
        pulses.add_updater(pul_upd)
        self.add(pulses)

        # ---- energy banner
        b_txt = T("Up to 80% Energy Reduction", 28, INK, BOLD)
        b_sub = T("vs. producing from virgin ore", 15, MUTED)
        b_row = VGroup(bolt(GREEN, 0.42), b_txt, b_sub).arrange(RIGHT, buff=0.25)
        b_box = RoundedRectangle(width=b_row.width + 0.6, height=0.72, corner_radius=0.18) \
            .set_fill(GREEN, 0.14).set_stroke(GREEN, 2)
        banner = VGroup(b_box, b_row).move_to([0.6, -3.38, 0])
        b_row.move_to(b_box)
        self.play(FadeIn(banner, scale=0.85, shift=UP * 0.2), run_time=0.6)
        self.play(Indicate(b_txt, color=GREEN, scale_factor=1.05), run_time=0.8)
        self.finish()


# ================================================================ SCENE 5
class S5_ClosedLoop(Base):
    """1:10-1:25  Linear pipeline morphs into an infinity loop."""
    DUR = 15

    def construct(self):
        self.add_chrome(5, "Closed-Loop Supply Chain", GREEN)
        A, CY = 4.6, -0.4
        lem = lambda t: np.array([A * np.cos(t) / (1 + np.sin(t) ** 2),
                                  A * np.sin(t) * np.cos(t) / (1 + np.sin(t) ** 2) + CY, 0])

        def grad(m, width=5):
            m.set_stroke([ORANGE, BLUE, GREEN], width).set_sheen_direction(RIGHT)
            return m

        line = grad(Line([-A, CY, 0], [A, CY, 0]))
        heads = VGroup(*[Triangle().set_fill(interpolate_color(ManimColor(c0), ManimColor(c1), 0.5), 1)
                         .set_stroke(width=0).scale_to_fit_width(0.26).rotate(-PI / 2).move_to([x, CY, 0])
                         for x, c0, c1 in [(-2.3, ORANGE, BLUE), (2.3, BLUE, GREEN)]])

        def node(icon, col):
            return VGroup(Circle(radius=0.58).set_fill(BG, 1).set_stroke(width=0),
                          Circle(radius=0.58).set_fill(col, 0.15).set_stroke(col, 4), icon)
        mine = node(pickaxe(ORANGE).scale(0.68), ORANGE).move_to([-A, CY, 0])
        market = node(battery_icon(0.6, 0.32, BLUE), BLUE).move_to([0, CY, 0])
        recyc = node(recycle_icon(GREEN, 0.28), GREEN).move_to([A, CY, 0])
        nl = VGroup(T("Mine", 22, ORANGE, SEMIBOLD).next_to(mine, DOWN, buff=0.2),
                    T("Market", 22, BLUE, SEMIBOLD).next_to(market, DOWN, buff=0.5),
                    T("Recycle", 22, GREEN, SEMIBOLD).next_to(recyc, DOWN, buff=0.2))
        cap = T("Linear:  extract  →  use  →  discard", 20, MUTED).move_to([0, 1.75, 0])

        self.play(Create(line), LaggedStart(*[FadeIn(n, scale=0.6) for n in (mine, market, recyc)], lag_ratio=0.25),
                  LaggedStart(*[FadeIn(l, shift=UP * 0.1) for l in nl], lag_ratio=0.25), run_time=1.5)
        self.play(FadeIn(heads), FadeIn(cap, shift=DOWN * 0.1), run_time=0.5)
        self.wait(0.6)

        # ---- morph into infinity
        curve = grad(ParametricFunction(lem, t_range=[0, TAU, 0.01]))
        curve_glow = VGroup(*[curve.copy().set_stroke(width=5 + i * 5, opacity=0.06) for i in range(4, 0, -1)])
        for g in curve_glow:
            g.set_sheen_direction(RIGHT)
        cap2 = T("Closed loop:  materials stay in circulation", 20, GREEN).move_to(cap)
        self.play(ReplacementTransform(line, curve), FadeOut(heads),
                  nl[0].animate.next_to(mine, LEFT, buff=0.2), nl[2].animate.next_to(recyc, RIGHT, buff=0.2),
                  ReplacementTransform(cap, cap2), run_time=2.0)
        self.bring_to_back(curve_glow)
        self.add_foreground_mobjects(mine, market, recyc)

        def lem_col(x):
            u = x / A
            return interpolate_color(ManimColor(ORANGE), ManimColor(BLUE), u + 1) if u < 0 else \
                interpolate_color(ManimColor(BLUE), ManimColor(GREEN), u)

        arrows = VGroup()
        for t in (1.25 * PI, 1.75 * PI, 0.25 * PI, 0.75 * PI):
            p, q = lem(t), lem(t + 0.02)
            tri = Triangle().set_fill(lem_col(p[0]), 1).set_stroke(width=0).scale_to_fit_width(0.26)
            arrows.add(tri.rotate(angle_of_vector(q - p) - PI / 2).move_to(p))

        ph = ValueTracker(0).add_updater(lambda m, dt: m.increment_value(dt * 0.7))
        self.add(ph)
        pulses = VGroup(*[glow_dot(BLUE, 0.06) for _ in range(16)])

        def pul_upd(m):
            for i, d in enumerate(m):
                p = lem(PI + ph.get_value() + i * TAU / len(m))
                d.move_to(p).set_color(lem_col(p[0]))
        pulses.add_updater(pul_upd)
        pul_upd(pulses)
        self.play(FadeIn(curve_glow), FadeIn(arrows, scale=0.5), FadeIn(pulses), run_time=0.8)
        self.wait(1.4)

        title = VGroup(T("From Mine to Market:", 46, INK, BOLD), T("Sustainable & Resilient.", 46, GREEN, BOLD)) \
            .arrange(DOWN, buff=0.18).move_to([0, 2.55, 0])
        self.play(FadeOut(self.hdr), FadeOut(self.tag), FadeOut(cap2),
                  LaggedStart(*[FadeIn(l, shift=UP * 0.2) for l in title], lag_ratio=0.35), run_time=1.5)
        self.finish(fade=False)
