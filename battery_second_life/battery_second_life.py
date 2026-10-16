"""From Mine to Market: A Battery's Second Life — 5-scene Manim explainer.

Earth-science styling on light paper with hand drawn trees. Animation plus a presenter on camera.
Render everything with ./render.sh, or one scene:
    manim -r 1920,1080 --fps 30 battery_second_life.py S1_Hook
Set PRESENTER_PLACEHOLDER=1 to draw a silhouette inside the presenter windows.
Render S5_Closing with --transparent to get an overlay for full-screen footage.
"""
import os

from manim import *
import numpy as np

# ---------------------------------------------------------------- palette
BG = "#F1EEE4"        # warm paper
BG_D = "#E5E5D7"      # panels
PAPER = "#FBF9F3"     # light fills
INK = "#2F3B33"       # deep forest text
MUTED = "#6F7B70"     # sage grey
STONE = "#9A9E92"
STONE_M = "#CDD0C2"
STONE_D = "#DFE0D2"
STONE_L = "#5E665C"   # outlines and line icons
WATER = "#4E8C93"     # lithium and water
WATER_D = "#C7DCD9"
LEAF = "#6B9455"      # recycling and renewal
OCHRE = "#C27A3F"     # raw minerals
RUST = "#A65A3A"      # pegmatite vein
SUN = "#D6A43C"       # sun and caution
CLAY = "#B5523B"      # bottleneck and landfill
COBALT = "#5F72B0"
NICKEL = "#7F887D"
SOIL = ["#CDB89D", "#BCA488", "#AA9175"]
LAND = "#C3CFB0"
OCEAN = "#DCE6E1"

SERIF = "Fraunces"
SANS = "Nunito Sans"
PLACEHOLDER = os.environ.get("PRESENTER_PLACEHOLDER", "0") == "1"

config.background_color = BG


# ---------------------------------------------------------------- helpers
def T(s, size=24, color=INK, weight=NORMAL):
    return Text(s, font=SANS, font_size=size, color=color, weight=weight)


def H(s, size=30, color=INK, weight=SEMIBOLD, italic=False):
    return Text(s, font=SERIF, font_size=size, color=color, weight=weight, slant=ITALIC if italic else NORMAL)


TREE_GREENS = ["#A9BA8F", "#9EB287", "#B3C199", "#94A97F"]
TREE_LINE = "#7C8E67"
TRUNK = "#A08A6E"


def blob(center, r, rng, k=0.07, squash=0.92):
    """Closed wobbly outline that reads as drawn by hand."""
    ph = rng.uniform(0, TAU, 3)
    pts = [center + r * (1 + k * np.sin(3 * t + ph[0]) + 0.6 * k * np.sin(5 * t + ph[1]) + 0.3 * k * np.sin(11 * t + ph[2]))
           * np.array([np.cos(t), squash * np.sin(t), 0]) for t in np.linspace(0, TAU, 90)]
    return VMobject().set_points_smoothly(pts)


def tree(x, y, h, rng, col=None, line=TREE_LINE):
    col = col or rng.choice(TREE_GREENS)
    th, tw = h * 0.36, h * 0.06
    trunk = Polygon([x - tw, y, 0], [x + tw, y, 0], [x + tw * 0.5, y + th, 0], [x - tw * 0.5, y + th, 0]) \
        .set_fill(TRUNK, 1).set_stroke(line, 1.2)
    twig = Line([x, y + th * 0.6, 0], [x + h * 0.12, y + th * 1.05, 0], stroke_width=2, color=TRUNK)
    r = h * 0.32
    c = np.array([x, y + th + r * 0.7, 0])
    canopy = Union(blob(c, r, rng), blob(c + [-r * 0.6, -r * 0.3, 0], r * 0.68, rng),
                   blob(c + [r * 0.6, -r * 0.28, 0], r * 0.7, rng)).set_fill(col, 1).set_stroke(line, 1.6)
    shine = blob(c + [-r * 0.3, r * 0.3, 0], r * 0.35, rng).set_fill(WHITE, 0.18).set_stroke(width=0)
    sketch = canopy.copy().set_fill(opacity=0).set_stroke(line, 1, opacity=0.35).shift([0.02, -0.015, 0]).rotate(0.02)
    return VGroup(trunk, twig, canopy, shine, sketch)


def pine(x, y, h, rng, col="#93A886", line=TREE_LINE):
    trunk = Rectangle(width=h * 0.07, height=h * 0.18).set_fill(TRUNK, 1).set_stroke(line, 1.2) \
        .move_to([x, y + h * 0.09, 0])
    g = VGroup(trunk)
    for w, y0, y1 in [(0.62, 0.14, 0.55), (0.5, 0.36, 0.78), (0.36, 0.58, 1.0)]:
        j = rng.uniform(-0.02, 0.02, 3) * h
        tri = Polygon([x - w * h / 2 + j[0], y + y0 * h, 0], [x + w * h / 2 + j[1], y + y0 * h, 0], [x + j[2], y + y1 * h, 0])
        g.add(tri.round_corners(h * 0.04).set_fill(col, 1).set_stroke(line, 1.5))
    return g


def landscape():
    """Soft hills dotted with hand drawn trees under a pale sky."""
    rng = np.random.default_rng(21)
    far = lambda x: -3.3 + 0.2 * np.sin(0.55 * x + 0.8) + 0.08 * np.sin(1.3 * x)
    near = lambda x: -3.78 + 0.12 * np.sin(0.7 * x + 2.0) + 0.05 * np.sin(1.9 * x + 1.0)

    def hill(f, col):
        xs = np.linspace(-7.4, 7.4, 120)
        body = Polygon(*[[x, f(x), 0] for x in xs], [7.4, -4.3, 0], [-7.4, -4.3, 0]).set_fill(col, 1).set_stroke(width=0)
        rim = VMobject().set_points_smoothly([[x, f(x), 0] for x in xs]).set_stroke(TREE_LINE, 1.3, opacity=0.35)
        return VGroup(body, rim)

    clouds = VGroup(*[Union(blob(np.array([cx, cy, 0]), 0.42, rng, squash=0.55),
                            blob(np.array([cx + 0.45, cy - 0.06, 0]), 0.32, rng, squash=0.6),
                            blob(np.array([cx - 0.42, cy - 0.08, 0]), 0.28, rng, squash=0.6))
                      .set_fill(WHITE, 0.55).set_stroke(TREE_LINE, 1, opacity=0.18) for cx, cy in [(-5.6, 3.3), (5.9, 2.9)]])
    birds = VGroup(*[VGroup(Arc(radius=0.09, start_angle=PI / 6, angle=2 * PI / 3),
                            Arc(radius=0.09, start_angle=PI / 6, angle=2 * PI / 3).shift(RIGHT * 0.15))
                     .set_stroke(INK, 1.4, opacity=0.35).move_to([bx, by, 0])
                     for bx, by in [(4.6, 3.45), (4.95, 3.6), (5.25, 3.4)]])
    g = VGroup(clouds, birds, hill(far, "#E2E4D0"))
    for x in (-5.4, -4.8, 5.5):
        g.add(tree(x, far(x) - 0.05, rng.uniform(0.5, 0.7), rng, col="#C9D3B6", line="#A3B08E"))
    g.add(hill(near, "#D5DCC2"))
    for x, h, kind in [(-6.85, 2.0, tree), (-6.32, 1.45, pine), (5.75, 1.0, pine), (6.25, 1.6, tree), (6.85, 2.1, pine)]:
        g.add(kind(x, near(x) - 0.04, h, rng))
    for x in np.linspace(-6.6, 6.6, 16):
        x += rng.uniform(-0.2, 0.2)
        y = near(x) - rng.uniform(0.08, 0.25)
        g.add(VGroup(*[Line([x + d, y, 0], [x + d * 1.8, y + 0.1, 0], stroke_width=1.3, color=TREE_LINE)
                       .set_stroke(opacity=0.5) for d in (-0.04, 0, 0.04)]))
    return g


def presenter_window(center, w, h):
    frame = RoundedRectangle(width=w, height=h, corner_radius=0.25).set_fill(BG_D, 1).set_stroke(INK, 2, opacity=0.3)
    g = VGroup(frame).move_to(center)
    if PLACEHOLDER:
        c = frame.get_center()
        head = Circle(radius=min(w, h) * 0.16).set_fill(STONE_M, 1).set_stroke(width=0).move_to(c + UP * h * 0.12)
        body = Intersection(Ellipse(width=w * 0.75, height=h * 0.5).move_to(c + DOWN * h * 0.36), frame)
        body.set_fill(STONE_M, 1).set_stroke(width=0)
        g.add(head, body, T("Presenter on camera", 15, MUTED).move_to(c + UP * h * 0.4))
    return g


def chem(parts, size=72, colors=None):
    """Chemical formula with real subscripts e.g. [("Li","2"),("CO","3")]."""
    g, x = VGroup(), 0.0
    for i, (base, sub) in enumerate(parts):
        c = colors[i] if colors else INK
        b = H(base, size, c, BOLD).move_to([x, 0, 0], aligned_edge=LEFT)
        g.add(b)
        x = b.get_right()[0] + 0.03
        if sub:
            s = H(sub, size * 0.58, c, BOLD)
            s.move_to([x, b.get_bottom()[1] - 0.08, 0], aligned_edge=DL)
            g.add(s)
            x = s.get_right()[0] + 0.05
    return g


def soft(path, color, width=4):
    """A line with a gentle halo."""
    halo = path.copy().set_fill(opacity=0).set_stroke(color, width + 8, opacity=0.12)
    return VGroup(halo, path.copy().set_fill(opacity=0).set_stroke(color, width, opacity=1))


def glow_dot(color, r=0.07):
    g = VGroup(*[Circle(radius=r * k, stroke_width=0).set_fill(color, o) for k, o in [(2.0, 0.1), (1.4, 0.22)]])
    g.add(Circle(radius=r, stroke_width=0).set_fill(color, 1))
    return g


def tip_at(path, color, size=0.24):
    """Filled arrowhead at the end of an arbitrary path."""
    p1 = path.point_from_proportion(1)
    d = normalize(p1 - path.point_from_proportion(0.98))
    tri = Triangle().set_fill(color, 1).set_stroke(width=0).scale_to_fit_width(size)
    tri.rotate(angle_of_vector(d) - PI / 2)
    return tri.move_to(p1 - d * size * 0.3)


def flow(dots, paths, c0, c1, run_time, lag):
    """Move each dot along its path while blending its color from c0 to c1."""
    anims = []
    for d, p in zip(dots, paths):
        def upd(m, a, p=p):
            m.move_to(p.point_from_proportion(smooth(a)))
            m.set_color(interpolate_color(ManimColor(c0), ManimColor(c1), a))
        anims.append(UpdateFromAlphaFunc(d, upd))
    return LaggedStart(*anims, lag_ratio=lag, run_time=run_time)


# ---------------------------------------------------------------- icons
def battery_icon(w=1.0, h=0.5, color=WATER, level=1.0, stroke=3):
    body = RoundedRectangle(width=w, height=h, corner_radius=h * 0.18, stroke_color=color, stroke_width=stroke)
    nub = RoundedRectangle(width=w * 0.07, height=h * 0.4, corner_radius=0.02,
                           stroke_width=0).set_fill(color, 1).next_to(body, RIGHT, buff=0.02)
    pad = h * 0.15
    fw = max((w - 2 * pad) * level, 0.001)
    fill = Rectangle(width=fw, height=h - 2 * pad, stroke_width=0).set_fill(color, 0.9)
    fill.move_to(body.get_left() + RIGHT * (pad + fw / 2))
    return VGroup(body, nub, fill)


def caution(size=0.6):
    tri = Triangle().set_stroke(width=0).set_fill(SUN, 1).scale_to_fit_width(size)
    tri.round_corners(size * 0.08)
    mark = T("!", 40, INK, BOLD).scale_to_fit_height(size * 0.42)
    mark.move_to(tri.get_center() + DOWN * size * 0.08)
    return VGroup(tri, mark)


def pickaxe(color=OCHRE):
    head = Arc(radius=0.45, start_angle=PI / 2 - 1.0, angle=2.0, stroke_width=7, color=color)
    handle = Line(UP * 0.42, DOWN * 0.5, stroke_width=7, color=STONE_L)
    return VGroup(handle, head).rotate(-PI / 5).move_to(ORIGIN)


def flask(color=STONE_L, liquid=WATER):
    body = Polygon([-0.12, 0.45, 0], [0.12, 0.45, 0], [0.12, 0.1, 0], [0.4, -0.4, 0],
                   [-0.4, -0.4, 0], [-0.12, 0.1, 0], stroke_color=color, stroke_width=4)
    liq = Polygon([-0.22, -0.08, 0], [0.22, -0.08, 0], [0.37, -0.36, 0], [-0.37, -0.36, 0],
                  stroke_width=0).set_fill(liquid, 0.85)
    rim = Line([-0.2, 0.45, 0], [0.2, 0.45, 0], stroke_width=4, color=color)
    return VGroup(liq, body, rim).move_to(ORIGIN)


def lattice(color=WATER):
    hexa = RegularPolygon(6, radius=0.42, stroke_color=color, stroke_width=4)
    pts = hexa.get_vertices()
    spokes = VGroup(*[Line(ORIGIN, p, stroke_width=2, color=color, stroke_opacity=0.6) for p in pts])
    dots = VGroup(*[Dot(p, radius=0.06, color=color) for p in pts], Dot(ORIGIN, radius=0.08, color=INK))
    return VGroup(spokes, hexa, dots)


def truck(color=STONE_L):
    cargo = RoundedRectangle(width=0.7, height=0.38, corner_radius=0.04, stroke_width=0).set_fill(color, 1)
    cab = RoundedRectangle(width=0.26, height=0.28, corner_radius=0.05, stroke_width=0).set_fill(color, 1)
    cab.next_to(cargo, RIGHT, buff=0.03).align_to(cargo, DOWN)
    win = Rectangle(width=0.1, height=0.09, stroke_width=0).set_fill(BG, 1).move_to(cab.get_center() + UR * 0.05)
    wheels = VGroup(*[Circle(radius=0.07).set_fill(BG, 1).set_stroke(color, 2.5).move_to([x, cargo.get_bottom()[1], 0])
                      for x in (cargo.get_left()[0] + 0.15, cargo.get_right()[0] - 0.12, cab.get_center()[0])])
    return VGroup(cargo, cab, win, wheels).move_to(ORIGIN)


def ship(color=STONE_L):
    hull = Polygon([-0.5, 0, 0], [0.5, 0, 0], [0.38, -0.2, 0], [-0.42, -0.2, 0],
                   stroke_width=0).set_fill(color, 1)
    boxes = VGroup(*[Rectangle(width=0.16, height=0.12).set_fill(c, 1).set_stroke(BG, 1.5)
                     for c in (OCHRE, WATER, LEAF, STONE)]).arrange(RIGHT, buff=0.02)
    boxes.next_to(hull, UP, buff=0).shift(RIGHT * 0.08)
    bridge = Rectangle(width=0.12, height=0.24, stroke_width=0).set_fill(color, 1)
    bridge.next_to(hull, UP, buff=0).align_to(hull, LEFT).shift(RIGHT * 0.06)
    return VGroup(hull, boxes, bridge).move_to(ORIGIN)


def sun(color=SUN):
    c = Circle(radius=0.28, stroke_width=0).set_fill(color, 1)
    rays = VGroup(*[Line(RIGHT * 0.4, RIGHT * 0.56, stroke_width=4, color=color).rotate(a, about_point=ORIGIN)
                    for a in np.arange(0, TAU, TAU / 8)])
    return VGroup(c, rays)


def factory(color=STONE_L):
    pts = [(-0.5, -0.3), (0.5, -0.3), (0.5, 0.12), (0.2, -0.04), (0.2, 0.12), (-0.1, -0.04),
           (-0.1, 0.12), (-0.32, 0.0), (-0.32, 0.42), (-0.46, 0.42), (-0.46, -0.0), (-0.5, 0.0)]
    body = Polygon(*[[x, y, 0] for x, y in pts], stroke_width=0).set_fill(color, 1)
    wins = VGroup(*[Square(0.09, stroke_width=0).set_fill(BG, 1).move_to([x, -0.15, 0]) for x in (-0.25, 0.0, 0.25)])
    return VGroup(body, wins).move_to(ORIGIN)


def recycle_icon(color=LEAF, r=0.3):
    g = VGroup()
    for k in range(3):
        a = Arc(radius=r, start_angle=k * TAU / 3 + 0.3, angle=TAU / 3 - 0.6, stroke_width=5, color=color)
        a.add_tip(tip_length=0.14, tip_width=0.14)
        g.add(a)
    return g


def leaf(color=LEAF, h=0.42):
    a = ArcBetweenPoints(DOWN * 0.5, UP * 0.5, angle=-1.4)
    b = ArcBetweenPoints(UP * 0.5, DOWN * 0.5, angle=-1.4)
    shape = VMobject()
    shape.append_points(a.points)
    shape.append_points(b.points)
    shape.set_fill(color, 1).set_stroke(width=0)
    vein = Line(DOWN * 0.6, UP * 0.38, stroke_width=2.5, color=BG_D)
    return VGroup(shape, vein).rotate(-0.6).scale_to_fit_height(h)


def landfill(w=2.4):
    hw = w / 2
    mound = Polygon(*[[x, 0.6 * np.sqrt(max(1 - (x / hw) ** 2, 0)), 0] for x in np.linspace(-hw, hw, 40)],
                    stroke_width=0).set_fill(SOIL[0], 1)
    rng = np.random.default_rng(11)
    junk = VGroup()
    for _ in range(9):
        x = rng.uniform(-hw * 0.65, hw * 0.65)
        y = rng.uniform(0.06, 0.6 * np.sqrt(1 - (x / hw) ** 2) - 0.12)
        junk.add(battery_icon(0.26, 0.13, rng.choice([STONE, MUTED, SUN]), rng.uniform(0.1, 0.5), 1.5)
                 .rotate(rng.uniform(-0.6, 0.6)).move_to([x, y, 0]))
    ground = Line([-hw - 0.3, 0, 0], [hw + 0.3, 0, 0], stroke_width=3, color=SOIL[0])
    return VGroup(mound, junk, ground)


def pin(color):
    ring = Circle(radius=0.11).set_fill(BG_D, 1).set_stroke(color, 3)
    return VGroup(ring, Dot(radius=0.045, color=color))


def tint(color, a):
    """Opaque mix of paper and a color."""
    return interpolate_color(ManimColor(PAPER), ManimColor(color), a)


def pill(text, color, size=15, fill_op=0.16):
    t = T(text, size, color, BOLD)
    box = RoundedRectangle(width=t.width + 0.4, height=t.height + 0.24, corner_radius=0.14)
    box.set_fill(tint(color, fill_op), 1).set_stroke(color, 1.5).move_to(t)
    return VGroup(box, t)


# ---------------------------------------------------------------- base scene
class Base(Scene):
    DUR = 10
    FADE = 0.6

    def setup(self):
        self.bg = VGroup()
        if not config.transparent:
            self.bg = landscape()
            self.add(self.bg)

    def intro(self, *extra):
        if extra:
            self.play(*[FadeIn(e) for e in extra], run_time=0.8)
        else:
            self.wait(0.8)

    def finish(self, fade=True):
        rest = self.DUR - self.renderer.time - (self.FADE if fade else 0)
        if rest > 0:
            self.wait(rest)
        if fade:
            self.play(*[FadeOut(m) for m in self.mobjects if m is not self.bg], run_time=self.FADE)


# ================================================================ SCENE 1
class S1_Hook(Base):
    """0:00-0:14  Presenter holds a phone while the graphic opens it up."""
    DUR = 14

    def construct(self):
        win = presenter_window([-4.05, -0.35, 0], 4.9, 6.2)
        self.intro(win)

        W, H_ = 2.0, 4.0

        def slab(fill, op, stroke, sw=2, w=W, h=H_, cr=0.3):
            return RoundedRectangle(width=w, height=h, corner_radius=cr).set_fill(fill, op).set_stroke(stroke, sw)

        housing = VGroup(slab(STONE_M, 1, STONE_L), slab(STONE_D, 1, STONE, 1, W - 0.18, H_ - 0.18, 0.22))
        battery = VGroup(
            RoundedRectangle(width=1.55, height=2.3, corner_radius=0.12).set_fill(WATER_D, 1).set_stroke(WATER, 2.5),
            T("Li-ion", 22, WATER, BOLD),
        ).shift(DOWN * 0.6)
        board_main = RoundedRectangle(width=1.75, height=1.05, corner_radius=0.1).set_fill(BG_D, 1) \
            .set_stroke(STONE_L, 1.5).move_to(UP * 1.3)
        chip = Square(0.38).set_fill(STONE_M, 1).set_stroke(STONE_L, 1).move_to([0.35, 1.2, 0])
        mem = Rectangle(width=0.32, height=0.18).set_fill(STONE_M, 1).set_stroke(STONE_L, 1).move_to([0.35, 1.6, 0])
        cams = VGroup(*[Circle(radius=0.13).set_fill(BG, 1).set_stroke(STONE_L, 1.5).move_to([-0.5, y, 0])
                        for y in (1.52, 1.1)])
        sub = RoundedRectangle(width=1.1, height=0.25, corner_radius=0.06).set_fill(BG_D, 1) \
            .set_stroke(STONE_L, 1.5).move_to(DOWN * 1.82)
        flex = Line([0.84, 0.78, 0], [0.84, -1.7, 0], stroke_width=3, color=STONE)
        board = VGroup(board_main, chip, mem, cams, sub, flex)

        app_cols = [WATER, LEAF, OCHRE, STONE, SUN]
        icons = VGroup(*[RoundedRectangle(width=0.28, height=0.28, corner_radius=0.07, stroke_width=0)
                         .set_fill(app_cols[(i * 3) % 5], 0.85) for i in range(20)]) \
            .arrange_in_grid(5, 4, buff=0.14).move_to(UP * 0.15)
        dock = VGroup(*[RoundedRectangle(width=0.28, height=0.28, corner_radius=0.07, stroke_width=0)
                        .set_fill(c, 0.85) for c in (LEAF, WATER, STONE, OCHRE)]) \
            .arrange(RIGHT, buff=0.14).move_to(DOWN * 1.55)
        clock = T("9:41", 15, PAPER, BOLD).move_to(UP * 1.5)
        display = VGroup(slab("#3E4A43", 1, STONE_L, 1.5), icons, dock, clock)
        glass = VGroup(slab(WHITE, 0.3, STONE_L, 2),
                       RoundedRectangle(width=0.55, height=0.14, corner_radius=0.07, stroke_width=0)
                       .set_fill("#3E4A43", 1).move_to(UP * 1.78))

        layers = [housing, battery, board, display, glass]
        orig = [L.copy() for L in layers]
        phone = VGroup(*layers).scale(0.85).move_to([2.7, -0.3, 0])
        self.play(FadeIn(phone, shift=UP * 0.3), run_time=1.0)

        # flat to isometric (rotate 50 deg then squash vertically)
        th = 50 * DEGREES
        M = 0.8 * np.diag([1, 0.55]) @ np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        base_c = np.array([1.3, -0.35, 0])
        targets = [o.apply_matrix(M, about_point=ORIGIN).shift(base_c + UP * 0.05 * i) for i, o in enumerate(orig)]
        self.play(*[Transform(L, t) for L, t in zip(layers, targets)], run_time=1.2)

        edge_cols = ["#AEB2A4", "#9EBFBB", "#C2C5B6", "#2A332D", "#C7CBBC"]
        edges = [L[0].copy().set_fill(c, 1).set_stroke(width=0).shift(DOWN * 0.06) for L, c in zip(layers, edge_cols)]
        self.remove(*layers)
        stacks = [VGroup(e, L) for e, L in zip(edges, layers)]
        self.add(*stacks)

        gap = 0.62
        shifts = [UP * ((i - 2) * gap - 0.05 * i) for i in range(5)]
        self.play(LaggedStart(*[s.animate.shift(v) for s, v in zip(stacks, shifts)], lag_ratio=0.08), run_time=1.2)

        names = [("Housing", "Al  Mg"), ("Battery", "Li  Co  Ni  Mn  C"), ("Circuit board", "Au  Ag  Cu  Ta  Pd"),
                 ("Display", "In  Sn  Ga"), ("Glass", "Si  Al  K")]
        labels = VGroup()
        for i, (s, (name, sub_t)) in enumerate(zip(stacks, names)):
            pts = s[1][0].get_points()
            anchor = pts[np.argmax(pts[:, 0])] + RIGHT * 0.08
            label_y = base_c[1] + (i - 2) * 0.95 - 0.2
            col = WATER if i == 1 else STONE_L
            dot = Dot(anchor, radius=0.04, color=col)
            ln = VMobject().set_points_as_corners([anchor, [3.15, anchor[1], 0], [3.35, label_y, 0], [3.45, label_y, 0]]) \
                .set_stroke(col, 1.5)
            txt = VGroup(T(name, 19, WATER if i == 1 else INK, BOLD), T(sub_t, 14, MUTED)) \
                .arrange(DOWN, aligned_edge=LEFT, buff=0.05)
            txt.next_to(ln.get_end(), RIGHT, buff=0.12)
            labels.add(VGroup(dot, ln, txt))
        self.play(LaggedStart(*[Create(l[1]) for l in labels], *[FadeIn(l[0]) for l in labels],
                              *[FadeIn(l[2], shift=RIGHT * 0.15) for l in labels], lag_ratio=0.04), run_time=1.3)
        self.wait(0.6)

        # zoom into the battery
        body = RoundedRectangle(width=1.35, height=1.9, corner_radius=0.14).set_fill(WATER_D, 1).set_stroke(WATER, 3)
        tab = RoundedRectangle(width=0.4, height=0.13, corner_radius=0.04, stroke_width=0).set_fill(WATER, 1) \
            .next_to(body, UP, buff=0)
        stripes = VGroup(*[RoundedRectangle(width=1.05, height=0.07, corner_radius=0.03, stroke_width=0)
                           .set_fill(c, 0.9) for c in [WATER, STONE_L, STONE] * 4]).arrange(DOWN, buff=0.045)
        stripes.move_to(body).shift(DOWN * 0.13)
        cell_name = T("Li-ion cell", 14, INK, BOLD).next_to(stripes, UP, buff=0.1)
        big = VGroup(body, tab, stripes, cell_name).move_to([0.5, 0.9, 0])
        others = [s for i, s in enumerate(stacks) if i != 1]
        self.play(*[FadeOut(o, shift=DOWN * 0.2) for o in others], FadeOut(labels),
                  ReplacementTransform(stacks[1], big), run_time=1.1)

        # counter 1 to 30+
        v = ValueTracker(1)
        num_anchor = np.array([4.75, 1.5, 0])
        num = always_redraw(lambda: H(str(int(round(v.get_value()))), 100, WATER, BOLD).move_to(num_anchor, aligned_edge=RIGHT))
        plus = H("+", 100, WATER, BOLD).next_to(num_anchor, RIGHT, buff=0.06)
        caption = H("minerals dug from the Earth", 28, INK, NORMAL, italic=True).move_to([4.4, 0.45, 0])

        syms = ["Li", "Co", "Ni", "Mn", "C", "Al", "Cu", "Au", "Ag", "Pd", "Sn", "In", "Ga", "Ta", "W",
                "Nd", "Pr", "Dy", "Tb", "Gd", "La", "Ce", "Eu", "Y", "Sb", "Bi", "Mg", "Si", "Ge", "Zn"]
        tiles = VGroup()
        for i, s in enumerate(syms):
            c = WATER if i < 5 else STONE_L
            sq = RoundedRectangle(width=0.46, height=0.46, corner_radius=0.08).set_fill(c, 0.12).set_stroke(c, 1.5)
            tiles.add(VGroup(sq, T(s, 16, c, BOLD).move_to(sq)))
        tiles.arrange_in_grid(3, 10, buff=0.08).move_to([3.0, -1.3, 0])
        note = T("Five of them power the battery alone", 18, WATER, BOLD).next_to(tiles, DOWN, buff=0.28)
        hook = H("So what happens when it dies?", 32, INK, SEMIBOLD, italic=True).move_to([2.6, -3.2, 0])

        self.add(num)
        self.play(FadeIn(caption, shift=UP * 0.1), run_time=0.4)
        self.play(v.animate.set_value(30), LaggedStart(*[FadeIn(t, scale=0.6) for t in tiles], lag_ratio=1.0),
                  run_time=2.6, rate_func=linear)
        self.play(FadeIn(plus, scale=1.5), run_time=0.4)
        self.play(FadeIn(note, shift=UP * 0.1), *[t[0].animate.set_fill(WATER, 0.32) for t in tiles[:5]], run_time=0.5)
        self.play(FadeIn(hook, shift=UP * 0.15), run_time=0.8)
        self.finish()


# ================================================================ SCENE 2
class S2_Extraction(Base):
    """0:14-0:30  Hard rock pegmatite versus salar brine ending in Li2CO3."""
    DUR = 16

    def construct(self):
        self.intro()
        divider = DashedLine([0, 2.85, 0], [0, -1.45, 0], dash_length=0.08, stroke_width=1.5, color=STONE)
        lt = VGroup(H("Hard rock mining", 26), T("Lithium locked inside pegmatite rock", 16, OCHRE, BOLD)) \
            .arrange(DOWN, buff=0.1).move_to([-3.6, 2.4, 0])
        rt = VGroup(H("Brine evaporation", 26), T("Salty water pumped up from under the desert", 16, WATER, BOLD)) \
            .arrange(DOWN, buff=0.1).move_to([3.0, 2.4, 0])

        # left: geological cross-section
        x0, x1 = -6.6, -0.6

        def band(yl0, yr0, yl1, yr1, col):
            return Polygon([x0, yl0, 0], [x1, yr0, 0], [x1, yr1, 0], [x0, yl1, 0], stroke_width=0).set_fill(col, 1)

        strata = VGroup(band(1.2, 1.2, 0.55, 0.7, SOIL[0]), band(0.55, 0.7, -0.35, -0.1, SOIL[1]),
                        band(-0.35, -0.1, -1.4, -1.4, SOIL[2]))
        vtop = lambda x: -0.9 + (x - x0) / 6 * 1.8
        vein = Polygon([x0, vtop(x0), 0], [x1, vtop(x1), 0], [x1, vtop(x1) - 0.5, 0], [x0, vtop(x0) - 0.5, 0],
                       stroke_width=0).set_fill(RUST, 1)
        rng = np.random.default_rng(3)
        crystals = VGroup()
        for _ in range(26):
            x = rng.uniform(x0 + 0.2, x1 - 0.2)
            y = vtop(x) - rng.uniform(0.1, 0.4)
            crystals.add(Rectangle(width=rng.uniform(0.12, 0.26), height=0.07, stroke_width=0)
                         .set_fill(OCHRE, 0.95).rotate(rng.uniform(0, PI)).move_to([x, y, 0]))
        pit = Polygon(*[[x, y, 0] for x, y in [(-5.4, 1.2), (-5.05, 0.75), (-4.75, 0.75), (-4.4, 0.3), (-4.1, 0.3),
                                               (-3.85, -0.1), (-2.85, -0.1), (-2.6, 0.3), (-2.3, 0.3), (-1.95, 0.75),
                                               (-1.65, 0.75), (-1.3, 1.2)]]).set_fill(BG, 1).set_stroke(STONE_L, 1.5)
        surface = Line([x0, 1.2, 0], [x1, 1.2, 0], stroke_width=2, color=STONE_L)
        pit_lbl = T("Open pit mine", 14, MUTED, BOLD).move_to([-3.35, 0.95, 0])
        vein_lbl = T("A vein full of spodumene crystals", 14, PAPER, BOLD).move_to([-4.6, -1.2, 0])
        left = VGroup(strata, vein, crystals, pit, surface, pit_lbl, vein_lbl)

        # right: salar with evaporation ponds
        gx0, gx1 = 0.6, 6.6
        crust = Polygon([gx0, 0.6, 0], [gx1, 0.6, 0], [gx1, -0.3, 0], [gx0, -0.3, 0], stroke_width=0).set_fill("#E6DCC4", 1)
        aquifer = Polygon([gx0, -0.3, 0], [gx1, -0.3, 0], [gx1, -1.4, 0], [gx0, -1.4, 0], stroke_width=0).set_fill(WATER_D, 1)
        aq_lbl = T("Lithium rich brine lies deep underground", 14, INK, BOLD).move_to([3.8, -0.95, 0])
        surface2 = Line([gx0, 0.6, 0], [gx1, 0.6, 0], stroke_width=2, color=STONE_L)
        PX = [1.95, 3.6, 5.25]
        ponds = VGroup(*[Polygon([cx - 0.7, 0.6, 0], [cx - 0.55, 0.1, 0], [cx + 0.55, 0.1, 0], [cx + 0.7, 0.6, 0])
                         .set_fill(BG_D, 1).set_stroke(STONE_L, 1.5) for cx in PX])
        pipe = VMobject().set_points_as_corners([[0.85, -0.85, 0], [0.85, 0.95, 0], [1.45, 0.95, 0], [1.45, 0.55, 0]]) \
            .set_stroke(STONE_L, 4)
        right = VGroup(crust, aquifer, aq_lbl, surface2, ponds)

        levels = [ValueTracker(0.001) for _ in PX]
        pond_cols = ["#8DBDC2", "#A3C9AE", "#D3D49A"]   # real ponds turn yellow green as lithium concentrates

        def liquid(cx, lv, col):
            def make():
                y = 0.1 + 0.42 * lv.get_value()
                hw = 0.55 + (y - 0.1) / 0.5 * 0.15
                return Polygon([cx - 0.55, 0.1, 0], [cx + 0.55, 0.1, 0], [cx + hw, y, 0], [cx - hw, y, 0],
                               stroke_width=0).set_fill(col, 0.9)
            return always_redraw(make)

        liquids = VGroup(*[liquid(cx, lv, c) for cx, lv, c in zip(PX, levels, pond_cols)])
        transfer = VGroup(*[Arrow([a, 0.82, 0], [b, 0.82, 0], buff=0, stroke_width=2.5, color=STONE_L,
                                  max_tip_length_to_length_ratio=0.35) for a, b in [(2.5, 3.05), (4.15, 4.7)]])
        sunny = sun().move_to([6.25, 1.7, 0])
        conc_arrow = Arrow([1.2, 1.3, 0], [5.3, 1.3, 0], buff=0, stroke_width=2, color=MUTED,
                           max_tip_length_to_length_ratio=0.04)
        conc_lbl = T("The sun concentrates the lithium over many months", 14, MUTED, BOLD).move_to([3.15, 1.56, 0])

        self.play(Create(divider), FadeIn(lt, shift=DOWN * 0.1), FadeIn(rt, shift=DOWN * 0.1), run_time=0.9)
        self.play(FadeIn(left, lag_ratio=0.05), FadeIn(right, lag_ratio=0.05), run_time=1.3)
        self.add(liquids)
        self.play(Create(pipe), FadeIn(sunny, scale=0.6), run_time=0.7)
        self.play(LaggedStart(*[lv.animate.set_value(1) for lv in levels], lag_ratio=0.45),
                  LaggedStart(*[GrowArrow(a) for a in transfer], lag_ratio=0.5),
                  Rotate(sunny[1], PI / 4), run_time=2.2)

        def wisp(x):
            return ParametricFunction(lambda t: np.array([0.05 * np.sin(t * 18), t, 0]), t_range=[0, 0.35],
                                      stroke_width=2, color=STONE_L).move_to([x, 0.75, 0]).set_stroke(opacity=0.7)
        wisps = VGroup(*[wisp(cx + d) for cx in PX for d in (-0.25, 0.25)])
        self.play(FadeIn(wisps, shift=UP * 0.2), GrowArrow(conc_arrow), FadeIn(conc_lbl), run_time=0.6)
        self.play(wisps.animate.shift(UP * 0.45).set_opacity(0), *[lv.animate.set_value(0.82) for lv in levels], run_time=0.9)

        # processing then particles converge on the formula
        def box(text, x):
            b = RoundedRectangle(width=2.5, height=0.55, corner_radius=0.14).set_fill(BG_D, 1).set_stroke(STONE_L, 1.5)
            return VGroup(b, T(text, 15, STONE_L, BOLD)).move_to([x, -2.05, 0])
        box_l = box("Crushed and heated", -2.3)
        box_r = box("Purified and dried", 2.3)
        formula = chem([("Li", "2"), ("CO", "3")], 64, [WATER, INK]).move_to([0, -2.75, 0])
        f_sub = T("Lithium carbonate ready to become a battery", 17, MUTED, BOLD).next_to(formula, DOWN, buff=0.16)
        f_msg = H("Every charge begins with tonnes of rock and months of sunlight", 21, INK, NORMAL, italic=True) \
            .move_to(f_sub)
        self.play(FadeIn(box_l, shift=UP * 0.1), FadeIn(box_r, shift=UP * 0.1), run_time=0.5)

        N = 28
        fl, fr = formula.get_left()[0], formula.get_right()[0]
        fy0, fy1 = formula.get_bottom()[1], formula.get_top()[1]

        def target():
            return np.array([rng.uniform(fl + 0.1, fr - 0.1), rng.uniform(fy0 + 0.1, fy1 - 0.1), 0])

        dl, pl, dr, pr = VGroup(), [], VGroup(), []
        for _ in range(N):
            s = np.array([rng.uniform(-3.7, -3.0), rng.uniform(-0.05, 0.2), 0])
            pl.append(VMobject().set_points_smoothly([s, [s[0] + 0.3, 1.65, 0], [-1.0, 0.9, 0], [-1.4, -1.6, 0],
                                                      [-1.7, -2.05, 0], target()]))
            dl.add(glow_dot(OCHRE, 0.045).move_to(s))
            s = np.array([rng.uniform(4.85, 5.65), rng.uniform(0.15, 0.4), 0])
            pr.append(VMobject().set_points_smoothly([s, [4.9, 1.1, 0], [1.0, 0.95, 0], [1.4, -1.6, 0],
                                                      [1.7, -2.05, 0], target()]))
            dr.add(glow_dot("#A9AC5E", 0.045).move_to(s))
        self.add(dl, dr)
        self.play(flow(dl, pl, OCHRE, WATER, 3.6, 0.03), flow(dr, pr, "#A9AC5E", WATER, 3.6, 0.03))
        self.play(FadeIn(formula, scale=0.85), FadeOut(dl), FadeOut(dr), run_time=0.8)
        self.play(FadeIn(f_sub, shift=UP * 0.1), run_time=0.5)
        self.wait(0.6)
        self.play(FadeOut(f_sub, shift=UP * 0.25), run_time=0.35)
        self.play(FadeIn(f_msg, shift=UP * 0.15), run_time=0.6)
        self.finish()


# ================================================================ SCENE 3
class S3_Journey(Base):
    """0:30-0:48  Split screen: presenter beside the supply chain and its fragile end."""
    DUR = 18

    def construct(self):
        win = presenter_window([-4.75, -0.3, 0], 3.6, 5.6)
        self.intro(win)
        CX = 2.15
        BY = 1.0
        top = BY + 0.16
        belt = RoundedRectangle(width=8.9, height=0.32, corner_radius=0.16).set_fill(BG_D, 1) \
            .set_stroke(STONE_L, 2).move_to([CX, BY, 0])
        bl, br = CX - 4.3, CX + 4.3

        def chevron(x):
            return VMobject().set_points_as_corners([[x - 0.04, BY + 0.06, 0], [x + 0.04, BY, 0], [x - 0.04, BY - 0.06, 0]]) \
                .set_stroke(STONE, 2)
        chevrons = VGroup(*[chevron(x) for x in np.arange(bl, br, 0.43)])

        def chev_upd(m, dt):
            for c in m:
                c.shift(RIGHT * dt * 0.9)
                if c.get_x() > br:
                    c.shift(LEFT * 0.43 * len(m))
        chevrons.add_updater(chev_upd)

        SX = [-1.2, 1.15, 3.5, 5.85]
        names = ["Mining", "Refining", "Cathode making", "Cell assembly"]
        cols = [OCHRE, STONE_L, WATER, WATER]
        icons = [pickaxe(OCHRE).scale(0.58), flask(STONE_L, WATER).scale(0.58), lattice(WATER).scale(0.58),
                 battery_icon(0.54, 0.28, WATER)]
        badges, slabels = VGroup(), VGroup()
        for x, c, ic, n in zip(SX, cols, icons, names):
            circ = Circle(radius=0.4).set_fill(c, 0.1).set_stroke(c, 3).move_to([x, 2.15, 0])
            badges.add(VGroup(circ, ic.move_to(circ)))
            slabels.add(T(n, 17, INK, BOLD).move_to([x, 0.55, 0]))
        stage_arrows = VGroup(*[Arrow([a + 0.5, 2.15, 0], [b - 0.5, 2.15, 0], buff=0, stroke_width=2,
                                      color=STONE, max_tip_length_to_length_ratio=0.12) for a, b in zip(SX[:-1], SX[1:])])

        self.play(FadeIn(belt), FadeIn(chevrons),
                  LaggedStart(*[FadeIn(b, scale=0.7) for b in badges], lag_ratio=0.15),
                  LaggedStart(*[GrowArrow(a) for a in stage_arrows], lag_ratio=0.2),
                  LaggedStart(*[FadeIn(l, shift=UP * 0.1) for l in slabels], lag_ratio=0.15), run_time=1.2)

        rock = Polygon(*[[x, y, 0] for x, y in [(-0.2, -0.12), (0.05, -0.15), (0.22, -0.05), (0.15, 0.12),
                                                 (-0.05, 0.17), (-0.2, 0.05)]]).set_fill(OCHRE, 1).set_stroke(RUST, 2)
        powder = Polygon(*[[x, 0.2 * (1 - (x / 0.26) ** 2), 0] for x in np.linspace(-0.26, 0.26, 24)],
                         stroke_width=0).set_fill(PAPER, 1).set_stroke(STONE_L, 1.5)
        cathode = RegularPolygon(6, radius=0.19).set_fill(WATER, 0.9).set_stroke(INK, 1.5)
        cell = battery_icon(0.5, 0.28, WATER)

        def at(m, x):
            return m.copy().move_to([x, top + m.height / 2 + 0.02, 0])

        item = at(rock, SX[0])
        self.play(FadeIn(item, shift=DOWN * 0.3), Indicate(badges[0], color=OCHRE, scale_factor=1.12), run_time=0.5)
        for k, shape in zip(range(1, 4), [powder, cathode, cell]):
            self.play(item.animate.move_to([SX[k], item.get_y(), 0]), run_time=0.65)
            self.play(Transform(item, at(shape, SX[k])),
                      Indicate(badges[k], color=cols[k] if k != 1 else INK, scale_factor=1.12), run_time=0.45)
            if k == 1:
                ring = Circle(radius=0.52).set_stroke(CLAY, 3).move_to(badges[1])
                tag = pill("Bottleneck", CLAY, 14).move_to([SX[1], 2.95, 0])
                bn_note = T("Most refining happens in just one country", 14, CLAY, BOLD).move_to([SX[1], 0.2, 0])
                self.play(Create(ring), FadeIn(tag, shift=DOWN * 0.1), FadeIn(bn_note),
                          badges[1][0].animate.set_stroke(CLAY), run_time=0.6)
                ring.add_updater(lambda m, dt: m.set_stroke(opacity=0.55 + 0.45 * np.sin(self.renderer.time * 6)))
        self.play(item.animate.move_to([br + 0.2, item.get_y(), 0]).set_opacity(0), run_time=0.5)
        self.remove(item)

        # world map with real continent shapes
        frame = RoundedRectangle(width=8.9, height=2.9, corner_radius=0.2).set_fill(OCEAN, 1).set_stroke(STONE, 1.2) \
            .move_to([CX, -1.95, 0])

        def proj(lon, lat):
            return np.array([CX + lon / 180 * 4.2, -1.95 + lat / 90 * 1.35, 0])

        continents = [
            [(-165, 65), (-140, 70), (-95, 72), (-75, 62), (-55, 52), (-80, 25), (-97, 18), (-105, 22), (-125, 40), (-130, 55)],
            [(-80, 10), (-60, 8), (-35, -7), (-40, -22), (-58, -38), (-68, -55), (-75, -45), (-71, -18), (-81, -5)],
            [(-10, 36), (-10, 44), (0, 50), (5, 58), (25, 71), (40, 68), (40, 45), (28, 41), (15, 38)],
            [(-17, 15), (-10, 35), (10, 37), (32, 31), (43, 12), (51, 11), (40, -15), (32, -28), (20, -35), (12, -18), (9, 4), (-8, 5)],
            [(40, 45), (40, 68), (70, 73), (110, 76), (140, 72), (178, 68), (160, 60), (140, 50), (130, 40), (122, 30),
             (120, 22), (108, 10), (100, 1), (95, 15), (80, 8), (72, 20), (60, 25), (48, 30), (35, 33)],
            [(114, -22), (130, -12), (142, -11), (153, -28), (147, -38), (135, -34), (115, -34)],
        ]
        land = VGroup(*[Polygon(*[proj(*q) for q in c]).round_corners(0.04).set_fill(LAND, 1).set_stroke(width=0)
                        for c in continents])
        P = {"chile": proj(-69, -24), "aus": proj(121, -28), "china": proj(110, 33),
             "korea": proj(129, 37), "eu": proj(9, 50), "mkt": proj(24, 44)}
        pin_cols = {"chile": OCHRE, "aus": OCHRE, "china": CLAY, "korea": WATER, "eu": WATER, "mkt": LEAF}
        pin_txt = {"chile": ("Chile", DOWN), "aus": ("Australia", DOWN), "china": ("China", UL),
                   "korea": ("Korea", RIGHT), "eu": ("Europe", UP), "mkt": ("You", DR)}
        pins, plabels = VGroup(), VGroup()
        for k, p in P.items():
            pins.add(pin(pin_cols[k]).move_to(p))
            txt, d = pin_txt[k]
            plabels.add(T(txt, 13, pin_cols[k] if k != "mkt" else INK, BOLD).next_to(p, d, buff=0.14))
        map_title = T("Thousands of kilometres by sea", 13, MUTED, BOLD).move_to(frame.get_corner(UL) + [1.45, -0.25, 0])
        self.play(FadeIn(frame), FadeIn(land, lag_ratio=0.1), FadeIn(map_title), run_time=0.8)
        self.play(LaggedStart(*[FadeIn(p, scale=0.5) for p in pins], lag_ratio=0.1),
                  LaggedStart(*[FadeIn(l) for l in plabels], lag_ratio=0.1), run_time=0.6)

        legs = [("chile", "china", 0.6, OCHRE, ship().scale(0.36)),
                ("aus", "china", 0.5, OCHRE, ship().scale(0.36)),
                ("china", "korea", -0.9, STONE_L, ship().scale(0.36)),
                ("korea", "eu", -0.5, WATER, ship().scale(0.36).flip(UP)),
                ("eu", "mkt", 0, WATER, truck().scale(0.34))]
        built = []
        for a, b, ang, col, veh in legs:
            path = (ArcBetweenPoints(P[a], P[b], angle=ang) if ang else Line(P[a], P[b])).set_stroke(col, 2.5)
            built.append((path, veh.move_to(P[a])))
        for idx, rt in [((0, 1), 1.2), ((2,), 0.5), ((3,), 0.9), ((4,), 0.5)]:
            self.add(*[built[i][1] for i in idx])
            self.play(*[a for i in idx for a in (Create(built[i][0]), MoveAlongPath(built[i][1], built[i][0]))],
                      run_time=rt, rate_func=smooth)
            self.play(*[FadeOut(built[i][1]) for i in idx], run_time=0.15)

        # the battery wears out
        map_all = VGroup(frame, land, pins, plabels, map_title, *[b[0] for b in built])
        soh = ValueTracker(1.0)

        def soh_col():
            return interpolate_color(ManimColor(WATER), ManimColor(SUN), np.clip((1 - soh.get_value()) / 0.28, 0, 1))
        bat = always_redraw(lambda: battery_icon(1.6, 0.78, soh_col(), soh.get_value(), 4).move_to([0.0, -1.8, 0]))
        pct = always_redraw(lambda: H(f"{int(round(soh.get_value() * 100))}%", 44, soh_col(), BOLD).move_to([2.2, -1.65, 0]))
        pct_lbl = T("Battery health", 15, MUTED, BOLD).move_to([2.2, -2.25, 0])
        warn = caution(0.62).move_to([4.2, -1.55, 0])
        worn = T("Worn out after a few years", 15, SUN, BOLD).next_to(warn, DOWN, buff=0.15)
        self.play(map_all.animate.fade(0.85), FadeIn(pct_lbl), run_time=0.5)
        self.add(bat, pct)
        self.play(soh.animate.set_value(0.72), run_time=1.2, rate_func=linear)
        self.play(FadeIn(warn, scale=1.8), FadeIn(worn, shift=UP * 0.1), run_time=0.4)

        # and too often it is buried
        dump = landfill(2.4).move_to([0.4, -2.55, 0])
        dead = battery_icon(0.5, 0.25, SUN, 0.2, 2).move_to([0.0, -1.8, 0])
        self.play(FadeOut(pct), FadeOut(pct_lbl), FadeOut(warn), FadeOut(worn), FadeOut(bat), FadeIn(dead),
                  FadeIn(dump, shift=UP * 0.15), run_time=0.6)
        self.play(dead.animate.move_to([0.3, -2.3, 0]).rotate(0.4).set_opacity(0.6), run_time=0.5)
        buried = T("Too many end up buried in landfills", 18, CLAY, BOLD).move_to([3.95, -2.05, 0])
        hope = H("but they don't have to", 22, LEAF, NORMAL, italic=True).next_to(buried, DOWN, buff=0.18)
        self.play(FadeIn(buried, shift=UP * 0.1), run_time=0.5)
        self.play(FadeIn(hope, shift=UP * 0.1), run_time=0.6)
        self.finish()


# ================================================================ SCENE 4
class S4_Recycling(Base):
    """0:48-1:10  Hydrometallurgy returns most of the metal to manufacturing."""
    DUR = 22

    def construct(self):
        self.intro()
        clock = ValueTracker(0).add_updater(lambda m, dt: m.increment_value(dt))
        self.add(clock)

        RX = -4.6
        body = RoundedRectangle(width=1.9, height=2.6, corner_radius=0.35).set_fill(STONE_D, 1) \
            .set_stroke(STONE_L, 2.5).move_to([RX, -0.4, 0])
        hopper = Polygon([RX - 0.6, 1.6, 0], [RX + 0.6, 1.6, 0], [RX + 0.25, 1.0, 0], [RX - 0.25, 1.0, 0]) \
            .set_fill(STONE_M, 1).set_stroke(STONE_L, 2)
        neck = Rectangle(width=0.5, height=0.12).set_fill(STONE_M, 1).set_stroke(STONE_L, 2).move_to([RX, 0.94, 0])
        WC = np.array([RX, -0.45, 0])
        window = RoundedRectangle(width=1.4, height=1.6, corner_radius=0.2).set_fill(BG_D, 1).set_stroke(STONE, 1.5).move_to(WC)
        liq = Intersection(window, Rectangle(width=2, height=1.2).move_to(WC + DOWN * 0.5)).set_fill("#8FBDB2", 0.9).set_stroke(width=0)
        shaft = Line([RX, 0.9, 0], [RX, -1.0, 0], stroke_width=3, color=STONE_L)
        blades = always_redraw(lambda: Line(LEFT, RIGHT, stroke_width=4, color=STONE_L)
                               .scale(0.4 * max(abs(np.cos(clock.get_value() * 4)), 0.05)).move_to([RX, -1.0, 0]))
        legs = VGroup(*[Line([RX + d, -1.7, 0], [RX + d * 1.15, -2.05, 0], stroke_width=4, color=STONE) for d in (-0.65, 0.65)])
        rng = np.random.default_rng(7)
        bubbles = VGroup(*[Circle(radius=rng.uniform(0.025, 0.05)).set_stroke(WHITE, 1.2, opacity=0.8)
                           .move_to([RX + rng.uniform(-0.55, 0.55), rng.uniform(-1.15, 0.0), 0]) for _ in range(9)])

        def bub_upd(m, dt):
            for b in m:
                b.shift(UP * dt * 0.45)
                if b.get_y() > 0.05:
                    b.set_y(-1.15)
        bubbles.add_updater(bub_upd)
        pipe = Line([RX + 0.95, -0.2, 0], [-2.8, -0.2, 0], stroke_width=6, color=STONE)
        r_lbl = T("Recycling reactor", 17, INK, BOLD).move_to([RX, -2.35, 0])
        reactor = VGroup(legs, pipe, body, hopper, neck, window, liq, shaft)

        self.play(FadeIn(reactor, lag_ratio=0.05), FadeIn(r_lbl), run_time=1.0)
        self.add(blades, bubbles)

        spent = VGroup(battery_icon(0.9, 0.45, SUN, 0.25))
        spent.add(caution(0.34).next_to(spent[0], UR, buff=-0.08))
        spent.move_to([-6.6, 2.4, 0])
        self.play(FadeIn(spent, shift=RIGHT * 0.3), run_time=0.4)
        self.play(spent.animate.move_to([RX, 2.15, 0]), run_time=0.9)
        self.play(spent.animate.move_to([RX, 1.15, 0]).scale(0.3).set_opacity(0), run_time=0.6)
        self.remove(spent)

        shred_cols = [WATER, COBALT, NICKEL, OCHRE, "#B87333", SUN]
        shreds = VGroup()
        for _ in range(34):
            a, r = rng.uniform(0, TAU), rng.uniform(0.1, 0.55)
            shreds.add(Square(rng.uniform(0.05, 0.1), stroke_width=0).set_fill(rng.choice(shred_cols), 0.95)
                       .rotate(rng.uniform(0, PI)).move_to(WC + r * np.array([np.cos(a), np.sin(a), 0])))
        shreds.add_updater(lambda m, dt: m.rotate(dt * 3.2, about_point=WC))
        steps = VGroup(T("Shredded and dissolved", 14, MUTED, BOLD), T("so each metal can come back out", 14, MUTED, BOLD)) \
            .arrange(DOWN, buff=0.06).next_to(r_lbl, DOWN, buff=0.12)
        self.play(FadeIn(shreds, scale=0.5), FadeIn(steps, shift=UP * 0.1), run_time=0.5)
        self.wait(1.5)

        # element nodes sit on the loop
        E = lambda th: np.array([1.0 + 3.4 * np.cos(th), -0.2 + 2.2 * np.sin(th), 0])
        elems = [("Li", "Lithium", WATER, 150), ("Co", "Cobalt", COBALT, 180), ("Ni", "Nickel", NICKEL, 210)]
        nodes = VGroup()
        for sym, name, col, deg in elems:
            c = Circle(radius=0.4).set_fill(BG, 1).set_stroke(col, 3)
            c2 = Circle(radius=0.4, stroke_width=0).set_fill(col, 0.18)
            n = VGroup(c, c2, H(sym, 24, col, BOLD), T(name, 14, MUTED, BOLD).next_to(c, RIGHT, buff=0.12))
            n[2].move_to(c)
            nodes.add(n.shift(E(deg * DEGREES) - c.get_center()))
        pipe_end = np.array([-2.8, -0.2, 0])
        for n in nodes:
            n.save_state()
            n.scale(0.1).move_to(pipe_end).set_opacity(0)
        self.play(LaggedStart(*[Restore(n) for n in nodes], lag_ratio=0.25), shreds.animate.scale(0.2).set_opacity(0), run_time=1.2)
        self.remove(shreds)

        MX = E(0)
        mfg = VGroup(Circle(radius=0.72).set_fill(STONE_D, 1).set_stroke(LEAF, 3), factory(STONE_L).scale(0.9)).move_to(MX)
        mfg_lbl = VGroup(T("New battery", 18, INK, BOLD), T("factory", 18, INK, BOLD)) \
            .arrange(DOWN, aligned_edge=LEFT, buff=0.04).next_to(mfg, RIGHT, buff=0.2)
        self.play(FadeIn(mfg, scale=0.7), FadeIn(mfg_lbl, shift=LEFT * 0.1), run_time=0.6)

        units, unit_paths, losses, loss_paths = VGroup(), [], VGroup(), []
        bin_c = np.array([-1.5, -2.75, 0])
        for (sym, name, col, deg), n in zip(elems, nodes):
            cen = n[0].get_center()
            for j in range(10):
                d = glow_dot(col, 0.045).move_to(cen)
                if j < 9:
                    ths = np.linspace(deg, 26, 50) * DEGREES
                    units.add(d)
                    unit_paths.append(VMobject().set_points_smoothly([cen] + [E(t) for t in ths]))
                else:
                    losses.add(d)
                    loss_paths.append(VMobject().set_points_smoothly([cen, cen + DOWN * 0.55 + RIGHT * 0.3,
                                                                      bin_c + UP * 0.25 + RIGHT * rng.uniform(-0.2, 0.2)]))
        self.add(units, losses)
        self.bring_to_back(units, losses)
        self.bring_to_back(self.bg)
        self.play(*[Indicate(n[0], color=n[0].get_stroke_color(), scale_factor=1.12) for n in nodes], run_time=0.5)

        top = ParametricFunction(E, t_range=[218 * DEGREES, 22 * DEGREES, -0.01])
        top_line = soft(top, LEAF, 5)
        top_tip = tip_at(top, LEAF, 0.3)
        lead = glow_dot(LEAF, 0.08).move_to(top.get_start())
        self.add(lead)
        self.bring_to_back(top_line)
        self.bring_to_back(self.bg)
        self.play(Create(top_line), MoveAlongPath(lead, top), run_time=1.5)
        rec = pill("About 90% of the metal comes back", LEAF, 15).move_to(E(PI / 2) + UP * 0.42)
        self.play(FadeIn(top_tip, scale=0.5), FadeOut(lead), FadeIn(rec, shift=DOWN * 0.1), run_time=0.5)

        bin_box = VGroup(RoundedRectangle(width=1.4, height=0.36, corner_radius=0.1).set_fill(BG_D, 1)
                         .set_stroke(STONE, 1.2), T("A little is lost", 12, MUTED, BOLD)).move_to(bin_c)
        bin_box[1].move_to(bin_box[0])
        self.play(FadeIn(bin_box), run_time=0.3)
        self.play(LaggedStart(*[MoveAlongPath(u, p, rate_func=smooth) for u, p in zip(units, unit_paths)],
                              lag_ratio=0.025, run_time=3.2),
                  LaggedStart(*[MoveAlongPath(u, p) for u, p in zip(losses, loss_paths)], lag_ratio=0.3, run_time=1.6))
        self.play(FadeOut(units, scale=0.3), FadeOut(losses), Indicate(mfg[0], color=LEAF, scale_factor=1.15), run_time=0.6)

        bot = ParametricFunction(E, t_range=[338 * DEGREES, 236 * DEGREES, -0.01])
        bot_line = soft(bot, STONE_L, 3)
        bot_tip = tip_at(bot, STONE_L, 0.26)
        newcell = battery_icon(0.42, 0.22, WATER).move_to(bot.get_start())
        bot_lbl = T("New batteries head back out into the world", 15, MUTED, BOLD).move_to([1.7, -2.78, 0])
        self.bring_to_back(bot_line)
        self.bring_to_back(self.bg)
        self.play(Create(bot_line), MoveAlongPath(newcell, bot), FadeIn(bot_lbl), run_time=1.3)
        self.play(FadeIn(bot_tip, scale=0.5), FadeOut(newcell), run_time=0.3)

        ph = ValueTracker(0).add_updater(lambda m, dt: m.increment_value(dt * 0.6))
        self.add(ph)
        pulses = VGroup(*[glow_dot(LEAF, 0.05) for _ in range(8)])

        def pul_upd(m):
            for i, d in enumerate(m):
                th = (-ph.get_value() + i * TAU / len(m)) % TAU
                d.move_to(E(th)).set_color(LEAF if np.sin(th) > 0 else STONE_L)
        pulses.add_updater(pul_upd)
        self.add(pulses)

        b_row = VGroup(leaf(LEAF, 0.42), T("Up to 80% less energy than mining new metal", 24, INK, BOLD)).arrange(RIGHT, buff=0.25)
        b_box = RoundedRectangle(width=b_row.width + 0.6, height=0.72, corner_radius=0.2).set_fill(tint(LEAF, 0.14), 1).set_stroke(LEAF, 2)
        banner = VGroup(b_box, b_row).move_to([0.6, -3.38, 0])
        b_row.move_to(b_box)
        self.play(FadeIn(banner, scale=0.85, shift=UP * 0.2), run_time=0.6)
        self.wait(0.8)

        statement = H("We can mine our past instead of the planet", 36, INK, SEMIBOLD, italic=True).move_to([0, 3.3, 0])
        self.play(FadeIn(statement, shift=UP * 0.15), run_time=1.0)
        self.finish()


# ================================================================ SCENE 5
class S5_Closing(Base):
    """1:10-1:25  Presenter full screen with the infinity loop as a lower third overlay.

    Render with --transparent for an overlay that sits on top of the presenter footage.
    """
    DUR = 15

    def construct(self):
        if PLACEHOLDER and not config.transparent:
            head = Circle(radius=1.0).set_fill(STONE_M, 1).set_stroke(width=0).move_to([0, 1.6, 0])
            body = Ellipse(width=6.0, height=3.6).set_fill(STONE_M, 1).set_stroke(width=0).move_to([0, -1.3, 0])
            self.add(head, body, T("Presenter full screen", 15, MUTED, BOLD).move_to([0, 3.5, 0]))
        self.placeholder = [m for m in self.mobjects if m is not self.bg]

        card = RoundedRectangle(width=13.6, height=2.9, corner_radius=0.3).set_fill(BG, 0.92).set_stroke(INK, 1.5, opacity=0.2) \
            .move_to([0, -2.45, 0])
        self.play(FadeIn(card, shift=UP * 0.2), run_time=0.6)

        A, CX, CY = 2.6, -3.6, -2.45
        lem = lambda t: np.array([CX + A * np.cos(t) / (1 + np.sin(t) ** 2),
                                  CY + A * np.sin(t) * np.cos(t) / (1 + np.sin(t) ** 2), 0])

        def grad(m, width=4):
            m.set_stroke([OCHRE, WATER, LEAF], width).set_sheen_direction(RIGHT)
            return m

        line = grad(Line([CX - A, CY, 0], [CX + A, CY, 0]))
        heads = VGroup(*[Triangle().set_fill(interpolate_color(ManimColor(c0), ManimColor(c1), 0.5), 1)
                         .set_stroke(width=0).scale_to_fit_width(0.22).rotate(-PI / 2).move_to([x, CY, 0])
                         for x, c0, c1 in [(CX - A / 2, OCHRE, WATER), (CX + A / 2, WATER, LEAF)]])

        def node(icon, col):
            return VGroup(Circle(radius=0.42).set_fill(PAPER, 1).set_stroke(width=0),
                          Circle(radius=0.42).set_fill(col, 0.15).set_stroke(col, 3.5), icon)
        mine = node(pickaxe(OCHRE).scale(0.5), OCHRE).move_to([CX - A, CY, 0])
        market = node(battery_icon(0.44, 0.24, WATER), WATER).move_to([CX, CY, 0])
        recyc = node(recycle_icon(LEAF, 0.2), LEAF).move_to([CX + A, CY, 0])
        nl = VGroup(T("Mine", 17, OCHRE, BOLD).next_to(mine, DOWN, buff=0.15),
                    T("Market", 17, WATER, BOLD).next_to(market, DOWN, buff=0.15),
                    T("Recycle", 17, LEAF, BOLD).next_to(recyc, DOWN, buff=0.15))
        TX = 3.2
        cap = H("Today most batteries travel a one way road", 22, INK, NORMAL, italic=True).move_to([TX, CY, 0])

        self.play(Create(line), LaggedStart(*[FadeIn(n, scale=0.6) for n in (mine, market, recyc)], lag_ratio=0.25),
                  LaggedStart(*[FadeIn(l, shift=UP * 0.1) for l in nl], lag_ratio=0.25), run_time=1.4)
        self.play(FadeIn(heads), FadeIn(cap, shift=UP * 0.1), run_time=0.5)
        self.wait(0.8)

        curve = grad(ParametricFunction(lem, t_range=[0, TAU, 0.01]))
        halo = curve.copy().set_stroke(width=12, opacity=0.1).set_sheen_direction(RIGHT)
        cap2 = H("Recycling bends that road into a loop", 22, LEAF, NORMAL, italic=True).move_to(cap)
        self.play(ReplacementTransform(line, curve), FadeOut(heads), ReplacementTransform(cap, cap2),
                  nl[0].animate.move_to([CX - A + 0.95, CY, 0]), nl[1].animate.move_to([CX, CY - 0.72, 0]),
                  nl[2].animate.move_to([CX + A - 1.05, CY, 0]), run_time=1.8)
        self.add(halo)
        self.bring_to_back(halo)
        self.bring_to_back(card)
        self.bring_to_back(*self.placeholder)
        self.bring_to_back(self.bg)
        self.add_foreground_mobjects(mine, market, recyc)

        def lem_col(x):
            u = (x - CX) / A
            return interpolate_color(ManimColor(OCHRE), ManimColor(WATER), u + 1) if u < 0 else \
                interpolate_color(ManimColor(WATER), ManimColor(LEAF), u)

        arrows = VGroup()
        for t in (1.25 * PI, 1.75 * PI, 0.25 * PI, 0.75 * PI):
            p, q = lem(t), lem(t + 0.02)
            tri = Triangle().set_fill(lem_col(p[0]), 1).set_stroke(width=0).scale_to_fit_width(0.2)
            arrows.add(tri.rotate(angle_of_vector(q - p) - PI / 2).move_to(p))

        ph = ValueTracker(0).add_updater(lambda m, dt: m.increment_value(dt * 0.7))
        self.add(ph)
        pulses = VGroup(*[glow_dot(WATER, 0.05) for _ in range(14)])

        def pul_upd(m):
            for i, d in enumerate(m):
                p = lem(PI + ph.get_value() + i * TAU / len(m))
                d.move_to(p).set_color(lem_col(p[0]))
        pulses.add_updater(pul_upd)
        pul_upd(pulses)
        self.play(FadeIn(arrows, scale=0.5), FadeIn(pulses), run_time=0.7)
        self.wait(1.5)

        title = VGroup(H("From Mine to Market", 44, INK, SEMIBOLD),
                       H("A sustainable and resilient future for our planet", 20, LEAF, NORMAL, italic=True)) \
            .arrange(DOWN, buff=0.16).move_to([TX, CY, 0])
        self.play(FadeOut(cap2, shift=UP * 0.1), LaggedStart(*[FadeIn(l, shift=UP * 0.15) for l in title], lag_ratio=0.4),
                  run_time=1.4)
        self.finish(fade=False)
