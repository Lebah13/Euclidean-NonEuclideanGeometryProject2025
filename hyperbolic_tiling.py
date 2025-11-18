import cmath
import math
import matplotlib.pyplot as plt
from collections import deque


def upper_to_disc(z: complex) -> complex:
 
    if z.imag <= 0:
        raise ValueError("upper_to_disc expects Im(z) > 0")
    return (z - 1j) / (z + 1j)


def disc_to_upper(w: complex) -> complex:
   
    if abs(w) >= 1:
        raise ValueError("disc_to_upper expects |w| < 1")
    return 1j * (1 + w) / (1 - w)

def _geodesic_circle_center(z1: complex, z2: complex):
   
    eps = 1e-9
    if abs(z1) < eps or abs(z2) < eps:
        return None
    if abs((z1 / z2).imag) < eps:
        return None

    a, b = z1.real, z1.imag
    c, d = z2.real, z2.imag
    u, v = (a - c), (b - d)

    A1, B1, C1 = 2 * a, 2 * b, 1 + a * a + b * b
    A2, B2, C2 = 2 * u, 2 * v, a * a + b * b - c * c - d * d

    det = A1 * B2 - A2 * B1
    if abs(det) < 1e-12:
        return None

    x = (C1 * B2 - C2 * B1) / det
    y = (A1 * C2 - A2 * C1) / det
    center = complex(x, y)
    R = abs(z1 - center)
    return center, R


def geodesic_points(z1: complex, z2: complex, n: int = 100):
  
    eps = 1e-9

    if abs(z1) < eps or abs(z2) < eps or abs((z1 / z2).imag) < eps:
        return [z1 + (z2 - z1) * (i / (n - 1)) for i in range(n)]

    res = _geodesic_circle_center(z1, z2)
    if res is None:
        return [z1 + (z2 - z1) * (i / (n - 1)) for i in range(n)]

    center, R = res
    t1 = cmath.phase(z1 - center)
    t2 = cmath.phase(z2 - center)

    delta = (t2 - t1 + math.pi) % (2 * math.pi) - math.pi

    points = [
        center + R * cmath.exp(1j * (t1 + delta * (i / (n - 1))))
        for i in range(n)
    ]
    return points

def reflect_in_geodesic(z: complex, z1: complex, z2: complex) -> complex:
   
    res = _geodesic_circle_center(z1, z2)
    if res is None:

        angle = cmath.phase(z1)
        rot = cmath.exp(-1j * angle)
        z_rot = z * rot
        z_ref = z_rot.conjugate()
        return z_ref * cmath.exp(1j * angle)

    center, R = res
    w = z - center
    if abs(w) < 1e-15:
        return z  
    w_inv = (R * R) / w.conjugate()
    return center + w_inv


def reflect_polygon_across_side(vertices, k):
    
    n = len(vertices)
    z1 = vertices[k]
    z2 = vertices[(k + 1) % n]
    return [reflect_in_geodesic(z, z1, z2) for z in vertices]

def regular_hyperbolic_polygon(p: int, q: int):
   
    if p < 3 or q < 3:
        raise ValueError("p and q must be >= 3")
    if (p - 2) * (q - 2) <= 4:
        raise ValueError("This {p,q} is not hyperbolic (needs (p-2)(q-2) > 4)")

    A = math.tan(math.pi / 2 - math.pi / q)
    B = math.tan(math.pi / p)
    d_sq = (A - B) / (A + B)
    if d_sq <= 0:
        raise ValueError("Radius squared non-positive; check p,q")
    d = math.sqrt(d_sq)

    vertices = [
        d * cmath.exp(2j * math.pi * k / p)
        for k in range(p)
    ]
    return vertices


def polygon_centroid(vertices):
    return sum(vertices) / len(vertices)

def generate_tiling(p: int, q: int, max_depth: int = 2):
    
    root = regular_hyperbolic_polygon(p, q)
    polys = []
    queue = deque()

    def key(verts):
        c = polygon_centroid(verts)
        return (round(c.real, 3), round(c.imag, 3))

    seen = set()
    seen.add(key(root))
    queue.append((root, 0))

    while queue:
        poly, depth = queue.popleft()
        polys.append(poly)

        if depth >= max_depth:
            continue

        n = len(poly)
        for k in range(n):
            new_poly = reflect_polygon_across_side(poly, k)
            c = polygon_centroid(new_poly)

            if abs(c) > 0.999:
                continue

            k_new = key(new_poly)
            if k_new not in seen:
                seen.add(k_new)
                queue.append((new_poly, depth + 1))

    return polys

def plot_tiling(polygons, p, q):
    fig, ax = plt.subplots()

    t_vals = [2 * math.pi * i / 400 for i in range(401)]
    circle_x = [math.cos(t) for t in t_vals]
    circle_y = [math.sin(t) for t in t_vals]
    ax.plot(circle_x, circle_y)

    for poly in polygons:
        n = len(poly)
        for k in range(n):
            z1 = poly[k]
            z2 = poly[(k + 1) % n]
            pts = geodesic_points(z1, z2, n=80)
            xs = [pnt.real for pnt in pts]
            ys = [pnt.imag for pnt in pts]
            ax.plot(xs, ys, linewidth=0.8)

    ax.set_aspect('equal', 'box')
    ax.set_title(f"Hyperbolic {{{p},{q}}} tiling in the Poincaré disc")

    plt.show()

if __name__ == "__main__":
    # Choose your (p, q): must satisfy (p-2)(q-2) > 4 for hyperbolic tiling
    p = 4
    q = 8
    depth = 3   # 2 or 3 is a good starting point

    print(f"Generating hyperbolic tiling for {{{p},{q}}} with depth {depth}...")
    polygons = generate_tiling(p, q, max_depth=depth)
    print(f"Generated {len(polygons)} polygons.")

    # In CI (GitHub Actions sets CI=true), save to file instead of opening a window
    running_in_ci = os.getenv("CI") == "true"

    if running_in_ci:
        fig, ax = plt.subplots()
        # quick reuse of the plotting logic
        t_vals = [2 * math.pi * i / 400 for i in range(401)]
        circle_x = [math.cos(t) for t in t_vals]
        circle_y = [math.sin(t) for t in t_vals]
        ax.plot(circle_x, circle_y)

        for poly in polygons:
            n = len(poly)
            for k in range(n):
                z1 = poly[k]
                z2 = poly[(k + 1) % n]
                pts = geodesic_points(z1, z2, n=80)
                xs = [pnt.real for pnt in pts]
                ys = [pnt.imag for pnt in pts]
                ax.plot(xs, ys, linewidth=0.8)

        ax.set_aspect('equal', 'box')
        ax.set_title(f"Hyperbolic {{{p},{q}}} tiling in the Poincaré disc")
        plt.savefig("tiling.png", dpi=200, bbox_inches="tight")
        print("Saved tiling.png (for CI artifact).")
    else:
        plot_tiling(polygons, p, q)
