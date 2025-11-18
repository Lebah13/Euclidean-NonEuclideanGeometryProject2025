# Hyperbolic Tiling in the Poincaré Disc

This project visualizes **regular hyperbolic tessellations** $$\({p, q}\)$$ in the **Poincaré disc model** of the hyperbolic plane.

- Each face is a **regular hyperbolic polygon** with $$\(p\)$$ sides.
- Exactly $$\(q\)$$ polygons meet at each vertex.
- The condition $$\((p - 2)(q - 2) > 4\)$$ ensures the geometry is hyperbolic.

The output is a picture of a finite patch of the hyperbolic plane drawn inside the unit disc, where the boundary circle represents “points at infinity”.

---

## Mathematical idea

We work in the **Poincaré disc model**:

$$
\[
\mathbb{D} = \{ z \in \mathbb{C} : |z| < 1 \}
\]
$$
with hyperbolic metric
$$
\[
ds = \frac{2|dz|}{1 - |z|^2}.
\]
$$

Key facts used in the code:

- **Geodesics** (hyperbolic straight lines) in the disc are:
  - circular arcs orthogonal to the unit circle, or  
  - diameters through the origin.
- A regular $$\(\{p, q\}\)$$ tiling is built from a **regular hyperbolic $$\(p\)$$-gon** where $$\(q\)$$such polygons meet at each vertex.
- Reflections in geodesics generate a discrete group of **hyperbolic isometries** that tile the plane.

The script:

1. Constructs a regular hyperbolic $$\(\{p, q\}\)$$ polygon centered at the origin.
2. Uses **reflections across its sides** (hyperbolic isometries) to generate neighboring tiles.
3. Repeats this process up to a chosen depth (a small breadth-first search on the tiling graph).
4. Draws all polygon edges as **true hyperbolic geodesics** (orthogonal circle arcs) inside the unit disc using Matplotlib.

---

## Repository contents

- `hyperbolic_tiling.py`  
  Main script that:
  - Defines the hyperbolic geometry utilities,
  - Generates the $$\(\{p, q\}\)$$ tiling,
  - Plots the result.

If GitHub Actions is enabled (see below), a generated image
(e.g. `tiling.png`) can be produced automatically on each push.

---

## How to run locally

### Requirements

- Python 3.10+ (3.11 is fine)
- `matplotlib`

Install the dependency:

```bash
pip install matplotlib
