"""
make_example_gds.py -- Genere un GDS de demonstration respectant la
convention de layers attendue par le solveur.

    layer 1 : metal supraconducteur (plan de masse fendu, ligne CPW)
    layer 2 : vue de dessus des airbridges
    layer 3 : polygones d'injection du courant, un calcul par polygone,
              nommes par un label GDS (ici "IN_A", "IN_B")
    layer 4 : disques de retour vers la masse
    layer 5 : surfaces ou l'on veut le flux magnetique, nommees par un label
    layer 6 : rectangles des jonctions Josephson, avec un label "EJ=20GHz"

Usage :  python -m scdc.make_example_gds example.gds
"""

import sys
import numpy as np
import gdstk


def disk(cx, cy, r, layer, n=96):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return gdstk.Polygon(np.column_stack([cx + r * np.cos(t),
                                          cy + r * np.sin(t)]), layer=layer)


def rect(x0, y0, x1, y1, layer):
    return gdstk.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)], layer=layer)


def build(path="example.gds"):
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = lib.new_cell("TOP")

    # --- layer 1 : 200 x 120 um ground plane with a CPW slot ---
    plane = gdstk.Polygon([(0, 0), (200, 0), (200, 120), (0, 120)], layer=1)
    gaps = []
    for sgn in (-1, +1):
        gaps.append(rect(30, 60 + sgn * 6 - 4, 170, 60 + sgn * 6 + 4, 1))
    metal = gdstk.boolean(plane, gaps, "not", layer=1)

    # two circular holes in the ground plane, and a cut of the CPW centre
    # conductor that the junction closes
    holes = [disk(70, 100, 8, 1), disk(140, 25, 10, 1),
             rect(100, 57.9, 101, 62.1, 1)]
    metal = gdstk.boolean(metal, holes, "not", layer=1)
    for p in metal:
        cell.add(p)

    # --- layer 2 : two airbridges across the lower slot ---
    for x in (80, 120):
        cell.add(rect(x - 5, 40, x + 5, 72, 2))

    # --- layer 3 : two injection polygons, one computation each ---
    cell.add(disk(15, 60, 5, 3))
    cell.add(gdstk.Label("IN_A", (15, 60), layer=3))
    cell.add(disk(40, 105, 5, 3))
    cell.add(gdstk.Label("IN_B", (40, 105), layer=3))

    # --- layer 4 : ground returns ---
    for cx, cy in [(190, 20), (190, 100), (105, 15)]:
        cell.add(disk(cx, cy, 5, 4))

    # --- layer 6 : junction closing the cut, E_J given by a label ---
    cell.add(rect(100, 58, 101, 62, 6))
    cell.add(gdstk.Label("EJ=20GHz", (100.5, 60), layer=6))

    # --- layer 5 : flux surfaces, named by labels ---
    cell.add(disk(70, 100, 8, 5))            # ground-plane hole
    cell.add(gdstk.Label("HOLE", (70, 100), layer=5))
    cell.add(rect(95, 45, 115, 75, 5))       # under an airbridge
    cell.add(gdstk.Label("BRIDGE_AREA", (105, 60), layer=5))
    # the third surface is left unnamed on purpose, it gets the name L<k>

    cell.add(disk(40, 30, 12, 5))

    lib.write_gds(path)
    print("written:", path)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "example.gds")
