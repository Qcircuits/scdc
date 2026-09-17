"""
near_pairs.py -- Recherche des paires proches sans listes Python.

query_ball_point renvoie une liste de listes d'entiers Python, soit environ
36 octets par paire (objet int + pointeur). Sur un maillage gradue ou un
grand triangle voisin d'une zone fine a plusieurs dizaines de milliers de
voisins, cela atteint des dizaines de Go. Ici :

  1. comptage des paires (return_length=True, aucune liste) et estimation
     memoire affichee avant toute allocation ;
  2. extraction par classes de rayon avec sparse_distance_matrix en sortie
     ndarray (24 octets par paire, tableaux numpy contigus), classe par classe.
"""
import numpy as np
from scipy.spatial import cKDTree


def count_near_pairs(tree, C3, rad):
    """Nombre de paires (i, j) avec |C_i - C_j| <= rad_i, sans symetrisation."""
    try:
        lens = tree.query_ball_point(C3, r=rad, workers=-1, return_length=True)
    except TypeError:
        lens = tree.query_ball_point(C3, r=rad, return_length=True)
    return np.asarray(lens, dtype=np.int64)


def near_pairs(C3, rad, verbose=True, nbins=8):
    """Paires symetrisees : (i, j) retenu si d <= rad_i ou d <= rad_j.

    Retourne (I, J) int64, incluant la diagonale, ordonnes.
    """
    n = len(C3)
    tree = cKDTree(C3)
    lens = count_near_pairs(tree, C3, rad)
    n_dir = int(lens.sum())
    if verbose:
        print(f"  near pairs: {n_dir} before symmetrisation "
              f"({n_dir / n:.0f}/triangle, max {lens.max()} for one triangle), "
              f"working memory ~{n_dir * 40 / 1e9:.1f} GB")
    # classes de rayon en progression geometrique
    r_lo, r_hi = float(rad.min()), float(rad.max())
    if r_hi / max(r_lo, 1e-300) < 1.001:
        edges = np.array([r_lo, r_hi * (1 + 1e-9)])
    else:
        edges = np.geomspace(r_lo, r_hi * (1 + 1e-9), nbins + 1)
    cls = np.clip(np.searchsorted(edges, rad, side='right') - 1, 0, len(edges) - 2)
    keys = []
    for c in range(len(edges) - 1):
        sel = np.where(cls == c)[0]
        if len(sel) == 0:
            continue
        sub = cKDTree(C3[sel])
        rmax = float(edges[c + 1])
        m = sub.sparse_distance_matrix(tree, rmax, output_type='ndarray')
        i = sel[m['i']]
        j = m['j'].astype(np.int64)
        keep = m['v'] <= rad[i]
        i, j = i[keep], j[keep]
        del m, keep
        # cle non ordonnee (min, max) -> symetrie gratuite
        a = np.minimum(i, j)
        b = np.maximum(i, j)
        keys.append(a * n + b)
        del i, j, a, b
    key = np.unique(np.concatenate(keys))
    del keys
    # diagonale garantie
    diag = np.arange(n, dtype=np.int64) * (n + 1)
    key = np.union1d(key, diag)
    a = key // n
    b = key % n
    off = a != b
    I = np.concatenate([a, b[off]])
    J = np.concatenate([b, a[off]])
    del key, a, b, off
    order = np.lexsort((J, I))
    return I[order], J[order]


if __name__ == "__main__":
    # test contre la version originale (listes) sur un nuage gradue
    rng = np.random.default_rng(1)
    n = 20000
    C3 = np.column_stack([rng.uniform(0, 500, n), rng.uniform(0, 500, n),
                          rng.choice([0.0, 3.0], n)])
    h = np.where(C3[:, 0] < 100, 0.5, 8.0) * rng.uniform(0.5, 1.5, n)
    rad = 4 * np.maximum(2.0, h)
    tree = cKDTree(C3)
    lists = tree.query_ball_point(C3, r=rad)
    I0 = np.repeat(np.arange(n), [len(l) for l in lists])
    J0 = np.concatenate([np.asarray(l) for l in lists])
    key = np.unique(np.concatenate([I0 * n + J0, J0 * n + I0]))
    Iref, Jref = key // n, key % n
    I, J = near_pairs(C3, rad)
    assert len(I) == len(Iref), (len(I), len(Iref))
    assert np.array_equal(I, Iref) and np.array_equal(J, Jref)
    print("near_pairs identical to the list-based version:", len(I), "pairs")
