# scdc — densité de courant DC et inductances mutuelles dans un circuit supraconducteur

Entrée : un GDS. Sortie : une carte de courant par polygone d'injection du
layer 3, la matrice des inductances mutuelles `M_ij = Φ_i / I_j` reliant le
flux à travers chaque surface du layer 5 au courant injecté dans chaque
polygone du layer 3, et un rapport Markdown. Tous les fichiers produits et
messages du programme sont en anglais.

Toutes les longueurs sont en micron, les courants en ampère, les flux en weber.

---

## 1. Convention de layers

| layer | contenu | rôle |
|---|---|---|
| 1 | toutes les structures métalliques | film supraconducteur, plan `z = 0` |
| 2 | vue de dessus des airbridges | feuilles surélevées à `z = h` |
| 3 | un ou plusieurs polygones, disques ou rectangles | injection du courant `I`, un calcul par polygone |
| 4 | des polygones | retours vers la masse, chacun via `R = 0.01 Ω` |
| 5 | des surfaces | contours où l'on calcule `Φ` puis `M = Φ/I` |

Les polygones des layers 3 et 5 peuvent recevoir un label GDS (objet TEXT
placé à l'intérieur du polygone). Ce texte sert de nom dans le rapport, dans
la matrice `M`, dans les noms de fichiers et dans le titre des cartes. À
défaut le nom est `S0`, `S1`, … pour les sources et `L0`, `L1`, … pour les
surfaces, dans l'ordre des centroïdes croissants en `x` puis `y`.

### Plusieurs sources et matrice `M`

Chaque polygone `j` du layer 3 donne lieu à un calcul indépendant, comme
s'il était seul, avec le courant `I` injecté dans ce polygone et retourné
par les polygones du layer 4 (répartition `I_k = I (1/R_k)/Σ 1/R`). Les
autres polygones du layer 3 ne sont alors que du métal. Le noyau pFFT ne
dépend que de la géométrie et n'est construit qu'une fois. On obtient

```
M_ij = Φ_i^{(j)} / I ,     i = surface du layer 5,  j = polygone du layer 3
```

Si `M` est carrée et régulière, le rapport donne aussi les valeurs propres
`λ_k` de `M⁻¹`, exprimées en `mA/Φ0` par `λ_k Φ0 × 10³` avec `Φ0 = h/2e`. La
valeur propre `λ_k Φ0` est le courant, le long de la `k`-ième direction
propre, qui produit un quantum de flux dans la combinaison de surfaces
correspondante. `M` n'étant pas symétrique en général (surfaces et sources
sont des objets distincts), des valeurs propres complexes sont possibles et
sont alors imprimées avec leur partie imaginaire. Si `M` n'est pas carrée le
rapport le signale et donne à la place les inverses des valeurs singulières
de `M`, qui coïncident avec les modules des valeurs propres quand `M` est
carrée et normale.
| 6 | des rectangles | jonctions Josephson, patchs à forte inductance cinétique |

---

## 2. Modèle physique

### 2.1 Énergie

Chaque film est traité comme une nappe de courant `K` (A/µm) à l'altitude `z_s`.
L'énergie totale est la somme d'un terme cinétique et d'un terme magnétique

```
E[K] = (µ0 Λ / 2) ∫ |K|² dA  +  (µ0 / 8π) ∫∫ K(r)·K(r') / |r − r'| dA dA'
```

avec la longueur de pénétration effective

```
Λ = λ coth(d/λ)      →  λ²/d  quand d ≪ λ,   →  λ  quand d ≫ λ
```

La longueur de Pearl correspondante vaut `2Λ`. Les deux paramètres demandés,
épaisseur `d` et longueur de London `λ`, entrent uniquement par `Λ`, ce qui
est la signature du modèle de London bidimensionnel.

L'inductance cinétique par carré est `L_□ = µ0 Λ`, soit `µ0 λ²/d` pour un film
mince, ou encore `ħ R_□ / (π Δ)` en limite sale. On peut la donner directement
avec `--Lsq` (en pH par carré) à la place de `d` et `λ`, ce qui est en général
la grandeur mesurée.

### 2.2 Pourquoi une minimisation d'énergie

En régime DC établi, `∂(Λ J)/∂t = E = 0` dans le supraconducteur. Le film est
donc **équipotentiel** et le réseau est purement inductif. Deux conséquences.

- Les disques du layer 4 sont tous au même potentiel. Avec des résistances
  identiques `R_k = 0.01 Ω` vers la masse, le courant se répartit à parts
  **exactement égales**, `I_k = I/N`. Plus généralement le code utilise
  `I_k = I (1/R_k) / Σ_j (1/R_j)`, si bien que la valeur commune de `R`
  n'influence pas le résultat. Seul un déséquilibre entre résistances le ferait.
- La distribution interne n'est pas fixée par le régime DC seul, elle dépend de
  l'histoire. Pour un refroidissement en champ nul suivi d'une montée du courant
  de 0 à `I`, la distribution atteinte est l'unique minimum de `E` sous les
  contraintes de conservation et de quantification du fluxoïde avec `n = 0`
  dans chaque trou. C'est cette solution qui est calculée.

La stationnarité de `E` redonne l'équation de London intégrée sur l'épaisseur

```
µ0 Λ K(r) + A(r) = − (ħ/2e) ∇θ
```

et, autour de chaque trou `C_k`,

```
∮_{C_k} (µ0 Λ K + A)·dl = n_k Φ0 ,   n_k = 0 par défaut
```

Le paramètre `fluxoid={indice_de_trou: n}` permet d'imposer du flux piégé.

### 2.3 Discrétisation

`K` est constant par triangle, et décomposé en

```
K = K_seed  +  Σ_b α_b F_b  +  ẑ × ∇g
```

- `K_seed` est une distribution admissible quelconque portant l'injection aux
  terminaux, obtenue par un problème résistif auxiliaire `∇·(σ∇V) = 0`.
- `F_b` sont les champs de boucle associés aux airbridges. Leur nombre est le
  rang cyclique du graphe dont les sommets sont les composantes connexes de
  métal et de ponts, et les arêtes les pieds de ponts. Un pont qui court-circuite
  deux points d'un même plan crée un degré de liberté. Un pont qui est le seul
  chemin entre deux plans disjoints n'en crée pas.
- `g` est une fonction de courant P1 nodale, constante sur chaque contour, avec
  une constante libre par trou et la jauge `g = 0` sur le contour extérieur de
  chaque composante connexe.

L'énergie devient une forme quadratique et le minimum est obtenu par un système
linéaire dense de taille `ndof(g) + nombre de boucles`.

Le noyau `G_{TT'} = ∫_T ∫_{T'} dA dA' / |r − r'|` est évalué par une quadrature
à trois points par triangle hors diagonale. Le terme diagonal utilise la valeur
exacte pour le disque équivalent d'aire `A`,

```
∫∫ dA dA' / |r − r'| = (16π/3) a³ ,   a = √(A/π)
```

qui découle de la distance inverse moyenne `⟨1/r⟩ = 16/(3πa)` entre deux points
d'un disque.

### 2.4 Jonctions Josephson

Pour `I ≪ I_c` une jonction est une inductance linéaire

```
L_J = Φ0 / (2π I_c) = (Φ0/2π)² / E_J ,      E_J = Φ0 I_c / (2π)
```

Un rectangle de longueur `ℓ` selon le courant et de largeur `w`, de sheet
inductance `L_□^J`, a `L = L_□^J ℓ / w` si le courant y est uniforme, ce qui est
assuré dès que le terme cinétique domine. Chaque rectangle du layer 6 reçoit
donc

```
Λ_J = L_J w / (µ0 ℓ)
```

En pratique le patch n'est pas toujours un rectangle propre, ses extrémités
peuvent aboutir sur des côtés obliques ou déborder sur les électrodes. Le code
ne suppose donc rien sur la forme. La zone active est la partie du polygone qui
dépasse des électrodes (tout le polygone s'il est entièrement sur du métal).
Après maillage, les nœuds du patch partagés avec le métal extérieur sont
regroupés en zones de contact, les deux plus grandes sont les électrodes A et B,
et l'on résout Laplace dans le patch avec `V = 0` sur A et `V = 1` sur B. Alors

```
G = ∫ |∇V|² dA ,     F = L / L_□ = 1/G ,     Λ_J = L_J / (µ0 F) ,
I_J = ∫_J K · ∇V dA
```

Pour un rectangle `ℓ × w` on retrouve `F = ℓ/w` exactement, et `I_J` est exact
pour tout `K` conservatif dans le patch. Si le patch touche plus de deux zones
de métal, un avertissement l'indique.

Le rectangle doit être la seule connexion entre ses deux électrodes. Un patch
plus étroit que la piste qu'il traverse est court-circuité par le film
alentour, c'est le comportement physique d'une inductance en parallèle d'un
supraconducteur, et le test de validation le montre.

Les paramètres sont fournis, par ordre de priorité, par un label GDS dans le
rectangle (`EJ=20GHz`, `Ic=40nA`, `LJ=8nH`, `Rn=8kOhm Delta=180ueV`). Il faut
un objet TEXT du GDS, dans KLayout un objet Texte, et non un texte converti en
polygones, qui n'est qu'un dessin. À défaut, par un
fichier JSON `--junction-file` de la forme
`[{"x": 100.5, "y": 60, "EJ_GHz": 20}, ...]`, ou par une valeur par défaut
`--EJ`, `--Ic`, `--LJ`. La relation d'Ambegaokar-Baratoff `I_c R_n = πΔ/(2e)`
sert pour `Rn`.

`--nonlinear` itère `L_J(I) = L_J0 / √(1 − (I/I_c)²)` jusqu'à convergence sur le
courant de chaque jonction. Le noyau pFFT ne dépend que de la géométrie, il est
construit une fois et réutilisé. Le résultat donne `I_J`, `I_J/I_c` et la phase
`φ = arcsin(I_J/I_c)`.

Le régime de validité est celui que vous décrivez, `I_J` petit devant `I_c`,
pas de tension continue sur la jonction. Le modèle ne contient ni la
capacité de la jonction ni sa dynamique, ce qui est cohérent avec un calcul DC.

### 2.5 Ponts, pieds et piliers

Le layer 2 est la vue de dessus du pont. Ses pieds sont les composantes
connexes de `pont ∩ métal`. Un pont qui enjambe une CPW en recouvre trois,
les deux masses et le conducteur central. Par défaut (`--bridge-feet ends`)
seules les deux composantes les plus éloignées l'une de l'autre sont des
pieds, les autres sont survolées sans contact. `--bridge-feet all` rétablit
le comportement où toute zone de recouvrement est un contact (vias multiples).

Chaque pied est relié au plan par un pilier vertical, parallélépipède de base
carrée de côté `a` (`--pillar-side`, défaut 30 µm) et de hauteur `h`
(`--pillar-height`, défaut `--bridge-height`). Un courant vertical ne produit
pas de `B_z` et ne se couple pas aux nappes horizontales, le produit scalaire
`J·J'` du noyau de Neumann étant nul. Les piliers n'entrent donc que par leur
matrice d'inductances propres et mutuelles `L_p`, et l'énergie reçoit le terme

```
E_p = (1/2) Iᵀ L_p I ,     I_k = ∮ courant montant dans le pilier k
```

Inductance propre d'un conducteur droit à courant de surface, `R = 0.5823 a`
étant la distance moyenne géométrique du contour carré avec lui-même
(Rosa 1908, Grover 1946)

```
L = (µ0 h / 2π) [ ln(2h/R) − 1 + R/h ]  +  µ0 λ h / (4a)
```

Mutuelle entre deux piliers parallèles d'axes distants de `d`

```
M = (µ0 h / 2π) [ ln(h/d + √(1 + h²/d²)) − √(1 + d²/h²) + d/h ]
```

Pour `a = 30 µm`, `h = 200 µm`, `L = 88.9 pH` et `M(100 µm) = 33.0 pH`, soit
`L₁ + L₂ − 2M = 111.8 pH` pour l'aller-retour d'un pont. Ces formules sont
valables pour `h ≫ a` et `d ≫ a`, l'erreur relative étant `O((a/h)²)` et
`O((a/d)²)`.

Le courant `I_k` est la divergence faible de `K` sommée sur les nœuds du pied
(`solver.contact_current_operator`), fonctionnelle linéaire de `K` qui ne
dépend que du seed et des champs de boucle. Le terme `Dᵀ L_p D` est de rang
`n_contacts`, il est ajouté à l'opérateur du gradient conjugué et au bloc des
boucles du préconditionneur. `run.py` imprime le courant et l'inductance
propre de chaque pilier. Le test avec deux piliers artificiels retrouve
`L₁ + L₂ − 2M` à 10⁻⁴ par identification sur un modèle à constantes
localisées. `--no-pillars` revient aux contacts parfaits.

### 2.6 Flux

`Φ` est obtenu par le théorème de Stokes, `Φ = ∮_{∂S} A·dl`, plutôt que par
`∫_S B_z dA`. Le potentiel vecteur d'une nappe de courant est continu partout,
y compris dans le plan du film, alors que `B_z` y est discontinu. Cette
formulation reste donc valable même pour une surface du layer 5 posée dans le
plan du métal.

```
A(r) = (µ0/4π) Σ_T A_T K_T / |r − r_T|
```

avec raffinement par quadrature pour les triangles proches. Les trous éventuels
des polygones du layer 5 sont traités par orientation des contours intérieurs.
La normale est `+ẑ` et le contour extérieur est parcouru dans le sens direct.

---

## 3. Installation

```bash
pip install numpy scipy matplotlib shapely gdstk triangle
```

`triangle` est l'interface Python du mailleur de Shewchuk. Le noyau physique
(`core.py`, `solver.py`) ne dépend que de numpy et scipy, ce qui permet de le
tester sans la chaîne GDS.

## 4. Utilisation

Deux arborescences sont acceptees.

**A plat**, tous les `.py` dans le meme dossier, on lance les scripts
directement.

```bash
cd dossier_contenant_les_py
python make_example_gds.py exemple.gds
python run.py exemple.gds --thickness 0.1 --lambda-london 0.09 --seg 1.5 \
              --current 1e-3 --bridge-height 3.0 --out resultat
python validation.py
```

**En package**, les `.py` dans un sous-dossier `scdc/` contenant `__init__.py`,
et on lance depuis le dossier parent.

```
projet/
  scdc/
    __init__.py  core.py  solver.py  geometry.py  run.py  make_example_gds.py
  validation.py
```

```bash
cd projet
python -m scdc.make_example_gds exemple.gds
python -m scdc.run exemple.gds --thickness 0.1 --lambda-london 0.09 --seg 1.5 \
       --current 1e-3 --bridge-height 3.0 --out resultat
python -m scdc.run exemple.gds --Lsq 0.14 --seg 1.5 --EJ 20 --nonlinear
python validation.py
```

L'erreur `No module named 'scdc'` signifie que le dossier `scdc/` n'existe pas,
ou que la commande n'est pas lancee depuis son parent. Dans ce cas utiliser la
forme a plat.

Sorties (préfixe `--out`, ici `resultat`).

- `resultat_report.md`, rapport Markdown commençant par la matrice `M` en pH
  et les valeurs propres de `M⁻¹` en mA/Φ0, puis paramètres, maillage,
  énergie et inductance `2E/I²` par source, flux `Φ/Φ0` par surface et par
  source, courants des piliers et des jonctions.
- `resultat_current_<k>_<nom>.png`, une carte de `|K|` avec lignes de courant
  par polygone du layer 3.
- `resultat_mutual.csv`, une ligne par surface du layer 5, une colonne
  `M_pH_<source>` par polygone du layer 3.
- `resultat_solution.npz`, maillage complet, `K` de forme
  `(n_sources, n_tri, 2)`, `M`, flux et noms des sources et surfaces.

En bibliothèque.

```python
from scdc.run import run
model, sols, res = run("exemple.gds", thickness=0.1, lambda_L=0.09,
                       seg_len=1.5, current=1e-3)
print(res["M"] * 1e12)            # matrice (n_surfaces, n_sources) en pH
print(res["inv_eigenvalues"])     # valeurs propres de M^-1 en mA/Phi0
J = sols[0].current_density(0.1)  # A/µm² pour la première source
```

`run.solve_model(model, ...)` accepte un dictionnaire construit à la main
(maillage, `terminal_sets`, `rings`, …) et permet de tester la boucle sur les
sources et le rapport sans la chaîne GDS.

Options utiles.

| option | effet |
|---|---|
| `--seg` | longueur des segments du maillage, contrôle coût et précision |
| `--bridge-height` | altitude des airbridges |
| `--bridge-thickness`, `--bridge-lambda` | propriétés distinctes pour les ponts |
| `--bridge-feet` | `ends` (défaut, deux pieds extrêmes) ou `all` |
| `--pillar-side`, `--pillar-height`, `--pillar-lambda`, `--no-pillars` | piliers des ponts |
| `--loop-z` | altitude des surfaces du layer 5 |
| `--ground-resistances` | liste, utile seulement si les `R_k` diffèrent |
| `--scale` | conversion vers le micron si le GDS n'est pas en µm |
| `--Lsq` | inductance par carré du layer 1 en pH, remplace `--thickness --lambda-london` |
| `--EJ`, `--Ic`, `--LJ`, `--junction-file` | paramètres des jonctions du layer 6 |
| `--nonlinear` | itération `L_J(I)` |
| `--show`, `--log` | fenêtre interactive, échelle logarithmique |
| `--dense` | assemblage dense de référence |
| `--near-cells` | rayon de précorrection pFFT, défaut 4 |
| `--backend` | `numpy`, `mlx`, `torch` ou `auto` |
| `--seg-far`, `--fine-radius`, `--seg-max` | gradation du maillage |

### Maillage gradué

Sur une puce de 10 mm, le layer 1 peut totaliser plus d'un mètre de périmètre,
et un maillage uniforme à 1,5 µm dépasserait `10^8` triangles. Le maillage est
donc gradué. Les bords sont discrétisés à `--seg` à moins de `--fine-radius` des
surfaces du layer 5 et des jonctions, à `--seg-far` au-delà (défaut `5*seg`),
et Triangle laisse les éléments grossir vers l'intérieur jusqu'à `--seg-max`
(défaut `30*seg`). Le courant étant concentré sur les bords, c'est là que la
résolution compte. Les triangles intérieurs, où `K` est quasi nul par
écrantage, peuvent être grands.

Le champ proche de la pFFT s'adapte à cette gradation. Une paire de triangles
est traitée exactement si sa distance est inférieure à `near_cells` fois le
maximum de la maille de grille et des tailles des deux triangles. Sur un
maillage de rapport 21 entre plus grand et plus petit élément, l'écart avec
l'assemblage dense est de 0,04 % en énergie.

### Coût et méthode rapide

Par défaut le calcul utilise `fast.py`, une méthode des moments accélérée du
type de celle de FastHenry. FastHenry accélère le produit matrice-vecteur par
un développement multipolaire. Ici, les films étant plans, on utilise la FFT
précorrigée de Phillips et White, mieux adaptée à des nappes. Le principe.

1. Les sources `A_T K_T` de chaque triangle sont projetées sur une grille
   régulière par leurs trois points de quadrature et des poids bilinéaires,
   opérateur creux `S`. Le noyau `1/r` est une convolution sur la grille,
   donc une FFT, et le retour aux triangles utilise `S^T`. L'opérateur
   `S K S^T` est symétrique.
2. Pour les paires de triangles à moins de `near_cells` mailles, on retranche
   la valeur grille et on ajoute la valeur exacte. Cette correction est une
   matrice creuse.
3. L'énergie est minimisée par gradient conjugué. Le préconditionneur est la
   factorisation LU creuse de la partie locale de l'opérateur (terme cinétique
   et champ proche), analogue de la matrice creusifiée de FastHenry.

Coût par itération `O(N log N)`, mémoire `O(N)`, typiquement 10 à 40
itérations. Sur un portable, 30 000 triangles prennent 7 s et 120 000
triangles 30 s, là où l'assemblage dense demanderait 115 Go. L'écart avec la
référence dense est de l'ordre de 0,1 % sur l'énergie, réglable par
`--near-cells`.

L'option `--dense` conserve l'assemblage direct de `solver.py` pour vérifier
sur de petits maillages.

Raffiner près des bords reste ce qui change le plus le résultat, puisque `K`
diverge en `1/√s` à distance `s` d'un bord quand `Λ ≪ W`.

### Accélération sur Apple Silicon

Le profil à 120 000 triangles est, en ordre de grandeur, 20 s de construction
du noyau (rassemblements de paires proches), 5 s de factorisation du
préconditionneur et 8 s de gradient conjugué dominé par les FFT. Trois niveaux,
par confiance décroissante.

1. **numpy multithread, par défaut.** Les FFT passent par `scipy.fft` avec
   `workers` égal au nombre de cœurs, et les calculs de paires sont découpés en
   blocs traités par un pool de threads, numpy libérant le GIL sur ces
   opérations. Les produits denses utilisent Accelerate, le BLAS d'Apple, qui
   est celui des roues numpy et scipy pour macOS arm64 depuis numpy 2.0.
   Vérification avec `python -c "import numpy; numpy.show_config()"`, on doit
   lire `accelerate`.
2. **CHOLMOD pour le préconditionneur.** La matrice locale est symétrique
   définie positive, une factorisation de Cholesky multithread remplace
   SuperLU si `scikit-sparse` est présent.
   ```bash
   brew install suite-sparse
   pip install scikit-sparse
   ```
3. **GPU Metal, expérimental.** `accel.py` contient deux backends, MLX
   (`pip install mlx`) et PyTorch sur le device `mps`, pour les paires proches
   et les FFT. Ils fonctionnent en simple précision avec des coordonnées
   relatives au centre de chaque triangle, ce qui préserve la précision sur des
   puces de plusieurs millimètres. Ces backends n'ont pas pu être exécutés sur
   la machine qui a écrit ce code. Un autotest compare leurs résultats à numpy
   au démarrage sur un échantillon, paires exactes, paires grille et FFT, et le
   calcul retombe sur numpy en cas d'écart ou d'erreur. La tolérance du gradient
   conjugué passe à `3e-6` en simple précision.

```bash
python -m scdc.run exemple.gds --Lsq 0.14 --seg 1.5 --backend mlx
SCDC_BACKEND=torch python -m scdc.run ...
```

Ce qu'il faut attendre. La construction du noyau et les FFT peuvent gagner un
ordre de grandeur sur GPU. Les produits creux et la résolution du
préconditionneur restent sur CPU, donc le gain global sur une grosse structure
sera plutôt d'un facteur 2 à 4, à mesurer avec `--backend numpy` contre
`--backend mlx` sur le même GDS. Si l'autotest échoue, le message indique
quelle partie diffère, et je pourrai corriger avec ce retour.

### Visualisation

`run.py` écrit une carte de `|K|` avec lignes de courant. Pour explorer, la
fenêtre interactive s'ouvre avec `--show`, ou plus tard depuis la sauvegarde.

```bash
python -m scdc.viewer resultat_solution.npz            # interactif
python -m scdc.viewer resultat_solution.npz --png zoom.png --log
```

Touches dans la fenêtre. `i` passe à la carte de courant du polygone
d'injection suivant (`I` pour le précédent), le titre indiquant le nom de la
source et son rang, `l` bascule linéaire/logarithmique, `S` (majuscule, `s` étant
le raccourci matplotlib de sauvegarde) recalcule les lignes de courant sur la
zone visible après un zoom, `q` affiche des flèches,
`c` efface. Un clic gauche imprime `K`, `|K|` et `J = K/d` au point cliqué.
L'échelle de couleur est commune à toutes les cartes, le courant injecté
étant le même pour chaque source. Avec `--png`, une image par source est
écrite (`--map k` pour n'en écrire qu'une).

---

## 5. Validation

`validation.py` compare le code à des résultats analytiques.

| test | résultat |
|---|---|
| inductance mutuelle nappe annulaire / boucle coaxiale, formule de Maxwell | écart 0.001 % à 0.006 % |
| conservation du courant à travers une coupe | 1.000000 mA pour 1 mA imposé |
| profil `K(y)` d'une piste avec `Λ → 0` comparé à `(I/π)/√((W/2)² − y²)` | écart RMS 1.6 % |
| inductance cinétique d'une piste, `Λ` grand, comparée à `µ0 Λ ℓ/W` | écart 2 % |
| fluxoïde imposé `n = 1` dans un anneau, cohérence `2E = Φ0 I` | exacte à 10⁻⁶ |
| patch jonction `L_J = 0,5` et `5 nH` dans une piste, `L_tot − L_piste` | écart 2×10⁻⁵ |
| solveur pFFT contre assemblage dense | 0,1 % sur `E`, 10⁻³ sur `K` |
| anneau 6–10 µm, `L = Φ0/I = 23.1 pH` contre `µ0 R [ln(8R/a) − 2] ≈ 21.7 pH` | cohérent |

```bash
python validation.py
```

---

## 6. Limites, et pourquoi l'accord expérimental reste imparfait

Ce code résout proprement le modèle de London 2D. Il ne résout pas les causes
usuelles des écarts entre simulation et mesure, qu'il vaut mieux connaître.

1. **Modèle de nappe.** La formulation suppose `d` petit devant les dimensions
   latérales. Pour `d ≳ λ` la répartition du courant dans l'épaisseur n'est plus
   uniforme, `Λ = λ coth(d/λ)` n'est qu'une correction au premier ordre, et une
   extraction volumique du type FastHenry ou TetraHenry est plus fidèle.
2. **`λ` effectif.** La valeur pertinente est `λ_eff = λ_L √(1 + ξ0/ℓ)` pour un
   film sale, dépend du dépôt, de la température et du champ. C'est en pratique
   la première source d'écart sur l'inductance cinétique, souvent 10 % à 30 %.
3. **Bords et géométrie réelle.** Le GDS n'est pas le circuit fabriqué. Une
   sous-gravure ou une sur-gravure de quelques dizaines de nanomètres modifie
   fortement le courant de bord dès que `Λ ≪ W`, et donc l'inductance.
4. **Airbridges.** Ils sont modélisés comme des nappes planes à `z = h` reliées
   au layer 1 par des piliers droits à inductance localisée (§2.5). La courbure
   réelle du pont, la résistance de contact et la répartition du courant sur
   la section du pilier au voisinage des pieds sont négligées.
5. **Contacts d'injection.** Le profil de courant vertical dans les disques des
   layers 3 et 4 est imposé par le problème résistif auxiliaire et n'est pas
   relaxé. Cela affecte la concentration de courant à l'intérieur du plot, très
   peu le champ lointain. La discussion de ce point se trouve chez Khapaev et
   Kupriyanov (arXiv:1412.3231).
6. **Vortex piégés.** Un vortex d'Abrikosov piégé modifie localement `K` et
   ajoute un quantum de flux. Le code accepte `n_k ≠ 0` pour les trous du
   layout, pas pour des vortex arbitrairement placés dans le métal.
7. **Non-linéarité.** `Λ(J)` augmente quand `J` approche `J_c`, ce qui n'est pas
   pris en compte. À 1 mA dans une piste de quelques microns, la correction
   reste en général sous le pour cent, mais elle devient dominante dans les
   constrictions.
8. **Effets de bord de puce.** Plan de masse fini, boîtier, et couplages hors
   du layout ne sont pas décrits.

---

## 7. Références

Modèle de London bidimensionnel et courants d'écran.

- J. Pearl, *Current distribution in superconducting films carrying quantized
  fluxoids*, Appl. Phys. Lett. **5**, 65 (1964).
  <https://doi.org/10.1063/1.1754056>
- E. H. Brandt et M. Indenbom, *Type-II-superconductor strip with current in a
  perpendicular magnetic field*, Phys. Rev. B **48**, 12893 (1993).
  <https://doi.org/10.1103/PhysRevB.48.12893>
- R. Meservey et P. M. Tedrow, *Measurements of the kinetic inductance of
  superconducting linear structures*, J. Appl. Phys. **40**, 2028 (1969).
  <https://doi.org/10.1063/1.1657905>
- J. R. Clem et K. K. Berggren, *Geometry-dependent critical currents in
  superconducting nanocircuits*, Phys. Rev. B **84**, 174510 (2011).
  <https://arxiv.org/abs/1109.4881>

Extraction d'inductance et codes existants.

- L. Bishop-Van Horn et K. A. Moler, *SuperScreen*, Comput. Phys. Commun.
  **280**, 108464 (2022). <https://arxiv.org/abs/2203.13388>
- M. M. Khapaev, *Inductance extraction of multilayer finite-thickness
  superconductor circuits*, IEEE Trans. Microw. Theory Tech. **49**, 217 (2001).
- M. M. Khapaev, A. Yu. Kidiyarova-Shevchenko, P. Magnelind, M. Yu. Kupriyanov,
  *3D-MLSI*, IEEE Trans. Appl. Supercond. **11**, 1090 (2001).
  <https://ieeexplore.ieee.org/document/919537/>
- M. M. Khapaev et M. Yu. Kupriyanov, *Inductance extraction of superconductor
  structures with internal current sources*, <https://arxiv.org/abs/1412.3231>
- M. Kamon, M. J. Tsuk, J. K. White, *FastHenry, a multipole-accelerated 3-D
  inductance extraction program*, IEEE Trans. Microw. Theory Tech. **42**,
  1750 (1994). <https://doi.org/10.1109/22.310584>
- J. R. Phillips et J. K. White, *A precorrected-FFT method for electrostatic
  analysis of complicated 3-D structures*, IEEE Trans. Comput.-Aided Des.
  Integr. Circuits Syst. **16**, 1059 (1997).
  <https://doi.org/10.1109/43.662670>
