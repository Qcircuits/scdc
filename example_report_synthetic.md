# scdc report: synthetic

Generated 2026-09-13 05:23. Injected current I = 1 mA in each layer-3 polygon in turn, returned through the layer-4 polygons. All lengths in um.

## Mutual-inductance matrix M = Phi / I (pH)

Rows are the flux surfaces of layer 5, columns the injection polygons of layer 3. M[i, j] is the flux through surface i when the current I is injected in polygon j alone, divided by I.

| surface \ source | IN_A | Port B |
|---|---|---|
| HOLE | -0.00683205 | -0.00479399 |
| L1 | 0.00439255 | 0.00737709 |

## Eigenvalues of M^-1 (mA / Phi0)

Phi0 = h / 2e = 2.067833848e-15 Wb. Each eigenvalue lambda_k of M^-1 (A/Wb) is expressed as lambda_k Phi0 in mA, the current along the k-th eigen-direction producing one flux quantum through the corresponding combination of surfaces. Sorted by increasing modulus.

Condition number of M: 4.672.

| k | eigenvalue (mA/Phi0) | modulus (mA/Phi0) |
|---|---|---|
| 0 | 363.015 | 363.015 |
| 1 | -401.425 | 401.425 |

## Model parameters

- GDS file: `synthetic`
- Effective penetration depth Lambda = 0.1119 um, L_square = mu0 Lambda = 0.1406 pH/square
- Injected current I = 1 mA
- Ground current shares (layer 4): [0.5, 0.5]
- Flux surfaces evaluated at z = 0 um
- Edge segment length: 1 um near the region of interest
- Solver: pFFT + conjugate gradient, near_cells = 4

## Mesh

- Sheets: 1 (metal)
- Triangles: 4800, nodes: 2501
- Triangle size sqrt(area): min 0.707, median 0.707, max 0.707 um

## Injection polygons (layer 3)

| source | x_c | y_c | area (um^2) | mesh nodes | energy E (J) | L = 2E/I^2 (pH) | CG iterations | residual |
|---|---|---|---|---|---|---|---|---|
| IN_A | 5.000 | 5.000 | 28.27 | 25 | 9.486999e-18 | 18.9740 | 21 | 3.70e-09 |
| Port B | 5.000 | 35.000 | 28.27 | 25 | 9.486159e-18 | 18.9723 | 21 | 3.77e-09 |

## Flux surfaces (layer 5)

| surface | x_c | y_c | area (um^2) | Phi/Phi0 (IN_A) | Phi/Phi0 (Port B) |
|---|---|---|---|---|---|
| HOLE | 30.000 | 10.000 | 50.27 | -0.00330397 | -0.00231836 |
| L1 | 30.000 | 30.000 | 50.27 | 0.00212423 | 0.00356755 |

## Output files

- `test_solution.npz`
- `test_mutual.csv`
- `test_current_0_IN_A.png`
- `test_current_1_Port_B.png`
- `test_report.md`
