# Convergence study (frozen scenarios S1-S4; records `records/runs/*conv*`, summary `records/convergence_summary.json`)

Idealised controlled problems, not a forecast. All values are generated from immutable per-run records by `convergence_analysis.py` / `convergence_table.py`.

## Definitions
* **R1, R2**: the two routes out of the lake. R1: 100 m long, 0.4 m wide, bed slope 0.005, lower sill. R2: 60 m long, 0.4 m wide, bed slope 0.02, sill above R1 by 0.05 m (S1, S4, S5), 0.02 m (S2) or 0.10 m (S3).
* **Incision** = lowering of the route's sill (bed elevation at the lake end of the route) after 8 h, in mm. **Incision error** = reduced minus reference; **relative incision error** = that difference divided by the reference sill lowering of the same route at the same time (the 800-cell reference in this table; the 400-cell reference in `FROZEN_S1_S5_ANALYSIS.md`).
* **Hydrograph RMSE/I**: root-mean-square difference between reduced and reference outflow of both routes on the reference's 60 s grid, combined as sqrt(R1^2 + R2^2), divided by the (mean) inflow.
* **Delivered-volume error**: difference of the 8 h volume through each route divided by the reference volume of that route. **Route (outlet choice)**: the route with the larger delivered volume. **t_act2**: first time Q_R2 > 1 % of the inflow. **t_dom2 (switching time)**: first time Q_R2 > Q_R1; "-" when it does not occur in either run (no switch).
* Baseline configuration: 17 stations per route, geometry step 450 s, 16 hydraulic substeps per segment, ODE rtol 1e-7 / atol 1e-9. One factor is varied at a time: stations 5/9/17/33/65 (spacing 25/12.5/6.25/3.1/1.6 m for R1, 15/7.5/3.75/1.9/0.94 m for R2); geometry step 112.5-1800 s; substeps 2-32; rtol 1e-9 to 1e-5. Reference: 100/200/400/800 cells at CFL 0.4, and 200 cells at CFL 0.2 and 0.1 (time-step test); target = 800 cells.

## How the reference uncertainty was estimated
The reference is a first-order upwind scheme (observed order 0.5-1.06 in the 100/200/400 frozen runs). The uncertainty of the 800-cell target is estimated as the change between the 400- and 800-cell references, |Q_800 - Q_400|: for an order-1 scheme this equals the 800-cell discretisation error; for an order-0.5 scheme the error would be 2.4 times larger, so the values below are a lower-range estimate (both are given for the hydrograph). It is an estimate from two meshes, not a bound; the 200-cell CFL test shows the reference is mesh-limited, not time-step-limited (hydrograph error vs 800 cells 0.0129 / 0.0125 / 0.0123 at CFL 0.4 / 0.2 / 0.1 in S1).

## S1

Target (800 cells): route R1; t_act2 4820.8 s; t_dom2 22763 s; sill lowering R1/R2 45.98 / 55.21 mm; delivered R1/R2 546 / 490 m3.
**Reference uncertainty (|400 - 800|)**: hydrograph RMSE/I 0.0044 (x2.4 if order 0.5: 0.011); t_act2 2.3 s; t_dom2 680 s; sill lowering R1/R2 0.11 / 0.28 mm.

| stations | route ok | hydrograph RMSE/I | order vs prev. | t_act2 err (s) | t_dom2 err (s) | delivered err R1/R2 (%) | incision err R1/R2 (mm) | rel. incision err R1/R2 (%) | runtime* (s) | law evals | s/eval (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | NO | 0.0544 | - | -304.9 | -10177 | -4.1 / +12.4 | -0.99 / +11.97 | -2.2 / +21.7 | 27.7 | 49704 | 0.56 |
| 9 | NO | 0.047 | 0.25 | -296.8 | -9320.3 | -3.7 / +10.7 | +0.22 / +11.81 | +0.5 / +21.4 | 15.8 | 14840 | 1.07 |
| 17 | yes | 0.0259 | 0.94 | -256.4 | -6127 | -1.5 / +6.0 | +0.83 / +9.21 | +1.8 / +16.7 | 10.5 | 5852 | 1.79 |
| 33 | yes | 0.00579 | 2.3 | 13.65 | -1178.4 | +1.1 / +1.0 | +1.11 / +5.64 | +2.4 / +10.2 | 71.2 | 8716 | 8.17 |
| 65 | yes | 0.00504 | 0.21 | 8.666 | -726.12 | +1.1 / +0.8 | +0.86 / +2.12 | +1.9 / +3.8 | 48.4 | 5392 | 8.97 |

| factor | value | route ok | hydrograph RMSE/I | t_act2 err (s) | incision err R1/R2 (mm) | rel. R2 incision err (%) | runtime* (s) | law evals |
|---|---|---|---|---|---|---|---|---|
| dt_geom | 112.5 | yes | 0.0264 | -256.2 | +0.84 / +9.25 | +16.8 | 17.4 | 10432 |
| dt_geom | 225 | yes | 0.0262 | -255.9 | +0.83 / +9.24 | +16.7 | 13.7 | 7760 |
| dt_geom | 450 | yes | 0.0259 | -256.4 | +0.83 / +9.21 | +16.7 | 10.5 | 5852 |
| dt_geom | 900 | yes | 0.0253 | -257.4 | +0.83 / +9.17 | +16.6 | 9.31 | 5122 |
| dt_geom | 1800 | yes | 0.0244 | -246.3 | +0.82 / +9.08 | +16.5 | 8.67 | 4686 |
| msub | 2 | yes | 0.0236 | -277 | +0.83 / +8.92 | +16.2 | 17 | 36492 |
| msub | 4 | yes | 0.0254 | -256.8 | +0.83 / +9.15 | +16.6 | 28.2 | 42892 |
| msub | 8 | yes | 0.0258 | -256.4 | +0.83 / +9.20 | +16.7 | 18.7 | 16500 |
| msub | 16 | yes | 0.0259 | -256.4 | +0.83 / +9.21 | +16.7 | 10.5 | 5852 |
| msub | 32 | yes | 0.0259 | -256.4 | +0.83 / +9.22 | +16.7 | 14.9 | 4834 |
| rtol | 1e-09 | yes | 0.0259 | -256.4 | +0.83 / +9.21 | +16.7 | 13.8 | 8234 |
| rtol | 1e-08 | yes | 0.0259 | -256.4 | +0.83 / +9.21 | +16.7 | 12 | 7058 |
| rtol | 1e-07 | yes | 0.0259 | -256.4 | +0.83 / +9.21 | +16.7 | 10.5 | 5852 |
| rtol | 1e-06 | yes | 0.0259 | -256.4 | +0.83 / +9.21 | +16.7 | 9.85 | 5340 |
| rtol | 1e-05 | yes | 0.0259 | -256.4 | +0.83 / +9.21 | +16.7 | 8.16 | 4282 |

| reference mesh / CFL | hydrograph RMSE/I vs 800 | runtime* (s) | steps | t_act2 (s) | R2 sill lowering (mm) |
|---|---|---|---|---|---|
| 100 cells, CFL 0.4 | 0.0293 | 19.6 | 147967 | 4817 | 53.57 |
| 200 cells, CFL 0.1 | 0.0123 | 213 | 1204434 | 4806.1 | 54.72 |
| 200 cells, CFL 0.2 | 0.0125 | 103 | 602226 | 4808.7 | 54.62 |
| 200 cells, CFL 0.4 | 0.0129 | 53.6 | 301130 | 4813.7 | 54.41 |
| 400 cells, CFL 0.4 | 0.00438 | 173 | 611478 | 4818.6 | 54.93 |
| 800 cells, CFL 0.4 | - | 614 | 1237808 | 4820.8 | 55.21 |

## S2

Target (800 cells): route R2; t_act2 3530.4 s; t_dom2 7735.2 s; sill lowering R1/R2 45.00 / 60.47 mm; delivered R1/R2 474 / 604 m3.
**Reference uncertainty (|400 - 800|)**: hydrograph RMSE/I 0.0047 (x2.4 if order 0.5: 0.011); t_act2 2 s; t_dom2 39.3 s; sill lowering R1/R2 0.11 / 0.29 mm.

| stations | route ok | hydrograph RMSE/I | order vs prev. | t_act2 err (s) | t_dom2 err (s) | delivered err R1/R2 (%) | incision err R1/R2 (mm) | rel. incision err R1/R2 (%) | runtime* (s) | law evals | s/eval (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | yes | 0.0569 | - | -296.8 | -2403.8 | -5.6 / +10.8 | -1.13 / +12.12 | -2.5 / +20.0 | 28 | 49914 | 0.56 |
| 9 | yes | 0.0479 | 0.29 | -285.7 | -2384.8 | -4.7 / +9.1 | +0.16 / +11.81 | +0.4 / +19.5 | 15.6 | 14318 | 1.09 |
| 17 | yes | 0.0258 | 0.97 | -243.3 | -1738.9 | -2.0 / +4.9 | +0.83 / +9.20 | +1.9 / +15.2 | 10.3 | 5616 | 1.84 |
| 33 | yes | 0.00612 | 2.2 | 12.07 | -92.737 | +1.1 / +0.9 | +1.15 / +5.79 | +2.6 / +9.6 | 66.4 | 8122 | 8.17 |
| 65 | yes | 0.00537 | 0.19 | 10.27 | -96.298 | +1.1 / +0.7 | +0.90 / +2.21 | +2.0 / +3.6 | 48.9 | 5208 | 9.38 |

| factor | value | route ok | hydrograph RMSE/I | t_act2 err (s) | incision err R1/R2 (mm) | rel. R2 incision err (%) | runtime* (s) | law evals |
|---|---|---|---|---|---|---|---|---|
| dt_geom | 112.5 | yes | 0.0264 | -248 | +0.84 / +9.24 | +15.3 | 17.4 | 9824 |
| dt_geom | 225 | yes | 0.0262 | -243.2 | +0.84 / +9.23 | +15.3 | 13.3 | 7356 |
| dt_geom | 450 | yes | 0.0258 | -243.3 | +0.83 / +9.20 | +15.2 | 10.3 | 5616 |
| dt_geom | 900 | yes | 0.025 | -235.9 | +0.83 / +9.15 | +15.1 | 9.59 | 4946 |
| dt_geom | 1800 | yes | 0.0239 | -236.2 | +0.84 / +9.07 | +15.0 | 9.13 | 4518 |
| msub | 2 | yes | 0.0234 | -261.8 | +0.84 / +8.92 | +14.8 | 22.4 | 44304 |
| msub | 4 | yes | 0.0252 | -243.6 | +0.84 / +9.14 | +15.1 | 31.3 | 45038 |
| msub | 8 | yes | 0.0257 | -243.3 | +0.84 / +9.19 | +15.2 | 19.4 | 16584 |
| msub | 16 | yes | 0.0258 | -243.3 | +0.83 / +9.20 | +15.2 | 10.3 | 5616 |
| msub | 32 | yes | 0.0258 | -243.3 | +0.83 / +9.21 | +15.2 | 13.8 | 4584 |
| rtol | 1e-09 | yes | 0.0258 | -243.3 | +0.83 / +9.20 | +15.2 | 14 | 7908 |
| rtol | 1e-08 | yes | 0.0258 | -243.3 | +0.83 / +9.20 | +15.2 | 12 | 6810 |
| rtol | 1e-07 | yes | 0.0258 | -243.3 | +0.83 / +9.20 | +15.2 | 10.3 | 5616 |
| rtol | 1e-06 | yes | 0.0258 | -243.3 | +0.83 / +9.20 | +15.2 | 9.44 | 5040 |
| rtol | 1e-05 | yes | 0.0258 | -243.3 | +0.83 / +9.20 | +15.2 | 7.64 | 4074 |

| reference mesh / CFL | hydrograph RMSE/I vs 800 | runtime* (s) | steps | t_act2 (s) | R2 sill lowering (mm) |
|---|---|---|---|---|---|
| 100 cells, CFL 0.4 | 0.0316 | 20.4 | 160933 | 3525.5 | 58.73 |
| 200 cells, CFL 0.1 | 0.0133 | 234 | 1311487 | 3512.5 | 59.97 |
| 200 cells, CFL 0.2 | 0.0135 | 121 | 655710 | 3516.3 | 59.85 |
| 200 cells, CFL 0.4 | 0.014 | 58.5 | 327820 | 3523.8 | 59.64 |
| 400 cells, CFL 0.4 | 0.00474 | 193 | 665739 | 3528.4 | 60.18 |
| 800 cells, CFL 0.4 | - | 669 | 1347659 | 3530.4 | 60.47 |

## S3

Target (800 cells): route R1; t_act2 4535.5 s; t_dom2 - s; sill lowering R1/R2 52.26 / 61.75 mm; delivered R1/R2 995 / 741 m3.
**Reference uncertainty (|400 - 800|)**: hydrograph RMSE/I 0.0034 (x2.4 if order 0.5: 0.0083); t_act2 0.96 s; t_dom2 - s; sill lowering R1/R2 0.10 / 0.26 mm.

| stations | route ok | hydrograph RMSE/I | order vs prev. | t_act2 err (s) | t_dom2 err (s) | delivered err R1/R2 (%) | incision err R1/R2 (mm) | rel. incision err R1/R2 (%) | runtime* (s) | law evals | s/eval (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | yes | 0.0399 | - | -226.1 | - | -2.9 / +9.5 | -0.90 / +11.10 | -1.7 / +18.0 | 26.7 | 48080 | 0.56 |
| 9 | yes | 0.0324 | 0.35 | -219.9 | - | -2.3 / +7.8 | +0.10 / +10.62 | +0.2 / +17.2 | 12.5 | 12052 | 1.04 |
| 17 | yes | 0.0168 | 1 | -186.5 | - | -0.7 / +4.1 | +0.63 / +8.31 | +1.2 / +13.5 | 9.38 | 5754 | 1.63 |
| 33 | yes | 0.00443 | 2 | -1.005 | - | +0.8 / +0.9 | +0.87 / +5.49 | +1.7 / +8.9 | 99.4 | 10526 | 9.44 |
| 65 | yes | 0.00398 | 0.16 | -4.006 | - | +0.9 / +0.7 | +0.70 / +2.13 | +1.3 / +3.5 | 45.8 | 5538 | 8.26 |

| factor | value | route ok | hydrograph RMSE/I | t_act2 err (s) | incision err R1/R2 (mm) | rel. R2 incision err (%) | runtime* (s) | law evals |
|---|---|---|---|---|---|---|---|---|
| dt_geom | 112.5 | yes | 0.017 | -192.2 | +0.64 / +8.34 | +13.5 | 15.7 | 10194 |
| dt_geom | 225 | yes | 0.0169 | -192.6 | +0.64 / +8.33 | +13.5 | 12.1 | 7550 |
| dt_geom | 450 | yes | 0.0168 | -186.5 | +0.63 / +8.31 | +13.5 | 9.38 | 5754 |
| dt_geom | 900 | yes | 0.0165 | -186.1 | +0.63 / +8.28 | +13.4 | 8.6 | 5028 |
| dt_geom | 1800 | yes | 0.0164 | -187.7 | +0.63 / +8.22 | +13.3 | 6.93 | 4116 |
| msub | 2 | yes | 0.0155 | -196.9 | +0.61 / +8.08 | +13.1 | 11.6 | 27412 |
| msub | 4 | yes | 0.0165 | -186.6 | +0.63 / +8.26 | +13.4 | 26.2 | 41546 |
| msub | 8 | yes | 0.0167 | -186.5 | +0.63 / +8.30 | +13.4 | 18.2 | 15674 |
| msub | 16 | yes | 0.0168 | -186.5 | +0.63 / +8.31 | +13.5 | 9.38 | 5754 |
| msub | 32 | yes | 0.0168 | -186.5 | +0.63 / +8.31 | +13.5 | 13.2 | 4736 |
| rtol | 1e-09 | yes | 0.0168 | -186.5 | +0.63 / +8.31 | +13.5 | 13.2 | 8142 |
| rtol | 1e-08 | yes | 0.0168 | -186.5 | +0.63 / +8.31 | +13.5 | 11.4 | 6906 |
| rtol | 1e-07 | yes | 0.0168 | -186.5 | +0.63 / +8.31 | +13.5 | 9.38 | 5754 |
| rtol | 1e-06 | yes | 0.0168 | -186.5 | +0.63 / +8.31 | +13.5 | 8.25 | 4866 |
| rtol | 1e-05 | yes | 0.0168 | -186.5 | +0.63 / +8.31 | +13.5 | 7.27 | 4188 |

| reference mesh / CFL | hydrograph RMSE/I vs 800 | runtime* (s) | steps | t_act2 (s) | R2 sill lowering (mm) |
|---|---|---|---|---|---|
| 100 cells, CFL 0.4 | 0.0232 | 24.3 | 174523 | 4531.6 | 60.17 |
| 200 cells, CFL 0.1 | 0.00988 | 250 | 1422113 | 4529.5 | 61.27 |
| 200 cells, CFL 0.2 | 0.00998 | 125 | 711059 | 4530.6 | 61.18 |
| 200 cells, CFL 0.4 | 0.0102 | 65.8 | 355524 | 4532.8 | 61.01 |
| 400 cells, CFL 0.4 | 0.00344 | 213 | 721485 | 4534.6 | 61.49 |
| 800 cells, CFL 0.4 | - | 720 | 1459315 | 4535.5 | 61.75 |

## S4

Target (800 cells): route R1; t_act2 4827.9 s; t_dom2 - s; sill lowering R1/R2 13.75 / 15.75 mm; delivered R1/R2 525 / 440 m3.
**Reference uncertainty (|400 - 800|)**: hydrograph RMSE/I 0.0036 (x2.4 if order 0.5: 0.0086); t_act2 2.3 s; t_dom2 - s; sill lowering R1/R2 0.04 / 0.05 mm.

| stations | route ok | hydrograph RMSE/I | order vs prev. | t_act2 err (s) | t_dom2 err (s) | delivered err R1/R2 (%) | incision err R1/R2 (mm) | rel. incision err R1/R2 (%) | runtime* (s) | law evals | s/eval (ms) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | yes | 0.0455 | - | -309.6 | - | -3.0 / +11.8 | -0.26 / +4.10 | -1.9 / +26.0 | 19.4 | 35504 | 0.55 |
| 9 | yes | 0.0396 | 0.24 | -299.7 | - | -2.9 / +10.2 | +0.12 / +4.06 | +0.9 / +25.8 | 10.6 | 10598 | 1.00 |
| 17 | yes | 0.0218 | 0.94 | -258.3 | - | -1.1 / +5.7 | +0.30 / +3.20 | +2.2 / +20.3 | 8.46 | 5162 | 1.64 |
| 33 | yes | 0.00419 | 2.5 | 6.613 | - | +1.2 / +0.9 | +0.38 / +1.96 | +2.7 / +12.4 | 30.6 | 5772 | 5.30 |
| 65 | yes | 0.00392 | 0.099 | 1.867 | - | +1.1 / +0.8 | +0.29 / +0.70 | +2.1 / +4.5 | 44.6 | 5162 | 8.64 |

| factor | value | route ok | hydrograph RMSE/I | t_act2 err (s) | incision err R1/R2 (mm) | rel. R2 incision err (%) | runtime* (s) | law evals |
|---|---|---|---|---|---|---|---|---|
| dt_geom | 112.5 | yes | 0.0219 | -258.4 | +0.30 / +3.21 | +20.4 | 14.4 | 9490 |
| dt_geom | 225 | yes | 0.0219 | -258.2 | +0.30 / +3.20 | +20.4 | 11.1 | 6890 |
| dt_geom | 450 | yes | 0.0218 | -258.3 | +0.30 / +3.20 | +20.3 | 8.46 | 5162 |
| dt_geom | 900 | yes | 0.0216 | -258.6 | +0.30 / +3.20 | +20.3 | 7.8 | 4468 |
| dt_geom | 1800 | yes | 0.0212 | -255.3 | +0.30 / +3.19 | +20.3 | 6.73 | 3902 |
| msub | 2 | yes | 0.0197 | -277.9 | +0.30 / +3.11 | +19.7 | 11.9 | 27168 |
| msub | 4 | yes | 0.0213 | -258.5 | +0.30 / +3.18 | +20.2 | 21.7 | 32628 |
| msub | 8 | yes | 0.0217 | -258.3 | +0.30 / +3.20 | +20.3 | 12.2 | 11228 |
| msub | 16 | yes | 0.0218 | -258.3 | +0.30 / +3.20 | +20.3 | 8.46 | 5162 |
| msub | 32 | yes | 0.0218 | -258.3 | +0.30 / +3.20 | +20.3 | 13.5 | 4600 |
| rtol | 1e-09 | yes | 0.0218 | -258.3 | +0.30 / +3.20 | +20.3 | 11.7 | 7268 |
| rtol | 1e-08 | yes | 0.0218 | -258.3 | +0.30 / +3.20 | +20.3 | 9.66 | 6030 |
| rtol | 1e-07 | yes | 0.0218 | -258.3 | +0.30 / +3.20 | +20.3 | 8.46 | 5162 |
| rtol | 1e-06 | yes | 0.0218 | -258.3 | +0.30 / +3.20 | +20.3 | 7.53 | 4502 |
| rtol | 1e-05 | yes | 0.0218 | -258.3 | +0.30 / +3.20 | +20.3 | 6.21 | 3728 |

| reference mesh / CFL | hydrograph RMSE/I vs 800 | runtime* (s) | steps | t_act2 (s) | R2 sill lowering (mm) |
|---|---|---|---|---|---|
| 100 cells, CFL 0.4 | 0.0245 | 19.3 | 144115 | 4824 | 15.47 |
| 200 cells, CFL 0.1 | 0.0103 | 204 | 1165775 | 4812.8 | 15.69 |
| 200 cells, CFL 0.2 | 0.0104 | 99.8 | 582992 | 4815.5 | 15.66 |
| 200 cells, CFL 0.4 | 0.0106 | 50.9 | 291603 | 4820.6 | 15.61 |
| 400 cells, CFL 0.4 | 0.00356 | 164 | 590632 | 4825.6 | 15.70 |
| 800 cells, CFL 0.4 | - | 581 | 1193539 | 4827.9 | 15.75 |

*Runtimes were measured with two campaign workers running concurrently on a four-core machine (and, for part of the study, one further job); they are indicative, not matched-resource timings. The serial timing set is reported separately.

## Which errors refine and which persist (numbers from the tables above)
* **Outlet choice**: wrong only at 5 and 9 stations in S1 (a near-tie: reference delivered 546 against 490 m3, 5.6 % margin); correct at every setting in S2-S4 and from 17 stations in S1. Not affected by geometry step, substeps or rtol.
* **Switching time (S1, S2 only)**: error shrinks with stations in S1 from -10 177 s (5) to -726 s (65), against a reference uncertainty of 680 s: resolved to the reference uncertainty only at 65 stations. S2: -2404 s to -96 s (reference uncertainty 39 s: not yet resolved at 65). Always early (reduced model overtakes sooner).
* **Activation time of R2**: -226 to -310 s at 5-17 stations, +14 to -4 s at 33-65 stations; reference uncertainty 1-2 s. Refines, and flips sign between 17 and 33 stations.
* **Discharge (hydrograph)**: 0.054 -> 0.0050 (S1) with stations; second order between 17 and 33 stations, then flat (0.0058 -> 0.0050) at roughly the size of the reference uncertainty (0.0044; 0.011 if the reference converges at order 0.5). Independent of geometry step, substeps and rtol (changes below 0.001 absolute between the extremes).
* **Delivered volume**: R2 over-delivered by +12.4 % (5 stations) -> +0.8 % (65) in S1, +9-12 % -> +0.7-0.9 % in S2-S4; R1 within -5.6 % to +1.2 %. Settles to +0.7-1.1 % (R1 and R2), above any estimated reference uncertainty of volume but small.
* **Incision**: the R2 over-prediction **decreases monotonically with station refinement**, so it is a resolution error and not a structural offset in the tested range: S1 21.7 -> 16.7 -> 10.2 -> 3.8 % (5, 17, 33, 65 stations), S2 20.0 -> 3.6 %, S3 18.0 -> 3.5 %, S4 26.0 -> 4.5 % (error in mm: 11.97 -> 2.12 for S1). It does **not** depend on the geometry step (16.5-16.8 % over 112-1800 s), substeps (16.2-16.7 %) or rtol (S1), so the first-order geometry update is not limiting and no higher-order coupling was implemented. Residual at 65 stations: +2.1 mm (R2) against a reference uncertainty of 0.28 mm, i.e. still 7x the estimated reference uncertainty; whether it reaches zero cannot be told from five stations values (observed order 33 -> 65: 1.4).
* **A possibly structural remainder is visible in R1 incision**: +0.83, +1.11, +0.86 mm (S1; 1.8-2.4 %) at 17, 33, 65 stations does not shrink, against a reference uncertainty of 0.11 mm. Candidate causes (not tested): the sill value in the reduced model is a station value averaged from adjacent segments while the reference uses a cell value; different wall-retreat idealisation (D' = D + dz) from the reference's per-cell melt; unequal thermal-entrance treatment. This is a ~2 % difference in R1 and is reported as unexplained.
* **Cost**: not monotone in stations. Time per discharge-law evaluation rises from 0.56 ms (5 stations) to 1.8 ms (17) to 8.2 ms (33) and 9.0 ms (65), while the number of evaluations varies 49 704 / 14 840 / 5 852 / 8 716 / 5 392; the 33-station run is the most expensive (71 s) and the 65-station run costs less (48 s). The profile (`records/profile/`) attributes 90-97 % of the time to the discharge-law solve, with mean Brent iterations rising from 4.9 (17 stations) to 7.6 (33) and max 18 -> 36 (bracket search at near-dry states). No optimisation was implemented.
