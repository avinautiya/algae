## Interpolation: parameter-evaluation set (64 Sobol points inside the development ranges)

### r17 (primary configuration)

* **Execution success**: 64 of 64 runs completed (0 failed/stalled/timed out). **Acceptable predictive performance** (class ok): 59 of 64 (92%); 3 unexplained, 2 outside regime, 0 numerical (execution failure).
* Within both models' assumptions (reduced envelope + shared checks): 34 completed cases; classes there {'ok': 31, 'failure_unexplained': 3, 'failure_outside_regime': 0}; outside: {'ok': 28, 'failure_unexplained': 0, 'failure_outside_regime': 2}.
* **Degenerate cases** (the reference delivers no water through either route in 8 h; agreement is trivial and carries no information about routing): 6 of 64 completed cases; the reduced model also delivers nothing in 6.
* Route (outlet) choice, clear cases (reference margin > 10 %): 56 of 56 correct; genuine near-ties (reference margin <= 10 %, not scored): 2; of those the reduced model picks the reference route in 2.
* Switching (R2 overtakes R1): reference 6 cases, reduced 8; agreement on existence 6 of 8 cases where either has it; time error where both: |err| median / max 2019 / 3245 s (n = 6)
* Error distributions over completed cases (median / 90th percentile / maximum):
  * hydrograph RMSE / inflow: 0.0265 / 0.0952 / 0.209; cumulative delivered-volume error (relative to total): 0.0683 / 0.223 / 0.662
  * |t_act2 error| (s): 284.7 / 662.1 / 6724 (n = 43)
  * |sill-lowering error| R1 / R2 (mm): 1.44 / 9.43 / 24  /  0.968 / 17.5 / 62
  * reduced solver time (s; two concurrent workers): 11.6 / 27.4 / 42.5

### r9 (secondary configuration)

* **Execution success**: 64 of 64 runs completed (0 failed/stalled/timed out). **Acceptable predictive performance** (class ok): 59 of 64 (92%); 3 unexplained, 2 outside regime, 0 numerical (execution failure).
* Within both models' assumptions (reduced envelope + shared checks): 34 completed cases; classes there {'ok': 31, 'failure_unexplained': 3, 'failure_outside_regime': 0}; outside: {'ok': 28, 'failure_unexplained': 0, 'failure_outside_regime': 2}.
* **Degenerate cases** (the reference delivers no water through either route in 8 h; agreement is trivial and carries no information about routing): 6 of 64 completed cases; the reduced model also delivers nothing in 6.
* Route (outlet) choice, clear cases (reference margin > 10 %): 56 of 56 correct; genuine near-ties (reference margin <= 10 %, not scored): 2; of those the reduced model picks the reference route in 2.
* Switching (R2 overtakes R1): reference 6 cases, reduced 9; agreement on existence 6 of 9 cases where either has it; time error where both: |err| median / max 2321 / 4019 s (n = 6)
* Error distributions over completed cases (median / 90th percentile / maximum):
  * hydrograph RMSE / inflow: 0.0346 / 0.0968 / 0.202; cumulative delivered-volume error (relative to total): 0.0764 / 0.295 / 0.675
  * |t_act2 error| (s): 321.4 / 757.1 / 7258 (n = 43)
  * |sill-lowering error| R1 / R2 (mm): 1.31 / 6.42 / 24  /  1.2 / 18.7 / 56.7
  * reduced solver time (s; two concurrent workers): 13.3 / 42.2 / 183

Non-ok reduced runs (all failures kept):

| point | config | class | hydrograph RMSE/I | route ok | near-tie | violated assumptions | note |
|---|---|---|---|---|---|---|---|
| E014 | r9 | failure_outside_regime | 0.104 | True | False | subcritical |  |
| E014 | r17 | failure_outside_regime | 0.105 | True | False | subcritical |  |
| E029 | r9 | failure_outside_regime | 0.107 | True | False | subcritical |  |
| E029 | r17 | failure_outside_regime | 0.122 | True | False | subcritical |  |
| E031 | r9 | failure_unexplained | 0.196 | True | False | none |  |
| E031 | r17 | failure_unexplained | 0.156 | True | False | none |  |
| E047 | r9 | failure_unexplained | 0.115 | True | False | none |  |
| E047 | r17 | failure_unexplained | 0.122 | True | False | none |  |
| E063 | r9 | failure_unexplained | 0.202 | True | False | none |  |
| E063 | r17 | failure_unexplained | 0.209 | True | False | none |  |

## Extrapolation: challenge set (32 Sobol points outside the development ranges on every axis)

### r17 (primary configuration)

* **Execution success**: 30 of 32 runs completed (2 failed/stalled/timed out). **Acceptable predictive performance** (class ok): 26 of 32 (81%); 0 unexplained, 4 outside regime, 2 numerical (execution failure).
* Within both models' assumptions (reduced envelope + shared checks): 6 completed cases; classes there {'ok': 6, 'failure_unexplained': 0, 'failure_outside_regime': 0}; outside: {'ok': 20, 'failure_unexplained': 0, 'failure_outside_regime': 4}.
* **Degenerate cases** (the reference delivers no water through either route in 8 h; agreement is trivial and carries no information about routing): 7 of 30 completed cases; the reduced model also delivers nothing in 6.
* Route (outlet) choice, clear cases (reference margin > 10 %): 23 of 23 correct; genuine near-ties (reference margin <= 10 %, not scored): 0; of those the reduced model picks the reference route in 0.
* Switching (R2 overtakes R1): reference 12 cases, reduced 9; agreement on existence 8 of 13 cases where either has it; time error where both: |err| median / max 15 / 1358 s (n = 8)
* Error distributions over completed cases (median / 90th percentile / maximum):
  * hydrograph RMSE / inflow: 0.0194 / 0.145 / 0.215; cumulative delivered-volume error (relative to total): 0.0418 / 0.266 / 1.74
  * |t_act2 error| (s): 146.9 / 797 / 1180 (n = 15)
  * |sill-lowering error| R1 / R2 (mm): 1.78 / 13.2 / 136  /  1.72 / 43.3 / 72.5
  * reduced solver time (s; two concurrent workers): 13.4 / 34.8 / 72.3

### r9 (secondary configuration)

* **Execution success**: 29 of 32 runs completed (3 failed/stalled/timed out). **Acceptable predictive performance** (class ok): 24 of 32 (75%); 0 unexplained, 5 outside regime, 3 numerical (execution failure).
* Within both models' assumptions (reduced envelope + shared checks): 8 completed cases; classes there {'ok': 8, 'failure_unexplained': 0, 'failure_outside_regime': 0}; outside: {'ok': 16, 'failure_unexplained': 0, 'failure_outside_regime': 5}.
* **Degenerate cases** (the reference delivers no water through either route in 8 h; agreement is trivial and carries no information about routing): 7 of 29 completed cases; the reduced model also delivers nothing in 7.
* Route (outlet) choice, clear cases (reference margin > 10 %): 21 of 22 correct; genuine near-ties (reference margin <= 10 %, not scored): 0; of those the reduced model picks the reference route in 0.
* Switching (R2 overtakes R1): reference 12 cases, reduced 9; agreement on existence 8 of 13 cases where either has it; time error where both: |err| median / max 0 / 3901 s (n = 8)
* Error distributions over completed cases (median / 90th percentile / maximum):
  * hydrograph RMSE / inflow: 0.0278 / 0.175 / 0.518; cumulative delivered-volume error (relative to total): 0.0423 / 0.28 / 1.69
  * |t_act2 error| (s): 93.39 / 816.5 / 1180 (n = 14)
  * |sill-lowering error| R1 / R2 (mm): 1.31 / 20.3 / 147  /  1.12 / 40.6 / 57.5
  * reduced solver time (s; two concurrent workers): 6.94 / 53.8 / 64.9

Non-ok reduced runs (all failures kept):

| point | config | class | hydrograph RMSE/I | route ok | near-tie | violated assumptions | note |
|---|---|---|---|---|---|---|---|
| C003 | r9 | failure_outside_regime | 0.518 | False | False | subcritical |  |
| C003 | r17 | failure_outside_regime | 0.170 | True | False | subcritical |  |
| C006 | r9 | failure_outside_regime | 0.201 | True | False | subcritical |  |
| C006 | r17 | failure_outside_regime | 0.215 | True | False | subcritical |  |
| C009 | r9 | failure_numerical | - | - | - | not diagnosed (no completed run) | stalled at simulated t = 900.0 s |
| C009 | r17 | failure_outside_regime | 0.142 | True | False | subcritical, storage_neglected |  |
| C010 | r9 | failure_outside_regime | 0.169 | True | False | subcritical |  |
| C010 | r17 | failure_outside_regime | 0.171 | True | False | subcritical |  |
| C016 | r9 | failure_outside_regime | 0.209 | True | False | subcritical |  |
| C019 | r9 | failure_numerical | - | - | - | not diagnosed (no completed run) | stalled at simulated t = 6300.0 s |
| C025 | r9 | failure_numerical | - | - | - | not diagnosed (no completed run) | AssertionError: sill incised below lake bottom: lake would drain completely (beyond model  |
| C025 | r17 | failure_numerical | - | - | - | not diagnosed (no completed run) | AssertionError: sill incised below lake bottom: lake would drain completely (beyond model  |
| C031 | r9 | failure_outside_regime | 0.125 | True | False | subcritical |  |
| C031 | r17 | failure_numerical | - | - | - | not diagnosed (no completed run) | RuntimeError: ODE failure at t=14850.0: Unexpected istate in LSODA. |

### Every challenge case (r17): is it within the physical assumptions of both models?

Reduced-model envelope: quasi-steady (p90 r_qs <= 1), subcritical (<= 5 % active snapshots with Froude >= 1), channel storage neglected (p90 <= 5 % of node storage). Shared: banks not overtopped (< 1 m), Manning plausible (median max Re >= 2000). The reference additionally represents supercritical flow and jumps.

| point | inflow m3/s | pulse amp / period h | deficit m | R2-sill offset m | slopes R1/R2 | ref outlet (margin) | r17 class | hydrograph RMSE/I | switching ref/red | violated assumptions | within both models |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C000 | 0.007 | 0.45 / 5.4 | +0.17 | +0.016 | 0.0020 / 0.097 | none (0.00) | ok | 0.000 | False/False | Re<2000 (Manning) | False |
| C001 | 0.070 | 0.05 / 1.1 | +0.39 | +0.126 | 0.0433 / 0.013 | R1 (0.66) | ok | 0.023 | False/False | subcritical | False |
| C002 | 0.182 | 0.70 / 2.5 | +0.29 | -0.076 | 0.0115 / 0.003 | R2 (0.57) | ok | 0.014 | True/True | subcritical, bank_overtopped | False |
| C003 | 0.020 | 0.25 / 0.4 | -0.09 | +0.032 | 0.0016 / 0.027 | R1 (0.15) | failure_outside_regime | 0.170 | True/True | subcritical | False |
| C004 | 0.027 | 0.81 / 0.8 | -0.02 | -0.144 | 0.0287 / 0.004 | R2 (0.79) | ok | 0.030 | True/False | subcritical | False |
| C005 | 0.173 | 0.37 / 4.9 | +0.22 | +0.097 | 0.0042 / 0.042 | R1 (0.61) | ok | 0.045 | False/False | subcritical, storage_neglected | False |
| C006 | 0.058 | 0.57 / 0.4 | +0.47 | -0.052 | 0.0008 / 0.060 | R2 (0.63) | failure_outside_regime | 0.215 | True/True | subcritical | False |
| C007 | 0.010 | 0.16 / 1.8 | +0.10 | +0.190 | 0.0175 / 0.010 | none (0.00) | ok | 0.000 | False/False | Re<2000 (Manning) | False |
| C008 | 0.013 | 0.77 / 0.3 | +0.38 | -0.121 | 0.0011 / 0.007 | none (0.00) | ok | 0.000 | False/False | Re<2000 (Manning) | False |
| C009 | 0.045 | 0.31 / 1.6 | +0.15 | +0.075 | 0.0086 / 0.067 | R1 (0.92) | failure_outside_regime | 0.142 | False/False | subcritical, storage_neglected | False |
| C010 | 0.132 | 0.56 / 0.6 | -0.03 | -0.029 | 0.0594 / 0.037 | R1 (0.14) | failure_outside_regime | 0.171 | True/False | subcritical | False |
| C011 | 0.035 | 0.06 / 4.2 | +0.34 | +0.168 | 0.0028 / 0.005 | none (0.00) | ok | 0.000 | False/False | Re<2000 (Manning) | False |
| C012 | 0.015 | 0.67 / 3.3 | +0.26 | -0.006 | 0.0126 / 0.018 | none (0.00) | ok | 0.000 | False/False | Re<2000 (Manning) | False |
| C013 | 0.237 | 0.18 / 0.5 | +0.05 | +0.148 | 0.0006 / 0.003 | R1 (0.68) | ok | 0.015 | False/False | subcritical | False |
| C014 | 0.091 | 0.88 / 7.0 | +0.07 | -0.098 | 0.0057 / 0.012 | R2 (0.82) | ok | 0.044 | True/False | subcritical | False |
| C015 | 0.005 | 0.42 / 1.4 | +0.45 | +0.054 | 0.0405 / 0.150 | none (0.00) | ok | 0.000 | False/False | Re<2000 (Manning) | False |
| C016 | 0.006 | 0.87 / 1.7 | -0.06 | +0.114 | 0.0038 / 0.126 | R1 (0.29) | ok | 0.073 | False/True | subcritical | False |
| C017 | 0.107 | 0.43 / 0.3 | +0.30 | +0.003 | 0.0230 / 0.011 | R1 (0.74) | ok | 0.048 | False/False | subcritical | False |
| C018 | 0.280 | 0.63 / 3.6 | +0.42 | +0.042 | 0.0217 / 0.003 | R1 (0.56) | ok | 0.011 | False/False | subcritical | False |
| C019 | 0.017 | 0.22 / 0.7 | +0.18 | -0.067 | 0.0008 / 0.020 | R2 (0.95) | ok | 0.044 | True/True | none | True |
| C020 | 0.032 | 0.51 / 0.6 | +0.11 | +0.109 | 0.0541 / 0.005 | R1 (1.00) | ok | 0.009 | False/False | subcritical | False |
| C021 | 0.119 | 0.11 / 2.9 | +0.49 | -0.131 | 0.0022 / 0.034 | R2 (1.00) | ok | 0.004 | True/True | none | True |
| C022 | 0.040 | 0.76 / 1.2 | +0.22 | +0.180 | 0.0014 / 0.080 | R1 (1.00) | ok | 0.000 | False/False | none | True |
| C023 | 0.012 | 0.32 / 7.7 | +0.01 | -0.060 | 0.0093 / 0.007 | R2 (0.37) | ok | 0.063 | True/False | none | True |
| C024 | 0.009 | 0.61 / 0.9 | +0.33 | +0.087 | 0.0006 / 0.008 | none (0.00) | ok | 0.002 | False/False | Re<2000 (Manning) | False |
| C025 | 0.053 | 0.11 / 5.8 | -0.06 | -0.110 | 0.0163 / 0.054 | R2 (0.92) | failure_numerical | - | - | n/a (run did not complete) | n/a |
| C026 | 0.158 | 0.82 / 0.4 | +0.14 | +0.159 | 0.0315 / 0.050 | R2 (0.13) | ok | 0.053 | True/True | subcritical | False |
| C027 | 0.024 | 0.36 / 2.2 | +0.35 | -0.039 | 0.0053 / 0.004 | R1 (0.12) | ok | 0.035 | True/True | subcritical | False |
| C028 | 0.023 | 0.70 / 4.3 | +0.43 | +0.135 | 0.0066 / 0.023 | R1 (1.00) | ok | 0.002 | False/False | none | True |
| C029 | 0.211 | 0.25 / 0.8 | +0.06 | -0.018 | 0.0011 / 0.002 | R1 (0.63) | ok | 0.005 | True/True | subcritical | False |
| C030 | 0.081 | 0.50 / 2.0 | +0.02 | +0.063 | 0.0030 / 0.016 | R1 (0.48) | ok | 0.024 | False/False | none | True |
| C031 | 0.008 | 0.00 / 0.3 | +0.26 | -0.089 | 0.0765 / 0.107 | R2 (1.00) | failure_numerical | - | - | n/a (run did not complete) | n/a |
