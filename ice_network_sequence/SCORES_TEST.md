## Frozen test: sequence effect AB - BA (predicted by RTM_seq before the runs; observed in the network and the resolved reference)

| id | config | pred share_end | net r17 | net r33 | ref 200 | ref 400 | fixed-geometry net | pred heat rel | net heat rel | ref heat rel | pred meaningful routing | observed (net / ref) | dominant AB/BA/AVG (net) | tie | flip |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| X01 | scan config 2 | -0.0019 | -0.0005 | - | - | - | +0.0000 | +0.037 | +0.037 | - | False | False / - | R1/R1/R1 | [False, False] | False |
| X02 | scan config 13 | +0.0000 | +0.0000 | - | +0.0000 | - | +0.0000 | +0.000 | +0.000 | +0.000 | False | False / False | R2/R2/R2 | [False, False] | False |
| X03 | scan config 23 | -0.0001 | +0.0001 | - | - | - | +0.0001 | +0.013 | +0.012 | - | False | False / - | R2/R2/R2 | [False, False] | False |
| X04 | scan config 42 | +0.0240 | +0.0242 | +0.0233 | +0.0216 | +0.0223 | +0.0008 | +0.128 | +0.126 | +0.127 | False | False / False | R2/R2/R2 | [False, False] | False |
| X05 | scan config 77 | -0.0223 | -0.0091 | - | - | - | +0.0000 | +0.108 | +0.111 | - | False | False / - | R1/R1/R1 | [False, False] | False |
| X06 | scan config 74 | -0.0156 | -0.0149 | - | - | - | +0.0000 | +0.029 | +0.027 | - | False | False / - | R1/R1/R1 | [False, False] | False |
| X07 | scan config 51 | -0.0180 | -0.0181 | - | -0.0169 | -0.0175 | +0.0034 | -0.037 | -0.038 | -0.039 | False | False / False | R2/R2/R2 | [False, False] | False |
| X08 | scan config 43 | +0.0093 | +0.0084 | - | +0.0067 | - | +0.0149 | +0.185 | +0.179 | +0.179 | False | False / False | R2/R2/R2 | [False, True] | False |
| X09 | scan config 64 | -0.0471 | -0.0354 | - | -0.0283 | - | +0.0001 | +0.129 | +0.120 | +0.131 | True | True / True | R1/R1/R1 | [False, False] | False |
| X10 | scan config 96 | +0.0451 | +0.0400 | +0.0402 | +0.0400 | - | -0.0000 | +0.093 | +0.093 | +0.095 | False | False / False | R2/R2/R2 | [False, False] | False |

## Cumulative delivered-share effect, persistence, incision, averaging error (network r17)

| id | delivered-share effect | cum2 effect (m3) | persistence (end/start-of-continuation share effect) | incision rel | pred avg error (AVG-AB, AVG-BA) | observed avg error | averaging adequate | |AVG - mean(AB,BA)| / |AB - BA| (share_end) |
|---|---|---|---|---|---|---|---|---|
| X01 | +0.0024 | +7 | -0.14 | 0.036 | [0.001, -0.001] | [0.0, -0.0] | True | 0.1 |
| X02 | +0.0000 | +0 | - | 0.000 | [0.0, 0.0] | [0.0, 0.0] | True | 0.0 |
| X03 | +0.0017 | +8 | 0.00 | 0.021 | [0.001, 0.001] | [0.001, 0.001] | True | 6.2 |
| X04 | +0.0252 | +143 | 0.63 | 0.120 | [-0.012, 0.012] | [-0.012, 0.012] | True | 0.0 |
| X05 | +0.0010 | +10 | 0.85 | 0.103 | [0.011, -0.012] | [0.004, -0.005] | True | 0.1 |
| X06 | -0.0096 | -32 | -0.14 | 0.044 | [0.012, -0.003] | [0.012, -0.003] | True | 0.3 |
| X07 | -0.0065 | -30 | -0.11 | 0.132 | [0.017, -0.001] | [0.017, -0.001] | True | 0.5 |
| X08 | +0.0070 | +86 | 0.15 | 0.163 | [0.001, 0.01] | [0.001, 0.01] | False | 0.6 |
| X09 | -0.0712 | -517 | -0.13 | 0.111 | [0.156, 0.109] | [0.154, 0.119] | False | 3.9 |
| X10 | +0.0285 | +415 | 0.87 | 0.090 | [-0.023, 0.022] | [-0.02, 0.02] | True | 0.0 |

## Closure multipliers and tolerance (network r17)

| id | base share effect | nu x0.5 | nu x2 | base heat rel | nu x0.5 | nu x2 | destination flip (any) | legacy vs corrected tolerance (share effect) |
|---|---|---|---|---|---|---|---|---|
| X04 | +0.0242 | +0.0100 | +0.0456 | +0.126 | +0.136 | +0.098 | False | +0.0242 vs +0.0242 |
| X08 | +0.0084 | +0.0101 | +0.0078 | +0.179 | +0.184 | +0.166 | False | - vs +0.0084 |
| X09 | -0.0354 | -0.0192 | -0.0459 | +0.120 | +0.123 | +0.100 | False | -0.0354 vs -0.0354 |

## Numerical uncertainty of the effect (3 x rule)

* X04: net_share_end 0.0008, net_heat_rel 0.0001, ref_share_end 0.0007, ref_heat_rel 0.0003
* X07: ref_share_end 0.0006, ref_heat_rel 0.0012
* X10: net_share_end 0.0002, net_heat_rel 0.0004

## Summary
* sign of the share effect agrees in 7 of 7 configurations with |observed effect| > 0.003 (network); against the resolved reference 5 of 5
* predicted/observed magnitude (network): median 1.10, range 0.99-2.45; (reference): median 1.13, range 1.03-1.67
* heat-difference prediction (RTM_seq) vs network: max |error| 0.009, median 0.0013
* destination flips observed (network r17): 0 of 10; resolved reference (informative cases): 0
* meaningful routing effect predicted for ['X09'], observed (network) for ['X09'], reference ['X09']
* meaningful incision effect predicted for ['X04', 'X05', 'X08'], observed (network) for ['X04', 'X05', 'X08', 'X09']