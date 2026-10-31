# Scores (dev) — frozen rules

Criterion A: tp 3 fp 1 fn 0 tn 3; recall 1.0, precision 0.75; established: **False**

| cfg | category | R | unsafe | AVG misrepresents | class TR/AVG/PW6/PW3/EA | share err (pp) AVG / PW6 / PW3 / PWn / EA | wet err (h) AVG / PW6 / PW3 / EA | RTM-TR share err |
|---|---|---|---|---|---|---|---|---|
| D01 | wet | 0.01 | False | False | wet/wet/wet/-/wet | 0.3 / 0.0 / nan / 0.1 / 0.0 | 0.0 / 0.1 / nan / 0.1 | 0.10 |
| D02 | intermittent* | 25.23 | True | True | int/int/int/-/int | 8.5 / 1.5 / nan / 4.5 / 1.8 | 0.0 / 0.0 / nan / 0.2 | 0.07 |
| D03 | intermittent* | 23.15 | True | True | int/dry/dry/-/dry | 0.3 / 0.3 / nan / 0.3 / 0.2 | 3.8 / 3.8 / nan / 2.4 | 0.03 |
| D04 | intermittent | 4.85 | True | True | int/int/int/-/int | 6.5 / 0.1 / nan / 1.3 / 1.4 | 2.1 / 0.0 / nan / 0.8 | 0.08 |
| D05 | intermittent | 23.24 | True | False | int/int/int/-/int | 3.9 / 0.9 / nan / 2.9 / 0.4 | 0.8 / 1.9 / nan / 0.5 | 0.16 |
| D06 | dry | 0.00 | False | False | dry/dry/dry/-/dry | 0.0 / 0.0 / nan / 0.0 / 0.0 | 0.0 / 0.0 / nan / 0.0 | 0.00 |
| D07 | dry* | 0.00 | False | False | dry/dry/dry/-/dry | 0.0 / 0.0 / nan / 0.0 / 0.0 | 0.0 / 0.0 / nan / 0.0 | 0.00 |

Truth-positive configs: 3. Median |share error| (pp): AVG 6.47, PW6 0.32, PW3 nan, EA 1.43, PWn 1.28. EA worse than AVG on 0. **EA useful improvement (rule 2): False**. RTM-TR median error 0.07176335969268782 pp → success True.

## Robustness
```
{}
```