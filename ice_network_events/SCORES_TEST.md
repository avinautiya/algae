# Scores (test) — frozen rules

Criterion A: tp 3 fp 2 fn 0 tn 5; recall 1.0, precision 0.6; established: **False**

| cfg | category | R | unsafe | AVG misrepresents | class TR/AVG/PW6/PW3/EA | share err (pp) AVG / PW6 / PW3 / PWn / EA | wet err (h) AVG / PW6 / PW3 / EA | RTM-TR share err |
|---|---|---|---|---|---|---|---|---|
| E01 | wet | 0.00 | False | False | wet/wet/int/wet/wet | 1.0 / 0.5 / 0.1 / 1.0 / 1.0 | 0.5 / 0.8 / 0.2 / 0.5 | 1.13 |
| E02 | wet | 0.00 | False | False | wet/wet/wet/wet/wet | 0.3 / 0.0 / 0.0 / 0.3 / 0.3 | 0.0 / 0.0 / 0.0 / 0.0 | 0.14 |
| E03 | intermittent* | 8.41 | True | True | int/dry/dry/int/dry | 0.7 / 0.6 / 0.3 / 0.5 / 0.7 | 4.2 / 4.2 / 0.4 / 4.2 | 0.01 |
| E04 | intermittent* | 38.57 | True | True | int/int/int/int/int | 5.9 / 0.7 / 0.3 / 2.3 / 1.4 | 5.3 / 0.1 / 0.0 / 1.3 | 0.14 |
| E05 | intermittent* | 516.50 | True | False | dry/dry/dry/dry/dry | 0.0 / 0.0 / 0.0 / 0.0 / 0.0 | 0.0 / 0.0 / 0.0 / 0.0 | 0.00 |
| E06 | intermittent | 7.36 | True | False | dry/dry/dry/dry/dry | 0.0 / 0.0 / 0.0 / 0.0 / 0.0 | 0.0 / 0.0 / 0.0 / 0.0 | 0.00 |
| E07 | intermittent | 17.51 | True | True | int/wet/int/int/int | 4.8 / 0.7 / 0.4 / 3.5 / 0.7 | 12.1 / 0.3 / 0.1 / 0.8 | 0.22 |
| E08 | intermittent | 0.00 | False | False | int/int/int/int/int | 0.9 / 0.1 / 0.0 / 0.9 / 0.9 | 1.0 / 0.0 / 0.0 / 1.0 | 6.13 |
| E09 | dry | 0.00 | False | False | dry/dry/dry/dry/dry | 0.0 / 0.0 / 0.0 / 0.0 / 0.0 | 0.0 / 0.0 / 0.0 / 0.0 | 0.00 |
| E10 | dry* | 0.63 | False | False | dry/dry/dry/dry/dry | 0.0 / 0.0 / 0.0 / 0.0 / 0.0 | 0.0 / 0.0 / 0.0 / 0.0 | 0.00 |

Truth-positive configs: 3. Median |share error| (pp): AVG 4.75, PW6 0.66, PW3 0.32, EA 0.74, PWn 2.27. EA worse than AVG on 0. **EA useful improvement (rule 2): False**. RTM-TR median error 0.13784026510650985 pp → success True.

## Robustness
```
{
 "E03": {
  "ref": {
   "class_preserved": true,
   "EA_beats_AVG": false,
   "EA_beats_AVG_r17": false,
   "share_deliv_pp": {
    "TR": 0.3512853571187762,
    "AVG": 0.0,
    "EA": 0.0
   },
   "share_deliv_r17_pp": {
    "TR": 0.6604561656780253,
    "AVG": 0.0,
    "EA": 0.0
   }
  },
  "net_nu0.5": {
   "class_preserved": false,
   "EA_beats_AVG": true,
   "EA_beats_AVG_r17": false,
   "share_deliv_pp": {
    "TR": 2.1662825841526105,
    "AVG": 0.0,
    "EA": 0.5815264695816864
   },
   "share_deliv_r17_pp": {
    "TR": 0.6604561656780253,
    "AVG": 0.0,
    "EA": 0.0
   }
  },
  "net_nu2": {
   "class_preserved": false,
   "EA_beats_AVG": false,
   "EA_beats_AVG_r17": false,
   "share_deliv_pp": {
    "TR": 0.0,
    "AVG": 0.0,
    "EA": 0.0
   },
   "share_deliv_r17_pp": {
    "TR": 0.6604561656780253,
    "AVG": 0.0,
    "EA": 0.0
   }
  },
  "legacy_tol_TR_share_pp": [
   0.6604561039868965,
   0.6604561656780253
  ]
 },
 "E04": {
  "ref": {
   "class_preserved": true,
   "EA_beats_AVG": true,
   "EA_beats_AVG_r17": true,
   "share_deliv_pp": {
    "TR": 4.56812510274298,
    "AVG": 0.10598436887657217,
    "EA": 2.9959634229114545
   },
   "share_deliv_r17_pp": {
    "TR": 6.201728567036223,
    "AVG": 0.31297848435574954,
    "EA": 4.837303953873214
   }
  },
  "net_nu0.5": {
   "class_preserved": true,
   "EA_beats_AVG": true,
   "EA_beats_AVG_r17": true,
   "share_deliv_pp": {
    "TR": 8.734233944396582,
    "AVG": 0.8100283927826978,
    "EA": 7.437802667938789
   },
   "share_deliv_r17_pp": {
    "TR": 6.201728567036223,
    "AVG": 0.31297848435574954,
    "EA": 4.837303953873214
   }
  },
  "net_nu2": {
   "class_preserved": true,
   "EA_beats_AVG": true,
   "EA_beats_AVG_r17": true,
   "share_deliv_pp": {
    "TR": 3.6106221078558716,
    "AVG": 0.05956424598881933,
    "EA": 1.5467576218313293
   },
   "share_deliv_r17_pp": {
    "TR": 6.201728567036223,
    "AVG": 0.31297848435574954,
    "EA": 4.837303953873214
   }
  },
  "legacy_tol_TR_share_pp": [
   6.201728610324238,
   6.201728567036223
  ]
 },
 "E05": {
  "net_nu0.5": {
   "class_preserved": false,
   "EA_beats_AVG": false,
   "EA_beats_AVG_r17": false,
   "share_deliv_pp": {
    "TR": 0.3053524289025002,
    "AVG": 0.0,
    "EA": 0.0
   },
   "share_deliv_r17_pp": {
    "TR": 0.0,
    "AVG": 0.0,
    "EA": 0.0
   }
  },
  "net_nu2": {
   "class_preserved": true,
   "EA_beats_AVG": false,
   "EA_beats_AVG_r17": false,
   "share_deliv_pp": {
    "TR": 0.0,
    "AVG": 0.0,
    "EA": 0.0
   },
   "share_deliv_r17_pp": {
    "TR": 0.0,
    "AVG": 0.0,
    "EA": 0.0
   }
  }
 },
 "E06": {
  "ref": {
   "class_preserved": true,
   "EA_beats_AVG": false,
   "EA_beats_AVG_r17": false,
   "share_deliv_pp": {
    "TR": 0.0,
    "AVG": 0.0,
    "EA": 0.0
   },
   "share_deliv_r17_pp": {
    "TR": 0.0,
    "AVG": 0.0,
    "EA": 0.0
   }
  }
 },
 "E10": {
  "ref": {
   "class_preserved": true,
   "EA_beats_AVG": false,
   "EA_beats_AVG_r17": false,
   "share_deliv_pp": {
    "TR": 0.0,
    "AVG": 0.0,
    "EA": 0.0
   },
   "share_deliv_r17_pp": {
    "TR": 0.0,
    "AVG": 0.0,
    "EA": 0.0
   }
  }
 }
}
```