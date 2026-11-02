# Conditional field-ablation transfer result

**No consistent field-ablation superiority is established.** All twelve directed transfers were executed for all four families and all36declared metadata scenarios (1728rows). The first2020→2023 year transfer is the primary direction; the other nine edges are retrospective spatial/source sensitivity, not nine independent glaciers.

|Thickness|Model|Primary2020→2023 MAPE|All12edge MAPE|
|---|---|---:|---:|
|profile|pde|23.04%|29.60%|
|profile|quasi_steady|23.04%|29.60%|
|profile|one_node|23.03%|29.59%|
|profile|constant_daily|17.54%|24.03%|
|stake|pde|21.08%|29.29%|
|stake|quasi_steady|21.08%|29.28%|
|stake|one_node|21.07%|29.28%|
|stake|constant_daily|17.54%|24.03%|

Primary default clock/current early geometry is the predeclared representative, not the best metadata scenario. Profile-thickness PDE predictions are695.43/608.57/729.38kg/m² against1035/873/688 observed; default errors−32.81/−30.29/+6.01%. Conditional conductivity transfer fails to reproduce the observed site contrast. A source-origin constant daily rate outperforms the three thermal families on aggregate. More distributed thermal physics does not fix heterogeneous effective conductivity.

Numerical64→128mode comparison changes any PDE full-season prediction by at most0.003338kg/m² across432PDE rows; all eight fitted diffusivities remain at the same grid winners. Signed energy is separately checked against the exact thermal-moment identity in SIGNED_ENERGY_AUDIT.json. Numerical convergence does not establish field validity.

## Protocol integrity

Initialv1 scoring mistakenly included the first24h and endpoint in the training objective. Allv1code/results are preserved in legacy_objective_v1. v2 corrected the objective to the original frozen[24h,14d) window; no grid/model/target/scenario changed. The corrected grid winners and resulting predictions happen to be identical. This is a disclosed postscore implementation correction, not a new blind protocol.

## Observation and physics scope

Temperatures identify alpha=k/C, not absolute conductivity. Every thermal family receives its OWN source-ablation-fitted conductivity; no published conductivity or target ablation enters its fit. First14days of target temperatures calibrate target diffusivity as predeclared, so this is conditional thermal transfer, not parameter-free prediction. Dates have no survey hour; all nine clock combinations remain reported. All four author-linked locations are from Kennicott, not Himalayan validation. Only4stake observations underlie1728scenario rows.

2023 leading2+days are unobserved: full-period point scores use phase-matched repetition of the first complete observed24h unit flux. Covered-period loss is reported separately and compared with the FULL mass only as a one-sided necessary check under nonnegative omitted-period melt. It is not a covered-period accuracy score. The maximum-first14day daily rate scaled to gap duration is an engineering scenario, not a proved partial-day upper bound; conductivity scenario endpoints are not statistical intervals.

2020 raw two-probe support covers the stake interval but installation transition is unknown. Earlier/current leading geometries are separately modeled, then explicitly reinitialized at the confirmed processed geometry timestamp. This is a structural initial-state scenario. Negative signed basal heat remains separately reported; positive-flux melt is a declared approximation and does not solve freezing/phase enthalpy. P12 stake-depth geometry excludes31.5cmprobe because it is outside the assumed29cmdebris domain.

The result tests a dry homogeneous conductive sub-debris component with an assumed0°Cbasal boundary. It does not validate wet debris, cliff contact recovery, stream coupling, collapse or glacier-wide melt. No new environmental superiority claim follows.

## Reproduce

`PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python3 run.py`

`ABLATION_MODES=128 PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python3 run.py`

`python3 summarize.py; python3 audit_energy.py; python3 finalize.py`

Raw archives and thermal solver live in sibling directories; no downloads occur during scoring.
