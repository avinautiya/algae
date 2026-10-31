# Failures and deviations

F1. First launch of the 12 resolved-reference test runs (version events_v1) failed after ~1 s: `run_reference()` has no `atol_energy` option (the reference solver is not an LSODA energy-state integrator). Records are kept (status failed). Relaunched as version events_v1b with the option removed for kind=reference; frozen configurations, forcings and metrics unchanged (cases_test_ref.json).
