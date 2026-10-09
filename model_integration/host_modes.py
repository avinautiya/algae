"""Host-coupling modes with double-counting prevention (Phase 5 item 6).

Mode A, forward parameterisation
  The host computes its own bare-ice albedo from its own state. The algae optics supply the algal
  albedo change for a PRESCRIBED abundance (no bloom-growth model exists in this repository, so
  abundance is an input, not a prognostic variable).
  Allowed:
    * host algae treatment 'absent' -> 'anomaly' (host + a_with - a_without) or 'replace';
    * host algae treatment 'explicit' -> 'replace' only, and only if the host's own algae scheme
      is declared switched off.
  Refused: a host whose albedo is observed ('implicit_observed'): observed albedo already contains
  the algae.

Mode B, satellite-informed state estimation
  The albedo comes from observation and is used unchanged. The algae optics only ATTRIBUTE part of
  the observed darkening (counterfactual no-algae albedo = observed + delta). Scores of Mode B
  against the same or co-registered observations are not independent validation, and every output
  carries that flag.

The host's algae treatment is a REQUIRED declaration; it is never inferred.

The named host for real tests is the reference point SEB `phase4/seb.py` (KAN_M forcing). There is
no E3SM/MAR integration here.

Output quantities are labelled distinctly:
  * modelled_increment_mwe: difference of two SEB runs;
  * potential_melt_mwe: sum of SW x delta alpha / L_f;
  * observed_ablation: never produced here;
  * runoff: not modelled.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .coupling import couple_host_albedo, _array

TREATMENTS = ("absent", "explicit", "implicit_observed")
ALBEDO_SOURCES = ("prognostic", "observed")
NOT_INDEPENDENT = "state estimation: constrained by the same observations; not an independent validation"


@dataclass(frozen=True)
class HostContract:
    host: str
    algae_treatment: str                 # absent | explicit | implicit_observed  (required)
    albedo_source: str                   # prognostic | observed
    own_algae_scheme_disabled: bool = False
    note: str = ""

    def __post_init__(self):
        if self.algae_treatment not in TREATMENTS:
            raise ValueError(f"declare algae_treatment as one of {TREATMENTS}")
        if self.albedo_source not in ALBEDO_SOURCES:
            raise ValueError(f"declare albedo_source as one of {ALBEDO_SOURCES}")
        if self.albedo_source == "observed" and self.algae_treatment != "implicit_observed":
            raise ValueError("an observed host albedo already contains any algae: algae_treatment must be "
                             "'implicit_observed'")


def check_mode(contract: HostContract, mode: str, method: str | None = None):
    """Raise unless the combination cannot double count. Returns the validation-status label."""
    if mode == "A":
        if contract.albedo_source != "prognostic":
            raise ValueError("Mode A needs a host that computes its own albedo; an observed albedo -> Mode B")
        if method not in ("anomaly", "replace"):
            raise ValueError("Mode A method must be 'anomaly' or 'replace'")
        if contract.algae_treatment == "explicit":
            if method != "replace" or not contract.own_algae_scheme_disabled:
                raise ValueError("host with an explicit algae scheme: only 'replace' with its own scheme "
                                 "declared switched off (otherwise algae are counted twice)")
        return "forward parameterisation (prescribed abundance); conditional on host and optics"
    if mode == "B":
        if contract.albedo_source != "observed":
            raise ValueError("Mode B uses an observed albedo")
        return NOT_INDEPENDENT
    raise ValueError("mode must be 'A' or 'B'")


def mode_a_albedo(contract, host_albedo, modeled_without, modeled_with, method):
    status = check_mode(contract, "A", method)
    treatment = contract.algae_treatment
    if treatment == "explicit":                      # host scheme switched off -> behaves as algae-free
        treatment = "absent"
    a = couple_host_albedo(host_albedo, modeled_without, modeled_with, mode=method,
                           host_algae_treatment=treatment)
    return a, dict(mode="A", method=method, validation_status=status, contract=asdict(contract))


def mode_b_attribution(contract, observed_albedo, dalpha_algae):
    """Observed albedo unchanged; counterfactual algae-free albedo = observed + delta (signed: a
    negative delta, i.e. brightening by scattering, is kept)."""
    status = check_mode(contract, "B")
    obs = _array(observed_albedo, "observed_albedo", 0, 1)
    d = _array(dalpha_algae, "dalpha_algae")
    if d.shape != obs.shape:
        raise ValueError("delta and observed albedo must have identical shapes")
    cf = obs + d
    flag = (cf < 0) | (cf > 1)
    return dict(albedo=obs.copy(), counterfactual_no_algae=np.where(flag, np.nan, cf),
                flag_counterfactual_out_of_range=flag, mode="B", validation_status=status,
                contract=asdict(contract))


def reference_seb_paired(contract, forcing, alpha_with, alpha_without, *, mode, method="anomaly",
                         host_albedo=None, **seb_kw):
    """Paired runs of the reference SEB (phase4/seb.py) on identical forcing.

    Mode A: alpha_with/alpha_without are the optics model's own pair; the host albedo is
    host_albedo (anomaly) or the optics albedo (replace).
    Mode B: alpha_with is the observed albedo (used as is); alpha_without = observed + delta.
    """
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase4"))
    import seb
    if mode == "A":
        h = alpha_without if host_albedo is None else host_albedo
        on, meta = mode_a_albedo(contract, h, alpha_without, alpha_with, method)
        off = _array(h if method == "anomaly" else alpha_without, "host albedo without algae", 0, 1)
    else:
        b = mode_b_attribution(contract, alpha_with, np.asarray(alpha_without, float) - np.asarray(alpha_with, float))
        if b["flag_counterfactual_out_of_range"].any():
            raise ValueError("counterfactual albedo outside [0, 1]; the attribution is not physical here")
        on, off = b["albedo"], b["counterfactual_no_algae"]
        meta = {k: b[k] for k in ("mode", "validation_status", "contract")}
    r = seb.paired_algae(None, on, off, forcing=forcing, **seb_kw)
    out = dict(modelled_increment_mwe=r["modelled_algal_melt_increment_mwe"],
               potential_melt_mwe=r["potential_algal_melt_mwe"],
               modelled_melt_with_mwe=r["melt_on_mwe"], modelled_melt_without_mwe=r["melt_off_mwe"],
               valid_hours=r["valid_hours"], missing_hours=r["missing_hours"],
               observed_ablation="not produced here (see phase4/seb.validate)", runoff="not modelled",
               host="reference point SEB phase4/seb.py", **meta)
    return out
