"""Model-validity contract. `python -m ice_integration.release.validity_contract METADATA.json` -> exit 0 only if an ENVIRONMENTAL prediction may be issued.
Anything else may be issued only as a labelled conditional/component hindcast. This checker never fills a missing field."""
import json, sys
REQUIRED = {
 "geometry": ["floor_elevation_survey_dates", "floor_survey_sigma_m", "cross_section_width_source", "width_measured_at_sensor_section", "sensor_clearance_m_or_on_bed_verified", "depth_basis", "section_to_station_binding"],
 "temperature": ["bulk_water_temperature_or_validated_floor_to_bulk_relation", "sensor_calibration_offset_applied_once", "clock_convention_resolved", "observed_support_fraction", "reconstructed_support_labelled"],
 "discharge": ["discharge_basis", "discharge_at_target_section_or_conservation_justified", "discharge_uncertainty"],
 "contact": ["wetted_contact_observed_or_modelled", "missing_record_not_treated_as_zero", "ice_temperature_known_or_bounded"],
 "measurement": ["target_definition_matches_observation_operator", "target_uncertainty_is_source_sigma", "independent_of_calibration_data", "heat_transfer_closure_not_fitted_to_target_event"],
}
def check(meta):
    miss = {g: [k for k in ks if meta.get(g, {}).get(k) in (None, False, "", "unknown")] for g, ks in REQUIRED.items()}
    miss = {g: v for g, v in miss.items() if v}
    return dict(may_issue_environmental_prediction=not miss, unmet=miss, label_if_not="conditional component hindcast; not an environmental prediction" if miss else None)
if __name__ == "__main__":
    r = check(json.load(open(sys.argv[1]))); print(json.dumps(r, indent=1)); sys.exit(0 if r["may_issue_environmental_prediction"] else 2)
