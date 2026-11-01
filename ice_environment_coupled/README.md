# Coupled environmental cliff forecast — integration milestone

This additive package connects hydraulic routing, finite water heat, atmospheric
energy ports, cold/debris storage, fresh-ice warming and separate face/floor
retreat. It changes no active Claude solver, checkpoint or frozen campaign.

**16 new checks pass.** Five of six synthetic development configurations complete;
the sixth is explicitly refused because the route/contact closures require more
head than the upstream pool supplies. No run is silently dropped. These checks
establish implementation properties, not environmental accuracy or novelty.

## Execute

```sh
python3 verify.py
python3 model.py example.json example_result.json
python3 run_development.py --plot
python3 data/fetch_meteorology.py
python3 data/audit_meteorology.py
python3 data/plot_weather.py
```

Core model/tests use only the standard library. Plotting requires Matplotlib.
Meteorological acquisition requires network access; tests do not. Vendored
components are byte-pinned copies from the earlier environmental package and
verified before tests.

## Coupling and conservation

1. Prescribed inflow is allocated across a toe route and a bypass using monotone
   effective ratings at a common upstream stage. Obstruction height and an
   under-rock aperture/Darcy option are explicit scenario inputs.
2. Downstream face wetting is evaluated independently using rectangular Manning
   resistance with an explicitly imposed friction slope. It is not set equal to
   upstream pond level. A necessary energy-head check refuses incompatible
   combinations. This does not resolve the backwater profile or prove that two
   effective ratings represent a real boulder deposit.
3. One sensible-water heat budget is depleted serially across the bare floor
   and wetted bare face. Bypass water has its own allocated share. Atmosphere
   acts only on dry surfaces in this model.
4. The upper face may carry dry debris, with finite heat storage and conduction
   to a zero-C reference interface. Water reaching that covered face is refused
   until a compatible submerged-debris thermal model exists. The obstruction
   itself is a thin effective upstream control, not a resolved thermally active
   boulder blanket covering the downstream floor.
5. Surface cold debt is an extensive J reservoir. Additional fresh-ice warming
   costs `rho_i cp_i DeltaT` per removed m³. Only remaining latent energy enters
   geometry. The finite surface debt is an explicitly separate extra reservoir;
   do not initialize it as the same cold material already charged volumetrically.
6. New area receives declared incoming debris enthalpy/cold debt; storage is
   never manufactured by keeping the old J/m² state on a larger area. Face
   fractions remain lumped prescribed reservoirs, not advected material masks.
7. Face retreat and floor lowering remove disjoint geometric volumes. Meltwater
   is accounted downstream at the melting point. It does not feed back into
   the steady hydraulic solve; a declared maximum phase-water/inflow ratio is
   enforced. Atmospheric melting at zero inflow is currently refused by this
   scope check rather than silently ignoring meltwater-generated flow.

The global check is

`E_atm + E_water,in - E_water,out + H_area_transfer
 = Delta(H_debris - C_surface) + E_fresh_ice_warming + E_latent`.

Upstream/contact hydraulics are sampled at the beginning of a step; geometry
and cumulative melt are at its end. Both timestamps are labelled. Forcing and
removal boundaries are hit exactly. Time refinement remains necessary.

## Actual environmental input acquired

The public [CliffEBM](https://github.com/pburi/CliffEBM) Lirung AWS example is
pinned at commit `ea5ca8bae9b3c67ec3cf2a8b56c724d87bb1909c`; original Git blob
and SHA256 are verified. It has **4,969 hourly records, 8 May–1 December 2013**.
The audit preserves 1,435 negative-shortwave readings and flags four humidity
values outside 0–100%. Precipitation has 926 missing values and its other values
are zero; it must not be treated as a validated rainfall record.

This is genuine station meteorology, not measured face radiation or stream
forcing. Timezone is not inferred. The model currently requires independently
prepared net atmospheric surface fluxes. Direct/diffuse radiation, terrain
shading, face aspect, pressure and bulk-transfer assumptions must be supplied
before deriving those fluxes from station measurements. Air temperature is
never substituted for stream temperature or runoff.

The 2013 station interval overlaps the earlier acquired Lirung geometry survey
season. This supports assembling a realistic seasonal forcing experiment;
cliff-ID correspondence, exposure masks and colocated stream forcing still
require verification. A common glacier and year do not establish colocated
hydrological validation.

## Development result and interpretation

The synthetic base predicts about 0.03044 m³ face melt over 24 hours at a
30-minute coupling step. Reducing 15 to 7.5 minutes changes face melt by about
0.064%. All geometry, weather-flux, runoff and closure values in that example
are assumed. The results are a coupling demonstration, not a glacier forecast,
independent skill test or proof of superior prediction.

The professional-model comparison protocol is in `benchmark/`. It separates
agreement with shared closures, reproduction of public professional code, and
independent predictive skill. The specific candidate contribution remains a
consequential face-melt error from prescribing hydraulic contact after clearing.
That result must survive defensible environmental and transmission uncertainty.

Collapse, nonhydrostatic flow, full face-shape evolution and glacier advection
are not implemented here. The independent audit and review are included so the
coupled model's scope stays reviewable. Adding PINNs would not remove any of
these physical or evidential requirements.
