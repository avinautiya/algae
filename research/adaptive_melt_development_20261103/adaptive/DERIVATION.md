# Causal calibrated weak-energy algorithm and its limits

Let T_j be temperatures on a one-dimensional P1 thermal graph. With fixed geometry/capacity, consistent finite-element nodal conservation over a window [a,b] gives

    C M (T(b)−T(a)) + K(G) integral_a^b T dt = 0

at interior nodes; boundary temperatures are supplied. Each edge contributes equal and opposite heat to its two nodes. G_e>0 is a passive edge conductance. The code estimates G from multiple six-hour integrated balances with bounded least squares. Source sensors enter both regressors and storage: this is errors-in-variables inference, not Gaussian parameter identification. Independent synthetic tests of this kernel are in the sibling inverse-foundations package.

Initially G is inferred from the full local14-day TRAIN temperature interval. At each later day, only the preceding3,7or14days of temperatures are used;7days is primary. The next day's signed lower-edge heat is G_b (T_deep−0) dt, conditional on a fixed melting boundary20mm below that probe. This is a contemporaneous state-conditioned inference, not a future weather forecast. Neither material conductivity nor absence of unresolved heat advection is established by a good residual.

The phase accumulator repays negative heat (cold debt) before assigning positive heat to latent release. Its exact discrete identity is sum(q)=sum(latent_release)+cold_final, cold_final<=0. It is an accounting closure, not a resolved cold-ice temperature field. It does not feed interface cold temperature back into the observed profile or fitted G. Calling it a fully coupled Stefan solver would be false.

## Heat-scale quotient

If every heat/capacity coefficient shares an unknown constant positive multiplier s, the temperature-only system is unchanged. The phase accumulator is positively homogeneous: release(s*q)=s*release(q). A single TRAIN latent amount Y0 therefore fixes a common scale via g=Y0/F0, and the predicted amount is Yi=Y0*Fi/F0. The common unknown multiplier cancels in this quotient.

This elementary invariance does NOT cancel an arbitrary time-varying multiplier, geometry change, unmodeled heat source, error in the known TRAIN amount, or sensor error. It also does not give the actual bounded inference exact scale invariance when conductivity constraints are held fixed while C is changed. Bounds must be scaled together, or remain inactive, for that estimator property. No new general mathematical theorem is claimed.

One diagnostic interpretation of g is C_effective=g*C_nominal. On these development data the inferred effectiveC is about0.998MJ/m3/K atTower1 and1.159MJ/m3/K atTower2. A calibration gain is not independent measurement of capacity; those values must be compared to external material evidence. Temperature fit plus local output calibration cannot certify physical transfer.

## Actual development result

Using corrected source thermal depths, primary7day MAE is25.805mm water equivalent atTower1 and5.339atTower2. Air-only gives7.906/7.265 and calendar gives2.253/3.435. Thus the algorithm improves over air atTower2 but fails atTower1 and against calendar atboth. This is a useful conditional gain at one already inspected stake, not validated superiority. All3window lengths remain in RESULTS.json. No parameter is chosen using Tower3 outcomes.

CORRECTIONS.md records an initial Tower2 depth error and retains that diagnostic output. Only corrected source-depth runs belong in comparisons. Neither error nor correction touched prior frozen evidence.
