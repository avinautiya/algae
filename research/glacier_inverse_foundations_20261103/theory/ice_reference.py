"""Conditional cold-ice reference inversion; established thermal-admittance math.

No historical novelty or field-accuracy assertion. The formula requires a
complete linear no-melt experiment with known initial/bottom ice conditions.
Cutting cold snippets out of a melt history does not satisfy that requirement.
"""
import math

def positive(x,name):
    if type(x) not in (int,float) or not math.isfinite(x) or x<=0:raise ValueError(name+' must be finite and positive')
    return float(x)

def impedance(s,*,ice_conductivity_W_mK,ice_capacity_J_m3K,ice_depth_m):
    s=positive(s,'Laplace rate')
    k=positive(ice_conductivity_W_mK,'ice conductivity');C=positive(ice_capacity_J_m3K,'ice volumetric capacity')
    if ice_depth_m==math.inf:return 1/math.sqrt(k*C*s)
    H=positive(ice_depth_m,'ice depth')
    return math.tanh(H*math.sqrt(s*C/k))/math.sqrt(k*C*s)

def check_context(context):
    required=['entire_transformed_history_nonmelting','initial_ice_temperature_anomaly_zero',
        'bottom_temperature_anomaly_zero','common_clock_and_transform_convention',
        'debris_boundary_ratio_at_known_ice_interface','ice_properties_independent_of_melt_target']
    if not isinstance(context,dict) or any(context.get(k)is not True for k in required):
        raise ValueError('Reference assumptions must be explicit; selected cold snippets/unknown initial ice are unsupported')

def identify(A1_m,A2_m,s1,s2,*,diffusivity_m2_s,ice_conductivity_W_mK,
             ice_capacity_J_m3K,ice_depth_m,context):
    check_context(context)
    a1=positive(A1_m,'first boundary ratio');a2=positive(A2_m,'second boundary ratio')
    alpha=positive(diffusivity_m2_s,'debris diffusivity')
    kwargs=dict(ice_conductivity_W_mK=ice_conductivity_W_mK,ice_capacity_J_m3K=ice_capacity_J_m3K,ice_depth_m=ice_depth_m)
    z1,z2=impedance(s1,**kwargs),impedance(s2,**kwargs);dz=z1-z2
    if dz==0:raise ValueError('Distinct reference impedances required')
    k=(a1-a2)/dz
    if not math.isfinite(k) or k<=0:raise ValueError('Data incompatible with positive conductivity; do not clip')
    ell=a1-k*z1;R=ell/k
    tol=1e-12*max(a1,a2,abs(k*z1),1.)/k
    if R < -tol:raise ValueError('Data incompatible with nonnegative contact resistance')
    return dict(conductivity_W_mK=k,volumetric_capacity_J_m3K=k/alpha,
        contact_resistance_m2K_W=max(0.,R),robin_length_m=max(0.,ell),
        reference_impedances_m2K_W=[z1,z2],reference_impedance_contrast_m2K_W=abs(dz),
        assumptions=context,measurement_uncertainty_propagated=False,field_validation=False,
        historical_mathematical_novelty_claim=False)

def enclosure(A1_interval,A2_interval,s1,s2,*,diffusivity_interval_m2_s,
              ice_conductivity_W_mK,ice_capacity_J_m3K,ice_depth_m,context):
    """Conservative parameter box from conditional boundary-ratio intervals.

    Positive Cartesian ratio intervals ignore correlation conservatively. The
    contact interval is intersected with R>=0; no noisy value is made physical.
    """
    check_context(context)
    pairs=[]
    for pair,name in [(A1_interval,'A1'),(A2_interval,'A2'),(diffusivity_interval_m2_s,'diffusivity')]:
        if len(pair)!=2:raise ValueError('two interval endpoints required')
        lo,hi=[positive(x,name)for x in pair]
        if hi<lo:raise ValueError('ordered intervals required')
        pairs.append((lo,hi))
    kw=dict(ice_conductivity_W_mK=ice_conductivity_W_mK,ice_capacity_J_m3K=ice_capacity_J_m3K,ice_depth_m=ice_depth_m)
    z1,z2=impedance(s1,**kw),impedance(s2,**kw)
    (l1,u1),(l2,u2),(al,au)=pairs
    if z1<z2:z1,z2,l1,u1,l2,u2=z2,z1,l2,u2,l1,u1
    dz=z1-z2
    if dz==0:raise ValueError('Distinct reference impedances required')
    kl,ku=(l1-u2)/dz,(u1-l2)/dz
    if ku<=0:return dict(status='incompatible_positive_conductivity',field_validation=False)
    if kl<=0:return dict(status='positive_conductivity_not_resolved_by_bounds',conductivity_upper_W_mK=ku,
        capacity_upper_J_m3K=ku/al,finite_contact_upper_established=False,
        uncertainty='Conditional supplied intervals; not statistical confidence',field_validation=False)
    # With A1>A2>0, R=(Z1*A2-Z2*A1)/(A1-A2) decreases in A1 and increases in A2.
    rl=(z1*l2-z2*u1)/(u1-l2);ru=(z1*u2-z2*l1)/(l1-u2)
    if ru<0:return dict(status='incompatible_nonnegative_contact_resistance',field_validation=False)
    return dict(status='conditionally_enclosed',conductivity_W_mK=[kl,ku],capacity_J_m3K=[kl/au,ku/al],
        contact_resistance_m2K_W=[max(0.,rl),ru],exact_attainable_joint_region_claim=False,
        uncertainty='Conditional supplied intervals; not statistical confidence',field_validation=False)
