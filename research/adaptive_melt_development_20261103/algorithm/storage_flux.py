"""Signed basal control-volume energy closure; no hidden melt clipping.

Static P1 geometry, basal T=0C, known constant volumetric C, no other source.
This is a conditional conductive melt-equivalent, not full field phase physics.
"""
import math
L=333700.

def slab_sensible_J_m2(T4,T5,spacing_m,gap_m,capacity_J_m3K):
 if not all(math.isfinite(x)for x in [T4,T5,spacing_m,gap_m,capacity_J_m3K])or spacing_m<=0 or gap_m<0 or capacity_J_m3K<=0:raise ValueError('finite valid geometry/C required')
 return capacity_J_m3K*(spacing_m*T4/8+(3*spacing_m/8+gap_m/2)*T5)

def gradient_integral(times_s,T4,T5,spacing_m):
 if not(len(times_s)==len(T4)==len(T5))or len(times_s)<2 or spacing_m<=0:raise ValueError('shape and spacing')
 if not all(math.isfinite(x)for a in [times_s,T4,T5]for x in a):raise ValueError('finite data')
 if any(b<=a for a,b in zip(times_s[:-1],times_s[1:])):raise ValueError('ordered times')
 return math.fsum((b-a)*((x-y)+(u-v))/(2*spacing_m)for a,b,x,y,u,v in zip(times_s[:-1],times_s[1:],T4[:-1],T5[:-1],T4[1:],T5[1:]))

def predict(training_melt_mmwe,integrated_gradient,storage_changes_J_m2):
 if len(integrated_gradient)!=len(storage_changes_J_m2)or not integrated_gradient or integrated_gradient[0]<=0 or training_melt_mmwe<=0:raise ValueError('positive TRAIN and matching intervals')
 if not all(math.isfinite(x)for x in [training_melt_mmwe]+list(integrated_gradient)+list(storage_changes_J_m2)):raise ValueError('finite data')
 k=(L*training_melt_mmwe+storage_changes_J_m2[0])/integrated_gradient[0]
 if k<=0:raise ValueError('TRAIN implies nonpassive conductivity')
 energies=[k*j-d for j,d in zip(integrated_gradient,storage_changes_J_m2)]
 return {'calibrated_conductivity_W_mK':k,'interval_net_energy_J_m2':energies,'interval_net_melt_equivalent_mmwe':[q/L for q in energies],
         'negative_interval_phase_conversion_unsupported':any(q<0 for q in energies),'no_energy_clipping':True}
