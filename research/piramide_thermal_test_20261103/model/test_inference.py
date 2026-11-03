import math,unittest
from inference import fit_grid,derivative_crh

class InferenceTests(unittest.TestCase):
    def kwargs(self):
        return dict(sensor_ids=['inside'],boundary_sensor_ids=['top','bottom'],train_end_s=2.,evaluation_start_s=3.,alpha_grid_m2_s=[1e-7,2e-7,3e-7],sensor_absolute_bound_C=0.)
    def test_train_only_recovery(self):
        r=fit_grid([0.,1.,2.],[[0.],[2.],[4.]],predict_training=lambda a:[[0.],[a/1e-7],[2*a/1e-7]],**self.kwargs())
        self.assertEqual(r['selected_alpha_m2_s'],2e-7)
        self.assertFalse(r['selected_at_grid_boundary'])
    def test_intervals_do_not_invent_identification(self):
        k=self.kwargs();k['sensor_absolute_bound_C']=10.
        r=fit_grid([0.,1.,2.],[[0.],[2.],[4.]],predict_training=lambda a:[[0.],[a/1e-7],[2*a/1e-7]],**k)
        self.assertEqual(len(r['fully_interval_compatible_grid_m2_s']),3)
        self.assertFalse(r['loss_identifies_unique_grid_point'])
    def test_eval_and_boundary_refusals(self):
        for change in (dict(train_end_s=1.),dict(evaluation_start_s=2.),dict(boundary_sensor_ids=['inside']),dict(alpha_grid_m2_s=[3e-7,1e-7])):
            k=self.kwargs();k.update(change)
            with self.assertRaises(ValueError):fit_grid([0.,1.,2.],[[0.],[2.],[4.]],predict_training=lambda a:[[0.],[2.],[4.]],**k)
    def test_callback_shape_nan(self):
        for callback in (lambda a:[[0.]],lambda a:[[0.],[float('nan')],[4.]]):
            with self.assertRaises(ValueError):fit_grid([0.,1.,2.],[[0.],[2.],[4.]],predict_training=callback,**self.kwargs())
    def test_exact_quadratic_pde_comparator_nonuniform_grid(self):
        a=1.7e-7;t=[0.,200.,600.,1400.];z=[0.,.03,.08,.13]
        # T=z²+2 alpha t solves Tt=alpha Tzz exactly.
        y=[[zz*zz+2*a*tt for zz in z] for tt in t]
        r=derivative_crh(t,y,z,train_end_s=1400.,evaluation_start_s=1500.)
        self.assertAlmostEqual(r['alpha_m2_s'],a,places=18)
        self.assertFalse(r['thermal_conductivity_identified'])
    def test_linear_profile_does_not_identify_alpha(self):
        with self.assertRaises(ValueError):derivative_crh([0.,1.,2.],[[0.,1.,2.]]*3,[0.,1.,2.],train_end_s=2,evaluation_start_s=3)
if __name__=='__main__':unittest.main()
