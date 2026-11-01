import unittest
from model import Case,evolve
from fixed_stock_forecast import predict_loss

class FixedStockTests(unittest.TestCase):
    def test_quadrature_matches_time_integrator(self):
        c=Case(discharge=.07,initial_debris_depth=.16,source_fractions=(0.,0.,1.),allow_supply=False,allow_export=False)
        stock=[0.,0.,c.initial_debris_depth*(1-c.porosity)*c.width*c.length]
        pred=predict_loss(c,stock)
        r=evolve([c],days=12,dt=225,keep_history=False)
        self.assertLess(abs(pred['loss_time_days']-r['first_dry_day'][0]),225/86400)
        lo,hi=pred['mathematical_time_bracket_days']
        self.assertLessEqual(lo,pred['loss_time_days']);self.assertLessEqual(pred['loss_time_days'],hi)
    def test_refuses_active_source(self):
        with self.assertRaises(ValueError):predict_loss(Case(),[0.,0.,.1])
    def test_no_bed_melt_no_event(self):
        c=Case(allow_supply=False,allow_ice_incision=False)
        self.assertIsNone(predict_loss(c,[0.,0.,.01])['loss_time_days'])
if __name__=='__main__':unittest.main()
