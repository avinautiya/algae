import unittest
import pandas as pd
from secondary2023 import normalize_pair

class TransferTests(unittest.TestCase):
    def test_year_specific_column_collisions_do_not_drop_training_year(self):
        a=pd.DataFrame({'Date_Time':['a','b'],'T_3m':[1.,3.]})
        first=pd.DataFrame({'Date_Time':['a'],'T_2m':[2.]})
        second=pd.DataFrame({'Date_Time':['b'],'T_2m':[4.],'T_3m':[999.]})
        joined=pd.concat([normalize_pair(first,a,{'T_2m':'T_3m'}),normalize_pair(second,a,{'T_2m':'T_3m'})])
        self.assertEqual(joined['T_2m__x'].tolist(),[1.,3.])
        self.assertEqual(joined['T_2m__y'].tolist(),[2.,4.])
if __name__=='__main__':unittest.main()
