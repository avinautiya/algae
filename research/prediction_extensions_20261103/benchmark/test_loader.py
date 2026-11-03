import datetime as dt,unittest
from loader import parse,scalar,time_utc,FORBIDDEN,ALLOWED
class Tests(unittest.TestCase):
 def test_explicit_nepal_offset(self):self.assertEqual(time_utc(dict(YEAR='2016',MONTH='1',DAY='2',HOUR='0',MINUTE='0')),'2016-01-01T18:15:00+00:00')
 def test_source_end_of_day24_not_machine_timezone(self):self.assertEqual(time_utc(dict(YEAR='2016',MONTH='1',DAY='1',HOUR='24',MINUTE='0')),'2016-01-01T18:15:00+00:00')
 def test_missing_is_not_zero(self):self.assertIsNone(scalar('NaN'));self.assertEqual(scalar('0'),0)
 def test_source_decimal_comma_scientific_literal(self):self.assertEqual(scalar('-5,00E-06'),-5e-6)
 def test_infinity_refuses(self):
  with self.assertRaises(ValueError):scalar('inf')
 def test_no_target_derived_inputs(self):self.assertFalse(FORBIDDEN.intersection(ALLOWED))
 def test_duplicate_and_unsorted_refuse(self):
  for raw in [b'YEAR;MONTH;DAY;T;T\n2016;1;1;0;0',b'YEAR;MONTH;DAY\n2016;1;2\n2016;1;1']:
   with self.assertRaises(ValueError):parse(raw)
 def test_header_offset_and_no_timezone_guessing(self):
  raw=b'description\nYEAR;MONTH;DAY;HOUR;MINUTE;T\n2016;1;2;0;0;NaN'
  r=parse(raw,1);self.assertEqual(r[0][0],'2016-01-01T18:15:00+00:00');self.assertIsNone(scalar(r[0][1]['T']))
if __name__=='__main__':unittest.main()
