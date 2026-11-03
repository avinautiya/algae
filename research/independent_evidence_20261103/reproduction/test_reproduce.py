import unittest,tempfile,zipfile
from pathlib import Path
from reproduce import percent_error,xlsx_rows,calculate
def workbook(path,rows):
 ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main';strings=[];xml=[]
 for i,row in enumerate(rows,1):
  cells=[]
  for col,value in row.items():
   if isinstance(value,str):strings.append(value);cells.append(f'<c r="{col}{i}" t="s"><v>{len(strings)-1}</v></c>')
   else:cells.append(f'<c r="{col}{i}"><v>{value}</v></c>')
  xml.append('<row>'+''.join(cells)+'</row>')
 with zipfile.ZipFile(path,'w') as z:
  z.writestr('xl/sharedStrings.xml',f'<sst xmlns="{ns}">'+''.join('<si><t>'+s+'</t></si>' for s in strings)+'</sst>')
  z.writestr('xl/worksheets/sheet1.xml',f'<worksheet xmlns="{ns}"><sheetData>'+''.join(xml)+'</sheetData></worksheet>')
class Tests(unittest.TestCase):
 def test_signed_percent(self):self.assertAlmostEqual(percent_error(.75,1),-25.)
 def test_perfect(self):self.assertEqual(percent_error(1.,1.),0.)
 def test_refuse_missing_or_invalid(self):
  for a,b in [(0,0),(1,-1),(float('nan'),1),(True,1),(-1,1)]:
   with self.assertRaises(ValueError):percent_error(a,b)
 def test_synthetic_end_to_end_not_independent_validation(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);workbook(p/'MMIC_model_run_dates.xlsx',[dict(A='CN-5470_20152017',B=1,C=2,F=2016,G=4,H=16,I=1,J=2016,K=4,L=16,M=2)])
   workbook(p/'published_table3.xlsx',[dict(B='DEBCF',D='0.0')])
   out=p/'outputs';(out/'DEBCF/CN').mkdir(parents=True);(out/'Dates_Rows').mkdir()
   (out/'DEBCF/CN/DEBCF_CN_EXP1_EXP2.csv').write_text(','.join(['0.045']+['0']*12)+'\n'+','.join(['0.045']+['0']*12)+'\n')
   (out/'Dates_Rows/Dates_hourly_CN.csv').write_text('2016-04-16 01:00\n2016-04-16 02:00\n')
   src=p/'data.zip'
   with zipfile.ZipFile(src,'w') as z:
    z.writestr('CN-5470_20152017/CN-5470_20152017_60.csv','Y;M;D;H;MIN\n2016;4;16;1;0\n2016;4;16;2;0\n')
    z.writestr('CN-5470_20152017/CN-5470_20152017_VAL_dailysurfaceheight.csv','metadata\nY;M;D;DIST;ABL;ACC;HEIGHT\n2016;4;16;1;.1;0;0\n')
   r=calculate(p,src,out);self.assertTrue(r['published_score_reproduced']);self.assertFalse(r['observation_alignment_certified']);self.assertFalse(r['independent_validation'])
 def test_xlsx_preserves_unicode_minus(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'a.xlsx';workbook(p,[dict(B='DEBCF',D='−23.3')]);self.assertEqual(xlsx_rows(p)[0]['D'],'−23.3')
if __name__=='__main__':unittest.main()
