import importlib.util, tempfile, unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('digitization',HERE/'extract.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class PixelTests(unittest.TestCase):
    def test_axes_map_anchors_with_subpixel_residual(self):
        for anchors in [module.X,module.Y]:
            m,b=module.fit(anchors)
            self.assertLess(max(abs(m*x+b-y) for x,y in anchors),.5)

    def test_source_substitution_is_refused(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'wrong.png';p.write_bytes(b'wrong')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):module.extract(p)

    def test_separate_station_colors(self):
        self.assertTrue(module.color_match('P3',(44,160,44)))
        self.assertFalse(module.color_match('P5',(44,160,44)))
        self.assertTrue(module.color_match('P5',(214,39,40)))
        self.assertFalse(module.color_match('P3',(214,39,40)))

if __name__=='__main__':unittest.main()
