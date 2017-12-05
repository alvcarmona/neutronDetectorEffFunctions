
from neutron_detector_eff_functions import Detector
from neutron_detector_eff_functions.efftools import efficiency4boron
from neutron_detector_eff_functions.efftools import mg_same_thick
import json
import os


class Detector_test:

    def __init__(self):
        self.detector_multigrid_mono = Detector.Detector.build_detector(15,1,0,[[10,100]], 90, 100, False)
        self.detector_multigrid_poli = Detector.Detector.build_detector(15, 1, 0, [[10, 90],[1,10]], 90, 100, False)
        self.detector_single_mono = Detector.Detector.build_detector(1,1,0,[[10,100]], 90, 100, True)

    def build_detector_test(self):
        assert len(self.detector_multigrid_mono.blades) == 15
        assert len(self.detector_multigrid_poli.blades) == 15
        assert len(self.detector_single_mono.blades) == 1
        d = Detector.Detector.build_detector(4, 1, 0, [[10, 100]], 90, 100, True)
        assert len(d.blades) == 1
        assert self.detector_multigrid_mono.blades[0].backscatter == 1
        assert len(self.detector_multigrid_mono.wavelength) == 1
        assert len(self.detector_multigrid_poli.wavelength) == 2
        assert self.detector_single_mono.single
        assert self.detector_multigrid_mono.single == False

    def json_parser_test(self):
        filepath = ''
        try:
            filepath = os.path.dirname(os.path.abspath(__file__)) + '/test.json'
            with open(filepath, "w") as outfile:
                outfile.write(json.dumps(self.detector_multigrid_mono.to_json(), sort_keys=True, indent=4, ensure_ascii=False))
                outfile.close()
            print('Export')
        except IOError:
            print "Path error"
        d = Detector.Detector.json_parser(filepath)
        assert len(d.blades) == len(self.detector_multigrid_mono.blades)
        c = 0
        for b in d.blades:
            assert b.backscatter == self.detector_multigrid_mono.blades[c].backscatter
        c=0
        for l in d.wavelength:
            assert l == self.detector_multigrid_mono.wavelength[c]
        assert d.single == self.detector_multigrid_mono.single
        assert d.threshold == self.detector_multigrid_mono.threshold
        os.remove(filepath)

    def calculate_eff(self):
        assert 1

    def calculate_sigma(self):
        assert 1

    def calculate_ranges(self):
        assert 1

    def plot_blade_figure(self):
        assert 1

    def plot_blade_figure_single(self):
        assert 1

    def plot_thick_vs_eff(self):
        assert 1

    def plot_thick_vs_eff2(self):
        assert 1

    def plot_wave_vs_eff(self):
        assert 1

    def optimize_thickness_same_test(self):
        assert 1

    def optimize_thickness_diff_mono(self):
        assert 1