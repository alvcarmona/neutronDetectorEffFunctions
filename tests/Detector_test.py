
import json
import os
import pytest
from neutron_detector_eff_functions import Detector
import numpy as np


class Test_Detector:


     def test_build_detector(self):
        detector_multigrid_mono = Detector.Detector.build_detector(15, 1, 0, [[10, 100]], 90, 100, False, '10B4C 2.20g/cm3')
        detector_multigrid_poli = Detector.Detector.build_detector(15, 1, 0, [[10, 90], [1, 10]], 90, 100, False, '10B4C 2.20g/cm3')
        detector_single_mono = Detector.Detector.build_detector(1, 1, 0, [[10, 100]], 90, 100, True, '10B4C 2.20g/cm3')
        assert len(detector_multigrid_mono.blades) == 15
        assert len( detector_multigrid_poli.blades) == 15
        assert len( detector_single_mono.blades) == 1
        d = Detector.Detector.build_detector(4, 1, 0, [[10, 100]], 90, 100, True, '10B4C 2.20g/cm3')
        assert len(d.blades) == 1
        assert  detector_multigrid_mono.blades[0].backscatter == 1
        assert len( detector_multigrid_mono.wavelength) == 1
        assert len( detector_multigrid_poli.wavelength) == 2
        assert  detector_single_mono.single
        assert  detector_multigrid_mono.single == False

     def test_json_parser(self):
        detector_multigrid_mono = Detector.Detector.build_detector(15, 1, 0, [[10, 100]], 90, 100, False, '10B4C 2.20g/cm3')
        detector_multigrid_poli = Detector.Detector.build_detector(15, 1, 0, [[10, 90], [1, 10]], 90, 100, False, '10B4C 2.20g/cm3')
        detector_single_mono = Detector.Detector.build_detector(1, 1, 0, [[10, 100]], 90, 100, True, '10B4C 2.20g/cm3')
        filepath = ''
        try:
            filepath = os.path.dirname(os.path.abspath(__file__)) + '/test.json'
            with open(filepath, "w") as outfile:
                outfile.write(json.dumps( detector_multigrid_mono.to_json(), sort_keys=True, indent=4, ensure_ascii=False))
                outfile.close()
            print('Export')
        except IOError:
            print ("Path error")
        d = Detector.Detector.json_parser(filepath)
        assert len(d.blades) == len( detector_multigrid_mono.blades)
        c = 0
        for b in d.blades:
            assert b.backscatter ==  detector_multigrid_mono.blades[c].backscatter
        c=0
        for l in d.wavelength:
            assert l ==  detector_multigrid_mono.wavelength[c]
        assert d.single ==  detector_multigrid_mono.single
        assert d.threshold ==  detector_multigrid_mono.threshold
        os.remove(filepath)

     def test_calculate_eff(self):
        detector_multigrid_mono = Detector.Detector.build_detector(15, 1, 0, [[10, 100]], 90, 100, False,'10B4C 2.24g/cm3')
        detector_multigrid_poli = Detector.Detector.build_detector(15, 1, 0, [[10, 90], [1, 10]], 90, 100, False,'10B4C 2.24g/cm3')
        detector_single_mono = Detector.Detector.build_detector(1, 1, 0, [[10, 100]], 90, 100, True, '10B4C 2.24g/cm3')
        assert detector_multigrid_mono.calculate_eff()[1] == pytest.approx(73.71, 0.01)
        assert detector_multigrid_poli.calculate_eff()[1] == pytest.approx(69.9, 0.01)
        assert detector_single_mono.calculate_eff()[0][0] == pytest.approx(14.8, 0.01)
        assert detector_single_mono.calculate_eff()[0][1] == pytest.approx(14.4, 0.01)
        assert detector_single_mono.calculate_eff()[1] == pytest.approx(26.4, 0.01)

     def test_calculate_phs(self):
        detector_multigrid_poli = Detector.Detector.build_detector(15, 1, 0, [[10, 90], [1, 10]], 90, 100, False,'10B4C 2.24g/cm3')
        detector_multigrid_poli.calculate_eff()
        detector_multigrid_poli.calculate_phs()
        eff = detector_multigrid_poli.metadata.get('eff')[0]
        phs = detector_multigrid_poli.metadata.get('phs')[1]
        efftr = np.sum(eff[0],axis=0)
        effbs = np.sum(eff[1], axis=0)
        phsb = np.sum(phs[0],axis=0)
        phsc = np.sum(phs[1], axis=0)
        print( 'trans '+str(efftr)+'. BS: '+str(effbs))
        print('trans ' + str(phsb) +'  '+ str(phsc))
        assert 0

     def test_calculate_sigma(self):
        assert 1

     def test_calculate_ranges(self):
        assert 1

     def test_plot_blade_figure(self):
        assert 1

     def test_plot_blade_figure_single(self):
        assert 1

     def test_plot_thick_vs_eff(self):
        assert 1

     def test_plot_thick_vs_eff2(self):
        assert 1

     def test_plot_wave_vs_eff(self):
        assert 1

     def test_optimize_thickness_same_test(self):
        assert 1

     def test_optimize_thickness_diff_mono(self):
        assert 1