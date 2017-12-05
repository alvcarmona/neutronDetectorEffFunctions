import numpy
from neutron_detector_eff_functions import B10
from neutron_detector_eff_functions.efftools import mg_same_thick, efficiency4boron

class Test_B10:

    #def __init__(self):
     #   self.b = B10()
    def test_works(self):
        assert 1==1

    def test_configurations(self):
        b = B10.B10()
        assert len(b.configurations) == 2
        conf = b.configurations.get('10B4C 2.24g/cm3')
        assert len(conf) == 4
        conf = b.configurations.get('10B4C 2.20g/cm3')
        assert len(conf) == 4

    def test_ranges(self):
        b = B10.B10()
        r1 = b.ranges(200, '10B4C 2.24g/cm3')
        r2 = b.ranges(100, '10B4C 2.24g/cm3')
        assert r1[0] == 2.8
        assert r1[1] == 1.1
        assert r1[2] == 3.7
        assert r1[3] == 1.3
        assert r2[0] == 3.1
        assert r2[1] == 1.35
        assert r2[2] == 4.0
        assert r2[3] == 1.6

    def test_get_th(self):
        b = B10.B10()
        config = b.configurations.get('10B4C 2.24g/cm3')
        assert B10.find_th(config.get('alpha94'), 200000) == 28000.0
        assert B10.find_th(config.get('Li94'), 200000) == 11000.0
        assert B10.find_th(config.get('alpha06'), 200000) == 37000.0
        assert B10.find_th(config.get('Li06'), 200000) == 13000.0

    def test_read_cross_section(self):
        # to test floats I need an aproximation given by numpy :/
        b = B10.B10()
        assert numpy.isclose(b.read_cross_section([[1.8, 100]])[0], [3844.3852472], rtol=1e-05, atol=1e-08, equal_nan=False)

    def test_macro_sigma(self):
        b = B10.B10()
        assert numpy.isclose(b.macro_sigma([3844.3852472]), [0.0398457257908], rtol=1e-05, atol=1e-08, equal_nan=False)

    def test_sigma_eq(self):
        b = B10.B10()
        sigmaeq = b.sigma_eq([0.0398457257908], 90)
        assert numpy.isclose([sigmaeq],  [0.0398457257908], rtol=1e-05, atol=1e-08, equal_nan=False)
        sigmaeq = b.sigma_eq([0.0398457257908], 5)
        assert numpy.isclose([sigmaeq], [0.457178431789], rtol=1e-05, atol=1e-08, equal_nan=False)

    def test_full_sigma_calculation(self):
        b = B10.B10()
        sigma = b.full_sigma_calculation([[1.8,100]], 90)
        assert numpy.isclose([sigma], [0.0398457257908], rtol=1e-05, atol=1e-08, equal_nan=False)
        sigma = b.full_sigma_calculation([[2,100]], 3)
        assert numpy.isclose([sigma], [0.845952315849], rtol=1e-05, atol=1e-08, equal_nan=False)

    # TODO add aluminium consideration
    def test_mg_same_thick(self):
        b = B10.B10()
        sigma = b.full_sigma_calculation([[1.8,100]], 5)
        r1 = b.ranges(200, '10B4C 2.24g/cm3')
        thick = 1
        # check if mg with 1 blade and single blade gives the same value
        assert numpy.isclose(mg_same_thick(sigma[0], r1, thick, 1), [efficiency4boron(thick, r1[0], r1[1], r1[2], r1[3], sigma[0])[0]], rtol=1e-05, atol=1e-08, equal_nan=False)
        assert numpy.isclose( mg_same_thick(sigma[0], r1, thick, 20), [0.69218537268800717], rtol=1e-05, atol=1e-08, equal_nan=False)
