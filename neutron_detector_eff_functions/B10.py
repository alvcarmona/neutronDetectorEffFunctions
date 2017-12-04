
import os
import sys
import math
import numpy as np
from bisect import bisect_left
import matplotlib.pyplot as plt
from scipy import interpolate

class B10:
    configurations = {}

    def __init__(self, parent=None):
        self.configurations = {
            '10B4C 2.24g/cm3': {
            'alpha06': np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C224/IONIZ_Linkoping_Alpha06.py", unpack=True, skiprows=27),
            'alpha94': np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C224/IONIZ_Linkoping_Alpha94.py", unpack=True, skiprows=27),
            'Li06': np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C224/IONIZ_Linkoping_Li06.py", unpack=True, skiprows=27),
            'Li94': np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C224/IONIZ_Linkoping_Li94.py", unpack=True, skiprows=27),
        },
            '10B4C 2.20g/cm3': {
            'alpha06': np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C220/IONIZ_Alpha06.py", unpack=True, skiprows=28),
            'alpha94': np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C220/IONIZ_Alpha94.py", unpack=True, skiprows=28),
            'Li06': np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C220/IONIZ_Li06.py", unpack=True, skiprows=28),
            'Li94': np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/10B4C220/IONIZ_Li94.py", unpack=True, skiprows=28),
        }

        }

    #
    def ranges(self, threshold, name):
        """calculates range of 4 particles given a threshold and configuration of B10.

        Args:
            threshold (int):  Threshold in KeV.
            name (String):  Name of the configuration of B10. Example: '10B4C 2.24g/cm3'
        Returns:
            r[1,4]: r(1,1)=ralpha_94;
                 r(1,2)=rLi_94;
                 r(1,3)=ralpha_06;
                 r(1,4)=rLi_06;

            ranges in um

        ..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/B10tools/rangesMAT.m?at=default&fileviewer=file-view-default

        """
        if threshold < 2:
            threshold = 2
            print ('Threshold set to 2 keV, no lower Th. accepted!')

        if threshold > 779:
            threshold = 779
            print('Threshold set to 779 keV, no higher Th. accepted!')

        r = [0,0,0,0]
        threshold *= 1000
        config = self.configurations.get(name)
        Ra94 = find_th(config.get('alpha94'), threshold)
        r[0] = Ra94/10000
        Rli94 = find_th(config.get('Li94'), threshold)
        r[1] = Rli94/10000
        Ralpha06 = find_th(config.get('alpha06'), threshold)
        r[2] = Ralpha06/10000
        Rli06 = find_th(config.get('Li06'), threshold)
        r[3] = Rli06/10000
        return r

    @staticmethod
    def read_cross_section(lambdalist):
        """calculates cross section for a list of lambdas

        Args:
            lambdalist (list):  List of lambda values in Amstrong.
        Returns:
            sigma(List): cross section list for lambda list

        ..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/B10tools/readCrossSect.m?at=default&fileviewer=file-view-default

        """
        ht = 1.054e-34
        nmass = 1.67e-27
        # lambdaList transformed to ev
        lamen = []
        sigma = []
        c0 = 0
        for l in lambdalist:
            lamen.append(((ht ** 2) * 4 * math.pi ** 2) / (2 * nmass * (l[0] ** 2)) * 1e20 * 6.24e18)
        x, y = np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/B10/B10CrossSect_(n,a).py", delimiter=',', unpack=True)
        xlog = []
        ylog = []
        lamenlog =[]
        for v in x:
            xlog.append(math.log10(v))
        for v in y:
            ylog.append(math.log10(v))
        for v in lamen:
            lamenlog.append(math.log10(v))
        f = interpolate.interp1d(xlog, ylog, kind='cubic')
        for v in lamenlog:
            sigma.append(10**f(v))
        # check if interpolation works
        # xnew = np.linspace(xlog[0], xlog[100], num=41, endpoint=True)
        # plt.plot(xlog, ylog, 'o', xnew, f(xnew), '--')
        plt.show()
        return sigma

    @staticmethod
    def macro_sigma(sigmainfin):
        """calculates macro sigma for a sigma infinitesimal

        Args:
            sigmainfin (list): cross section list for lambda list
        Returns:
            sigma(List):

        ..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/B10tools/macroB10sigma.m?at=default&fileviewer=file-view-default

        """
        # 2.242 g/cm3
        density = 2.38 * (52 / 55.2)
        composition = [[0.793, 0.025, 0.17, 0.007, 0.004, 0.001], [10, 11, 12, 1, 16, 14]]
        # should be 1.0
        suma = math.fsum(composition[0])
        mmol = np.dot(composition[0], composition[1])
        b10perc = composition[0][0]
        sigma = []
        for s in sigmainfin:
            sigma.append((density*6.022e23/mmol)*b10perc*s*1e-24*1e-4)
        return sigma

    @staticmethod
    def sigma_eq(macrosigma, theta):
        """calculates sigma equivalent of Boron for a macro sigma and an angle

        Args:
            macrosigma (list): macro sigma of B10
            theta (int): incident angle of neutron in degrees
        Returns:

            sigma (float): sigma equivalent
        ..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/B10tools/macroB10sigma.m?at=default&fileviewer=file-view-default

        """
        sigmaeq = []
        for m in macrosigma:
            sigmaeq.append(m/math.sin(math.radians(theta)))
        return sigmaeq

    def full_sigma_calculation(self, lambd, theta):
        sigma = self.sigma_eq(self.macro_sigma(self.read_cross_section(lambd)), theta)
        return sigma





# TODO separate integral and threshold point calculation in different functions
def find_th(array, threshold):
    """calculates integral and x value for a y in the integral data matrix

    Args:
        threshold (int):  y point in the integral to find an x.
        array (array[2]):  Array of 2 arrays, [0]contains array of x values, [1] contains array of y values
    Returns:
        thValue: point where y in the integral have the closest value to threshold



    ..  Original source in Matlab: https: // bitbucket.org / europeanspallationsource / dg_matlabborontools / src / bcbac538ad10d074c5150a228847efc2e0269e0d / B10tools / rangesMAT.m?at = default & fileviewer = file - view - default

    """
    x= array[0]
    E1= array[1]
    E3 = []
    c1 = 0
    i = 0
    w = x[1] - x[0]
    c0 = 0
    thValue=0
    thFound= False
    #Total integral value 'i' calculation
    for y in E1:
        if c0 != 0:
            # interpolation of distance from previous point
            # all points seem to be at the same distance!
            w = x[c0] - x[c0 - 1]
        i = i + y * w
        c0 += 1

    for a in x:
        area = 0
        c2 = 0
        c0 = 0
        w = x[1] - x[0]
        for s in range(0, c1):
            if c0 != 0:
                # interpolation of distance from previous point
                w = x[c0] - x[c0 - 1]
            area = area + E1[c2] * w
            c2 += 1
            c0 += 1
        localArea = i - area
        if thFound == False:
            if localArea <= threshold:
                thFound = True
                thValue = x[c2]
        # check if this is the threshold point
        E3.append(localArea)
        c1 += 1
    return thValue


if __name__ == '__main__':
    b = B10()
    print (b.ranges(100, '10B4C 2.24g/cm3'))
    #print b.read_cross_section([1.8, 3, 6])
