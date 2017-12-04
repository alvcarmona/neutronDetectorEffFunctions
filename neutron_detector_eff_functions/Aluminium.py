import os
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import interpolate


def aluminium(substrate, wavelength, inclination):
    """calculates delta of aluminium

    Args:
        substrate (int):  thickness of aluminium layer in um
        wavelength (array[1]):  array of
    Returns:
        thValue: point where y in the integral have the closest value to threshold



    ..  Original source in Matlab: https: // bitbucket.org / europeanspallationsource / dg_matlabborontools / src / bcbac538ad10d074c5150a228847efc2e0269e0d / B10tools / rangesMAT.m?at = default & fileviewer = file - view - default

    """
    #g/cm3
    densityAl = 2.7
    denAl = (densityAl * 6.022e23 / 26.98) * 1e-4
    temp = read_cross_section(wavelength)
    sigmaal=[]
    for t in temp:
        sigmaal.append(denAl * t * 1e-24)
    delta = []
    for s in sigmaal:
        delta.append(np.exp(-1*(s * substrate) / math.sin(math.radians(inclination))))
    return delta


def read_cross_section(lambdalist):
    ht = 1.054e-34
    nmass = 1.67e-27
    # lambdaList transformed to ev
    lamen = []
    sigma = []
    c0 = 0
    for l in lambdalist:
        lamen.append(((ht ** 2) * 4 * math.pi ** 2) / (2 * nmass * (l[0] ** 2)) * 1e20 * 6.24e18)

    x, y = np.loadtxt(fname=os.path.dirname(os.path.abspath(__file__)) + "/data/Aluminium/AlCrossSect_(n,g).py",
                      delimiter=',', unpack=True)
    xlog = []
    ylog = []
    lamenlog = []
    for v in x:
        xlog.append(math.log10(v))
    for v in y:
        ylog.append(math.log10(v))
    for v in lamen:
        lamenlog.append(math.log10(v))
    f = interpolate.interp1d(xlog, ylog, kind='cubic')
    for v in lamenlog:
        sigma.append(10 ** f(v))
    # check if interpolation works
    # xnew = np.linspace(xlog[0], xlog[100], num=41, endpoint=True)
    # plt.plot(xlog, ylog, 'o', xnew, f(xnew), '--')
    print (lamen)
    print (sigma)
    plt.show()
    return sigma

