#!/usr/bin/env python

import json
import neutron_detector_eff_functions.Detector as Detector
import matplotlib.pyplot as plt

def calculate_eff_multiblade(nb,converterThickness, substrateThickness, wavelength, angle, threshold, single, converter):
    '''
    :param nb(int): Number of blades
    :param converterThickness(float): Thickness of converter in microns
    :param substrateThickness(float) : Thickness of substrate in microns
    :param wavelength(list): List with tuples [weigth in %, wavelength in angstoms]
    :param angle(int): Inclination of neutrom beam in degrees
    :param threshold(int): in KeV
    :param single (boolean): True if is a detector with a single coated blade
    :param converter(string):  Converter material, see list of available converters
    :return: list of values, first list contains the efficiency for each blade in depth order (the last is the deepest), last value si the total efficiency of the detector
    '''
    assert nb >= 1
    assert len(wavelength) >= 1
    detector = Detector.Detector.build_detector(nb, converterThickness, substrateThickness, wavelength, angle, threshold, single, converter)
    return detector.calculate_eff()

def calculate_eff_json(path):
    '''

    :param path: path to a json file containing a configuration
    :return:
    '''
    print ('Script eff json')
    detector = Detector.Detector.json_parser(path)
    return detector.calculate_eff()

#TODO fix this function (crashes after enable polichromatic calculation)
def plot_eff_vs_thick(path):
    detector = Detector.Detector.json_parser(path)
    detector.plot_thick_vs_eff2()
    plt.show()

def plot_eff_vs_wave(path):
    detector = Detector.Detector.json_parser(path)
    detector.plot_eff_vs_wave()
    plt.show()

def optimize_config_same_thick(originPath, destinyPath):
    detector = Detector.Detector.json_parser(originPath)
    detector.optimize_thickness_same()
    try:
        filepath = destinyPath
        with open(str(filepath) + '/' + detector.name + '_optimized_config_same.json', "w") as outfile:
            outfile.write(json.dumps(detector.to_json(), sort_keys=True, indent=4, ensure_ascii=False))
            outfile.close()
        print('Export')
    except IOError:
        print ("Path error")

def optimize_config_diff_thick(originPath, destinyPath):
    detector = Detector.Detector.json_parser(originPath)
    detector.optimize_thickness_diff_mono()
    detector.to_json()
    try:
        filepath = destinyPath
        with open(str(filepath) + '/' + detector.name + '_optimized_config_diff.json', "w") as outfile:
            outfile.write(json.dumps(detector.to_json(), sort_keys=True, indent=4, ensure_ascii=False))
            outfile.close()
        print('Export')
    except IOError:
        print ("Path error")



if __name__ == '__main__':
  #print (calculate_eff_multiblade(10,1,0,[[1.8,100]], 90, 100,False, '10B4C 2.24g/cm3'))
  #print (calculate_eff_json('/Users/alvarocbasanez/PycharmProjects/dg_efficiencycalculator/efficiencyCalculator/exports/detector1.json'))
  #plot_eff_vs_thick('/Users/alvarocbasanez/workspace/dg_efficiencyCalculator/efficiencyCalculator/exports/single.json')
  #optimize_config_same_thick('/Users/alvarocbasanez/PycharmProjects/Git dg_efficiencyCalculator/efficiencyCalculator/exports/detector1.json', '/Users/alvarocbasanez/PycharmProjects/Git dg_efficiencyCalculator/efficiencyCalculator/exports')
  #optimize_config_diff_thick('/Users/alvarocbasanez/PycharmProjects/Git dg_efficiencyCalculator/efficiencyCalculator/exports/detector1.json','/Users/alvarocbasanez/PycharmProjects/Git dg_efficiencyCalculator/efficiencyCalculator/exports')
  #plot_eff_vs_wave('/Users/alvarocbasanez/workspace/dg_efficiencyCalculator/efficiencyCalculator/exports/detector1.json')

  plot_eff_vs_thick('/Users/alvarocbasanez/workspace/dg_efficiencyCalculator/efficiencyCalculator/exports/detector1.json')
  #plot_eff_vs_wave('/Users/alvarocbasanez/workspace/dg_efficiencyCalculator/efficiencyCalculator/exports/detector1.json')