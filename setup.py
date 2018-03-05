from distutils.core import setup

setup(
    name='neutron_detector_eff_functions',
    packages=['neutron_detector_eff_functions','neutron_detector_eff_functions.data.B10','neutron_detector_eff_functions.data.B10.10B4C220','neutron_detector_eff_functions.data.B10.10B4C224'],
    version='v1.3.0',
    description='A library to calculate Neutron detector theoretical efficiency',
    author='acarmona',
    author_email='acarmona@opendeusto.es',
    url='https://github.com/alvcarmona/neutronDetectorEffFunctions',  # URL to the github repo
    download_url='https://github.com/alvcarmona/neutronDetectorEffFunctions/archive/1.3.0.tar.gz',
    keywords=['science'],
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
    ],
    classifiers=[],
    include_package_data=True
)
