import numpy
from Cython.Build import cythonize
from setuptools import setup, find_packages
from distutils.core import setup, Extension

include_path = [numpy.get_include()]
cythonize('TransPAC/cython/aligner_opts.pyx',include_path)

setup(
    name='TransPAC',
    description='',
    packages=find_packages(),
    include_package_data=True,
    entry_points = {
        'console_scripts': ['transpac = TransPAC.main:main'],
        },
    ext_modules = [Extension('aligner_opts', ['TransPAC/cython/aligner_opts.c'])],
    include_dirs=[numpy.get_include()],
    platforms='any'
)
