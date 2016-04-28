# setup.py
from setuptools import setup, find_packages
from distutils.util import convert_path

main_ns = {}
ver_path = convert_path('mgfunc_v2/_version.py')
with open(ver_path) as ver_file:_
    exec(ver_file.read(), main_ns)

setup(name='mgfunc_v2',
      version=main_ns['__version__'],
      description='MGFunc is a functional and taxonomical annotation tool for large metagenomic data sets',
      url='https://github.com/MG-group-tools/MGFunc',
      author='Asli Ismihan Ozen',
      author_email='asli@cbs.dtu.dk',
      license='GPL 3.0',
      packages=['mgfunc_v2'],
      zip_safe=False)




