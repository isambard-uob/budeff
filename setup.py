"""Setup script for BUDEFF."""

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize


def readme():
    """Loads the readme file for AMPAL."""
    with open('README.md', 'r') as inf:
        return inf.read()


setup(name='BUDEFF',
      version='1.0.0',
      description='A Python implementation of the BUDE force field.',
      long_description=readme(),
      long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
      url='https://github.com/isambard-uob/buff',
      author='Woolfson Group, University of Bristol',
      author_email='chris.wood@bristol.ac.uk',
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],
      packages=find_packages('src'),
      package_dir={'': 'src'},
      include_package_data=True,
      # This code automatically builds the Cython extensions.
      ext_modules=cythonize(
          [Extension(
              "budeff.calculate_energy",
              ["src/budeff/calculate_energy.pyx"],
              include_dir=["src/budeff/"],
              language='c++'),
           ]
      ),
      install_requires=[
          'Cython',
          'ampal'
      ],
      zip_safe=False,
      )
