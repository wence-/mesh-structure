# -*- coding: utf-8 -*-
#
import os

from setuptools import find_packages, setup

# https://packaging.python.org/single_source_version/
base_dir = os.path.abspath(os.path.dirname(__file__))


setup(name="meshstructure",
      packages=find_packages(),
      install_requires=["fenics-fiat",
                        "fenics-ufl",
                        "pymbolic",
                        "islpy",
                        "numpy"],
      description="Mesh structure for code generation",
      classifiers=[
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Programming Language :: Python :: 3",
          "Topic :: Scientific/Engineering",
          "Topic :: Scientific/Engineering :: Mathematics",
      ],
      python_requires=">=3.5")
