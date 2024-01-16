#!/usr/bin/env python

from setuptools import setup

setup(name="neet",
      version="0.0.1",
      description="Nanopore error pattern exploration toolkit",
      author="Vincent Dietrich",
      author_email="dietricv@uni-mainz.de",
      url="https://github.com/dietvin/neet",
      license="LICENSE",
      packages=["neet", "neet.summary_style"],
      entry_points={"console_scripts": ["neet = neet.__main__:main",],},
      python_requires=">3.9",
      install_requires=["plotly", "tqdm", "intervaltree", "scipy", "pyfiglet", "setuptools"],
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      package_data = {"neet.summary_style": ["plotly.js", "style.css"]},
      classifiers=["Environment :: Console",
                   "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                   "Programming Language :: Python :: 3",
                   "Topic :: Scientific/Engineering :: Bio-Informatics",
                   "Operating System :: OS Independent"])