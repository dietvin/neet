#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name="neet",
      version="0.1.0",
      description="Nanopore error pattern exploration toolkit",
      author="Vincent Dietrich",
      author_email="",
      url="https://github.com/dietvin/neet",
      license="LICENSE",
      packages=find_packages(),
      entry_points={"console_scripts": ["neet = neet.neet:main",],},
      install_requires=["plotly", "tqdm", "intervaltree", "scipy", "pyfiglet"],
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',)