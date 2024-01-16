#!/usr/bin/env python

from setuptools import setup
from neet.__version__ import __version__

setup(name="neet",
      version=__version__,
      description="Nanopore error pattern exploration toolkit",
      author="Vincent Dietrich",
      author_email="",
      url="https://github.com/dietvin/neet",
      license="LICENSE",
      packages=["neet", "neet.summary_style"],
      entry_points={"console_scripts": ["neet = neet.__main__:main",],},
      python_requires=">3.9",
      install_requires=["plotly", "tqdm", "intervaltree", "scipy", "pyfiglet", "setuptools"],
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      package_data = {"neet.summary_style": ["plotly.js", "style.css"]},)