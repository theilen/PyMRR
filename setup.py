# -*- coding: utf-8 -*-

from distutils.core import setup

from mrr import __version__

setup(name='mrr',
      version=__version__,
      description="Pure python package to analyze MRR-phase-data",
      long_description="""
      mrr is a pure python package for analyzing phasedata obtained with
      Magnetic Resonance Rheology. It provides simple data structures and
      functions to work with the phase data, vizualize it and save the data.
      """,
      author='Sebastian Theilenberg',
      author_email='theilenberg@hiskp.uni-bonn.de',
      url="https://github.com/theilen/PyMRR.git",
      packages=["mrr",
                "mrr.arc",
                "mrr.bvalue",
                "mrr.coordinates",
                "mrr.mrrcore",
                "mrr.plotting",
                "mrr.read",
                "mrr.strain",
                "mrr.timeline",
                "mrr.unwrapping",
                "mrr.unwrapping.py_unwrap",
                "mrr.unwrapping.py_unwrap.algorithms",
                "mrr.unwrapping.py_unwrap.array_manipulation",
                "mrr.unwrapping.py_unwrap.qualitymaps",
                "mrr.waveform",
                "mrr.waveform.tektronix",
                "mrr.write"
                ]
      )
