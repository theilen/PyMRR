# -*- coding: utf-8 -*-

from distutils.core import setup

desc = 'Packagae for analyzing phase data obtained using magnetic resonance rheology.'

setup(name='mrr',
      version='1.4.2',
      description='Analyze MRR-data',
      long_description=desc,
      author='Sebastian Theilenberg',
      author_email='theilenberg@hiskp.uni-bonn.de',
      url="\\jarvis\Experimente\Kopfklappe\Eigene Programme\Python\mrr-package",
      packages=['mrr',
                'mrr.arc',
                'mrr.bvalue',
                'mrr.mrrcore',
                'mrr.plotting',
                'mrr.read',
                'mrr.timeline',
                'mrr.unwrapping',
                'mrr.unwrapping.py_unwrap',
                'mrr.unwrapping.py_unwrap.algorithms',
                'mrr.unwrapping.py_unwrap.array_manipulation',
                'mrr.unwrapping.py_unwrap.qualitymaps',
                'mrr.waveform',
                'mrr.write']
      )
