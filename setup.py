#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='clip_analysis',
    version='0.0.1',
    packages=['clip_analysis', 'clip_analysis_legacy'],
    package_dir={
        'clip_analysis':'clip_analysis'
    },
    include_package_data=True,
    url='',
    license='',
    author='byee4',
    author_email='',
    description='clip analysis plotting scripts',
    entry_points = {
        'console_scripts': [
            'clip_analysis = clip_analysis.plot_clip_analysis_figures:main',
        ]
    }
)
