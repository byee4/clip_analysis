#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='clip_analysis',
    version='0.0.2',
    packages=['clip_analysis', 'clip_analysis_legacy', 'rnaseq_analysis'],
    package_dir={
        'clip_analysis':'clip_analysis',
        'rnaseq_analysis':'rnaseq_analysis',
        'integrated':'integrated'
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
            'plot_repetitive_elements_sunburst = clip_analysis.plot_repetitive_elements_sunburst:main',
            'plot_repetitive_elements_bar = clip_analysis.plot_repetitive_elements_bar:main',
            'plot_kmer_enrichment = clip_analysis.plot_kmer_enrichment:main',
            'plot_region_distribution = clip_analysis.plot_region_distribution:main',
            'plot_histogram_enriched_regions = clip_analysis.plot_histogram_enriched_regions:main',
            'plot_ip_foldchange_over_input_reads = clip_analysis.plot_ip_foldchange_over_input_reads:main',
            'plot_compare_rbp_enriched_regions = clip_analysis.plot_compare_rbp_enriched_regions:main',
        ]
    }
)
