#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : setup.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/TheJacksonLaboratory/nanome
"""
Build package command:
    find . -name '*.egg-info' -type d | parallel -j1 -v rm -r {}
    rm -rf dist/*  && python -m build
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nanome-jax",
    version="1.3.8",
    author="Yang Liu",
    author_email="yang.liu@jax.org",
    description="NANOME (Nanopore methylation) pipeline developed by Li Lab at The Jackson Laboratory",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TheJacksonLaboratory/nanome",
    project_urls={
        'Bug Tracker': 'https://github.com/TheJacksonLaboratory/nanome/issues'
    },
    packages=(
        setuptools.find_packages(where="src", exclude=("*.resource", "*.*.saniti_ecoli",))
    ),
    package_dir={'nanocompare': 'src/nanocompare'},
    scripts=[
        'src/plot_figure.py',
        'src/nanocompare/read_level_eval.py',
        'src/nanocompare/site_level_eval.py',
        'src/nanocompare/tss_eval.py',
        'src/nanocompare/computeRawReadsCoverage.py',
        'src/nanocompare/report/gen_html_report.py',
        'utils/FilesSeparator.py',
        'utils/clean_old_basecall_in_fast5.py',
        'utils/extract_methylation_fast5_support_dir.py',
        'utils/combination_model_prediction.py',
        'utils/gen_readme.py',
        'utils/hm_cluster_predict.py',
        'utils/sum_chr_mod.py',
        'utils/tombo_extract_per_read_stats.py',
        'utils/validate_nanome_container.sh',
    ],
    include_package_data=True,
    package_data={'': ['src/nanocompare/*.csv']},
    classifiers=[
        "Programming Language :: Python :: 3",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'biopython',
        'pybedtools',
        'pandas',
        'seaborn',
        'scipy',
        'numpy',
        'statsmodels',
        'scikit-learn',
        'matplotlib',
        'jinja2',
        'openpyxl',
        'h5py',
        'tqdm',
        'joblib',
        'psutil'
    ]
)
