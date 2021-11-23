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
    version="1.3.18",
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
        'src/nanocompare/pcc_region_eval.py',
        'src/nanocompare/nanome_consensus.py',
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
        'utils/unify_format_for_calls.sh',
        'src/nanocompare/xgboost/xgboost_train.py',
        'src/nanocompare/xgboost/xgboost_predict.py',
        'src/nanocompare/xgboost/xgboost_prepdata.py',
    ],
    include_package_data=True,
    package_data={'': ['src/nanocompare/*.csv', 'src/nanocompare/xgboost/trained_model/*.pkl']},
    classifiers=[
        "Programming Language :: Python :: 3",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        'biopython',
        'pybedtools >=0.8.2', # lower version, NA values have issues, ref: https://daler.github.io/pybedtools/changes.html#changes-in-v0-8-0
        'pandas',
        'seaborn',
        'scipy',
        'numpy',
        'statsmodels',
        'scikit-learn <=0.23.2', # upper version may not load model success, ref: https://github.com/EpistasisLab/tpot/issues/1171
        'matplotlib',
        'jinja2',
        'openpyxl',
        'h5py',
        'tqdm',
        'joblib',
        'psutil',
        'xgboost'
    ]
)
