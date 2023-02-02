#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : setup.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome
"""
Install:
    pip install build twine

Build package command:
    conda activate py39
    find . -name '*.egg-info' -type d -exec rm -rf {} \+ &&\
        rm -rf dist/*  &&\
        python -m build

    twine upload dist/*

Test package:
    conda activate py36
    pip install dist/nanome-jax-2.0.6.tar.gz
    pip show nanome-jax
    ls /pod/2/li-lab/yang/anaconda3/envs/py36/lib/python3.6/site-packages/nanome
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nanome-jax",
    version="2.0.11",
    author="Yang Liu",
    author_email="yang.liu@jax.org",
    description="NANOME (Nanopore methylation) pipeline developed by Li Lab at The Jackson Laboratory",
    license='MIT License',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LabShengLi/nanome",
    project_urls={
        'Bug Tracker': 'https://github.com/LabShengLi/nanome/issues'
    },
    packages=(
        setuptools.find_packages(where="src", exclude=("*.*.resource", "*.*.*.saniti_ecoli", "*.xgboost.sanity",))
    ),
    package_dir={'nanome': 'src/nanome'},
    scripts=[
        'src/nanome/nanocompare/plot_figure.py',
        'src/nanome/nanocompare/read_level_eval.py',
        'src/nanome/nanocompare/site_level_eval.py',
        'src/nanome/nanocompare/tss_eval.py',
        'src/nanome/nanocompare/pcc_region_eval.py',
        'src/nanome/nanocompare/join_preds_eval.py',
        'src/nanome/nanocompare/region_intersect.py',
        'src/nanome/nanocompare/newtool_parser.py',
        'src/nanome/nanocompare/computeRawReadsCoverage.py',
        'src/nanome/nanocompare/report/gen_html_report.py',
        'src/nanome/nanocompare/report/gen_txt_readme.py',
        'src/nanome/xgboost/xgboost_train.py',
        'src/nanome/xgboost/xgboost_predict.py',
        'src/nanome/xgboost/xgboost_prepdata.py',
        'src/nanome/xgboost/cs_agg_site.py',
        'src/nanome/xgboost/cs_eval_read.py',
        'src/nanome/xgboost/cs_eval_site.py',
        'src/nanome/xgboost/cs_predict.py',
        'src/nanome/xgboost/cs_train.py',
        'src/nanome/other/phasing/hp_split.py',
        'src/nanome/other/phasing/mega_parser.py',
        'src/nanome/other/phasing/methcall2bed.py',
        'src/nanome/other/phasing/nanomethphase.py',
        'utils/FilesSeparator.py',
        'utils/clean_old_basecall_in_fast5.py',
        'utils/extract_methylation_fast5_support_dir.py',
        'utils/combination_model_prediction.py',
        'utils/hm_cluster_predict.py',
        'utils/sum_chr_mod.py',
        'utils/tombo_extract_per_read_stats.py',
        'utils/validate_nanome_container.sh',
        'utils/unify_format_for_calls.sh',
        'utils/getGuppyVersion.py',
    ],
    include_package_data=True,
    package_data={'': ['src/nanome/common/*.csv', 'src/nanome/xgboost/trained_model/*.pkl']},
    classifiers=[
        "Programming Language :: Python :: 3",
        'Intended Audience :: Science/Research',
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        'biopython',
        'pybedtools >=0.8.2',
        # lower version, NA values have issues, ref: https://daler.github.io/pybedtools/changes.html#changes-in-v0-8-0
        'pandas',
        'seaborn',
        'scipy',
        'numpy',
        'statsmodels',
        'scikit-learn',
        # upper version may not load model success, ref: https://github.com/EpistasisLab/tpot/issues/1171
        'matplotlib',
        'jinja2',
        'openpyxl',
        'h5py',
        'tqdm',
        'joblib',
        'psutil',
        'xgboost',
        'pytabix',
        'pysam',
        'ont-fast5-api'
    ]
)
