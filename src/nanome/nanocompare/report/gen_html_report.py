#!/usr/bin/env python3
# @Author   : Yang Liu
# @FileName : gen_html_report.py
# @Software : NANOME project
# @Organization : JAX Li Lab
# @Website  : https://github.com/LabShengLi/nanome

import glob
import os
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
from jinja2 import Environment, FileSystemLoader


def density_plot(df, outfn, tool, col_index=6, bins=None, the_range=None, x_label=None):
    # setting font size
    plt.rcParams.update({'font.size': 6})

    # plots
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.hist(df.iloc[:, col_index], bins=bins, range=the_range, color='blue', alpha=0.4)

    # Add labels
    if the_range:
        plt.xlim(the_range)
    plt.title(f'{tool} = {len(df):,} CpGs')
    plt.xlabel(x_label)
    plt.ylabel('# CpGs')

    # save figures
    os.makedirs(os.path.dirname(outfn), exist_ok=True)
    plt.savefig(outfn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"save to {outfn}")


def get_methcall_report_df(baseDir, outDir):
    """
    Generate the methylation calling report dataframe
    :param baseDir:
    :return:
    """
    ## Pre-check the max-coverage
    max_cov = None
    for tool in tool_list:
        fnlist = glob.glob(os.path.join(baseDir, f'*_{tool}-perSite-cov1.sort.bed.gz'))
        if len(fnlist) < 1:
            print(f"Not found file in baseDir={baseDir}, pattern={f'*_{tool}-perSite-cov1.sort.bed.gz'}")
            continue
        try:
            df_site_level = pd.read_csv(fnlist[0], sep='\t', header=None, index_col=False)
            the_max_cov = max(df_site_level.iloc[:, 7])
            if max_cov is None or max_cov < the_max_cov:
                max_cov = the_max_cov
        except:
            pass
    print(f"max_cov={max_cov}")

    ret_dict = defaultdict(list)
    for tool in tool_list:
        fnlist = glob.glob(os.path.join(baseDir, f'*_{tool}-perSite-cov1.sort.bed.gz'))
        if len(fnlist) < 1:
            print(f"Not found file in baseDir={baseDir}, pattern={f'*_{tool}-perSite-cov1.sort.bed.gz'}")
            continue
        try:
            df_site_level = pd.read_csv(fnlist[0], sep='\t', header=None, index_col=False)
            ret_dict['Tool'].append(tool)
            ret_dict['CpGs'].append(f"{len(df_site_level):,}")

            outfn = os.path.join(outDir, 'images', f'dist_{tool}.png')
            density_plot(df_site_level, outfn, tool, bins=10, the_range=(0, 1), x_label='Methylation %')

            outfn = os.path.join(outDir, 'images', f'cov_{tool}.png')
            density_plot(df_site_level, outfn, tool, bins=10,
                         col_index=7, x_label='Coverage', the_range=(0, max_cov + 1))
        except:  # can not call any results
            ret_dict['Tool'].append(tool)
            ret_dict['CpGs'].append('0')
    return pd.DataFrame.from_dict(ret_dict)


if __name__ == '__main__':
    dsname = sys.argv[1]
    infn_running_info = sys.argv[2]
    infn_basecall_info = sys.argv[3]
    indir_methcall = sys.argv[4]
    outdir = sys.argv[5]
    basedir = sys.argv[6]
    version_file = sys.argv[7]

    ## Running information summary
    df_running_info = pd.read_csv(infn_running_info, sep='\t', )
    # print(df_running_info.to_html())

    ## Basecalling report summary
    dataset_basecall = defaultdict(list)
    if infn_basecall_info != 'None':
        with open(infn_basecall_info) as stats:
            for line in stats:
                linesplit = line.strip().split(':')
                linesplit = [e.strip() for e in linesplit]
                if line.startswith('Top 5'):
                    break
                if len(linesplit) > 1:
                    dataset_basecall['Title'].append(linesplit[0])
                    dataset_basecall['Information'].append(linesplit[1])
                else:
                    dataset_basecall['Title'].append(linesplit[0])
                    dataset_basecall['Information'].append("")
        df_basecall_info = pd.DataFrame.from_dict(dataset_basecall)
    else:
        df_basecall_info = pd.DataFrame()

    print(df_basecall_info)

    df_version = pd.read_csv(version_file, index_col=None, sep='\t', header=0)
    tool_list = list(df_version['Tool'])
    print(tool_list)
    df_methcall_info = get_methcall_report_df(indir_methcall, outdir).dropna()

    if 'Tool' in df_methcall_info:
        df_methcall_info = df_methcall_info.merge(df_version, on='Tool', how='left')
        df_methcall_info = df_methcall_info.fillna('1.0').iloc[:, [0, 2, 1]]
        df_fig_inf = df_methcall_info[df_methcall_info['CpGs'] != '0']
    else:
        df_methcall_info = pd.DataFrame()
        df_fig_inf = pd.DataFrame()

    print(f"df_methcall_info={df_methcall_info}")

    env = Environment(loader=FileSystemLoader(basedir))
    template = env.get_template('index.html')

    filename = os.path.join(outdir, 'nanome_report.html')
    with open(filename, 'w') as fh:
        fh.write(template.render(
            dsname=dsname,
            df_running_info=df_running_info,
            df_basecall_info=df_basecall_info,
            df_methcall_info=df_methcall_info,
            df_fig_inf=df_fig_inf
        ))
    print("### gen_html_report DONE")
