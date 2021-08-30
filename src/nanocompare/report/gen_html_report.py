import glob
import os
import sys
from collections import defaultdict

import pandas as pd
from jinja2 import Environment, FileSystemLoader

from nanocompare.global_settings import ToolNameList


def get_methcall_report_df(baseDir):
    """
    Generate the methylation calling report dataframe
    :param baseDir:
    :return:
    """
    ret_dict = defaultdict(list)
    for tool in ToolNameList:
        fnlist = glob.glob(os.path.join(baseDir, f'*_{tool}-perSite-cov1.sort.bed.gz'))
        if len(fnlist) < 1:
            raise Exception(f"Not found file in baseDir={baseDir}, pattern={f'*_{tool}-perSite-cov1.sort.bed.gz'}")
        df_site_level = pd.read_csv(fnlist[0], sep='\t', header=None, index_col=False)
        ret_dict['Tool'].append(tool)
        ret_dict['CpGs'].append(f"{len(df_site_level):,}")
    return pd.DataFrame.from_dict(ret_dict)


if __name__ == '__main__':
    dsname = sys.argv[1]
    infn_running_info = sys.argv[2]
    infn_basecall_info = sys.argv[3]
    indir_methcall = sys.argv[4]
    outdir = sys.argv[5]

    ## Running information summary
    df_running_info = pd.read_csv(infn_running_info, sep='\t', )
    # print(df_running_info.to_html())

    ## Basecalling report summary
    dataset_basecall = defaultdict(list)
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
    print(df_basecall_info)

    df_methcall_info = get_methcall_report_df(indir_methcall)
    print(df_methcall_info)

    env = Environment(loader=FileSystemLoader('src/nanocompare/report'))
    template = env.get_template('index.html')

    filename = os.path.join(outdir, 'nanome_report.html')
    with open(filename, 'w') as fh:
        fh.write(template.render(
            dsname=dsname,
            df_running_info=df_running_info,
            df_basecall_info=df_basecall_info,
            df_methcall_info=df_methcall_info,
        ))

    # df_basecall_info = pd.read_csv(infn_basecall_info, sep='\t')
    # print(df_basecall_info)

    pass
