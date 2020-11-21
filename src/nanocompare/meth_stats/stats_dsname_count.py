import os

from lilab.tcga.global_tcga import pic_base_dir, logger, set_log_debug_level
from nanocompare.meth_stats.report_generation import load_ontcall_by_tool
from nanocompare.nanocompare_global_settings import tools_abbr
import pandas as pd
from pathlib import Path


def count_calls(ontcall):
    cnt = 0
    for key in ontcall.keys():
        cnt += len(ontcall[key])
    return cnt


def ontcall_stats_for_dslist():
    dslist = ['NA19240', 'K562', 'HL60', "APL"]

    toollist = tools_abbr

    df = pd.DataFrame()
    for dsname in dslist:
        for tool in toollist:
            ontcall = load_ontcall_by_tool(dsname=dsname, tool=tool)
            numcpgs = len(ontcall)
            numcalls = count_calls(ontcall)

            ret = {'dsname'   : dsname,
                    'tool'    : tool,
                    'numcpgs' : numcpgs,
                    'numcalls': numcalls
                    }

            logger.debug(f'ret={ret}')

            df = df.append(ret, ignore_index=True)

    outfn = os.path.join(pic_base_dir, 'ds4_results.xlsx')
    df.to_excel(outfn)


def summarize_basecalls():
    dsname_list = ['NA19240', 'K562', 'HL60', 'APL']

    for dsname in dsname_list:
        basedir = f"/fastscratch/liuya/nanocompare/{dsname}_basecalled"
        dflist = []
        for filename in Path(basedir).rglob('sequencing_summary.txt'):
            print(filename)

            df = pd.read_csv(filename, sep='\t', header=0)
            dflist.append(df)
            pass
            # logger.debug(f"df = {df}")

        dfall = pd.concat(dflist)
        logger.debug(f'dfall={dfall}')

        outfn = os.path.join(pic_base_dir, f'{dsname}_albacore_sequencing_summary.txt')
        dfall.to_csv(outfn, sep='\t', index=False)

    pass


if __name__ == '__main__':
    set_log_debug_level()

    # ontcall_stats_for_dslist()
    summarize_basecalls()

    pass
