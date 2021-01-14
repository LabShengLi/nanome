import glob
import logging
import os
import subprocess

import pandas as pd

from nanocompare.global_config import set_log_debug_level

basedir = '/fastscratch/liuya/nanocompare'
tool_list = ['Tombo', 'DeepMod', 'DeepSignal', 'Nanopolish']
datasets = [{'dsname': 'HL60', 'ntarget': 50}]

if __name__ == '__main__':
    set_log_debug_level()
    logging.info("hello")
    retdata = []
    for dsdict in datasets:
        for tool in tool_list:
            # /fastscratch/liuya/nanocompare/HL60-Runs/HL60-Nanopolish-N50/HL60-Nanopolish-N50-methcall/log
            logdir = os.path.join(basedir, f'{dsdict["dsname"]}-Runs', f'{dsdict["dsname"]}-{tool}-N{dsdict["ntarget"]}', f'{dsdict["dsname"]}-{tool}-N{dsdict["ntarget"]}-methcall', 'log')
            logging.debug(logdir)

            pat_fns = os.path.join(logdir, '*.mcal.*.batch*.*.out')
            fnlist = glob.glob(pat_fns)
            logging.debug(len(fnlist))

            for fn in fnlist:
                last1_index = fn.rindex('.')
                last2_index = fn[:last1_index].rindex('.')
                taskid = fn[last2_index + 1: last1_index]
                logging.debug(taskid)
                command = f"""
seff {taskid}
"""
                ret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")

                logging.debug(ret)

                data = {'dsname': dsdict["dsname"], 'ntarget': dsdict["ntarget"], 'tool': tool, 'seff-ret': ret}
                retdata.append(data)
            break
        break

    df = pd.DataFrame(retdata)
    logging.info(df)
