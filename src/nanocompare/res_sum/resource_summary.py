import argparse
import glob
import logging
import os
import subprocess
from collections import defaultdict
from datetime import datetime

import pandas as pd
from tqdm import tqdm

from nanocompare.global_config import set_log_debug_level, pic_base_dir, logger

basedir = '/fastscratch/liuya/nanocompare'
fastdir = '/projects/li-lab/Nanopore_compare/result/running-logs'

tool_list_on_sumner = ['Tombo', 'DeepMod', 'DeepSignal', 'Nanopolish']
tool_list_on_winter = ['Guppy']

ntarget_dict = {'HL60': 50, 'K562': 50, 'APL': 50, 'NA19240': 300}  # , 'NA19240': 300

pkldir = '/projects/li-lab/yang/results/share_prj/result/running-logs'
sunmer_pkl = os.path.join(pkldir, 'sumner.task.resource.summary.pkl')
winter_pkl = os.path.join(pkldir, 'winter.task.resource.summary.pkl')
batch_fast5_pkl = os.path.join(pkldir, 'dsname.batch.fast5.summary.pkl')


def get_jobid_and_taskid(fn):
    """
    Sample input file name bascal.Guppy.K562.N50.batch49.22818.err
    :param fn:
    :return:
    """
    fn = os.path.basename(fn)
    last1_index = fn.rindex('.')
    last2_index = fn[:last1_index].rindex('.')
    taskid = fn[last2_index + 1: last1_index]

    last3_index = fn[:last2_index].rindex('.')
    batchstr = fn[last3_index + 1:last2_index].replace('batch', '')

    return taskid, batchstr


def winter_task_summary():
    dataset = defaultdict(list)
    for dsname in ntarget_dict:
        logging.info(dsname)
        ## Resquiggle collection
        #   HL60-Runs/HL60-N50-basecall/log
        basecall_logdir = os.path.join(basedir, f'{dsname}-Runs', f'{dsname}-N{ntarget_dict[dsname]}-basecall', 'log')
        pat_fns = os.path.join(basecall_logdir, 'bascal.*.out')
        fnlist = glob.glob(pat_fns)
        logging.info(f'Basecall collect: {len(fnlist)}')

        for fn in tqdm(fnlist):
            taskid, batchid = get_jobid_and_taskid(fn)
            command = f"""
                   seff {taskid}
                   """
            ret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")

            dataset['dsname'].append(dsname)
            dataset['batchid'].append(int(batchid))
            dataset['type'].append('basecall')
            dataset['job.results'].append(ret)

        ## Tool methcall collection
        for tool in tool_list_on_winter:
            methcall_logdir = os.path.join(basedir, f'{dsname}-Runs', f'{dsname}-{tool}-N{ntarget_dict[dsname]}', f'{dsname}-{tool}-N{ntarget_dict[dsname]}-methcall', 'log')
            pat_fns = os.path.join(methcall_logdir, '*.mcal.*.batch*.*.out')
            fnlist = glob.glob(pat_fns)
            logging.info(f'Methcall of {tool} collect: {len(fnlist)}')

            for fn in tqdm(fnlist):
                taskid, batchid = get_jobid_and_taskid(fn)
                # logging.debug(taskid)
                command = f"""
        seff {taskid}
        """
                ret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")
                dataset['dsname'].append(dsname)
                dataset['batchid'].append(int(batchid))
                dataset['type'].append(tool)
                dataset['job.results'].append(ret)
    df = pd.DataFrame.from_dict(dataset)
    outfn = os.path.join(pic_base_dir, 'winter.task.resource.summary.pkl')
    df.to_pickle(outfn)
    logging.info(f'save to {outfn}')

    outfn = os.path.join(pic_base_dir, 'winter.task.resource.summary.xlsx')
    df.to_excel(outfn)


def sunmer_task_summary():
    dataset = defaultdict(list)
    for dsname in ntarget_dict:
        logging.info(dsname)
        ## Resquiggle collection
        #   HL60-Runs/HL60-N50-basecall/log
        basecall_logdir = os.path.join(basedir, f'{dsname}-Runs', f'{dsname}-N{ntarget_dict[dsname]}-resquiggle', 'log')
        pat_fns = os.path.join(basecall_logdir, 'rsquigl.*.out')
        fnlist = glob.glob(pat_fns)
        logging.info(f'Resquiggle collect: {len(fnlist)}')

        for fn in tqdm(fnlist):
            taskid, batchid = get_jobid_and_taskid(fn)
            command = f"""
                seff {taskid}
                """
            ret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")

            dataset['dsname'].append(dsname)
            dataset['batchid'].append(int(batchid))
            dataset['type'].append('resquiggle')
            dataset['job.results'].append(ret)

        ## Tool methcall collection
        for tool in tool_list_on_sumner:
            # HL60-Runs/HL60-Nanopolish-N50/HL60-Nanopolish-N50-methcall/log
            methcall_logdir = os.path.join(basedir, f'{dsname}-Runs', f'{dsname}-{tool}-N{ntarget_dict[dsname]}', f'{dsname}-{tool}-N{ntarget_dict[dsname]}-methcall', 'log')
            # logging.debug(meth_logdir)

            pat_fns = os.path.join(methcall_logdir, '*.mcal.*.batch*.*.out')
            fnlist = glob.glob(pat_fns)
            logging.info(f'Methcall of {tool} collect: {len(fnlist)}')

            for fn in tqdm(fnlist):
                taskid, batchid = get_jobid_and_taskid(fn)
                # logging.debug(taskid)
                command = f"""
    seff {taskid}
    """
                ret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")
                dataset['dsname'].append(dsname)
                dataset['batchid'].append(int(batchid))
                dataset['type'].append(tool)
                dataset['job.results'].append(ret)

    df = pd.DataFrame.from_dict(dataset)
    # logging.info(df)
    outfn = os.path.join(pic_base_dir, 'sumner.task.resource.summary.pkl')
    df.to_pickle(outfn)
    logging.info(f'save to {outfn}')
    outfn = os.path.join(pic_base_dir, 'sumner.task.resource.summary.xlsx')
    df.to_excel(outfn)


def dataset_batch_summary():
    dataset = defaultdict(list)
    for dsname in ntarget_dict:
        logging.info(dsname)
        for batchid in range(1, ntarget_dict[dsname] + 1):
            septdir = os.path.join(basedir, f'{dsname}-Runs', f'{dsname}-N{ntarget_dict[dsname]}-sept', f'{batchid}')
            command = f"""
                           ls {septdir}/*.fast5 | wc -l
                            """
            ret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")
            logging.debug(f'batchid={batchid}, ret={ret}')
            dataset['dsname'].append(dsname)
            dataset['batchid'].append(batchid)
            dataset['fast5'].append(int(ret.strip()))

    df = pd.DataFrame.from_dict(dataset)
    outfn = os.path.join(pic_base_dir, 'dsname.batch.fast5.summary.pkl')
    df.to_pickle(outfn)
    logging.info(f'save to {outfn}')
    outfn = os.path.join(pic_base_dir, 'dsname.batch.fast5.summary.xlsx')
    df.to_excel(outfn)


def str_extract_time_mem(jobret):
    """
    Example

    Job ID: 6289162
    Cluster: slurm_cluster
    User/Group: root/jaxuser
    State: COMPLETED (exit code 0)
    Cores: 1
    CPU Utilized: 00:18:08
    CPU Efficiency: 17.31% of 01:44:46 core-walltime
    Job Wall-clock time: 01:44:46
    Memory Utilized: 71.49 MB
    Memory Efficiency: 0.05% of 150.00 GB

    :param jobret:
    :return:
    """
    cpu_use = None
    mem_use = None

    for line in jobret.splitlines():
        if line.strip().startswith('State: F') or line.strip().startswith('State: CANCELLED'):
            return None, None

        if line.strip().startswith('CPU Utilized:'):
            cpu_use = line.strip().replace('CPU Utilized:', '').strip()
        elif line.strip().startswith('Memory Utilized:'):
            mem_use = line.strip().replace('Memory Utilized:', '').strip()
            break

    return cpu_use, mem_use


def running_resouce_extraction(row):
    """
    :param row:
    :return:
    """
    jobret = row['job.results']
    runtime, mem = str_extract_time_mem(jobret)

    start_time = datetime.strptime("00:00:00", '%H:%M:%S')

    try:
        end_time = datetime.strptime(runtime, '%H:%M:%S')
    except:
        end_time = datetime.strptime(runtime, '%d-%H:%M:%S')

    duration_time = end_time - start_time

    # 71.49 MB
    if mem.endswith('MB'):
        mem_gb = float(mem.replace('MB', '').strip()) / 1000
    elif mem.endswith('GB'):
        mem_gb = float(mem.replace('GB', '').strip())
    else:
        raise Exception(f'Unrecognized mem={mem}')

    return pd.Series([runtime, mem, duration_time.total_seconds(), mem_gb], index=['running.time', 'mem.usage', 'running.time.seconds', 'mem.usage.gb'])


def unify_data_df():
    df1 = pd.read_pickle(winter_pkl)
    df2 = pd.read_pickle(sunmer_pkl)
    df3 = pd.read_pickle(batch_fast5_pkl)
    df = pd.concat([df1, df2])

    df = df.merge(df3, on=['dsname', 'batchid'], how='left')

    df[['running.time', 'mem.usage', 'running.time.seconds', 'mem.usage.gb']] = df.apply(running_resouce_extraction, axis=1)

    logger.info(df)

    logger.info(list(df.columns))

    run_report_columns = ['dsname', 'batchid', 'type', 'fast5', 'running.time', 'mem.usage', 'running.time.seconds', 'mem.usage.gb', 'job.results']
    outdf = df[run_report_columns]
    outfn = os.path.join(pic_base_dir, 'running.summary.table.xlsx')
    outdf.to_excel(outfn)
    logger.info(f'save to {outfn}')

    outdf = df[run_report_columns[:-1]]
    outfn = os.path.join(pic_base_dir, 'running.summary.table.csv')
    outdf.to_csv(outfn)
    logger.info(f'save to {outfn}')


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Resource-summary')
    parser.add_argument('--cpu-task', action='store_true')
    parser.add_argument('--gpu-task', action='store_true')
    parser.add_argument('--dataset-batch', action='store_true')
    parser.add_argument('--unify', action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    set_log_debug_level()
    args = parse_arguments()
    logger.debug(args)

    if args.cpu_task:
        sunmer_task_summary()
    if args.gpu_task:
        winter_task_summary()
    if args.dataset_batch:
        dataset_batch_summary()
    if args.unify:
        unify_data_df()
    logger.info("DONE")
