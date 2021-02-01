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

## TODO: not use anymore
basedir = '/fastscratch/liuya/nanocompare1'

run_log_dir = '/projects/li-lab/Nanopore_compare/result/running-logs'
basedir = run_log_dir

tool_list_on_sumner = ['Tombo', 'DeepMod', 'DeepSignal', 'Nanopolish']
tool_list_on_winter = ['Guppy', 'Megalodon']

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


def winter_task_summary_na19240():
    dsname = 'NA19240'
    dataset = defaultdict(list)
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

        rt_secs, mem_gb, jobid, array_jobid = str_extract_time_mem(ret)

        if not rt_secs:  # Failed tasks are not consider
            continue

        dataset['dsname'].append(dsname)
        dataset['tool'].append('basecall')
        dataset['batchid'].append(int(batchid))
        dataset['running.time.seconds'].append(rt_secs)
        dataset['mem.usage.gb'].append(mem_gb)
        dataset['jobid'].append(jobid)
        dataset['array.jobid'].append(array_jobid)
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
            rt_secs, mem_gb, jobid, array_jobid = str_extract_time_mem(ret)

            if not rt_secs:  # Failed tasks are not consider
                continue

            dataset['dsname'].append(dsname)
            dataset['tool'].append(tool)
            dataset['batchid'].append(int(batchid))
            dataset['running.time.seconds'].append(rt_secs)
            dataset['mem.usage.gb'].append(mem_gb)
            dataset['jobid'].append(jobid)
            dataset['array.jobid'].append(array_jobid)
            dataset['job.results'].append(ret)

    df = pd.DataFrame.from_dict(dataset)
    outfn = os.path.join(pic_base_dir, 'na19240.winter.task.resource.summary.pkl')
    df.to_pickle(outfn)
    logging.info(f'save to {outfn}')

    outfn = os.path.join(pic_base_dir, 'na19240.winter.task.resource.summary.xlsx')
    df.to_excel(outfn)


def winter_megalodon_task_summary():
    dataset = defaultdict(list)
    for dsname in ntarget_dict:
        logging.info(dsname)
        methdir = os.path.join(basedir, f'{dsname}-Runs', f'{dsname}-Megalodon-N{ntarget_dict[dsname]}', f'{dsname}-Megalodon-N{ntarget_dict[dsname]}-methcall', 'log')

        pat_fns = os.path.join(methdir, '*.mcal.*.out')
        fnlist = glob.glob(pat_fns)
        logging.info(f'Megalodon of {dsname} collect: {len(fnlist)}')

        for fn in tqdm(fnlist):
            taskid, batchid = get_jobid_and_taskid(fn)
            command = f"""
                               seff {taskid}
                               """
            jobret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")

            dataset['dsname'].append(dsname)
            dataset['batchid'].append(int(batchid))
            dataset['tool'].append('Megalodon')
            dataset['job.results'].append(jobret)

            runtime_seconds, mem_gb, jobid, array_jobid = str_extract_time_mem(jobret)
            dataset['running.time.seconds'].append(runtime_seconds)
            dataset['mem.usage.gb'].append(mem_gb)
            dataset['jobid'].append(jobid)
            dataset['array.jobid'].append(array_jobid)
        df1 = pd.DataFrame.from_dict(dataset)
        df2 = pd.read_pickle(batch_fast5_pkl)
        df = df1.merge(df2, on=['dsname', 'batchid'], how='left')

        dfout = df[['dsname', 'tool', 'batchid', 'jobid', 'array.jobid', 'fast5', 'running.time.seconds', 'mem.usage.gb', 'job.results']]
        outfn = os.path.join(pic_base_dir, 'recalculate.running.summary.Megalodon.xlsx')
        dfout.to_excel(outfn)

        dfout = df[['dsname', 'tool', 'batchid', 'jobid', 'array.jobid', 'fast5', 'running.time.seconds', 'mem.usage.gb']]
        outfn = os.path.join(pic_base_dir, 'recalculate.running.summary.Megalodon.csv')
        dfout.to_csv(outfn)


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


def sunmer_task_summary_na19240():
    dsname = 'NA19240'
    dataset = defaultdict(list)
    resquiggle_dir = os.path.join(basedir, f'{dsname}-Runs', f'{dsname}-N{ntarget_dict[dsname]}-resquiggle', 'log')
    pat_fns = os.path.join(resquiggle_dir, 'rsquigl.*.out')
    logger.debug(f'pat_fns={pat_fns}')
    fnlist = glob.glob(pat_fns)
    logging.info(f'Resquiggle collect: {len(fnlist)}')

    for fn in tqdm(fnlist):
        taskid, batchid = get_jobid_and_taskid(fn)
        command = f"""
            seff {taskid}
            """
        ret = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")

        rt_secs, mem_gb, jobid, array_jobid = str_extract_time_mem(ret)

        if not rt_secs:  # Failed tasks are not consider
            continue

        dataset['dsname'].append(dsname)
        dataset['tool'].append('resquiggle')
        dataset['batchid'].append(int(batchid))
        dataset['running.time.seconds'].append(rt_secs)
        dataset['mem.usage.gb'].append(mem_gb)
        dataset['jobid'].append(jobid)
        dataset['array.jobid'].append(array_jobid)
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

            rt_secs, mem_gb, jobid, array_jobid = str_extract_time_mem(ret)

            if not rt_secs:  # Failed tasks are not consider
                continue

            dataset['dsname'].append(dsname)
            dataset['tool'].append(tool)
            dataset['batchid'].append(int(batchid))
            dataset['running.time.seconds'].append(rt_secs)
            dataset['mem.usage.gb'].append(mem_gb)
            dataset['jobid'].append(jobid)
            dataset['array.jobid'].append(array_jobid)
            dataset['job.results'].append(ret)

    df = pd.DataFrame.from_dict(dataset)
    # logging.info(df)
    outfn = os.path.join(pic_base_dir, 'na19240.sumner.task.resource.summary.pkl')
    df.to_pickle(outfn)
    logging.info(f'save to {outfn}')
    outfn = os.path.join(pic_base_dir, 'na19240.sumner.task.resource.summary.xlsx')
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

    Job ID: 23443
    Array Job ID: 23398_23
    Cluster: winter
    User/Group: liuya/jaxuser
    State: COMPLETED (exit code 0)
    Nodes: 1
    Cores per node: 16
    CPU Utilized: 00:24:01
    CPU Efficiency: 9.26% of 04:19:28 core-walltime
    Job Wall-clock time: 00:16:13
    Memory Utilized: 2.52 GB
    Memory Efficiency: 0.16% of 1.56 TB

    :param jobret:
    :return:
    """
    cpu_use_time = None
    mem_use = None
    jobid = None
    array_jobid = None

    for line in jobret.splitlines():
        if line.strip().startswith('State: F') or line.strip().startswith('State: CANCELLED'):
            return None, None, None, None

        if line.strip().startswith('Job ID:'):
            jobid = line.strip().replace('Job ID:', '').strip()

        if line.strip().startswith('Array Job ID:'):
            array_jobid = line.strip().replace('Array Job ID:', '').strip()

        if line.strip().startswith('CPU Utilized:'):
            cpu_use_time = line.strip().replace('CPU Utilized:', '').strip()
        elif line.strip().startswith('Memory Utilized:'):
            mem_use = line.strip().replace('Memory Utilized:', '').strip()
            break
    start_time = datetime.strptime("00:00:00", '%H:%M:%S')

    try:
        end_time = datetime.strptime(cpu_use_time, '%H:%M:%S')
    except:
        end_time = datetime.strptime(cpu_use_time, '%d-%H:%M:%S')

    duration_time = end_time - start_time

    # 71.49 MB
    if mem_use.endswith('MB'):
        mem_gb = float(mem_use.replace('MB', '').strip()) / 1000
    elif mem_use.endswith('GB'):
        mem_gb = float(mem_use.replace('GB', '').strip())
    else:
        raise Exception(f'Unrecognized mem={mem_use} from jobret={jobret}')

    return duration_time.total_seconds(), mem_gb, jobid, array_jobid


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


def recalculate(fnlist=['na19240.sumner.task.resource.summary.pkl', 'na19240.winter.task.resource.summary.pkl']):
    dflist = []
    for fn in fnlist:
        dflist.append(pd.read_pickle(os.path.join(pkldir, fn)))
    df = pd.concat(dflist)

    df3 = pd.read_pickle(batch_fast5_pkl)
    df = df.merge(df3, on=['dsname', 'batchid'], how='left')
    logger.debug(df)

    dataset = defaultdict(list)
    for index, row in df.iterrows():
        if row['tool'] not in tool_list_on_sumner + tool_list_on_winter:
            continue
        dsname = row['dsname']
        batchid = row['batchid']
        runt = row['running.time.seconds']
        memg = row['mem.usage.gb']

        basecall_row = df[(df['dsname'] == dsname) & (df['batchid'] == batchid) & (df['tool'] == 'basecall')].iloc[0, :]
        resquiggle_row = df[(df['dsname'] == dsname) & (df['batchid'] == batchid) & (df['tool'] == 'resquiggle')].iloc[0, :]

        if row['tool'] in ['DeepSignal', 'Tombo']:
            runt += basecall_row['running.time.seconds'] + resquiggle_row['running.time.seconds']
            memg += basecall_row['mem.usage.gb'] + resquiggle_row['mem.usage.gb']
        elif row['tool'] in ['Nanopolish', 'DeepMod']:
            runt += basecall_row['running.time.seconds']
            memg += basecall_row['mem.usage.gb']
        dataset['dsname'].append(dsname)
        dataset['tool'].append(row['tool'])
        dataset['batchid'].append(row['batchid'])
        dataset['fast5'].append(row['fast5'])
        dataset['running.time.seconds'].append(runt)
        dataset['mem.usage.gb'].append(memg)
    outdf = pd.DataFrame.from_dict(dataset)
    logger.info(outdf)
    outfn = os.path.join(pic_base_dir, 'recalculate.running.summary.na19240.csv')
    outdf.to_csv(outfn)


def parse_arguments():
    """
    :return:
    """
    parser = argparse.ArgumentParser(description='Resource-summary')
    parser.add_argument('--cpu-task', action='store_true')
    parser.add_argument('--gpu-task', action='store_true')
    parser.add_argument('--dataset-batch', action='store_true')
    parser.add_argument('--unify', action='store_true')
    parser.add_argument('--recalculate', action='store_true')
    parser.add_argument('--megalodon', action='store_true')
    parser.add_argument('--na19240-winter', action='store_true')
    parser.add_argument('--na19240-sumner', action='store_true')

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
    if args.recalculate:
        recalculate()
    if args.megalodon:
        winter_megalodon_task_summary()
    if args.na19240_sumner:
        sunmer_task_summary_na19240()
    if args.na19240_winter:
        winter_task_summary_na19240()

    logger.info("DONE")
