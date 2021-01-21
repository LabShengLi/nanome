"""
Collect, pre-processing and save data to pkl

Collect and save all important experimental data and save them to pkl, etc. includes:

1. Performance data
2. Correlation data
3. Running time and resource usage data

etc.

"""
import glob
import os

import numpy as np
import pandas as pd

from nanocompare.global_config import logger, pic_base_dir
from nanocompare.global_settings import locations_category, locations_singleton


def collect_data_selected_locations(path, extension="tsv", prefix="APL_BSseq_cut10/APL_Bsseq_cut10.", sel_locations=["GW", "singletons", "nonsingletons", "cpgIslandExt", "promoters_500bp"]):
    """
    Collect performance raw results on specific path, there are tsv report files for performance results
    :param path:
    :param extension:
    :param prefix:
    :param sel_locations:
    :return:
    """
    # all required files should be placed in one directory
    # print(path)
    pattern = path + "/*." + extension
    # print(pattern)
    files = glob.glob(pattern)

    dfs = []

    for f in files:
        df = pd.read_csv(f, index_col=0, sep="\t")
        dfs.append(df)

    cobmined_data = pd.concat(dfs)

    tool_names = cobmined_data["prefix"].str.split(".", expand=True)

    #     tool_names = cobmined_data["prefix"].str.split(".",n=1,expand=True).str.split("_",n=1,expand=True)#,n=1,expand=True)
    #     print(tool_names)
    cobmined_data["Tool"] = tool_names[1]

    cobmined_data = cobmined_data.replace(to_replace="False", value="x.x.GW")
    #     cobmined_data=cobmined_data.replace(to_replace=prefix,value="x.")
    #     print(cobmined_data)
    cobmined_data = cobmined_data.replace(to_replace="hg38_singletons.bed", value="x.x.singletons")
    cobmined_data = cobmined_data.replace(to_replace="hg38_nonsingletons.bed", value="x.x.nonsingletons")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_100.bed", value="x.x.promoters_100bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_1000.bed", value="x.x.promoters_1000bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_200.bed", value="x.x.promoters_200bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_2000.bed", value="x.x.promoters_2000bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_500.bed", value="x.x.promoters_500bp")
    cobmined_data = cobmined_data.replace(to_replace="ONT.hg38.promoterFeature.flank_750.bed", value="x.x.promoters_750bp")

    #     if "hg38_singletons.mixed.bed

    #     location = cobmined_data["coord"].str.split(".",expand=True)
    location = cobmined_data["coord"].str.replace(prefix, "x.")
    location = location.str.split(".", expand=True)

    #     tool = cobmined_data["prefix"].str.split(".",n=4,expand=True)

    #     print(location)
    cobmined_data["Location"] = location[2]
    #     cobmined_data["BasecallTool"] = tool[1]

    ### Because I have toooooo many labels, I am removing some categories:
    #     cobmined_data = cobmined_data[[cobmined_data["Location"] == "GW" or cobmined_data["Location"] == "singletons"]]
    cobmined_data = cobmined_data[cobmined_data["Location"].isin(sel_locations)]

    # logger.debug(cobmined_data.info())
    logger.debug(cobmined_data)
    return cobmined_data


def statsFileParser(infileName):
    infile = open(infileName)

    switch = 0
    for row in infile:
        if "cput" in row and "mem" in row and "vmem" in row and "walltime" in row:
            tmp = row.replace(" ", "").replace("|", "").replace("=", "").replace("kb", "").replace("cput", "").replace("vmem", "").replace("mem", "").replace("walltime", "").strip().split(",")
            cput_tmp = tmp[0].split(":")
            #             print(cput_tmp, tmp)
            cput = int(cput_tmp[2]) + int(cput_tmp[1]) * 60 + int(cput_tmp[0]) * 3600
            mem = int(tmp[1])
            vmem = int(tmp[2])

            switch = 1

    infile.close()

    if switch == 1:
        return cput, mem, vmem
    else:
        print("#statsFileParser WARNING: no TORQUE stats detected in {} file".format(infileName))
        return 0, 0, 0


def parseFiles(path, filePrefix="*.o*"):
    os.chdir(path)
    files = []
    for file in glob.glob(filePrefix):
        files.append("{}/{}".format(path, file))
    #     print("#\tparseFiles: {} files parsed".format(len(files)))
    return files


def getStats(path, filePrefix):
    files = parseFiles(path, filePrefix)

    cput = mem = vmem = 0

    for file in files:
        a, b, c = statsFileParser(file)
        cput += a
        mem += b
        vmem += c
    #         print(a, b, c)
    return cput, mem, vmem


def convert2HumanReadable(cput, mem, vmem):
    """
    Converts mem and vmem to human readable format (i.e. KB to GB).
    Based on https://www.gbmb.org/kb-to-gb and https://blog.codinghorror.com/gigabyte-decimal-vs-binary/
    I am choosing to show results in binary format, but decimal option will be just commented in the next line.
    """

    # binary
    mem = mem / float(2 ** 20)
    vmem = vmem / float(2 ** 20)

    # decimal:
    #     mem = mem/float(10**6)
    #     vmem = vmem/float(10**6)

    """
    actually the above analyssi doesnt make sense. Its unfair comparison. Each tool needs to be run as one process,
    and for for example 10k reads. Repeated maybe a few times to report average time and resource consumption. 
    Moreover, this should preferably be run on a single cluster, with me alone, so that the running process would 
    not be disrupted by the external scripts.
    """

    #########################

    """
    time I will show in format of "XX days, XXh : XXm : XXs"
    """

    d = int(60 * 60 * 24)
    h = 60 * 60
    m = 60

    days = int(cput / d)
    cput = cput - days * d

    hours = int(cput / h)
    cput = cput - hours * h

    minutes = int(cput / m)
    seconds = cput - minutes * m

    output = "{} day(s), {}h : {}m : {}s".format(days, hours, minutes, seconds)
    print(output)

    return output


def convert2HumanReadable_average(cput, mem):
    """
    Converts mem and vmem to human readable format (i.e. KB to GB).
    Based on https://www.gbmb.org/kb-to-gb and https://blog.codinghorror.com/gigabyte-decimal-vs-binary/
    I am choosing to show results in binary format, but decimal option will be just commented in the next line.
    """

    # binary
    mem = mem / float(2 ** 20)
    #     vmem = vmem/float(2**20)

    # decimal:
    #     mem = mem/float(10**6)
    #     vmem = vmem/float(10**6)

    """
    actually the above analyssi doesnt make sense. Its unfair comparison. Each tool needs to be run as one process,
    and for for example 10k reads. Repeated maybe a few times to report average time and resource consumption. 
    Moreover, this should preferably be run on a single cluster, with me alone, so that the running process would 
    not be disrupted by the external scripts.
    """

    #########################

    """
    time I will show in format of "XX days, XXh : XXm : XXs"
    """

    cput = round(cput)

    d = int(60 * 60 * 24)
    h = 60 * 60
    m = 60

    days = int(cput / d)
    cput = cput - days * d

    hours = int(cput / h)
    cput = cput - hours * h

    minutes = int(cput / m)
    seconds = cput - minutes * m

    output = "{} day(s), {}h : {}m : {}s".format(days, hours, minutes, seconds)
    #     print(output)

    return output, mem


def save_nanocompare_running_pkl():
    """
    Save running time and resource usage data to DF
    :return:
    """
    finals_mem = []
    finals_time = []

    print("\n\n####### Tombo:")

    cputs = []
    mems = []

    filePrefix = "Tombo_automated.submit_01.sh.o*"

    for part in range(1, 101):
        cput, mem, vmem = getStats("/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/randomlySelected_configs/results/Tombo/Tombo_batch_{}".format(part), filePrefix)
        if cput > 0:
            cputs.append(cput)
            mems.append(mem)

    average_time, average_mem = convert2HumanReadable_average(np.mean(cputs), np.mean(mems))

    print("average cput:", average_time)
    print("average cput (s):", np.mean(cputs))
    print("average mem:", average_mem)
    print("average mem (kb):", np.mean(mems))
    print("points:", len(cputs))
    finals_mem.append(average_mem)
    finals_time.append(np.mean(cputs))

    print("\n\n####### Nanopolish:")
    # finals_mem.append(0)
    # finals_time.append(0)

    cputs = []
    mems = []

    filePrefix = "Nanopolish_automated.submit.sh.o*"

    for part in range(1, 101):
        cput, mem, vmem = getStats("/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/randomlySelected_configs/results/Nanopolish/Nanopolish_batch_{}".format(part), filePrefix)
        if cput > 0:
            cputs.append(cput)
            mems.append(mem)

    average_time, average_mem = convert2HumanReadable_average(np.mean(cputs), np.mean(mems))

    print("average cput:", average_time)
    print("average cput (s):", np.mean(cputs))
    print("average mem:", average_mem)
    print("average mem (kb):", np.mean(mems))
    print("points:", len(cputs))
    finals_mem.append(average_mem)
    finals_time.append(np.mean(cputs))

    print("\n\n####### DeepSignal:")

    cputs = []
    mems = []

    filePrefix = "DeepSignal_automated.submit_03.sh.o*"

    for part in range(1, 101):
        cput, mem, vmem = getStats("/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/randomlySelected_configs/results/DeepSignal/DeepSignal_batch_{}".format(part), filePrefix)
        if cput > 0:
            cputs.append(cput)
            mems.append(mem)

    average_time, average_mem = convert2HumanReadable_average(np.mean(cputs), np.mean(mems))

    print("average cput:", average_time)
    print("average cput (s):", np.mean(cputs))
    print("average mem:", average_mem)
    print("average mem (kb):", np.mean(mems))
    print("points:", len(cputs))
    finals_mem.append(average_mem)
    finals_time.append(np.mean(cputs))

    print("\n\n####### DeepMod:")

    cputs = []
    mems = []

    filePrefix = "DeepMod_automated.submit_01.sh.o*"

    for part in range(1, 101):
        cput, mem, vmem = getStats("/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/randomlySelected_configs/results/DeepMod/DeepMod_batch_{}".format(part), filePrefix)
        if cput > 0:
            cputs.append(cput)
            mems.append(mem)

    average_time, average_mem = convert2HumanReadable_average(np.mean(cputs), np.mean(mems))

    print("average cput:", average_time)
    print("average cput (s):", np.mean(cputs))
    print("average mem:", average_mem)
    print("average mem (kb):", np.mean(mems))
    print("points:", len(cputs))
    finals_mem.append(average_mem)
    finals_time.append(np.mean(cputs))

    # textOut = convert2HumanReadable(cput, mem, vmem)
    # textOut

    df = pd.DataFrame()
    df["mem"] = finals_mem
    df["time"] = finals_time
    df["tool"] = ["Tombo", "Nanopolish", "DeepSignal", "DeepMod"]

    print(df)

    outfn = os.path.join(pic_base_dir, "nanocompare_running_time.pkl")
    df.to_pickle(outfn)
    return df


def load_data(path, extension="tsv"):
    # all required files should be placed in one directory
    # print(path)
    pattern = path + "/*." + extension
    # print(pattern)
    files = glob.glob(pattern)

    dfs = []

    for f in files:
        df = pd.read_csv(f, index_col=0, sep="\t")
        dfs.append(df)

    cobmined_data = pd.concat(dfs)

    tool_names = cobmined_data["prefix"].str.split("_", n=1, expand=True)
    cobmined_data["Tool"] = tool_names[0]

    cobmined_data = cobmined_data.replace(to_replace="False", value="x.x.GW")

    cobmined_data = cobmined_data.replace(to_replace="hg38_singletons.bed", value="x.x.singletons")
    cobmined_data = cobmined_data.replace(to_replace="hg38_nonsingletons.bed", value="x.x.nonsingletons")

    # 	if "hg38_singletons.mixed.bed

    location = cobmined_data["coord"].str.split(".", n=4, expand=True)

    tool = cobmined_data["prefix"].str.split(".", n=4, expand=True)

    # 	print(location)
    cobmined_data["Location"] = location[2]
    cobmined_data["BasecallTool"] = tool[1]
    return cobmined_data


def select_locations_from_reportdf(df, sel_locations=locations_category + locations_singleton):
    """
    Select only interested locations
    :param df:
    :param sel_locations:
    :return:
    """
    return df[df['Location'].isin(sel_locations)]


def get_num_lines(fn):
    """
    Count how many lines for single sites start and end bed only
    :param fn:
    :return:
    """
    return sum(1 for line in open(fn))


def load_singleton_nonsingleton_sites():
    """
    Check all singleton and non-singleton CpG bed.sites file lines
    :return:
    """
    runPrefixList = ['K562_WGBS_joined_cut5', 'APL_Bsseq_cut5', 'HL60_AML_Bsseq_cut5', 'NA19240_RRBS_joined_cut5']

    basedir = '/projects/liuya/results/pkl/nanocompare'

    df = pd.DataFrame()
    for runPrefix in runPrefixList:
        rundir = os.path.join(basedir, runPrefix)
        singletonFn = os.path.join(rundir, f'{runPrefix}.hg38_singletons.absolute.bed')
        num_s = get_num_lines(singletonFn)
        concordantFn = os.path.join(rundir, f'{runPrefix}.hg38_nonsingletons.concordant.bed.sites')
        discordantFn = os.path.join(rundir, f'{runPrefix}.hg38_nonsingletons.discordant.bed.sites')
        num_con = get_num_lines(concordantFn)

        num_dis = get_num_lines(discordantFn)

        ret = {'dsname': runPrefix, 'Singletons': num_s, 'Non-singletons': num_con + num_dis, 'Concordant': num_con, 'Discordant': num_dis}
        df = df.append(ret, ignore_index=True)
    df = df[['dsname', 'Singletons', 'Non-singletons', 'Concordant', 'Discordant']]
    logger.info(f'df={df}')

    outfn = os.path.join(pic_base_dir, 'singleton.nonsingleton.sites.distribution.xlsx')
    df.to_excel(outfn, index=False)


def collect_singleton_vs_nonsingleton_df(runPrefixDict):
    dflist = []
    for runPrefix in runPrefixDict:
        pattern = os.path.join(runPrefixDict[runPrefix], "*.summary.singleton.nonsingleton.csv")
        flist = glob.glob(pattern)
        if len(flist) != 1:
            raise Exception(f"Too much/No summary of singleton vs non-singleton for {runPrefix} in folder {runPrefixDict[runPrefix]} with pattern={pattern}, len={len(flist)}")
        df = pd.read_csv(flist[0], index_col=0)
        dflist.append(df)
    retdf = pd.concat(dflist)
    retdf.index.name = 'Dataset'
    return retdf


def collect_performance_report_as_df(runPrefixDict):
    """
    create report from list of runPrefix, return specified columns
    :return:
    """
    dflist = []
    for runPrefix in runPrefixDict:
        logger.debug(f'runPrefix={runPrefix}')
        pattern = os.path.join(runPrefixDict[runPrefix], "performance?results", "*.performance.report.csv")
        files = glob.glob(pattern)
        for infn in files:
            df = pd.read_csv(infn, index_col=0, sep=",")
            dflist.append(df)
            logger.debug(f'Collect data from {runPrefix}:{os.path.basename(infn)} = {len(df)}')
    if len(dflist) == 0:
        raise Exception(f"Can not find report at {runPrefix} in folder {runPrefixDict[runPrefix]} with pattern={pattern}")
    combdf = pd.concat(dflist, ignore_index=True)
    logger.info(f'We collected total {len(dflist)} files with {len(combdf)} records')
    return combdf


def load_wide_format_performance_results(runPrefixDict, sel_locations=locations_category + locations_singleton):
    """
    Collect the currently new performance of exp results for paper
    :return:
    """
    df = collect_performance_report_as_df(runPrefixDict)
    seldf = select_locations_from_reportdf(df, sel_locations=sel_locations)

    logger.debug(f"collect_newly_exp_data, wide-format seldf={len(seldf)} using locations={sel_locations}")

    return seldf


def save_wide_format_performance_results(runPrefixDict, outdir, tagname):
    """
    Save all performance report results into a csv
    :return:
    """
    df = load_wide_format_performance_results(runPrefixDict)
    outfn = os.path.join(outdir, f'performance-results{f"-{tagname}" if tagname else ""}.csv')
    df.to_csv(outfn)
    logger.info(f'save to {outfn}')


if __name__ == '__main__':
    # df = load_singleton_nonsingleton_sites()
    df = load_wide_format_performance_results()
    logger.info(f'df={df}')

    # df = get_data_all()
    # save_alldata(df)
