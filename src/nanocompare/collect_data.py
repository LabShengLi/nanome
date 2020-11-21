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

from global_config import logger, pic_base_dir, current_time_str, data_base_dir
from nanocompare.legacy import performance_plots as pp
# from nanocompare.load_data import load_all_perf_data
from nanocompare.nanocompare_global_settings import nanocompare_basedir, locations_category2, locations_singleton2

ret_report_columns = ['Dataset', 'Tool', 'Location', 'Accuracy', 'ROC_AUC',
        'F1_5mC', 'F1_5C', 'Precision_5mC', 'Precision_5C', 'Recall_5mC', 'Recall_5C',
        'Corr_Mix', 'Corr_All', 'Corr_mixedSupport', 'Corr_allSupport',
        'mCsites_called', 'Csites_called', '5mCs', '5Cs', 'TotalSites', 'prefix', 'coord'
        ]


def collect_box_plot_all_data():
    """
    Collect all box plot data from original exp report tsv files from folder 'reports_dump'
    :return:
    """
    folder = "reports_dump"
    df = pp.load_data(f"{nanocompare_basedir}/reports/{folder}")
    logger.info(f"df={df}")
    return df


def save_box_plot_all_data(df=None):
    """
    Save box plot all data into files pkl and xlsx
    :param df:
    :return:
    """
    if df is None:
        df = collect_box_plot_all_data()

    outfn = os.path.join(pic_base_dir, 'box_plots_all_data.joinedReports.pkl')
    df.to_pickle(outfn)
    logger.info(f"save to {outfn}")

    outfn = os.path.join(pic_base_dir, 'box_plots_all_data.joinedReports.xlsx')
    df.to_excel(outfn)
    logger.info(f"save to {outfn}")
    return df


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


def get_data_all():
    """
    Get all performance row data to DF
    :return:
    """
    listOfSelected = ["GW", "cpgIslandExt", "promoters_500bp", "exonFeature", "intergenic", "intronFeature"]
    listOfSelected = listOfSelected + ["singletons", "nonsingletons", "discordant", "concordant"]

    # ["K562_WGBS_joined", "APL_BSseq_cut10", "APL_Bsseq", "K562_RRBS_joined", "HL60_RRBS_joined", "APL_oxBSseq_cut5"]

    dsname_list = ["K562_WGBS_joined", "K562_RRBS_joined", "APL_BSseq_cut10", "APL_Bsseq", "APL_oxBSseq_cut5", "HL60_RRBS_joined", 'HL60_AML_Bsseq_cut5', 'HL60_AML_oxBsseq_cut5']

    dsname_list = ["K562_WGBS_joined", "APL_BSseq_cut10", 'HL60_AML_Bsseq_cut5', 'NA19240_RRBS_joined']

    pdlist = []
    for prefix in dsname_list:
        data = collect_data_selected_locations(f"{nanocompare_basedir}/reports/{prefix}", sel_locations=listOfSelected)
        data['Dataset'] = prefix
        pdlist.append(data)
    alldata = pd.concat(pdlist)
    return alldata


def add_dataset_report_to_all(indir='/projects/liuya/results/pkl/nanocompare/NA19240_RRBS_joined', only_test=True):
    """
    add a new results to original datasets, used for NA19240 new coming results
    :param indir:
    :param only_test: True if not save, False will save results to pkl
    :return:
    """
    listOfSelected = ["GW", "cpgIslandExt", "promoters_500bp", "exonFeature", "intergenic"]
    listOfSelected = listOfSelected + ["singletons", "nonsingletons", "discordant", "concordant"]

    dstagname = os.path.basename(indir)

    dsreport = collect_data_selected_locations(indir, sel_locations=listOfSelected)
    dsreport['Dataset'] = dstagname
    logger.info(f"dsreport={dsreport}")

    if not only_test:
        outfn = os.path.join(pic_base_dir, f"{dstagname}_report.xlsx")
        dsreport.to_excel(outfn)

    dsall = load_all_perf_data()
    # logger.debug(f"dsall={dsall}")

    d2 = pd.concat([dsall, dsreport], ignore_index=True)
    logger.info(f"d2={d2}")

    if not only_test:
        save_alldata(d2)


def save_alldata(dfall=None):
    """
    Save all performance row data to DF
    :param dfall:
    :return:
    """
    if dfall is None:
        dfall = get_data_all()

    logger.info(f"dfall={dfall}")

    outfn = os.path.join(pic_base_dir, "Nanocompare_performance_alldata.xlsx")
    dfall.to_excel(outfn)
    logger.info(f"save to {outfn}")

    outfn = os.path.join(pic_base_dir, "Nanocompare_performance_alldata.pkl")
    dfall.to_pickle(outfn)
    logger.info(f"save to {outfn}")

    return dfall


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


def rename_reportdf(df):
    """
    Rename and change raw values of report df to more meaning full for display
    :param df:
    :return:
    """

    col_maps = {'Csites'   : 'Csites_called', 'mCsites': 'mCsites_called', 'Csites1': '5mCs', 'mCsites1': '5Cs',
            'precision_5mC': 'Precision_5mC', 'precision_5C': 'Precision_5C', 'recall_5mC': 'Recall_5mC', 'recall_5C': 'Recall_5C',
            'accuracy'     : 'Accuracy', 'roc_auc': 'ROC_AUC', 'referenceCpGs': 'TotalSites', 'corrMix': 'Corr_Mix', 'corrAll': 'Corr_All'}
    df = df.rename(columns=col_maps)

    df = df.replace(to_replace="False", value="x.x.Genome-wide")

    df = df.replace(to_replace="hg38_singletons.bed", value="x.x.Singletons")
    df = df.replace(to_replace="hg38_nonsingletons.bed", value="x.x.Non-singletons")

    df['coord'] = df['coord'].str.replace("promoterFeature.flank_", "promoterFeature")

    df["Location"] = df["coord"].str.split(".", n=3, expand=True)[2]
    df['Location'] = df['Location'].replace({'cpgIslandExt'       : 'CpG Island',
                                                    'discordant'  : 'Discordant',
                                                    'concordant'  : 'Concordant',
                                                    'cpgShoresExt': 'CpG Shores', 'cpgShelvesExt': 'CpG Shelves', 'exonFeature': 'Exons', 'intergenic': 'Intergenic', 'intronFeature': 'Introns', 'promoterFeature500': 'Promoters',
                                                    'absolute'    : 'Absolute',
                                                    })

    return df


def select_locations_from_reportdf(df, locations=locations_category2 + locations_singleton2):
    """
    Select only interested locations
    :param df:
    :param locations:
    :return:
    """
    return df[df['Location'].isin(locations)]


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


def save_newly_df():
    df = collect_newly_exp_data()
    outfn = os.path.join(pic_base_dir, f"combined_results_report_time_{current_time_str()}.xlsx")
    df.to_excel(outfn)


def create_report_datadf(runPrefixList=['APL_Bsseq_cut10', 'HL60_AML_Bsseq_cut5', 'K562_WGBS_joined_cut10', 'NA19240_RRBS_joined_cut10'], ret_col=ret_report_columns):
    """
    create report from list of runPrefix, return specified columns

    data are from: pkl/nanocompare/<runPrefix>

    APL_Bsseq_cut10/
    HL60_AML_Bsseq_cut5/
    K562_WGBS_joined_cut10/
    NA19240_RRBS_joined_cut5/

    :return:
    """
    dfs = []
    base_dir = os.path.join(data_base_dir, 'perf-plot-data')
    for runPrefix in runPrefixList:
        pattern = os.path.join(base_dir, runPrefix, "performance_results", "*.report.tsv")
        # print(pattern)
        files = glob.glob(pattern)
        for infn in files:
            df = pd.read_csv(infn, index_col=0, sep="\t")
            dfs.append(df)
    combdf = pd.concat(dfs, ignore_index=True)

    retdf = rename_reportdf(combdf)

    retdf = retdf[ret_col]
    return retdf


def collect_newly_exp_data():
    """
    Collect the currently new performance of exp results for paper
    :return:
    """
    # specify which runPrefix is the newly results you need
    runPrefixList = ['K562_WGBS_joined_cut5', 'APL_Bsseq_cut5', 'HL60_AML_Bsseq_cut5', 'NA19240_RRBS_joined_cut5']

    # specify your cared columns
    ret_col = ret_report_columns
    # ret_col = ret_report_columns[0:-2]

    df = create_report_datadf(runPrefixList=runPrefixList, ret_col=ret_col)
    retdf = select_locations_from_reportdf(df)

    logger.debug(f"retdf={retdf}")

    return retdf


if __name__ == '__main__':
    # df = load_singleton_nonsingleton_sites()
    df = collect_newly_exp_data()
    logger.info(f'df={df}')
    save_newly_df()

    # df = get_data_all()
    # save_alldata(df)
