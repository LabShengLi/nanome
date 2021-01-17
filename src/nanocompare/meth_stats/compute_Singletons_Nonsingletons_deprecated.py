import os

from pybedtools import BedTool

from lilab.tcga.global_tcga import logger, set_log_debug_level, pic_base_dir, pkl_base_dir, ensure_dir
from nanocompare.meth_stats.Universal_meth_stats_evaluation import importGroundTruth_bed_file_format, nonSingletonsPostprocessing, singletonsPostprocessing, singletonsPostprocessing2, nonSingletonsPostprocessing2, dict2txt
import pandas as pd

from nanocompare.global_settings import gbtruth_filedict, singletonsFile, nonsingletonsFile
from nanocompare.plot_figures_deprecated import plot_pie_chart

joined_bed_file = {
        'K562'   : '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/reports/K562_WGBS_joined/K562_WGBS_joined.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed',
        'APL'    : '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/reports/APL_BSseq_cut10/APL_Bsseq_cut10.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed',
        'HL60'   : '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/reports/HL60_AML_Bsseq_cut5/HL60_AML_Bsseq_cut5.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed',
        'NA19240': '/projects/liuya/results/pkl/nanocompare/NA19240_RRBS_joined/NA19240_RRBS_joined.DeepSignal_Tombo_Nanopolish_DeepMod_Bacground.bed'
        }

minCov_dict = {
        'K562'   : 10,
        'APL'    : 10,
        'HL60'   : 5,
        'NA19240': 10
        }


def gen_bed_files(minCov=10):
    """
    Generate bed files based on background truth, especially for Singletons and Non-Singletons
    :return:
    """

    reportdf = pd.DataFrame()
    for dsname, infn in gbtruth_filedict.items():
        logger.info(f"dsname={dsname}, minCov={minCov}, infn={infn}")

        RunPrefix = dsname

        bgTruth = importGroundTruth_bed_file_format(infn, covCutt=minCov)

        outdir = os.path.join(pic_base_dir, f'bed_sing_nonsing.{dsname}.minCov.{minCov}')
        ensure_dir(outdir)

        ret = singletonsPostprocessing2(bgTruth, singletonsFile, RunPrefix, outdir=outdir)
        ret1 = nonSingletonsPostprocessing2(bgTruth, nonsingletonsFile, RunPrefix, outdir=outdir)
        ret.update(ret1)
        ret.update({"dsname": dsname})

        reportdf = reportdf.append(ret, ignore_index=True)

    outfn = os.path.join(pic_base_dir, f'sing_nonsing_count.minCov.{minCov}.csv')
    reportdf.to_csv(outfn)

    logger.info(f"reportdf = {reportdf}")


def get_bed_by_type(dsname, minCov=10, type='absolute'):
    if type == 'absolute':
        infn = os.path.join(pkl_base_dir, 'nanocompare', 'CountDataset', f'bed_sing_nonsing.{dsname}.minCov.{minCov}', f'{dsname}.hg38_singletons.absolute.bed')
    elif type == 'concordant':
        infn = os.path.join(pkl_base_dir, 'nanocompare', 'CountDataset', f'bed_sing_nonsing.{dsname}.minCov.{minCov}', f'{dsname}.hg38_nonsingletons.concordant.bed')
    elif type == 'discordant':
        infn = os.path.join(pkl_base_dir, 'nanocompare', 'CountDataset', f'bed_sing_nonsing.{dsname}.minCov.{minCov}', f'{dsname}.hg38_nonsingletons.discordant.bed')
    else:
        raise Exception(f"unsupported type={type}")

    retbed = BedTool(infn)
    retbed = retbed.sort()
    return retbed


def joined_files_analyze(minCov=10):
    """
    Generate bed files based on background truth, especially for Singletons and Non-Singletons
    :return:
    """

    reportdf = pd.DataFrame()
    for dsname, gbtruthfn in gbtruth_filedict.items():

        minCov = minCov_dict[dsname]
        logger.info(f"dsname={dsname}, minCov={minCov}, infn={gbtruthfn}")

        bgTruth = importGroundTruth_bed_file_format(gbtruthfn, covCutt=minCov)
        # bgTruthBed = BedTool(dict2txt(bgTruth), from_string=True)
        # bgTruthBed = bgTruthBed.sort()

        joinedfn = joined_bed_file[dsname]
        joinedbed = BedTool(joinedfn)
        joinedbed = joinedbed.sort()

        absolutebed = get_bed_by_type(dsname, minCov=minCov, type='absolute')
        joined_absolute_bed = joinedbed.intersect(absolutebed, wa=True, u=True)

        conbed = get_bed_by_type(dsname, minCov=minCov, type='concordant')
        joined_con_bed = joinedbed.intersect(conbed, wa=True, u=True)

        discbed = get_bed_by_type(dsname, minCov=minCov, type='discordant')
        joined_disc_bed = joinedbed.intersect(discbed, wa=True, u=True)

        retdict = {
                'dsname'           : dsname,
                'minCov'           : minCov,
                'bgtruth'          : len(bgTruth),  # all bg
                'joined'           : len(joinedbed),  # joined with 4 tools and bg
                'Singletons'       : len(absolutebed),  # absolute in bg
                'Singletons.joined': len(joined_absolute_bed),  # absolute in bg join the joined
                'Concordant'       : len(conbed),
                'Concordant.joined': len(joined_con_bed),
                'Discordant'       : len(discbed),
                'Discordant.joined': len(joined_disc_bed),
                }

        reportdf = reportdf.append(retdict, ignore_index=True)

        logger.info(f"retdict={retdict}")
    logger.info(f"reportdf={reportdf}")
    outfn = os.path.join(pic_base_dir, f'sing_nonsing_joined.count.minCov.{minCov}.csv')
    reportdf.to_csv(outfn)

    #
    #
    #     outdir = os.path.join(pic_base_dir, f'bed_sing_nonsing.{dsname}.minCov.{minCov}')
    #     ensure_dir(outdir)
    #
    #     ret = singletonsPostprocessing2(bgTruth, singletonsFile, RunPrefix, outdir=outdir)
    #     ret1 = nonSingletonsPostprocessing2(bgTruth, nonsingletonsFile, RunPrefix, outdir=outdir)
    #     ret.update(ret1)
    #     ret.update({"dsname": dsname})
    #
    #
    # outfn = os.path.join(pic_base_dir, f'sing_nonsing_count.minCov.{minCov}.csv')
    # reportdf.to_csv(outfn)
    #
    # logger.info(f"reportdf = {reportdf}")


import glob


def count_results():
    basedir = "/projects/liuya/workspace/tcgajax/nanocompare/meth_stats/bed_singleton_nonsingletons"
    patstr = os.path.join(basedir, "*.bed")
    flist = glob.glob(patstr)

    col1 = []
    col2 = []
    for fn in flist:
        basefn = os.path.basename(fn)
        num_lines = sum(1 for line in open(fn))

        col1.append(basefn)
        col2.append(num_lines)

        logger.debug(f"basefn={basefn}, num_lines={num_lines}")

        pass
    data = {'ds_source': col1, 'count': col2}
    df = pd.DataFrame(data=data)
    logger.debug(f"df={df}")

    outfn = os.path.join(pic_base_dir, "singleton_nonsingleton_count.csv")
    df.to_csv(outfn)


def load_raw_count_df():
    infn = os.path.join(pic_base_dir, "singleton_nonsingleton_count.csv")
    df = pd.read_csv(infn, header=0, index_col=0)
    logger.debug(f"df={df}")
    return df


def getdsname(str):
    cut = str.find(".")
    return str[:cut]
    pass


def getlocation(str):
    str = str.replace(".bed", "")
    cut = str.rfind(".")
    return str[cut + 1:]


def convert_raw_df_to_new_df():
    df = load_raw_count_df()

    df['dsname'] = df['ds_source'].apply(getdsname)
    df['location'] = df['ds_source'].apply(getlocation)

    df = df[['dsname', 'location', 'count', 'ds_source']]

    df = df.sort_values(by=['dsname', 'location'])
    df.reset_index(drop=True, inplace=True)

    logger.debug(f"df={df}")

    outfn = os.path.join(pic_base_dir, "singleton_nonsingleton_count_with_locations.csv")
    df.to_csv(outfn)

    pass


if __name__ == '__main__':
    set_log_debug_level()

    joined_files_analyze()

    # for k in [1, 5, 10]:
    #     gen_bed_files(k)

    pass
