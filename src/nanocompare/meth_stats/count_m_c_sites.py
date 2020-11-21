import csv
import glob
import os

from pybedtools import BedTool

from lilab.tcga.global_tcga import logger, pic_base_dir, set_log_debug_level
from nanocompare.load_data import load_all_perf_data_for_dataset_list
import pandas as pd

from nanocompare.meth_stats.Universal_meth_stats_evaluation import NonSingletonsScanner, importPredictions_Nanopolish_2_nofilter, dict2txt
from nanocompare.meth_stats.report_generation import load_ontcall_by_tool
from nanocompare.nanocompare_global_settings import important_region_bed_fns


def main():
    dsname_dict = {'K562': 'K562_WGBS_joined',
            'APL'        : 'APL_BSseq_cut10',
            'HL60'       : 'HL60_AML_Bsseq_cut5',
            'NA19240'    : 'NA19240_RRBS_joined'}

    dsname_dict_inv = {v: k for k, v in dsname_dict.items()}

    df = load_all_perf_data_for_dataset_list()

    test = (df['Tool'] == 'Nanopolish') & df['Location'].isin(['singletons', 'nonsingletons', 'concordant', 'discordant'])

    df = df[test]

    df['total.sites'] = df['mCsites'] + df['Csites']

    df = df[['Dataset', 'Location', 'mCsites', 'Csites', 'total.sites', 'referenceCpGs']]

    df1 = df.pivot(index='Dataset', columns='Location')
    df1.columns = df1.columns.swaplevel()
    df1.sort_index(axis=1, inplace=True)
    df1 = df1.rename(index=dsname_dict_inv)

    outfn = os.path.join(pic_base_dir, "count_from_original_results_table.xlsx")
    df1.to_excel(outfn)

    outfn = os.path.join(pic_base_dir, "count_from_original_results_table.pkl")
    df1.to_pickle(outfn)

    outfn = os.path.join(pic_base_dir, "count_from_original_results_table_before_pivot.xlsx")
    df.to_excel(outfn)

    pass


def load_count_ds_original():
    infn = os.path.join(pic_base_dir, 'count_from_original_results_table.pkl')
    df = pd.read_pickle(infn)
    return df


def count_bilogogical_regions():
    basedir = '/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/reports'
    pattern = os.path.join(basedir, '*.bed')
    files = glob.glob(pattern)

    df = pd.DataFrame()
    for fn in files:  # TODO end - start to get range count
        num_lines = sum(1 for line in open(fn))
        logger.debug(f'fn={fn}, num_lines={num_lines}')

        fbed1 = BedTool(fn)
        fbed2 = BedTool(fn)

        # logger.debug(f'fbed={fbed.count()}')
        # logger.debug(f'len={len(fbed)}')

        fbed = fbed1.intersect(fbed2, wa=True, wb=True)

        for ovr in fbed:
            regionKey = "{}\t{}\t{}\n".format(ovr[0], ovr[1], ovr[2])
            methKey = "{}\t{}\t{}\n".format(ovr[3], ovr[4], ovr[5])
            logger.info(f'regionKey={regionKey}, methKey={methKey}, ovr={ovr}')

        return

        basename = os.path.basename(fn)
        ret = {'Location': basename, 'CpGs': num_lines}
        df = df.append(ret, ignore_index=True)

    df = df[['Location', 'CpGs']]

    rep_dict = {'ONT.hg38.promoterFeature.flank_500.bed': 'Promoters', 'ONT.hg38.exonFeature.bed': 'Exon', 'ONT.hg38.cpgIslandExt.bed': 'CpG Island', 'ONT.hg38.intergenic.bed': 'Intergenic', 'ONT.hg38.intronFeature.bed': 'Intron'}
    df = df.replace(rep_dict)

    df = df[df['Location'].isin(['Promoters', 'Exon', 'CpG Island', 'Intergenic', 'Intron'])]

    logger.info(f'df={df}')

    outfn = os.path.join(pic_base_dir, 'biological.region.count.xlsx')
    df.to_excel(outfn, index=False)


# def region_sites_count():
#     referenceGenomeFile = '/projects/li-lab/reference/hg38/hg38.fasta'
#     outfileName_s = os.path.join(pic_base_dir, 'singletons.bed')
#     outfileName_ns = os.path.join(pic_base_dir, 'nonsingletons.bed')
#
#     NonSingletonsScanner(referenceGenomeFile, outfileName_s, outfileName_ns)


def region_sites_count(tsvFilename):
    """
    Check all CpG sites in each region
    :return:
    """
    logger.debug(f"region_sites_count, load config file:{tsvFilename}")

    infile = open(tsvFilename, 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')

    bedregions = {}
    for bedfn in important_region_bed_fns:
        logger.debug(f'bedfn={bedfn}')
        basename = os.path.basename(bedfn)
        bedregions[basename] = BedTool(bedfn)

    df = pd.DataFrame()
    for row in csvfile:
        if row['status'] == "submit":
            logger.debug(f"read row={row}\n")
            runPrefix = row['RunPrefix']
            dsname = row['Dataset']
            encoder = row['parser']

            nanopolishCall = importPredictions_Nanopolish_2_nofilter(row['Nanopolish_calls'])
            ontCalls_bed = BedTool(dict2txt(nanopolishCall), from_string=True)
            ontCalls_bed = ontCalls_bed.sort()
            ret = {'dsname': dsname}
            for regionname, narrowBed in bedregions.items():
                # narrowBed = BedTool(bedFile)
                narrowBed = narrowBed.sort()
                ontCalls_intersect = ontCalls_bed.intersect(narrowBed, u=True, wa=True)
                ret[regionname] = ontCalls_intersect.count()
            df = df.append(ret, ignore_index=True)

    colrename = {'ONT.hg38.promoterFeature.flank_500.bed': 'Promoters', 'ONT.hg38.exonFeature.bed': 'Exon', 'ONT.hg38.intronFeature': 'Intron', 'ONT.hg38.intergenic.bed': 'Intergenic', 'ONT.hg38.cpgIslandExt.bed': 'CpG Island'}
    df = df.rename(columns=colrename)
    logger.info(f'df={df}')

    outfn = os.path.join(pic_base_dir, 'regions.sites.distribution.xlsx')
    df.to_excel(outfn, index=False)


if __name__ == '__main__':
    set_log_debug_level()
    tsvFilename = '/projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComarePerformance_paper2.tsv'
    region_sites_count(tsvFilename)
