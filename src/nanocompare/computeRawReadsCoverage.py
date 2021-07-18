#!/usr/bin/env python3
"""
Calculate Nanopore raw reads coverage at each CpG sites based on raw fast5 files
"""
import glob
import gzip
import os
import re
import sys
from collections import defaultdict
from multiprocessing import Pool

import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm

from nanocompare.eval_common import get_dna_seq_from_reference, open_file_gz_or_txt, find_bed_filename, get_ref_fasta
from nanocompare.global_config import set_log_debug_level, logger, pic_base_dir
from nanocompare.global_settings import humanChrSet, narrowCoordFileList, location_filename_to_abbvname

# used for convert region bed cov to base level cov
rawReadDir = '/pod/2/li-lab/Nanopore_compare/data/Nanopore_cov'
# rawReadDir = '/projects/li-lab/yang/results/2021-04-05'

# base level cov bed files (combined + - strand)
baseReadCovDir = '/projects/li-lab/yang/results/2021-04-05'
baseReadCovDir = '/projects/li-lab/Nanopore_compare/suppdata/raw.fast5.reads.coverage'

# modify this dir to newly concordant and discordant perf results base dir
datasetBedDir = '/projects/li-lab/yang/results/2021-06-26'


def convert_region_to_cpg_base(dsname):
    """
    Conver bed file of regions into CpG base bed files

    Assume infn is 0-based file, we print 1-based output. Because 1-based site level reports are easy to stats overlap with region files.
    :param infn:
    :return:
    """

    fnlist = glob.glob(os.path.join(rawReadDir, f'{dsname}*coverage.*strand.bed.gz'))
    logger.info(f'convert_region_to_cpg_base:{fnlist}')

    outfn = os.path.join(pic_base_dir, f'{dsname}.rawfast5.coverage.base.bed.gz')
    outfile = gzip.open(outfn, 'wt')

    print_first = True
    for infn in fnlist:
        logger.info(f'processing file={infn}')
        infile = open_file_gz_or_txt(infn)
        for row in tqdm(infile):
            tmp = row.strip().split(" ")
            if print_first:
                print_first = False
                logger.info(f'row={row}, tmp={tmp}')
            chr = tmp[0]
            start = int(tmp[1])
            end = int(tmp[2])
            cov = int(tmp[3])
            # TODO: get strand info, + strand link to CG's C, - strand link to CG's G
            strand = tmp[4]

            if chr not in humanChrSet:  # filter out non-human chrs
                continue

            # we want get seq at least two bases, evaluate 'CG' patterns
            # if end == start + 1:
            #     newend = end + 1
            # else:
            #     newend = end
            newend = end + 1

            dnastr = get_dna_seq_from_reference(chr, start, newend, ref_fasta=refFasta)
            for cpg in re.finditer("CG", dnastr):
                cpgstart = start + cpg.start()
                if strand == '+':  # point to CG's C, just report that position (0-based) into 1-based
                    outstart = cpgstart + 1
                elif strand == '-':  # negative strand, point to CG's G, based on seq of + strand
                    outstart = cpgstart + 2
                else:
                    raise Exception(f'strand={strand} is no suppport')

                out_txt = '\t'.join([chr, str(outstart), str(outstart), str(cov), '.', strand])
                outfile.write(f'{out_txt}\n')
                pass
        infile.close()
    outfile.close()


def preprocess_bed_to_cgp_base():
    dsname_list = ['HL60', 'K562', 'APL', 'NA19240']
    dsname_list = ['NA12878']

    for dsname in dsname_list:
        logger.info(f'dsname={dsname}')
        convert_region_to_cpg_base(dsname)


def count_sites_in_coord(readBed, coordfn, tagname, cutoff_list):
    """
    Counts how many CpG sites from readBed are also in range of coordinate file

    Sample readBed files:
    head HL60.rawfast5.coverage.base.bed
    chr1	10470	10470	8	.	-
    chr1	10472	10472	6	.	-
    chr1	10485	10485	4	.	-
    chr1	10490	10490	6	.	-

    :param readBed:
    :param coordfn:
    :param cutoff_list:
    :return:
    """
    ret = {}
    covList = []

    ##TODO: we can use to_dataframe() for BedTool
    for row in readBed:
        covList.append(int(row[3]))
    covSeries = pd.Series(covList)

    for i, cutoff in enumerate(cutoff_list):
        ret.update({f'total_sites_after_cutoff.cutoff{cutoff}': (covSeries >= cutoff).sum()})

    coordBed = BedTool(coordfn)
    coordBed = coordBed.sort()
    intersectBed = readBed.intersect(coordBed, u=True, wa=True)

    covList = []
    for row in intersectBed:
        cov = int(row[3])
        covList.append(cov)
    covSeries = pd.Series(covList)
    for i, cutoff in enumerate(cutoff_list):
        ret.update({f'{tagname}.cutoff{cutoff}': (covSeries >= cutoff).sum()})

    return ret


def get_len_of_bedtool(bed_obj, cov_col=3, cutoff=1):
    """
    Get the number of rows (sites) of BedTool object
    :param bed_obj:
    :param cutoff:
    :return:
    """
    # covList = []
    # for row in bed_obj:
    #     covList.append(int(row[cov_col]))
    # covSeries = pd.Series(covList)
    df = bed_obj.to_dataframe()
    return (df.iloc[:, cov_col] >= cutoff).sum()


def report_table():
    """
    Generate Figure 2 C and D data
    :return:
    """
    # Nanopore reads coverage cutoff
    cutoff_list = [1, cov_cutoff]

    dataset = []
    distribution_dataset = defaultdict(list)
    for dsname in dsname_list:
        fnlist = glob.glob(os.path.join(baseReadCovDir, f'{dsname}.rawfast5.coverage.base.bed'))
        fn = fnlist[0]
        logger.info(fn)

        rawReadBed = BedTool(fn)
        rawReadBed = rawReadBed.sort()
        logger.info(f'len(rawReadBed)={len(rawReadBed):,}')
        dataDict = {'dsname': dsname, 'filename': fn, 'total': len(rawReadBed)}
        retList = []
        with Pool(processes=8) as pool:
            for bedfn in narrowCoordFileList[1:]:
                tagname = location_filename_to_abbvname[os.path.basename(bedfn)]
                ret = pool.apply_async(count_sites_in_coord, (rawReadBed, bedfn, tagname,),
                                       kwds={'cutoff_list': cutoff_list})
                retList.append(ret)

            concordantFileName = find_bed_filename(basedir=datasetBedDir,
                                                   pattern=f'{dsname}*hg38_nonsingletons.concordant.bed')
            ret = pool.apply_async(count_sites_in_coord, (rawReadBed, concordantFileName, 'Concordant',),
                                   kwds={'cutoff_list': cutoff_list})
            retList.append(ret)

            discordantFileName = find_bed_filename(basedir=datasetBedDir,
                                                   pattern=f'{dsname}*hg38_nonsingletons.discordant.bed')
            ret = pool.apply_async(count_sites_in_coord, (rawReadBed, discordantFileName, 'Discordant',),
                                   kwds={'cutoff_list': cutoff_list})
            retList.append(ret)

            pool.close()
            pool.join()
        retList = [ret.get() for ret in retList]
        for ret in retList:
            dataDict.update(ret)

        dataset.append(dataDict)
        logger.debug(dataDict)

        ## Report singleton/nonsingleton at each genomic regions
        # rawReadBed - raw fast5 BedTool
        singletonFilename = narrowCoordFileList[1]
        singleton_bed = BedTool(singletonFilename).sort()

        nonsingletonFilename = narrowCoordFileList[2]
        nonsingleton_bed = BedTool(nonsingletonFilename).sort()

        concordantFileName = find_bed_filename(basedir=datasetBedDir,
                                               pattern=f'{dsname}*hg38_nonsingletons.concordant.bed')
        concordant_bed = BedTool(concordantFileName).sort()

        discordantFileName = find_bed_filename(basedir=datasetBedDir,
                                               pattern=f'{dsname}*hg38_nonsingletons.discordant.bed')
        discordant_bed = BedTool(discordantFileName).sort()

        for coordFilename in narrowCoordFileList[3:]:
            tagname = os.path.basename(coordFilename)
            coord_bed = BedTool(coordFilename).sort()
            intersect_coord_bed = rawReadBed.intersect(coord_bed, u=True, wa=True)
            num_total = get_len_of_bedtool(intersect_coord_bed, cutoff=cov_cutoff)

            intersect_singleton_bed = intersect_coord_bed.intersect(singleton_bed, u=True, wa=True)
            num_singleton = get_len_of_bedtool(intersect_singleton_bed, cutoff=cov_cutoff)

            intersect_nonsingleton_bed = intersect_coord_bed.intersect(nonsingleton_bed, u=True, wa=True)
            num_nonsingleton = get_len_of_bedtool(intersect_nonsingleton_bed, cutoff=cov_cutoff)

            intersect_concordant_bed = intersect_coord_bed.intersect(concordant_bed, u=True, wa=True)
            num_concordant = get_len_of_bedtool(intersect_concordant_bed, cutoff=cov_cutoff)

            intersect_discordant_bed = intersect_coord_bed.intersect(discordant_bed, u=True, wa=True)
            num_discordant = get_len_of_bedtool(intersect_discordant_bed, cutoff=cov_cutoff)

            distribution_dataset['Dataset'].append(dsname)
            distribution_dataset['Coord'].append(location_filename_to_abbvname[tagname])
            distribution_dataset['Total'].append(num_total)
            distribution_dataset['Singletons'].append(num_singleton)
            distribution_dataset['Non-singletons'].append(num_nonsingleton)
            distribution_dataset['Concordant'].append(num_concordant)
            distribution_dataset['Discordant'].append(num_discordant)

    ## dataframe for distribution of singleton and non-singletons
    dist_df = pd.DataFrame.from_dict(distribution_dataset)
    outfn = os.path.join(pic_base_dir, f'nanopore.raw.fast5.distribution.at.each.genomic.region.cov{cov_cutoff}.xlsx')
    dist_df.to_excel(outfn)

    ## dataframe for raw reads at each genomic regions
    df = pd.DataFrame(dataset)
    logger.info(df)

    for cutoff in cutoff_list:
        outdf = pd.concat([df.loc[:, ['dsname', 'filename', 'total']], df.filter(regex=f'.cutoff{cutoff}$', axis=1)],
                          axis=1)
        outfn = os.path.join(pic_base_dir, f'raw.fast5.reads.cpg.coverage.across.regions.cutoff{cutoff}.xlsx')
        outdf.to_excel(outfn)
        logger.info(f'save to {outfn}')


def combine_na12878_coverage_bed():
    baseDir = "/fastscratch/liuya/nanocompare/NA12878-coverage"
    outfn = os.path.join(pic_base_dir, "NA12878-allChrs.coverage.bothstrand.bed.gz")
    outf = gzip.open(outfn, 'wt')
    for chrName in humanChrSet:
        logger.info(f"Processing chr={chrName}")
        flist = glob.glob(os.path.join(baseDir, f"NA12878-{chrName.upper()}.coverage.*.bed.gz"))
        logger.info(flist)

        for fn in flist:
            logger.info(fn)
            with gzip.open(fn, 'rt') as inf:
                for row in tqdm(inf):
                    tmp = row.strip().split(" ")
                    chr = tmp[0]
                    if chr != chrName:
                        continue
                    outf.write(f"{row}")

    outf.close()
    logger.info(f"save to {outfn}")


if __name__ == '__main__':
    set_log_debug_level()

    # Combine NA12878 all chrs bed regions
    if False:
        combine_na12878_coverage_bed()

    # Generate bed sites file from region file using sequencing string
    if True:
        refFasta = get_ref_fasta()
        preprocess_bed_to_cgp_base()
    sys.exit(0)

    if True:
        dsname_list = ['HL60', 'K562', 'APL', 'NA19240']
        cov_cutoff = 3
        refFasta = None
        report_table()
