import sys

import pandas as pd

from nanocompare.global_config import set_log_debug_level, logger


def study_file(infn,  nrows=None):
    logger.info(f"infn={infn}")
    df = pd.read_csv(infn, sep='\t', nrows=nrows)
    logger.info(df)
    df.info()
    seldf = df[df['Pos']==cpg1]
    logger.info(seldf)


if __name__ == '__main__':
    set_log_debug_level()
    cpg1 = 159199943
    cpgs = [159199943, 159200094, 159200216, 159200386]

    logger.info("\n\nStudy METEORE results in paper")
    infn_meteore = '/projects/li-lab/Nanopore_compare/suppdata/METEORE_results/NA12878.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz'
    study_file(infn_meteore)

    logger.info("\n\nStudy METEORE new generated intermediate")
    infn_intermediate = '/projects/li-lab/yang/results/2021-08-09/METEORE-combine/NA12878.METEORE.megalodon_deepsignal-optimized-model-intermediate.combinedf.tsv.gz'
    study_file(infn_intermediate)

    logger.info("\n\nStudy METEORE new generated results")
    infn_meteore2 = '/projects/li-lab/yang/results/2021-08-09/METEORE-combine/NA12878.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz'
    study_file(infn_meteore2)

    logger.info("DONE")
