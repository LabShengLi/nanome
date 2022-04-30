import os.path
from collections import defaultdict

import pandas as pd

from nanocompare.eval_common import import_call
from nanocompare.global_config import set_log_debug_level, logger, pic_base_dir


def combine_both_strand_by_meteore(call_dict):
    """
    Combine both strand as METEORE paper
    :param call_dict:
    :return:
    """
    combine_dict = defaultdict(tuple)

    for key in call_dict:
        (chr, start, strand) = key
        if strand == '+':
            pos_key = (chr, start, '+')
            neg_key = (chr, start + 1, '-')
            combine_key = (chr, start)
        elif strand == '-':
            pos_key = (chr, start - 1, '+')
            neg_key = (chr, start, '-')
            combine_key = (chr, start - 1)
        else:
            raise Exception(f"Not correct key={key}")

        if combine_key in combine_dict:
            continue
        ## Combine both strand, and average the methylation freqency

        combine_list = list()
        if pos_key in call_dict and neg_key in call_dict:
            pos_cov = len(call_dict[pos_key])
            pos_freq = sum(call_dict[pos_key]) / len(call_dict[pos_key])
            neg_cov = len(call_dict[neg_key])
            neg_freq = sum(call_dict[neg_key]) / len(call_dict[neg_key])
            ## Ref: https://github.com/comprna/METEORE/blob/317484e6ae33c3bfc6d84f83dd73d1863965d483/script/run_deepmod.R#L46
            coverage = pos_cov + neg_cov
            meth_freq = 0.5 * (pos_freq + neg_freq)
        elif pos_key in call_dict:
            combine_list += list(call_dict[pos_key])
            coverage = len(call_dict[pos_key])
            meth_freq = sum(call_dict[pos_key]) / len(call_dict[pos_key])
        elif neg_key in call_dict:
            combine_list += list(call_dict[neg_key])
            coverage = len(call_dict[neg_key])
            meth_freq = sum(call_dict[neg_key]) / len(call_dict[neg_key])
        # meth_freq = sum(combine_list) / len(combine_list)
        # coverage = len(combine_list)
        combine_dict[combine_key] = (meth_freq, coverage)
        # logger.debug(f"combine_key={combine_key}, meth_freq={meth_freq}, coverage={coverage}")
    return combine_dict


def get_keys_from_meteore():
    """
    Using chr pos in METEORE paper for sanity check, we are sure it is 1-based in pos
    :param df:
    :return:
    """
    infn_meteore_supp_data11 = '/fastscratch/liuya/nanocompare/deepmod_problem/Supplementary_data_11.xlsx'
    df_meteore = pd.read_excel(infn_meteore_supp_data11, engine='openpyxl')
    logger.info(df_meteore)

    ret_keys = list()
    for index, row in df_meteore.iterrows():
        ret_keys += [(row['chr'], int(row['pos']))]
    return ret_keys


if __name__ == '__main__':
    set_log_debug_level()

    selected_keys = get_keys_from_meteore()

    infn_deepmod = '/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.deepmod.C.combine_allchrs.bed.gz'
    infn_nanopolish = '/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.nanopolish.methylation_calls.combine_allchrs.tsv.gz'
    infn_meteore = '/projects/li-lab/Nanopore_compare/suppdata/METEORE_results/NA12878.METEORE.megalodon_deepsignal-optimized-model-perRead.combine.tsv.gz'
    infn_megalodon = '/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.megalodon.per_read.combine_allchrs.bed.gz'
    infn_deepsignal = '/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.deepsignal.call_mods.combine_allchrs.tsv.gz'
    infn_tombo = '/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs/NA12878.tombo.perReadsStats.combine_allchrs.bed.gz'

    infn_bg1 = '/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF279HCL.bed.gz'
    infn_bg2 = '/projects/li-lab/Nanopore_compare/data/NA12878/ENCFF835NTC.bed.gz'

    baseFormat = 1
    using_cache = True
    enable_cache = True
    # using_cache = False
    # enable_cache = False

    ## Get input as key-> value: key=(chr, pos, strand), value = [0011000]
    dict_meteore = import_call(infn_meteore, "METEORE", baseFormat=baseFormat, include_score=False, siteLevel=False,
                               using_cache=using_cache, enable_cache=enable_cache)
    dict2_meteore = combine_both_strand_by_meteore(dict_meteore)

    dict_deepmod = import_call(infn_deepmod, "DeepMod.C", baseFormat=baseFormat, include_score=False, siteLevel=False,
                               using_cache=using_cache, enable_cache=enable_cache)
    dict2_deepmod = combine_both_strand_by_meteore(dict_deepmod)

    dict_nanopolish = import_call(infn_nanopolish, "Nanopolish", baseFormat=baseFormat, include_score=False,
                                  siteLevel=False,
                                  using_cache=using_cache, enable_cache=enable_cache)
    dict2_nanopolish = combine_both_strand_by_meteore(dict_nanopolish)

    dict_megalodon = import_call(infn_megalodon, "Megalodon", baseFormat=baseFormat, include_score=False,
                                 siteLevel=False,
                                 using_cache=using_cache, enable_cache=enable_cache)
    dict2_megalodon = combine_both_strand_by_meteore(dict_megalodon)

    dict_deepsignal = import_call(infn_deepsignal, "DeepSignal", baseFormat=baseFormat, include_score=False,
                                  siteLevel=False,
                                  using_cache=using_cache, enable_cache=enable_cache)
    dict2_deepsignal = combine_both_strand_by_meteore(dict_deepsignal)

    dict_tombo = import_call(infn_tombo, "Tombo", baseFormat=baseFormat, include_score=False,
                             siteLevel=False,
                             using_cache=using_cache, enable_cache=enable_cache)
    dict2_tombo = combine_both_strand_by_meteore(dict_tombo)

    dict2_list = [dict2_nanopolish, dict2_deepmod, dict2_meteore, dict2_megalodon, dict2_deepsignal, dict2_tombo]
    name_list = ['nanoplish', 'deepmod', 'meteore_rf', 'megalodon', 'deepsignal', 'tombo']

    logger.info("\n\nReport chr1 cpgs:")
    cpgs = [159199943, 159200094, 159200216, 159200386]
    for dict2, name in zip(dict2_list, name_list):
        for cpg in cpgs:
            key = ('chr1', cpg)
            if key in dict2:
                logger.info(f"key={key} in name={name}")
        pass

    dataset = defaultdict(list)
    for key in selected_keys:
        dataset['chr'].append(key[0])
        dataset['pos'].append(key[1])
        for k in range(len(name_list)):
            # if name_list[k] == 'meteore_rf' and (key not in dict2_deepsignal or key not in dict2_megalodon):
            #     dataset[f'{name_list[k]}_freq'].append(None)
            #     dataset[f'{name_list[k]}_cov'].append(None)
            #     continue
            if key in dict2_list[k]:
                dataset[f'{name_list[k]}_freq'].append(dict2_list[k][key][0])
                dataset[f'{name_list[k]}_cov'].append(dict2_list[k][key][1])
            else:
                dataset[f'{name_list[k]}_freq'].append(None)
                dataset[f'{name_list[k]}_cov'].append(None)

    df = pd.DataFrame.from_dict(dataset)
    outfn = os.path.join(pic_base_dir, 'sanity_meteore_deepmod.xlsx')
    df.to_excel(outfn)
    logger.info(f"save to {outfn}")

    pass
