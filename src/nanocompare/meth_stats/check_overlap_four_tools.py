import csv

from lilab.tcga.global_tcga import *
from nanocompare.meth_stats.Universal_meth_stats_evaluation import importPredictions_DeepSignal, save_ontcalls_to_pkl, importPredictions_Tombo, importPredictions_Tombo_nofilter, importPredictions_Nanopolish_2_nofilter, importPredictions_DeepMod, load_ontcalls_pkl, combine2programsCalls
import pandas as pd


def load_nofilter_ontcalls_by_tool(dsname, tool):
    """
    Load ont Calls from pkl
    :param runPrefix:
    :param tool:
    :return:
    """
    infn = os.path.join('/projects/li-lab/nmf_epihet/results/02-25', f'{dsname}_nofilter_ontcalls', f"Ontcalls.{tool}.{dsname}.pkl")
    ret = load_ontcalls_pkl(infn)
    logger.debug(f'load nofilter ontcall by {tool} from {dsname} ok. CpGs={len(ret)}')
    return ret


def compare_tools_ontcalls(tsvFilename):
    logger.debug(f"compare_tools_ontcalls, load config file:{tsvFilename}")

    infile = open(tsvFilename, 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')

    for row in csvfile:
        if row['status'] == "submit":
            logger.debug(f"read row={row}")
            dsname = row['Dataset']

            deepsignalCalls = load_nofilter_ontcalls_by_tool(dsname, 'DeepSignal')

            deepmodCalls = load_nofilter_ontcalls_by_tool(dsname, 'DeepMod')
            nanopolishCalls = load_nofilter_ontcalls_by_tool(dsname, 'Nanopolish')
            tomboCalls = load_nofilter_ontcalls_by_tool(dsname, 'Tombo')

            DeepSignal_Tombo = combine2programsCalls(deepsignalCalls, tomboCalls)
            DeepSignal_Tombo_Nanopolish = combine2programsCalls(DeepSignal_Tombo, nanopolishCalls)
            DeepSignal_Tombo_Nanopolish_DeepMod = combine2programsCalls(DeepSignal_Tombo_Nanopolish, deepmodCalls)

            logger.info(f'Joined CpGs:{len(DeepSignal_Tombo_Nanopolish_DeepMod)}')

            outdirONT = os.path.join(pic_base_dir, f'{dsname}_nofilter_ontcalls')
            ensure_dir(outdirONT)

            outfn = os.path.join(outdirONT, f"Ontcalls.Joined.{dsname}.pkl")
            save_ontcalls_to_pkl(DeepSignal_Tombo_Nanopolish_DeepMod, outfn)

            diffreportdf = pd.DataFrame()

            cnt = 0
            for ontCall in DeepSignal_Tombo_Nanopolish_DeepMod:
                if cnt >= 10000:
                    break
                ontsplit = ontCall[0:-1].split('\t')
                chr = ontsplit[0]
                start = ontsplit[1]
                end = ontsplit[2]
                ret = {'chr'        : chr, 'start': start, 'end': end,
                        'DeepSignal': len(deepsignalCalls[ontCall]),
                        'DeepMod'   : len(deepmodCalls[ontCall]),
                        'Nanopolish': len(nanopolishCalls[ontCall]),
                        'Tombo'     : len(tomboCalls[ontCall]), }
                diffreportdf = diffreportdf.append(ret, ignore_index=True)
                cnt += 1

            diffreportdf = diffreportdf[['chr', 'start', 'end', 'DeepSignal', 'DeepMod', 'Nanopolish', 'Tombo']]
            outfn = os.path.join(pic_base_dir, f'Call.diff.{dsname}.top100.tsv')
            diffreportdf.to_csv(outfn, index=False, sep='\t')


def get_and_save_ontcalls_nofilters(tsvFilename):
    logger.debug(f"get_and_save_ontcalls_no_filters, load config file:{tsvFilename}")

    infile = open(tsvFilename, 'r')
    csvfile = csv.DictReader(infile, delimiter='\t')

    for row in csvfile:
        if row['status'] == "submit":
            logger.debug(f"read row={row}")

            dsname = row['Dataset']
            runPrefix = row['RunPrefix']
            encoder = row['parser']
            bgtruth_fn = row['bgTruth']
            minCov = int(row['minCov'])

            outdirONT = os.path.join(pic_base_dir, f'{dsname}_nofilter_ontcalls')
            ensure_dir(outdirONT)

            DeepSignal_calls = importPredictions_DeepSignal(row['DeepSignal_calls'])
            outfn = os.path.join(outdirONT, f"Ontcalls.DeepSignal.{dsname}.pkl")
            save_ontcalls_to_pkl(DeepSignal_calls, outfn)

            Tombo_calls = importPredictions_Tombo_nofilter(row['Tombo_calls'])
            outfn = os.path.join(outdirONT, f"Ontcalls.Tombo.{dsname}.pkl")
            save_ontcalls_to_pkl(Tombo_calls, outfn)

            Nanopolish_calls = importPredictions_Nanopolish_2_nofilter(row['Nanopolish_calls'])
            outfn = os.path.join(outdirONT, f"Ontcalls.Nanopolish.{dsname}.pkl")
            save_ontcalls_to_pkl(Nanopolish_calls, outfn)

            DeepMod_calls = importPredictions_DeepMod(row['DeepMod_calls'])
            outfn = os.path.join(outdirONT, f"Ontcalls.DeepMod.{dsname}.pkl")
            save_ontcalls_to_pkl(DeepMod_calls, outfn)

        logger.info(f"####### get_and_save_ontcalls_no_filters finished on {dsname}!")

    pass


def main():
    pass


if __name__ == '__main__':
    set_log_debug_level()
    tsvFilename = '/projects/liuya/workspace/tcgajax/nanocompare/meth_stats/NanoComarePerformance_paper2.tsv'

    # get_and_save_ontcalls_nofilters(tsvFilename)
    compare_tools_ontcalls(tsvFilename)

    pass
