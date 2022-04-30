"""
Combine all chr methylation-calling results into a file for each tool
"""
import glob
import gzip
import os
import re
import sys

from nanome.common.global_config import set_log_debug_level, logger
from nanome.common.global_settings import HUMAN_CHR_SET

input_dir = "/projects/li-lab/Nanopore_compare/data/NA12878"
combine_dir = "/projects/li-lab/Nanopore_compare/data/NA12878/combine_allchrs"

## chr 1-22, chrX, chrY
# humanChrSet = [f'chr{k}' for k in range(1, 23)] + ['chrX', 'chrY']

toolSuffix = {
        'Nanopolish'     : 'nanopolish.methylation_calls',
        'Tombo'          : 'tombo.perReadsStats',
        'DeepSignal'     : 'deepsignal.call_mods',
        'DeepMod.C'      : 'deepmod.C',
        'DeepMod.Cluster': 'deepmod.C_clusterCpG',
        'Guppy'          : 'guppy.fast5mod_site_level',
        'Guppy.gcf52ref' : 'guppy.gcf52ref_read_level',
        'Megalodon'      : 'megalodon.per_read',
        }

sep_dict = {'Nanopolish' : '\t', 'Tombo': '\t', 'Megalodon': '\t', 'DeepSignal': '\t', 'DeepMod.C': ' ',
        'DeepMod.Cluster': ' ', 'Guppy': '\t', 'Guppy.gcf52ref': '\t'
        }
chr_col_dict = {'Nanopolish': 0, 'Tombo': 0, 'Megalodon': 1, 'DeepSignal': 0, 'DeepMod.C': 0,
        'DeepMod.Cluster'   : 0, 'Guppy': 0, 'Guppy.gcf52ref': 0}

if __name__ == '__main__':
    set_log_debug_level()

    # toolName = 'Nanopolish'
    toolName = sys.argv[1]

    os.makedirs(combine_dir, exist_ok=True)
    fnlist = glob.glob(os.path.join(input_dir, f'*/NA12878*.{toolSuffix[toolName]}.combine.*.gz'))
    logger.info(f'Find total files={len(fnlist)}')

    processDict = {}
    for fn in fnlist:
        ## such as: NA12878-CHR17.nanopolish.methylation_calls.combine.tsv.gz
        temp_name = os.path.basename(fn)
        matchObj = re.match(rf'NA12878.*(CHR.*)\.{toolSuffix[toolName]}\.combine\.(tsv|bed)\.gz', temp_name)
        if matchObj:
            chr = matchObj.group(1)
            bed_or_tsv = matchObj.group(2)
        else:
            raise Exception(f"regularEx match failed for: {temp_name}")
        processDict[chr] = fn

    logger.info(f'Construct total files={len(processDict)}')
    outfn = os.path.join(combine_dir, f'NA12878.{toolSuffix[toolName]}.combine_allchrs.{bed_or_tsv}.gz')
    with gzip.open(outfn, 'wt') as outf:
        for chr in HUMAN_CHR_SET:
            cntLines = 0
            if chr.upper() not in processDict:
                continue
            fn = processDict[chr.upper()]
            logger.info('\n\n==============================')
            logger.info(f'Start output: {fn}')
            with gzip.open(fn, 'rt')  as inf:
                for row in inf:
                    tmp = row.strip().split(sep_dict[toolName])
                    if tmp[chr_col_dict[toolName]] != chr:  ## filter out other chr outputs
                        continue
                    cntLines += 1
                    outf.write(row)
            logger.info(f'Finish {chr}, total lines={cntLines:,}')
    print(f"DONE combine all chrs for {toolName}")
