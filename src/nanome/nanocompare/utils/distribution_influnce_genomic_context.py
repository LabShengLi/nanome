import os
from collections import defaultdict

import pandas as pd
import scipy

from nanome.common.global_config import pic_base_dir
from nanome.common.global_settings import ToolNameList

perfInfileName = os.path.join('/projects/li-lab/yang/results/2021-08-07/',
                              'performance-results.csv')
distributionInfileName = os.path.join('/projects/li-lab/yang/results/2021-08-07/',
                                      'all.certain.sites.distribution.each.genomic.region.cov5.csv')

genomic_regions = \
    ["Promoters", "Exons", "Introns", "Intergenic", "CpG Island", "CpG Shores", "CpG Shelves"] + \
    ['CG_20', 'CG_40', 'CG_60', 'CG_80', 'CG_100'] + \
    ['rep_SINE', 'rep_LINE', 'rep_LTR', 'rep_DNA', 'rep_Others']

if False:
    perfInfileName = os.path.join('/projects/li-lab/yang/results/2021-08-15', 'performance-results-seven.csv')
    distributionInfileName = os.path.join('/projects/li-lab/yang/results/2021-08-15',
                                          'all.certain.sites.distribution.each.genomic.region.seven.cov5.table.s6.csv')

perfDf = pd.read_csv(perfInfileName).loc[:, ['Dataset', 'Tool', 'Location', 'Macro-F1']]
perfDf = perfDf.loc[perfDf['Location'].isin(genomic_regions)]
perfDf = perfDf[perfDf['Dataset'].isin(['NA12878', 'NA19240', 'APL', 'K562'])]
print(len(perfDf))

distDf = pd.read_csv(distributionInfileName).loc[:,
         ['Dataset', 'Coord', 'Singletons', 'Non-singletons', 'Concordant', 'Discordant']]
distDf = distDf.loc[distDf['Coord'].isin(genomic_regions)]
distDf = distDf[distDf['Dataset'].isin(['NA12878', 'NA19240', 'APL', 'K562'])]
distDf['%Non-singletons'] = distDf['Non-singletons'] / (distDf['Non-singletons'] + distDf['Singletons'])
distDf['%Discordant'] = distDf['Discordant'] / (distDf['Discordant'] + distDf['Concordant'])

print(len(distDf))

df = perfDf.merge(distDf, left_on=['Dataset', 'Location'], right_on=['Dataset', 'Coord'])
dataset = defaultdict(list)
for toolName in ToolNameList:
    ndf = df[df['Tool'] == toolName]
    print(len(ndf))
    xx = ndf['%Non-singletons']
    yy = ndf['Macro-F1']
    from scipy import stats

    coe, pvalue = scipy.stats.pearsonr(xx, yy)
    print(f"{toolName}:coe={coe}, pvalue={pvalue}")
    dataset['Method'].append(toolName)
    dataset['COE'].append(coe)
    dataset['P-value'].append(pvalue)
outdf = pd.DataFrame.from_dict(dataset)
outfn = os.path.join(pic_base_dir, 'coe.influence.genomic.regions.percent.nonsingleton.csv')
outdf.to_csv(outfn)

dataset = defaultdict(list)
for toolName in ToolNameList:
    ndf = df[df['Tool'] == toolName]
    print(len(ndf))
    xx = ndf['%Discordant']
    yy = ndf['Macro-F1']

    coe, pvalue = scipy.stats.pearsonr(xx, yy)
    print(f"{toolName}:coe={coe}, pvalue={pvalue}")
    dataset['Method'].append(toolName)
    dataset['COE'].append(coe)
    dataset['P-value'].append(pvalue)
outdf = pd.DataFrame.from_dict(dataset)
outfn = os.path.join(pic_base_dir, 'coe.influence.genomic.regions.percent.discordant.csv')
outdf.to_csv(outfn)
