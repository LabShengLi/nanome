#!/usr/bin/env python3

from sys import argv

import numpy as np
from tombo import tombo_helper, tombo_stats

# python ../Tombo_extract_per_read_stats.py /projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/hg38.chrom.sizes TomboTest.tombo.per_read_stats    TomboTest.tombo.per_read_stats.bed
"""
Ouput format:

chr14	19118837	19118837	18eacd99-d2b9-48d8-a16c-cf28b351c66b	2.610044548594933	+
chr14	19118838	19118838	18eacd99-d2b9-48d8-a16c-cf28b351c66b	2.7027532865530004	+
chr14	19118852	19118852	18eacd99-d2b9-48d8-a16c-cf28b351c66b	5.296365664897572	+
"""

chromSizesInfile = argv[1]  # "/projects/li-lab/NanoporeData/WR_ONT_analyses/NanoCompare/hg38.chrom.sizes"
perReadStatsInfile = argv[2]  # "TomboTest_001.batch_0.5mC.tombo.per_read_stats"
outfileName = argv[3]  # "test.bed"

verbose = False
if len(argv) >= 5:
    if argv[4] in ['debug', 'verbose']:
        verbose = True


def parseChromSizesFile(infileName):
    result = {}
    infile = open(infileName, 'r')
    for row in infile:
        tmp = row.strip().split("\t")
        result[tmp[0]] = int(tmp[1])
    infile.close()
    return result


chromSizes = parseChromSizesFile(chromSizesInfile)
outfile = open(outfileName, "w")

per_read_stats = tombo_stats.PerReadStats(perReadStatsInfile)

total_reads = 0
for chrm in chromSizes:
    ####save plus strand
    int_data = tombo_helper.intervalData(chrm=chrm, start=1, end=chromSizes[chrm], strand='+')
    reg_per_read_stats_plus = per_read_stats.get_region_per_read_stats(int_data)
    if verbose:
        print(reg_per_read_stats_plus)
    if isinstance(reg_per_read_stats_plus, np.ndarray):
        """
        Structure of each cpg is as:
         (3505579, -0.54981702, '21a26cb9-0be2-4670-8694-a3cee91d49b8')
        """
        for cpg in reg_per_read_stats_plus:
            if verbose:
                print(cpg)
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrm, cpg[0], cpg[0], cpg[2][:], cpg[1], "+"))
            total_reads += 1

    ####save minus strand
    int_data = tombo_helper.intervalData(chrm=chrm, start=1, end=chromSizes[chrm], strand='-')
    reg_per_read_stats_minus = per_read_stats.get_region_per_read_stats(int_data)
    if verbose:
        print(reg_per_read_stats_minus)
    if isinstance(reg_per_read_stats_minus, np.ndarray):
        for cpg in reg_per_read_stats_minus:
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrm, cpg[0] - 1, cpg[0] - 1, cpg[2][:], cpg[1], "-"))
            total_reads += 1

outfile.close()
print(f"### Total reads extracted={total_reads:,} from file={perReadStatsInfile}")
print("### Tombo extract bed file DONE")
