"""
remove duplicate reads of $1 fastq file, save it into $2 fastq file
"""

# example run command: python /home/rosikw/programs/scripts/nanopore_nanopolish.NA19240_pipeline.step_02.preindexing_checkDups.py NA19240_top400.dups.fastq NA19240_top400.noDups.fastq

from sys import argv

infile = open(argv[1], 'r')
outfile = open(argv[2], 'w')

line = 0
readID = -1
seenReadIDs = {}  # {readID : 1}
for row in infile:
    if (line % 4) == 0:  # it means we deal with a new read, so we need to process the previous one
        if (readID != -1) and (readID not in seenReadIDs):
            outfile.write(fqRead)
            seenReadIDs[readID] = 1
            fqRead = ""
            readID = row.split(" ")[0]
        else:
            readID = row.split(" ")[0]
            fqRead = ""
    fqRead += row
    line += 1

infile.close()
outfile.close()
