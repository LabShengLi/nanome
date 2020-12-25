# coding: utf-8

import argparse
import os
import subprocess
import sys

"""
ly usage:

HL60:

python BS_Seq_Pipeline.py -o /projects/li-lab/yang/results/2020-12-21/hl60-results -c /projects/li-lab/yang/results/2020-12-21/HL60/ENCSR000DDM_new_analysis/configFile.tsv -dt rrbs -g hg38 -gBis /projects/li-lab/yang/workspace/nano-compare/data/reference/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/ -gBwa /projects/li-lab/yang/workspace/nano-compare/data/reference/hg38/bwaMeth/hg38.fa --test

"""


# #### Coments, limitations, to do:
# ##### current version of the pipeline only supports Bismark and paired end reads analysis
# 
# Example command used to run this pipeline:
# python /home/rosikw/programs/pipelines/BS_Seq/BS_Seq_Pipeline_03.py -o /fastscratch/rosikw/bsseq_03_05 -c /home/rosikw/programs/pipelines/BS_Seq/exampleConfigFile.tsv -dt rrbs -g hg19 -gBis /projects/li-lab/reference/hg19/ --test
# python /home/rosikw/programs/pipelines/BS_Seq/BS_Seq_Pipeline_04.py -p 2 -o /fastscratch/rosikw/bsseq_04_02_hg38 -c /home/rosikw/programs/pipelines/BS_Seq/exampleConfigFile.tsv -dt rrbs -g hg38 -gBis /home/rosikw/projects/reference/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/ -gBwa /projects/li-lab/reference/hg19/hg19.fa --test
# python /home/rosikw/programs/pipelines/BS_Seq/BS_Seq_Pipeline_05.py -p 3 -o /fastscratch/rosikw/CasilioME_hg38 -c /home/rosikw/projects/CasilioME/configFile.BSseq_auto.tsv -dt rrbs -g hg38 -gBis /home/rosikw/projects/reference/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/ -gBwa /projects/li-lab/reference/hg19/hg19.fa --test
# python /home/rosikw/programs/pipelines/BS_Seq/BS_Seq_Pipeline_06.py -p 3 -o /fastscratch/rosikw/CasilioME_hg38_pip06_01 -c /home/rosikw/projects/CasilioME/configFile.BSseq_auto.tsv -dt rrbs -g hg38 -gBis /home/rosikw/projects/reference/hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/ -gBwa /home/rosikw/projects/reference/hg38/bwaMeth/hg38.fa --test


# python /home/rosikw/programs/pipelines/BS_Seq/BS_Seq_Pipeline_06.py -p 3 -o /fastscratch/rosikw/ACE_seq_mm10_hg38BWA_pip06_01 -c /home/rosikw/projects/ace-seq/configFile.tsv -dt wgbs -g mm10 -gBis /home/rosikw/li-lab/reference/mm10/ -gBwa /home/rosikw/projects/reference/hg38/bwaMeth/hg38.fa --test
# python /home/rosikw/programs/pipelines/BS_Seq/BS_Seq_Pipeline_06.py -p 3 -o /fastscratch/rosikw/ACE_seq_mm10_pip06_02 -c /home/rosikw/projects/ace-seq/configFile.tsv -dt wgbs -g mm10 -gBis /home/rosikw/li-lab/reference_BETA/Mus_musculus/UCSC/mm10/Sequence/BismarkIndex/0.16.1_bowtie2_2.3.1/ -gBwa /home/rosikw/li-lab/reference_BETA/Mus_musculus/UCSC/mm10/Sequence/BWAmethIndex/0.2.0/genome.fa --test


# ### Funtions:

def get_script_path():
    # based on: https://stackoverflow.com/questions/4934806/how-can-i-find-scripts-directory-with-python
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def str2bool(v):
    # based on https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def dataTypeParser(v):
    if v.upper() == "RRBS":
        return "RRBS"
    elif v.upper() == "WGBS":
        return "WGBS"
    else:
        raise argparse.ArgumentTypeError('### Input Parameters Error: Wrong data type specified. -df flag should be either "RRBS" or "WGBS"')


def prepareOutputsStructure(outDir, pipeline):
    #     This function will create the following output catalogs structure:
    #     outDir/
    #             ├── extractBismark
    #             ├── extractBWA_meth
    #             ├── mappingBismark
    #             ├── mappingBWA_meth
    #             ├── rawData
    #             │   └── QC
    #             ├── scripts
    #             └── trimData
    #                 └── QC
    #     Note that it will be affected by the pipeline type used (Bismark, BWA-meth or both)

    subprocess.Popen("mkdir {}".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    subprocess.Popen("mkdir {}/rawData".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    #     subprocess.Popen("mkdir {}/rawData/QC".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    subprocess.Popen("mkdir {}/rawData_BWA".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    #     subprocess.Popen("mkdir {}/rawData_BWA/QC".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    subprocess.Popen("mkdir {}/trimData".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    #     subprocess.Popen("mkdir {}/trimData/QC".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    subprocess.Popen("mkdir {}/trimData_BWA".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    #     subprocess.Popen("mkdir {}/trimData_BWA/QC".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    subprocess.Popen("mkdir {}/scripts".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    if pipeline == 3 or pipeline == 1:
        subprocess.Popen("mkdir {}/mappingBismark".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
        #         subprocess.Popen("mkdir {}/mappingBismark/QC".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
        subprocess.Popen("mkdir {}/extractBismark".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    #         subprocess.Popen("mkdir {}/extractBismark/QC".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    if pipeline == 3 or pipeline == 2:
        subprocess.Popen("mkdir {}/mappingBWA_meth".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    #         subprocess.Popen("mkdir {}/extractBWA_meth".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    print("### prepareOutputsStructure NOTE: The structure of the output was created.")


def parseConfigFile(configFile):
    #     Example config file should be a tab separated file that looks like this:
    #
    #     #fqName	prefix	readIdentity	replicate
    #     /fastscratch/rosikw/TF-1-RRBS/TF1-IDH2-Rep1-RRBS_R1.fastq.gz	TF1_IDH2	R1	1
    #     /fastscratch/rosikw/TF-1-RRBS/TF1-IDH2-Rep1-RRBS_R2.fastq.gz	TF1_IDH2	R2	1
    #     /fastscratch/rosikw/TF-1-RRBS/TF1-WT-Rep1-RRBS_R1.fastq.gz	TF1_WT	R1	1
    #     /fastscratch/rosikw/TF-1-RRBS/TF1-WT-Rep1-RRBS_R2.fastq.gz	TF1_WT	R2	1
    infile = open(configFile, 'r')
    config = {}  # will have the structure: {prefix.replicate : [R1.fastq, R2.fastq]}
    for row in infile:
        tmp = row.strip().split('\t')
        if row[:1] != "#":
            fastq = tmp[0]
            prefix = tmp[1]
            readIdentity = tmp[2]
            replicate = tmp[3]
            ID = "{}.Rep_{}".format(prefix, replicate)
            if ID not in config:
                config[ID] = [0, 0]
            if readIdentity == "R1":
                config[ID][0] = fastq
            elif readIdentity == "R2":
                config[ID][1] = fastq
            else:
                print("### parseConfigFile ERROR: wrong values have been used to describe reads identity. Should be either 'R1' or 'R2'. Please double check your configuration file.")
                return False
    infile.close()

    print("### parseConfigFile NOTE: Parser have detected {} samples for the processing: {}".format(len(config.keys()), config.keys()))
    return config


def runBismarkSingleOrPaired(sample, outDir, R1_oryginalName, R2_oryginalName, genomeBismark, referenceGenome, dataType, non_directional, zero_based, walltime, threads):
    if R2_oryginalName == 0:
        helixJobID = prepareBismarkPipeline_singleEnd(sample, outDir, R1_oryginalName, genomeBismark, referenceGenome, dataType, non_directional, zero_based, walltime, threads)
    else:
        helixJobID = prepareBismarkPipeline(sample, outDir, R1_oryginalName, R2_oryginalName, genomeBismark, referenceGenome, dataType, non_directional, zero_based, walltime, threads)
    return helixJobID


def prepareBismarkPipeline(sample, outDir, R1_oryginalName, R2_oryginalName, genomeBismark, referenceGenome, dataType, non_directional, zero_based, walltime, threads):
    script = "#!/bin/bash\n"

    if walltime > 72:
        qType = "long"
    else:
        qType = "batch"
    if walltime < 10:
        walltime = "0{}".format(walltime)

    script += """#PBS -q {} # batch - up to 72h, long - up to 504h. 
#PBS -l nodes=1:ppn={}
#PBS -l walltime={}:00:00
#PBS -o {}/scripts
#PBS -e {}/scripts\n\n
### Load modules:\n
module unload python
module load python/2.7.11
module load samtools/1.5
module load bamtools
module load bowtie2
module load bowtie/0.12.8
module load cutadapt/1.2.1
module load trim_galore/0.4.0
module load fastqc/0.11.5
module load bismark/0.16.1\n\n""".format(qType, threads, walltime, outDir, outDir)

    if 'q.gz' in R1_oryginalName:
        ext = "fastq.gz"
        extTrim = "fq.gz"
    else:
        ext = "fastq"
        extTrim = "fq"
    prefix = sample.split(".")[0]
    rep = sample.split(".")[1]
    R1 = "{}.Read_R1.{}.{}".format(prefix, rep, ext)
    R2 = "{}.Read_R2.{}.{}".format(prefix, rep, ext)
    R1_trimmed = "{}.Read_R1.{}_val_1.{}".format(prefix, rep, extTrim)
    R2_trimmed = "{}.Read_R2.{}_val_2.{}".format(prefix, rep, extTrim)
    bismarkDir = "{}/mappingBismark/".format(outDir)
    bismark_methylation_extractor_Dir = "{}/extractBismark/".format(outDir)
    bam = "{}/mappingBismark/{}.Read_R1.{}_val_1_bismark_bt2_pe.bam".format(outDir, prefix, rep)
    script += """### Input parameters:\n
outDir="{}"
R1_oryginalName="{}"
R2_oryginalName="{}"
R1="{}" #new name
R2="{}" #new name
R1_trimmed="{}" #new name
R2_trimmed="{}" #new name
genomeDir="{}"
bismarkDir="{}" # mappingBismark subcatalog
bismark_methylation_extractor_Dir="{}" # extractBismark subcatalog
bam="{}"
threads={}\n
#####################\n\n""".format(outDir, R1_oryginalName, R2_oryginalName, R1, R2, R1_trimmed, R2_trimmed, genomeBismark, bismarkDir, bismark_methylation_extractor_Dir, bam, threads)

    script += """### Analysis:\n
cd $outDir\n
# Copy fastq files:
cp $R1_oryginalName $outDir/rawData/$R1
cp $R2_oryginalName $outDir/rawData/$R2\n\n"""

    if non_directional:
        nd = " --non_directional"
    else:
        nd = ""
    if dataType == "RRBS":
        rrbs = " --rrbs"
    else:
        rrbs = ""
    #### I will need to add distinction beteween paired end and single end data here
    script += """# Trim reads:
trim_galore{}{} --stringency 3 -q 30 --output_dir $outDir/trimData --fastqc --paired $outDir/rawData/$R1 $outDir/rawData/$R2
\n""".format(nd, rrbs)

    script += """# Map reads:
bismark{} --output_dir $bismarkDir --multicore $threads -p 4 $genomeDir -1 $outDir/trimData/$R1_trimmed -2 $outDir/trimData/$R2_trimmed
\n""".format(nd)

    if dataType == "WGBS":
        dedup = "# Deduplicate reads (WGBS only)\n"
        dedup += "deduplicate_bismark -p --bam $bam\n"
        dedup += 'bam="{}/mappingBismark/{}.Read_R1.{}_val_1_bismark_bt2_pe.deduplicated.bam"\n\n'.format(outDir, prefix, rep)
    else:
        dedup = ""
    script += dedup

    if dataType == "RRBS":
        rrbs = " --ignore 10 --ignore_r2 10 --ignore_3prime 3 --ignore_3prime_r2 3"
    else:
        rrbs = ""
    if zero_based:
        zb = " --zero_based"
    else:
        zb = ""
    script += """# Extract methylation frequencies:
bismark_methylation_extractor -p {} --no_overlap --merge_non_CpG \
--report -o $bismark_methylation_extractor_Dir --gzip --multicore $threads --bedGraph {} \
--genome_folder $genomeDir --cytosine_report $bam""".format(rrbs, zb)

    ############
    outfile = open("{}/scripts/{}.bismarkPipeline.sh".format(outDir, sample), 'w')
    outfile.write(script)
    outfile.close()
    print("### prepareBismarkPipeline NOTE: Pipeline created for '{}' sample".format(sample))

    if test == False:
        helixJobID = subprocess.Popen("qsub {}/scripts/{}.bismarkPipeline.sh".format(outDir, sample), shell=True, stdout=subprocess.PIPE).stdout.read()
        print("### prepareBismarkPipeline NOTE: Pipeline '{}' sample have been submited to helix with ID {}".format(sample, helixJobID))
        return helixJobID.decode('ascii').strip()
    else:
        return "Test.Run"


def prepareBismarkPipeline_singleEnd(sample, outDir, R1_oryginalName, genomeBismark, referenceGenome, dataType, non_directional, zero_based, walltime, threads):
    script = "#!/bin/bash\n"

    if walltime > 72:
        qType = "long"
    else:
        qType = "batch"
    if walltime < 10:
        walltime = "0{}".format(walltime)

    script += """#PBS -q {} # batch - up to 72h, long - up to 504h. 
#PBS -l nodes=1:ppn={}
#PBS -l walltime={}:00:00
#PBS -o {}/scripts
#PBS -e {}/scripts\n\n

set -x

### Load modules:\n
module load python/2.7.11
module load samtools/1.5
module load bamtools
module load bowtie2
module load bowtie/0.12.8
module load cutadapt/1.2.1
module load trim_galore/0.4.0
module load fastqc/0.11.5
module load bismark/0.16.1\n\n""".format(qType, threads, walltime, outDir, outDir)

    if 'q.gz' in R1_oryginalName:
        ext = "fastq.gz"
        extTrim = "fq.gz"
    else:
        ext = "fastq"
        extTrim = "fq"
    prefix = sample.split(".")[0]
    rep = sample.split(".")[1]
    R1 = "{}.Read_R1.{}.{}".format(prefix, rep, ext)
    R1_trimmed = "{}.Read_R1.{}_trimmed.{}".format(prefix, rep, extTrim)
    bismarkDir = "{}/mappingBismark/".format(outDir)
    bismark_methylation_extractor_Dir = "{}/extractBismark/".format(outDir)
    bam = "{}/mappingBismark/{}.Read_R1.{}_trimmed_bismark_bt2.bam".format(outDir, prefix, rep)
    script += """### Input parameters:\n
outDir="{}"
R1_oryginalName="{}"
R1="{}" #new name
R1_trimmed="{}" #new name
genomeDir="{}"
bismarkDir="{}" # mappingBismark subcatalog
bismark_methylation_extractor_Dir="{}" # extractBismark subcatalog
bam="{}"
threads={}\n
#####################\n\n""".format(outDir, R1_oryginalName, R1, R1_trimmed, genomeBismark, bismarkDir, bismark_methylation_extractor_Dir, bam, threads)

    script += """### Analysis:\n
cd $outDir\n
# Copy fastq files:
cp $R1_oryginalName $outDir/rawData/$R1\n\n"""

    if non_directional:
        nd = " --non_directional"
    else:
        nd = ""
    if dataType == "RRBS":
        rrbs = " --rrbs"
    else:
        rrbs = ""
    #### I will need to add distinction beteween paired end and single end data here
    script += """# Trim reads:
trim_galore{}{} --stringency 3 -q 30 --output_dir $outDir/trimData --fastqc $outDir/rawData/$R1
module unload trim_galore/0.4.0\n
""".format(nd, rrbs)

    script += """# Map reads:
bismark{} --output_dir $bismarkDir --multicore $threads -p 4 $genomeDir $outDir/trimData/$R1_trimmed
\n""".format(nd)

    if dataType == "WGBS":
        dedup = "# Deduplicate reads (WGBS only)\n"
        dedup += "deduplicate_bismark -s --bam $bam\n"
        dedup += 'bam="{}/mappingBismark/{}.Read_R1.{}_trimmed_bismark_bt2.deduplicated.bam"\n\n'.format(outDir, prefix, rep)
    else:
        dedup = ""
    script += dedup

    if dataType == "RRBS":
        rrbs = " --ignore 10 --ignore_3prime 3 "
    else:
        rrbs = ""
    if zero_based:
        zb = " --zero_based"
    else:
        zb = ""
    script += """# Extract methylation frequencies:
bismark_methylation_extractor -s {} --no_overlap --merge_non_CpG \
--report -o $bismark_methylation_extractor_Dir --gzip --multicore $threads --bedGraph {} \
--genome_folder $genomeDir --cytosine_report $bam""".format(rrbs, zb)

    ############
    outfile = open("{}/scripts/{}.bismarkPipeline.sh".format(outDir, sample), 'w')
    outfile.write(script)
    outfile.close()
    print("### prepareBismarkPipeline NOTE: Pipeline created for '{}' sample".format(sample))

    if test == False:
        helixJobID = subprocess.Popen("qsub {}/scripts/{}.bismarkPipeline.sh".format(outDir, sample), shell=True, stdout=subprocess.PIPE).stdout.read()
        print("### prepareBismarkPipeline NOTE: Pipeline '{}' sample have been submited to helix with ID {}".format(sample, helixJobID))
        return helixJobID.decode('ascii').strip()
    else:
        return "Test.Run"


def QualityControl_PostBismark(dependentJobs, outDir, pipelinePath):
    dependentString = " -W depend=afterok:{} ".format(','.join(str(x) for x in dependentJobs))
    # for helixJobID in dependentJobs:
    #     dependentString += "','.join(str(x) for x in l) ".format(helixJobID)

    # QC for raw data:
    subfolder = "rawData"
    command = "qsub {0} -o {1}/{2}/QC -e {1}/{2}/QC -v outDir='{1}',pipelinePath='{3}',subfolder='{2}' {3}/fastqc_automated_masterScript.sh".format(dependentString, outDir, subfolder, pipelinePath)
    print(command)
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()

    # QC for trimmed data:
    subfolder = "trimData"
    command = "qsub {0} -o {1}/{2}/QC -e {1}/{2}/QC -v outDir='{1}',pipelinePath='{3}',subfolder='{2}' {3}/fastqc_automated_masterScript.sh".format(dependentString, outDir, subfolder, pipelinePath)
    print(command)
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()

    # QC for Bismark align folder:
    subfolder = "mappingBismark"
    command = "qsub {0} -o {1}/{2}/QC -e {1}/{2}/QC -v outDir='{1}',pipelinePath='{3}',subfolder='{2}' {3}/fastqc_automated_slave_multiqcComponent.bismark.sh".format(dependentString, outDir, subfolder, pipelinePath)
    print(command)
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()

    # QC for Bismark extracts:
    subfolder = "extractBismark"
    command = "qsub {0} -o {1}/{2}/QC -e {1}/{2}/QC -v outDir='{1}',pipelinePath='{3}',subfolder='{2}' {3}/fastqc_automated_slave_multiqcComponent.bismark.sh".format(dependentString, outDir, subfolder, pipelinePath)
    print(command)
    subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).stdout.read()


def runBWASingleOrPaired(sample, outDir, R1_oryginalName, R2_oryginalName, genomeBWA, referenceGenome, dataType, walltime, non_directional, threads):
    if R2_oryginalName == 0:
        helixJobID = prepareBWAmethPipeline_singleEnd(sample, outDir, R1_oryginalName, genomeBWA, referenceGenome, dataType, walltime, non_directional, threads)
    else:
        helixJobID = prepareBWAmethPipeline(sample, outDir, R1_oryginalName, R2_oryginalName, genomeBWA, referenceGenome, dataType, walltime, non_directional, threads)
    return helixJobID


def prepareBWAmethPipeline(sample, outDir, R1_oryginalName, R2_oryginalName, genomeBWA, referenceGenome, dataType, walltime, non_directional, threads):
    script = "#!/bin/bash\n"

    if walltime > 72:
        qType = "long"
    else:
        qType = "batch"
    if walltime < 10:
        walltime = "0{}".format(walltime)

    script += """#PBS -q {} # batch - up to 72h, long - up to 504h. 
#PBS -l nodes=1:ppn={}
#PBS -l walltime={}:00:00
#PBS -o {}/scripts
#PBS -e {}/scripts\n\n
### Load modules:\n
module load samtools/1.5
module load bwa/0.7.9a
module load cutadapt/1.2.1
module load trim_galore/0.4.0
module load fastqc/0.11.5\n\n
# create genome index. I assume it is already there, but if not, use commands like these to create them:
# cd /projects/li-lab/software/bwa-meth/bwa-meth-master/
# python bwameth.py index /projects/li-lab/reference/hg19\n\n\n""".format(qType, threads, walltime, outDir, outDir)

    if 'q.gz' in R1_oryginalName:
        ext = "fastq.gz"
        extTrim = "fq.gz"
    else:
        ext = "fastq"
        extTrim = "fq"
    prefix = sample.split(".")[0]
    rep = sample.split(".")[1]
    R1 = "{}.Read_R1.{}.{}".format(prefix, rep, ext)
    R2 = "{}.Read_R2.{}.{}".format(prefix, rep, ext)
    R1_trimmed = "{}.Read_R1.{}_val_1.{}".format(prefix, rep, extTrim)
    R2_trimmed = "{}.Read_R2.{}_val_2.{}".format(prefix, rep, extTrim)
    bwaDir = "{}/mappingBWA_meth/".format(outDir)
    bwa_methylation_extractor_Dir = "{}/extractBWA_meth/".format(outDir)
    bam = "{}.{}.bam".format(prefix, rep)
    sortedBam = "{}.{}.sorted.bam".format(prefix, rep)
    if dataType == "WGBS":
        biscuit = '\nbiscuit="/projects/li-lab/software/biscuit/biscuit-release/biscuit"'
        sortedDedupedBam = '\nsortedDedupedBam="{}.{}.sorted.dedup.bam"'.format(prefix, rep)
        sortedDedupedSortedBam = '\nsortedDedupedSortedBam="{}.{}.sorted.dedup.sorted.bam"'.format(prefix, rep)
    else:
        biscuit = ""
        sortedDedupedBam = ""
        sortedDedupedSortedBam = ""
    script += """### Input parameters:
reference="{}"
R1_oryginalName="{}"
R2_oryginalName="{}"
R1="{}" #new name
R2="{}" #new name
R1_trimmed="{}" #new name
R2_trimmed="{}" #new name
bam="{}"
sortedBam="{}"
outDir="{}"
threads={}

### These parameters are fixed, but may be modified if needed:
bwaPath="/projects/li-lab/software/bwa-meth/bwa-meth-master"
MethylDackel="/projects/li-lab/software/MethylDackel/MethylDackel"
bwameth="/projects/li-lab/software/bwa-meth/bwa-meth-master/bwameth.py"
{}
#####################\n\n""".format(genomeBWA, R1_oryginalName, R2_oryginalName, R1, R2, R1_trimmed, R2_trimmed, bam, sortedBam, outDir, threads, biscuit)

    script += """### Analysis:\n
cd $outDir\n
# Copy fastq files:
cp $R1_oryginalName $outDir/rawData_BWA/$R1
cp $R2_oryginalName $outDir/rawData_BWA/$R2\n\n"""

    if non_directional:
        nd = " --non_directional"
    else:
        nd = ""
    if dataType == "RRBS":
        rrbs = " --rrbs"
    else:
        rrbs = ""
    #### I will need to add distinction beteween paired end and single end data here
    script += """# Trim reads:
trim_galore{}{} --stringency 3 -q 30 --output_dir $outDir/trimData_BWA --fastqc --paired $outDir/rawData_BWA/$R1 $outDir/rawData_BWA/$R2
module unload trim_galore/0.4.0
# Trim Galore needs to be unloaded as together with it, python 2.7.3 is becoming the dafult option. As a result "toolshed" needed for bwa-meth is not loading, at least if it is installed via anaconda.
\n""".format(nd, rrbs)

    script += """# Map reads:
python $bwameth --reference $reference $outDir/trimData_BWA/$R1_trimmed $outDir/trimData_BWA/$R2_trimmed -t $threads | samtools view -b - > $outDir/mappingBWA_meth/$bam
\n
### Methylation extraction
samtools sort -@ $threads -o $outDir/mappingBWA_meth/$sortedBam $outDir/mappingBWA_meth/$bam
rm $outDir/mappingBWA_meth/$bam
samtools index -@ $threads $outDir/mappingBWA_meth/$sortedBam\n\n"""

    if dataType == "WGBS":
        script += """# Deduplication and sort:
$biscuit markdup -r $outDir/mappingBWA_meth/$sortedBam $outDir/mappingBWA_meth/$sortedDedupedBam
samtools sort -o $outDir/mappingBWA_meth/$sortedDedupedSortedBam $outDir/mappingBWA_meth/$sortedDedupedBam
samtools index $outDir/mappingBWA_meth/$sortedDedupedSortedBam
sortedBam=$sortedDedupedSortedBam\n\n"""

    script += """### Methylation extraction:
$MethylDackel extract -@ $threads --methylKit $reference $outDir/mappingBWA_meth/$sortedBam"""

    ############
    outfile = open("{}/scripts/{}.bwaPipeline.sh".format(outDir, sample), 'w')
    outfile.write(script)
    outfile.close()
    print("### prepareBWAmethPipeline NOTE: Pipeline created for '{}' sample".format(sample))

    if test == False:
        helixJobID = subprocess.Popen("qsub {}/scripts/{}.bwaPipeline.sh".format(outDir, sample), shell=True, stdout=subprocess.PIPE).stdout.read()
        print("### prepareBWAmethPipeline NOTE: Pipeline '{}' sample have been submited to helix with ID {}".format(sample, helixJobID))
        return helixJobID.decode('ascii').strip()
    else:
        return "Test.Run"


def prepareBWAmethPipeline_singleEnd(sample, outDir, R1_oryginalName, genomeBWA, referenceGenome, dataType, walltime, non_directional, threads):
    script = "#!/bin/bash\n"

    if walltime > 72:
        qType = "long"
    else:
        qType = "batch"
    if walltime < 10:
        walltime = "0{}".format(walltime)

    script += """#PBS -q {} # batch - up to 72h, long - up to 504h. 
#PBS -l nodes=1:ppn={}
#PBS -l walltime={}:00:00
#PBS -o {}/scripts
#PBS -e {}/scripts\n\n
### Load modules:\n
module load samtools/1.5
module load bwa/0.7.9a
module load cutadapt/1.2.1
module load trim_galore/0.4.0
module load fastqc/0.11.5\n\n
# create genome index. I assume it is already there, but if not, use commands like these to create them:
# module unload trim_galore/0.4.0
# cd /projects/li-lab/software/bwa-meth/bwa-meth-master/
# python bwameth.py index /projects/li-lab/reference/hg19.fa
# module load trim_galore/0.4.0\n\n\n""".format(qType, threads, walltime, outDir, outDir)

    if 'q.gz' in R1_oryginalName:
        ext = "fastq.gz"
        extTrim = "fq.gz"
    else:
        ext = "fastq"
        extTrim = "fq"
    prefix = sample.split(".")[0]
    rep = sample.split(".")[1]
    R1 = "{}.Read_R1.{}.{}".format(prefix, rep, ext)
    R1_trimmed = "{}.Read_R1.{}_trimmed.{}".format(prefix, rep, extTrim)
    bwaDir = "{}/mappingBWA_meth/".format(outDir)
    bwa_methylation_extractor_Dir = "{}/extractBWA_meth/".format(outDir)
    bam = "{}.{}.bam".format(prefix, rep)
    sam = "{}.{}.sam".format(prefix, rep)
    sortedBam = "{}.{}.sorted.bam".format(prefix, rep)
    if dataType == "WGBS":
        biscuit = '\nbiscuit="/projects/li-lab/software/biscuit/biscuit-release/biscuit"'
        sortedDedupedBam = '\nsortedDedupedBam="{}.{}.sorted.dedup.bam"'.format(prefix, rep)
        sortedDedupedSortedBam = '\nsortedDedupedSortedBam="{}.{}.sorted.dedup.sorted.bam"'.format(prefix, rep)
    else:
        biscuit = ""
        sortedDedupedBam = ""
        sortedDedupedSortedBam = ""
    script += """### Input parameters:
reference="{}"
R1_oryginalName="{}"
R1="{}" #new name
R1_trimmed="{}" #new name
sam="{}"
bam="{}"
sortedBam="{}"{}{}
outDir="{}"
threads={}

### These parameters are fixed, but may be modified if needed:
bwaPath="/projects/li-lab/software/bwa-meth/bwa-meth-master"
MethylDackel="/projects/li-lab/software/MethylDackel/MethylDackel"
bwameth="/projects/li-lab/software/bwa-meth/bwa-meth-master/bwameth.py"
{}
#####################\n\n""".format(genomeBWA, R1_oryginalName, R1, R1_trimmed, sam, bam, sortedBam, sortedDedupedBam, sortedDedupedSortedBam, outDir, threads, biscuit)

    script += """### Analysis:\n
cd $outDir\n
# Copy fastq files:
cp $R1_oryginalName $outDir/rawData_BWA/$R1\n\n"""

    if non_directional:
        nd = " --non_directional"
    else:
        nd = ""
    if dataType == "RRBS":
        rrbs = " --rrbs"
    else:
        rrbs = ""
    #### I will need to add distinction beteween paired end and single end data here
    script += """# Trim reads:
trim_galore{}{} --stringency 3 -q 30 --output_dir $outDir/trimData_BWA --fastqc $outDir/rawData_BWA/$R1
module unload trim_galore/0.4.0
# Trim Galore needs to be unloaded as together with it, python 2.7.3 is becoming the dafult option. As a result "toolshed" needed for bwa-meth is not loading, at least if it is installed via anaconda.
\n""".format(nd, rrbs)

    script += """# Map reads and sort:
python $bwameth --reference $reference $outDir/trimData_BWA/$R1_trimmed -t $threads > $outDir/mappingBWA_meth/$sam
samtools view -b $outDir/mappingBWA_meth/$sam > $outDir/mappingBWA_meth/$bam
rm $outDir/mappingBWA_meth/$sam
\n
samtools sort -@ $threads -o $outDir/mappingBWA_meth/$sortedBam $outDir/mappingBWA_meth/$bam
rm $outDir/mappingBWA_meth/$bam
samtools index -@ $threads $outDir/mappingBWA_meth/$sortedBam\n\n"""

    if dataType == "WGBS":
        script += """# Deduplication and sort:
$biscuit markdup -r $outDir/mappingBWA_meth/$sortedBam $outDir/mappingBWA_meth/$sortedDedupedBam
samtools sort -o $outDir/mappingBWA_meth/$sortedDedupedSortedBam $outDir/mappingBWA_meth/$sortedDedupedBam
samtools index $outDir/mappingBWA_meth/$sortedDedupedSortedBam
sortedBam=$sortedDedupedSortedBam\n\n"""

    script += """### Methylation extraction:
$MethylDackel extract -@ $threads --methylKit $reference $outDir/mappingBWA_meth/$sortedBam"""

    ############
    outfile = open("{}/scripts/{}.bwaPipeline.sh".format(outDir, sample), 'w')
    outfile.write(script)
    outfile.close()
    print("### prepareBWAmethPipeline NOTE: Pipeline created for '{}' sample".format(sample))

    if test == False:
        helixJobID = subprocess.Popen("qsub {}/scripts/{}.bwaPipeline.sh".format(outDir, sample), shell=True, stdout=subprocess.PIPE).stdout.read()
        print("### prepareBWAmethPipeline NOTE: Pipeline '{}' sample have been submited to helix with ID {}".format(sample, helixJobID))
        return helixJobID.decode('ascii').strip()
    else:
        return "Test.Run"


# ### Main:

WelcomeMessage = """
  ____   _____       _____            
 |  _ \ / ____|     / ____|           
 | |_) | (___ _____| (___   ___  __ _ 
 |  _ < \___ \______\___ \ / _ \/ _` |
 | |_) |____) |     ____) |  __/ (_| |
 |____/|_____/     |_____/ \___|\__, |
                                   | |
  PIPELINE v0.07                   |_|
                                   
 _____________________________________

Command used to run the pipeline: python {}
 _____________________________________
""".format(' '.join(str(x) for x in sys.argv))

print(WelcomeMessage)

parser = argparse.ArgumentParser(description="This is our super cool BS-Seq pipeline.")
## parameters for both pipelines (Bismark and BWA-meth):
parser = argparse.ArgumentParser(description='Parameters for both pipelines (Bismark and BWA-meth):')
parser.add_argument("-o", "--outDir", help="Full path to the output directory.", action="store", type=str, required=True, dest="outDir")
parser.add_argument("-c", "--configFile",
                    help="Full path to the configuration file. should have 4 columns: fastqc input file name, Prefix (like IDH_WT or IDH_Mut), Reads identity (either 'R1' or 'R2'; in case of paired end reads, 'R1' in everywhere in case of single ends), replicate id (eg. 1, 2 etc.). IMPORTANT NOTE! If the first line of the configuration file contains the header (which I highly recommend), it MUST be commented (i.e. begin with '#'). Also, you CANNOT use dots ('.') in the prefix or replicate IDs.",
                    action="store", type=str, required=True, dest="configFile")
parser.add_argument("-dt", "--dataType", help="Specify if the script should be optimized for RRBS or WGBS data analysis. If one selects 'RRBS' analysis type, two flags will be added to the mapping protocol: '--ignore_r2 1' and '--ignore_3prime_r2 2'. By default = 'WGBS'", default='WGBS', action="store", type=dataTypeParser, required=False, dest="dataType")
parser.add_argument("-p", "--pipeline", help="Which pipeline should be used (specify the option number): [1] Bismark only; [2] BWA-meth only; [3] Bismark and BWA-meth. By default = 3.", default=3, action="store", type=int, required=False, dest="pipeline")
parser.add_argument("-t", "--threads", help="Hom many processors do you want to use. By default = 10.", default=10, action="store", type=int, required=False, dest="threads")
parser.add_argument("-w", "--walltime", help="Hom much walltime would you like to reseve for your analysis (in hours). By default = 24.", default=24, action="store", type=int, required=False, dest="walltime")
parser.add_argument("-g", "--referenceGenome", help="Short name of the reference genome version (e.g. 'hg38')- will be a part of the bam file name to easily keep track. By default = ''.", action="store", type=str, required=False, dest="referenceGenome")

## Bismark specific parameters:
parser.add_argument("-gBis", "--genomeBismark", help="Full path to the reference genome for Bismark (e.g. '/projects/li-lab/reference/hg19/').", action="store", type=str, required=False, dest="genomeBismark")
parser.add_argument("--non_directional", help="Use this flag to activate 'Non Directional' mode in Bismark. By default Bismark is using directional mode.", action="store_true", dest="non_directional")
parser.add_argument("--zero_based", help="Use this flag to activate 'Zero Based' mode in Bismark. By default Bismark is using one based mode.", action="store_true", dest="zero_based")
parser.add_argument("--test", help="This flag will activate the test mode, which will create all output subcatalogs structure and scripts, but will not submit anything to helix server.", action="store_true", dest="test")

## BWA-meth specific parameters:
parser.add_argument("-gBwa", "--genomeBWA", help="Full path to the reference genome for BWA-meth (e.g. '/projects/li-lab/reference/hg19/hg19.fa').", action="store", type=str, required=False, dest="genomeBWA")

##########################
args = parser.parse_args()

outDir = args.outDir
configFile = args.configFile
dataType = args.dataType
pipeline = args.pipeline
threads = args.threads
walltime = args.walltime
referenceGenome = args.referenceGenome
genomeBismark = args.genomeBismark
non_directional = args.non_directional
zero_based = args.zero_based
test = args.test
genomeBWA = args.genomeBWA

# outDir = "/fastscratch/rosikw/bsseq"
# configFile = "/home/rosikw/programs/pipelines/BS_Seq/exampleConfigFile.tsv"
# dataType = "RRBS"
# pipeline = 3
# threads = 10
# walltime = 24
# referenceGenome = 'hg19'
# genomeBismark = "/projects/li-lab/reference/hg19/"
# non_directional = False
# zero_based = False
# test = True

pipelinePath = get_script_path()  # "/home/rosikw/programs/pipelines/BS_Seq/" # in the future I will have to automate this step

print("Parameters used in the analysis:")
print("pipelinePath:", pipelinePath)
print("outDir:", outDir)
print("configFile:", configFile)
print("dataType:", dataType)
print("pipeline:", pipeline)
print("threads:", threads)
print("walltime:", walltime)
print("referenceGenome:", referenceGenome)
print("genomeBismark:", genomeBismark)
print("non_directional:", non_directional)
print("zero_based:", zero_based)
print("test:", test)
print("genomeBWA:", genomeBWA)
print(" _____________________________________")

exitStatus = False
if (pipeline == 2 and genomeBWA == None) or (pipeline == 3 and genomeBWA == None):
    print("BS-seq pipeline ERROR: the reference genome for BWA-Meth was not provided")
    exitStatus = True
if (pipeline == 2 and genomeBismark == None) or (pipeline == 3 and genomeBismark == None):
    print("BS-seq pipeline ERROR: the reference genome for Bismark was not provided")
    exitStatus = True
if exitStatus == True:
    print("BS-seq pipeline ExitStatus=True - something went wrong, please trace back the error messages")
    exit()

# STEP 0: Provide full path to the pipeline scripts
# pipelinePath = "/home/rosikw/programs/pipelines/BS_Seq/" # in the future I will have to automate this step

# STEP 1: Prepare the output directory:
prepareOutputsStructure(outDir, pipeline)
print(" _____________________________________")

# STEP 2: Parse the configuration file:
config = parseConfigFile(configFile)
print(" _____________________________________")

if pipeline == 1 or pipeline == 3:
    # proceed with the analysis for Bismark:

    # STEP 3: Optimize pipeline per sample:
    dependentJobs = []
    for sample in config:
        helixJobID = runBismarkSingleOrPaired(sample, outDir, config[sample][0], config[sample][1], genomeBismark, referenceGenome, dataType, non_directional, zero_based, walltime, threads)
        #     tmp = helixJobID
        dependentJobs.append(helixJobID)
    #
    # 	# STEP 4: run QC dependent on previous jobs:
    # 	if test == False:
    # 		QualityControl_PostBismark(dependentJobs, outDir, pipelinePath)
    print(" _____________________________________")

if pipeline == 2 or pipeline == 3:

    # STEP 3: Optimize pipeline per sample:
    dependentJobs = []
    for sample in config:
        helixJobID = runBWASingleOrPaired(sample, outDir, config[sample][0], config[sample][1], genomeBWA, referenceGenome, dataType, walltime, non_directional, threads)
        #     tmp = helixJobID
        dependentJobs.append(helixJobID)
    print(" _____________________________________")
# 
# 	# STEP 4: run QC dependent on previous jobs:
# 	if test == False:
# 		QualityControl_PostBismark(dependentJobs, outDir, subfolder, pipelinePath)
# 	# proceed with the analysis for BWA-meth:


###### Dependencies:
# "toolshed" - python module, that needs to be installed for example using conda: "conda install -c bioconda toolshed" or with pip "pip install toolshed" (i was only testing with the former option). It is required for BWA-Meth mapping.
