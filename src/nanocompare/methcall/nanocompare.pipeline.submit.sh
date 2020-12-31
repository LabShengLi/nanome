#!/bin/bash

###################################################################################
### Run pipeline for each tool in ToolList                                      ###
###################################################################################

outbasedir=/fastscratch/liuya/nanocompare
mkdir -p ${outbasedir}

for Tool in ${ToolList[@]}; do
	### Building output folder configuration, such as
	#	/fastscratch/liuya/nanocompare/K562-Tombo-N50/
	#	├── K562-Tombo-N50-basecall
	#	├── K562-Tombo-N50-meth-call
	#	├── K562-Tombo-N50-sept
	#	└── K562-Tombo-N50-untar
	analysisPrefix=${dsname}-${Tool}-N${targetNum}
	untaredInputDir=${outbasedir}/${analysisPrefix}/${analysisPrefix}-untar
	septInputDir=${outbasedir}/${analysisPrefix}/${analysisPrefix}-sept
	basecallOutputDir=${outbasedir}/${analysisPrefix}/${analysisPrefix}-basecall
	methCallsDir=${outbasedir}/${analysisPrefix}/${analysisPrefix}-meth-call

	### Start script for submiting jobs
	echo "Start pipeline submit for tool with analysisPrefix=${analysisPrefix}"
	echo "Basedir=${outbasedir}/${analysisPrefix}"
	source /projects/li-lab/yang/workspace/nano-compare/src/nanocompare/methcall/nanocompare.pipeline.eachtool.submit.sh
done
###################################################################################
###################################################################################
###################################################################################