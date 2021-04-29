//params.list = 'NA19240_R_9_4_1.flist.txt'
//
//Channel.fromPath( params.list )
//    .splitCsv(header: false)
//    .map { file(it[0]) }
//    .toList()
//    .set{ online_filelist }
//
////online_filelist.flatten()
//
//process test{
//	input:
//	file x from online_filelist.flatten()
//
//	"""
//	ls -lh $x
//	"""
//}


log.info """\
	=================================
	dsname              :${params.dsname}
	input               :${params.input}
	benchmarking        :${params.benchmarking}
	eval                :${params.eval}
	debug               :${params.debug}
	online              :${params.online}
	reference_genome    :${params.referenceGenome}
	chromSizesFile      :${params.chromSizesFile}
	=================================
	"""
	.stripIndent()

// Check all tools work well on the platform
process EnvCheck {

	tag 'EnvCheck'
	// teminate all later processes if this process is not passed
	errorStrategy 'terminate'
	label 'with_gpus'

	when:
	! params.online

	"""
	set -x

	which tombo
    tombo -v

    which nanopolish
    nanopolish --version

    which megalodon
    megalodon -v

    which deepsignal
    deepsignal

    which DeepMod.py
    DeepMod.py

	which guppy_basecaller
    guppy_basecaller -v
    """
}



