params.input = 'test/benchmark.filelist.txt'

if (params.input.endsWith("filelist.txt")) { // filelist
	println("Detect filelist:");
	Channel.fromPath( params.input )
	    .splitCsv(header: false)
	    .map { file(it[0]) }
	    .toList()
	    .set{ online_filelist }
} else { // single file
	println("Detect one file:");
	Channel.fromPath( params.input ).set{online_filelist}
}

online_filelist.flatten().into{online_filelist1; online_filelist2}
online_filelist2.view()


process test{
	input:
	file x from online_filelist1.flatten()

	"""
	ls -lh $x
	"""
}


