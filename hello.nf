params.list = 'NA19240_R_9_4_1.flist.txt'

Channel.fromPath( params.list )
    .splitCsv(header: false)
    .map { file(it[0]) }
    .toList()
    .set{ online_filelist }

//online_filelist.flatten()

process test{
	input:
	file x from online_filelist.flatten()

	"""
	ls -lh $x
	"""
}