params.list = 'NA19240_R_9_4_1.flist.txt'

//Channel.fromPath(params.list)
//	.splitText()
//	.map { file(it) }
//	.set { file_list }
//
////file_list.view()


Channel.fromPath( params.list )
    .splitCsv(header: false)
    .map { file(it[0]) }
    .toList()
    .set{ fl2 }

fl2.flatten()

process test{
	input:
	file x from fl2.flatten()

	"""
	ls -lh $x
	"""
}