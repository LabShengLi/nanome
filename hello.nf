params.list = 'NA19240_R_9_4_1.flist.txt'

Channel.fromPath(params.list)
	.splitText()
	.map { file(it) }
	.set { file_list }

file_list.view()


