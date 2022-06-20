/*
=========================================================================================
  		NANOME(Nanopore methylation) pipeline for Oxford Nanopore sequencing
=========================================================================================
 NANOME Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/LabShengLi/nanome
 @Author   : Yang Liu
 @FileName : COMMONS.nf
 @Software : NANOME project
 @Organization : JAX Sheng Li Lab
----------------------------------------------------------------------------------------
*/
// check nextflow version, then declare DSL2 in two ways
def nextflowVersionCheck() {
	// We now support both latest and lower versions, due to Lifebit CloudOS is only support 20.04
	// Note: NXF_VER=20.04.1 nextflow run main.nf -profile test,singularity
	if( nextflow.version.matches(">= 20.07.1") ){
		nextflow.enable.dsl = 2
	} else {
		// Support lower version of nextflow
		nextflow.preview.dsl = 2
	}
}
