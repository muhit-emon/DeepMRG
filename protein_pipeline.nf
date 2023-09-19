#!/usr/bin/env nextflow
 
// Enable DSL 2 syntax
nextflow.enable.dsl = 2

params.prot = "36.faa" // 36.faa has been used just to initialize. It will be overridden when given as --prot
params.out_prefix = "output" // output has been used just to initialize. It will be overridden when given as --out_prefix
params.path_of_bacmet_exp_db_clusters_json = "${projectDir}/clusters_of_ref_exp_mrg.json"
params.path_of_deepmrg_model1 = "${projectDir}/models/model1.h5"
params.path_of_deepmrg_model2 = "${projectDir}/models/model2.h5"
params.path_of_deepmrg_model3 = "${projectDir}/models/model3.h5"
params.path_of_deepmrg_model4 = "${projectDir}/models/model4.h5"
params.path_of_deepmrg_model5 = "${projectDir}/models/model5.h5"

process run_diamond_with_ref_exp_db {

	input:
	path prot_fa

	output:
	path "diamond_with_ref_exp.tsv"

	"""
	bash $projectDir/run_diamond_against_ref_exp.sh ${prot_fa} ${projectDir}
	"""
}

process compute_bitscore_distribution_with_bacmet_exp_db {

	input:
	path x

	output:
	path "bitscore_distribution_with_ref_exp.json"

	"""
	python3 ${projectDir}/compute_bitscore_distribution.py $x
	"""
}


process run_deepmrg {

	publishDir "${projectDir}", mode: "copy"

	input:
	path bitscore_distribution

	output:
	path "${params.out_prefix}_DeepMRG_annotation.tsv"

	"""
	python3 ${projectDir}/deepmrg.py ${params.prot} $bitscore_distribution ${params.path_of_bacmet_exp_db_clusters_json} ${params.path_of_deepmrg_model1} ${params.path_of_deepmrg_model2} ${params.path_of_deepmrg_model3} ${params.path_of_deepmrg_model4} ${params.path_of_deepmrg_model5} ${params.out_prefix}
	"""

}

workflow {

	input_ch = Channel.from(params.prot) // this channel contains input seq fasta file
	
	diamond_against_bacmet_exp_db_ch = run_diamond_with_ref_exp_db(input_ch)

	bitscore_distribution_with_ref_exp_ch = compute_bitscore_distribution_with_bacmet_exp_db(diamond_against_bacmet_exp_db_ch)
	
	deepmrg_annotation_ch = run_deepmrg(bitscore_distribution_with_ref_exp_ch)
	
}
