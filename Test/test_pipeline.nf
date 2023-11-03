#!/usr/bin/env nextflow
 
// Enable DSL 2 syntax
nextflow.enable.dsl = 2

params.TEST_fasta = "${projectDir}/TEST.fasta" // TEST.fasta is the test dataset that was used to evaluate DeepMRG in the paper
params.LOW_fasta = "${projectDir}/BacMet_Predicted_MRG_DB_partition_2.fasta" // BacMet Predicted MRG DB partition 2 used to evaluate DeepMRG in the paper

params.fold1_fasta = "${projectDir}/Fold1_Validation.fasta" // validation dataset under 5-fold cross-validation (Fold 1)
params.fold2_fasta = "${projectDir}/Fold2_Validation.fasta" // validation dataset under 5-fold cross-validation (Fold 2)
params.fold3_fasta = "${projectDir}/Fold3_Validation.fasta" // validation dataset under 5-fold cross-validation (Fold 3)
params.fold4_fasta = "${projectDir}/Fold4_Validation.fasta" // validation dataset under 5-fold cross-validation (Fold 4)
params.fold5_fasta = "${projectDir}/Fold5_Validation.fasta" // validation dataset under 5-fold cross-validation (Fold 5)

params.test_set = "TEST" // TEST has been used as default. It will be overridden when given as --test_set

params.path_of_bacmet_exp_db_clusters_json = "${projectDir}/../clusters_of_ref_exp_mrg.json"
params.path_of_deepmrg_model1 = "${projectDir}/../models/model1.h5"
params.path_of_deepmrg_model2 = "${projectDir}/../models/model2.h5"
params.path_of_deepmrg_model3 = "${projectDir}/../models/model3.h5"
params.path_of_deepmrg_model4 = "${projectDir}/../models/model4.h5"
params.path_of_deepmrg_model5 = "${projectDir}/../models/model5.h5"

process run_diamond_with_ref_exp_db {

	input:
	path prot_fa

	output:
	path "diamond_with_ref_exp.tsv"

	"""
	bash $projectDir/../run_diamond_against_ref_exp.sh ${prot_fa} ${projectDir}/..
	"""
}

process compute_bitscore_distribution_with_bacmet_exp_db {

	input:
	path x

	output:
	path "bitscore_distribution_with_ref_exp.json"

	"""
	python3 ${projectDir}/../compute_bitscore_distribution.py $x
	"""
}


process run_deepmrg {

	publishDir "${projectDir}", mode: "copy"

	input:
	path bitscore_distribution
	
	output:
	path "*_classification_report.txt"
	

	"""
	python3 ${projectDir}/DeepMRG_Test.py ${params.test_set} ${params.test_seq} $bitscore_distribution ${params.path_of_bacmet_exp_db_clusters_json} ${params.path_of_deepmrg_model1} ${params.path_of_deepmrg_model2} ${params.path_of_deepmrg_model3} ${params.path_of_deepmrg_model4} ${params.path_of_deepmrg_model5}
	"""

}

workflow {
	
	if ( params.test_set == "TEST" ){
		input_ch = Channel.from(params.TEST_fasta) // this channel contains input seq fasta file
		params.test_seq = params.TEST_fasta
	}
	else if ( params.test_set == "LOW" ){
		input_ch = Channel.from(params.LOW_fasta) // this channel contains input seq fasta file
		params.test_seq = params.LOW_fasta
	}
	else if ( params.test_set == "fold1" ){
		input_ch = Channel.from(params.fold1_fasta) // this channel contains input seq fasta file
		params.test_seq = params.fold1_fasta
	}
	else if ( params.test_set == "fold2" ){
		input_ch = Channel.from(params.fold2_fasta) // this channel contains input seq fasta file
		params.test_seq = params.fold2_fasta
	}
	else if ( params.test_set == "fold3" ){
		input_ch = Channel.from(params.fold3_fasta) // this channel contains input seq fasta file
		params.test_seq = params.fold3_fasta
	}
	else if ( params.test_set == "fold4" ){
		input_ch = Channel.from(params.fold4_fasta) // this channel contains input seq fasta file
		params.test_seq = params.fold4_fasta
	}
	else if ( params.test_set == "fold5" ){
		input_ch = Channel.from(params.fold5_fasta) // this channel contains input seq fasta file
		params.test_seq = params.fold5_fasta
	}
		
	diamond_against_bacmet_exp_db_ch = run_diamond_with_ref_exp_db(input_ch)

	bitscore_distribution_with_ref_exp_ch = compute_bitscore_distribution_with_bacmet_exp_db(diamond_against_bacmet_exp_db_ch)
	
	deepmrg_annotation_ch = run_deepmrg(bitscore_distribution_with_ref_exp_ch)
	
}
