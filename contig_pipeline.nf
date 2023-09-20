#!/usr/bin/env nextflow
 
// Enable DSL 2 syntax
nextflow.enable.dsl = 2

params.contig = "36.fna" // 36 has been used just to initialize. It will be overridden when given as --contig
params.out_prefix = "output"
params.path_of_bacmet_exp_db_clusters_json = "${projectDir}/clusters_of_ref_exp_mrg.json" // change it while deploying to other machine
params.path_of_deepmrg_model1 = "${projectDir}/models/model1.h5"
params.path_of_deepmrg_model2 = "${projectDir}/models/model2.h5"
params.path_of_deepmrg_model3 = "${projectDir}/models/model3.h5"
params.path_of_deepmrg_model4 = "${projectDir}/models/model4.h5"
params.path_of_deepmrg_model5 = "${projectDir}/models/model5.h5"


process  make_single_line_seq_and_split {
	
	input:
	path contig_fna
	
	output:
	path 'x*.fna'
	
	"""
	g++ -o CPP.out $projectDir/make_single_line_seq_and_split_new_for_nf.cpp && ./CPP.out ${contig_fna}
	"""
}


process run_prodigal {

	input:
	path x
	
	output:
	path 'x*.faa'
	
	"""
	bash $projectDir/prodigal.sh '$x' ${projectDir}
	"""
}

process merge_predicted_aa_files {
	
	input:
	path x
	
	output:
	path "predicted_proteins.faa"
	
	"""
	cat "$x" > predicted_proteins.faa
	"""
}

process get_the_predicted_aa_file {
	
	publishDir "${projectDir}", mode: "copy"
	
	input:
	path x
	
	output:
	path "${params.out_prefix}_predicted_proteins.faa"
	
	"""
	cat $x >> ${params.out_prefix}_predicted_proteins.faa
	"""
}

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
	path prot_seq
	path bitscore_distribution

	output:
	path "${params.out_prefix}_DeepMRG_annotation.tsv"

	"""
	python3 ${projectDir}/deepmrg.py $prot_seq $bitscore_distribution ${params.path_of_bacmet_exp_db_clusters_json} ${params.path_of_deepmrg_model1} ${params.path_of_deepmrg_model2} ${params.path_of_deepmrg_model3} ${params.path_of_deepmrg_model4} ${params.path_of_deepmrg_model5} ${params.out_prefix}
	"""

}

workflow {

	input_ch = Channel.from(params.contig) // this channel contains input contig (nucleic acid sequences) fasta file
	
	splitted_contig_files_ch = make_single_line_seq_and_split(input_ch)
	
	predicted_aa_files_ch = run_prodigal(splitted_contig_files_ch.flatten())
	
	final_predicted_aa_file_ch = merge_predicted_aa_files(predicted_aa_files_ch) | collectFile
	
	protein_file_ch = get_the_predicted_aa_file(final_predicted_aa_file_ch)

	diamond_with_ref_exp_ch = run_diamond_with_ref_exp_db(protein_file_ch)

	bitscore_distribution_with_ref_exp_ch = compute_bitscore_distribution_with_bacmet_exp_db(diamond_with_ref_exp_ch)
	
	deepmrg_annotation_ch = run_deepmrg(protein_file_ch, bitscore_distribution_with_ref_exp_ch)
	
}
