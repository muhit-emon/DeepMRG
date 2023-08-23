# DeepMRG
A multi-label classifier for predicting bacterial metal resistance genes using deep learning 
# Installation
<div style="background-color: silver; padding: 10px;">
  <b>git clone https://<i></i>github.com/muhit-emon/DeepMRG.git </b><br>
  <b>cd DeepMRG </b><br>
  <b>bash install.sh </b><br>
</div>
# Conda environment activation
After installation of DeepMRG, a conda environment named <b>deepmrg</b> will be created<br>
To activate the environment, run the following command <br>
<b>conda activate deepmrg</b>
# Usage (demo)
Go inside DeepMRG directory, then execute the following command to run DeepMRG on a protein fasta file <br>
<b>nextflow run protein_pipeline.nf --prot /absolute/path/to/protein/fasta/file --out_prefix demo </b>
