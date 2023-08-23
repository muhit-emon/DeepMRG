# DeepMRG
A multi-label classifier for predicting bacterial metal resistance genes using deep learning 
# Requirements
<ol>
  <li>Linux operating system</li>
  <li>conda</li>
</ol>

# Installation
<pre>
git clone https://<i></i>github.com/muhit-emon/DeepMRG.git
cd DeepMRG
bash install.sh
</pre>
# Conda environment activation
After installation of DeepMRG, a conda environment named <b>deepmrg</b> will be created.<br>
To activate the environment, run the following command <br>
<pre>
conda activate deepmrg
</pre>
# Usage (demo)
Go inside DeepMRG directory, then execute the following command to run DeepMRG on a protein fasta file <br>
<b>nextflow run protein_pipeline.nf --prot /absolute/path/to/protein/fasta/file --out_prefix demo </b>
