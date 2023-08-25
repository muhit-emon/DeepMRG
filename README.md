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
# Usage (1) on protein sequences
Go inside DeepMRG directory. <br> <br>
<b>If you want to run DeepMRG on protein sequences (must be in fasta format) to predict bacterial metal resistance genes, use the following command</b> <br>
<pre>
nextflow run protein_pipeline.nf --prot /path/to/protein/fasta/file --out_prefix TEST
rm -r work
</pre>
An output tsv file named <b>TEST_DeepMRG_annotation.tsv</b> will be generated inside DeepMRG directory. <br><br>
You need to replace <b>/path/to/protein/fasta/file</b> with your protein fasta file path and the output prefix <b>TEST</b> with your own output prefix.
# Usage (2) on (meta)genomic contigs (nucleic acid sequences)
Go inside DeepMRG directory. <br> <br>
<b>If you want to run DeepMRG on (meta)genomic assembled contigs (must be in fasta format) to predict bacterial metal resistance genes, use the following command</b> <br>
<pre>
nextflow run contig_pipeline.nf --contig /path/to/contig/fasta/file --out_prefix TEST
rm -r work
</pre>
An output tsv file named <b>TEST_DeepMRG_annotation.tsv</b> will be generated inside DeepMRG directory. <br><br>
You need to replace <b>/path/to/contig/fasta/file</b> with your contig fasta file path and the output prefix <b>TEST</b> with your own output prefix.
