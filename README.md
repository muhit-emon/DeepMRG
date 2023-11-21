# DeepMRG
DeepMRG is a multi-label deep learning classifier for predicting bacterial metal resistance genes (MRGs). It can be used to predict MRGs from protein sequences provided in fasta file. It also can be applied on metagenomic or isolate assembled contigs (in fasta) to predict MRGs. <br><br>
This Github is the local version of DeepMRG. The web server of DeepMRG is available via <a href="https://deepmrg.cs.vt.edu/deepmrg">server version</a>. 

Preprint: <a href="https://doi.org/10.1101/2023.11.14.566903"> https://doi.org/10.1101/2023.11.14.566903 </a>
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
# conda environment activation
After installation of DeepMRG, a conda environment named <b>deepmrg</b> will be created.<br>
To activate the environment, run the following command <br>
<pre>
conda activate deepmrg
</pre>
# (1) Usage on protein sequences
Go inside DeepMRG directory. <br> <br>
<b>To run DeepMRG on protein sequences (must be in fasta format) to predict MRGs, use the following command</b> <br>
<pre>
nextflow run protein_pipeline.nf --prot &ltabsolute/path/to/protein/fasta/file&gt --out_prefix &ltprefix of output file name&gt
rm -r work
</pre>
The command line options for this script (<b>protein_pipeline.nf</b>) are: <br><br>
<b>--prot</b>: The absolute path of the fasta file containing protein sequences to be classified <br>
<b>--out_prefix</b>: The prefix of the output file name <br><br>
With <b>--out_prefix demo</b>, An output tsv file named <b>demo_DeepMRG_annotation.tsv</b> (contains MRG predictions by DeepMRG) will be generated inside DeepMRG directory. <br><br>
Replace <b>absolute/path/to/protein/fasta/file</b> with your protein fasta file absolute path and the output prefix <b>demo</b> with your own output prefix.
# (2) Usage on metagenomic or isolate assemblies (DNA sequences of assembled contigs)
Go inside DeepMRG directory. <br> <br>
<b>To run DeepMRG on metagenomic or isolate assembled contigs (must be in fasta format) to predict MRGs, use the following command</b> <br>
<pre>
nextflow run contig_pipeline.nf --contig &ltabsolute/path/to/contig/fasta/file&gt --out_prefix &ltprefix of output file name&gt
rm -r work
</pre>
The command line options for this script (<b>contig_pipeline.nf</b>) are: <br><br>
<b>--contig</b>: The absolute path of the fasta file containing contigs <br>
<b>--out_prefix</b>: The prefix of the output file name <br><br>
With <b>--out_prefix demo</b>, the following files will be generated inside DeepMRG directory.

<ul>
  <li><b>demo_DeepMRG_annotation.tsv</b> (contains MRG predictions by DeepMRG)</li>
  <li><b>demo_predicted_proteins.faa</b> (contains prodigal predicted proteins from contigs)</li>
</ul>

Replace <b>absolute/path/to/contig/fasta/file</b> with your contig fasta file absolute path and the output prefix <b>demo</b> with your own output prefix. <br><br>
Pipeline for predicting bacterial MRGs from assembled contigs using DeepMRG is shown below:

![Fig6](https://drive.google.com/uc?export=view&id=1CKf_NAzzfiuuRonf-eB6bWV7DTgCVQBS)

# Output
<b>&lt;prefix of output file name&gt;_DeepMRG_annotation.tsv</b> is the main output file that contains MRG prediction. The output file is a tab separated file with each line containing a protein sequence header and the corresponding MRG predictions. The sequences are in the same order as in the input fasta file. <br><br>

![demo_output](https://drive.google.com/uc?export=view&id=1xJeWZSlvlaNwTMnz609rAjJyDBtCKNgP)

The output file contains 2 columns:<br><br>
The 1st  column (Protein_ID) contains the header of the protein sequences in the input fasta file.<br><br>
The 2nd column (Prediction) contains the prediction results of DeepMRG. By default, the proteins with the prediction score less than 3.5 (out of 5) are considered as non-MRG.<br><br>

For example, protein 3 has been predicted by DeepMRG to confer resistances to both Cu and Zn. On the other hand, protein 5 has been predicted as a non-MRG as its prediction score falls below 3.5.
