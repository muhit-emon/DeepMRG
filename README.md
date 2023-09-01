# DeepMRG
DeepMRG is a multi-label deep learning classifier for predicting bacterial metal resistance genes (MRGs). It can be used to predict MRGs from protein sequences provided in fasta file. It also can be applied on (meta)genomic assembled contigs (in fasta) to predict MRGs.
The web server of DeepMRG is available via https://deepmrg.cs.vt.edu/deepmrg.
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
nextflow run protein_pipeline.nf --prot &lt/path/to/protein/fasta/file&gt --out_prefix &ltprefix of output file name&gt
rm -r work
</pre>
The command line options for this script (<b>protein_pipeline.nf</b>) are: <br><br>
<b>--prot</b>: The fasta file containing protein sequences to be classified <br>
<b>--out_prefix</b>: The prefix of the output file name <br><br>
With <b>--out_prefix demo</b>, An output tsv file named <b>demo_DeepMRG_annotation.tsv</b> (contains MRG predictions by DeepMRG) will be generated inside DeepMRG directory. <br><br>
Replace <b>/path/to/protein/fasta/file</b> with your protein fasta file path and the output prefix <b>demo</b> with your own output prefix.
# (2) Usage on (meta)genomic contigs (DNA sequences)
Go inside DeepMRG directory. <br> <br>
<b>To run DeepMRG on (meta)genomic assembled contigs (must be in fasta format) to predict MRGs, use the following command</b> <br>
<pre>
nextflow run contig_pipeline.nf --contig &lt/path/to/contig/fasta/file&gt --out_prefix &ltprefix of output file name&gt
rm -r work
</pre>
The command line options for this script (<b>contig_pipeline.nf</b>) are: <br><br>
<b>--contig</b>: The fasta file containing contigs <br>
<b>--out_prefix</b>: The prefix of the output file name <br><br>
With <b>--out_prefix demo</b>, the following files will be generated inside DeepMRG directory.

<ul>
  <li><b>demo_DeepMRG_annotation.tsv</b> (contains MRG predictions by DeepMRG)</li>
  <li><b>demo_predicted_proteins.faa</b> (contains prodigal predicted proteins from contigs)</li>
</ul>

Replace <b>/path/to/contig/fasta/file</b> with your contig fasta file path and the output prefix <b>demo</b> with your own output prefix. <br><br>
Pipeline for predicting bacterial MRGs from (meta)genomic assembled contigs using DeepMRG is shown below:

![Fig6](https://drive.google.com/uc?export=view&id=1Nph1cXD6rJN0VSrwdKKpTVfUisx0rB6H)

# Output
<b>&lt;prefix of output file name&gt;_DeepMRG_annotation.tsv</b> is the main output file that contains MRG prediction. The output file is a tab separated file with each line containing a protein sequence header and the corresponding MRG predictions. The sequences are in the same order as in the input fasta file. <br><br>

![demo_output](https://drive.google.com/uc?export=view&id=1-pw5s0s6-eZwOe8OZ-woVWkd6kLqmNMM)

The output file contains 2 columns:<br><br>
The 1st  column (Protein_ID) contains the header of the protein sequences in the input fasta file.<br><br>
The 2nd column (Prediction(probability %)) contains the prediction results with corresponding prediction probabilities calculated by DeepMRG. By default, the proteins with the prediction probability less than 70% are regarded as non-MRG.<br><br>

For example, protein 3 has been predicted to confer resistances to both Cu (with prediction probability of 99.84%) and Zn (with prediction probability of 86.56%). On the other hand, protein 5 has been predicted as a non-MRG as its prediction probability falls below 70%.
