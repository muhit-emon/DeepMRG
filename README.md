# DeepMRG
DeepMRG is a multi-label deep learning classifier for predicting bacterial metal resistance genes (MRGs). It can be used to predict MRGs from protein sequences provided in fasta file. It also can be applied on (meta)genomic assembled contigs (in fasta) to predict MRGs.
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
With <b>--out_prefix TEST</b>, An output tsv file named <b>TEST_DeepMRG_annotation.tsv</b> (contains MRG predictions by DeepMRG) will be generated inside DeepMRG directory. <br><br>
Replace <b>/path/to/protein/fasta/file</b> with your protein fasta file path and the output prefix <b>TEST</b> with your own output prefix.
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
With <b>--out_prefix TEST</b>, the following files will be generated inside DeepMRG directory.

<ul>
  <li><b>TEST_DeepMRG_annotation.tsv</b> (contains MRG predictions by DeepMRG)</li>
  <li><b>TEST_predicted_proteins.faa</b> (contains prodigal predicted proteins from contigs)</li>
</ul>

Replace <b>/path/to/contig/fasta/file</b> with your contig fasta file path and the output prefix <b>TEST</b> with your own output prefix. <br><br>
Pipeline for predicting bacterial MRGs from (meta)genomic assembled contigs using DeepMRG is shown below:

![Fig6](https://drive.google.com/uc?export=view&id=1Nph1cXD6rJN0VSrwdKKpTVfUisx0rB6H)
