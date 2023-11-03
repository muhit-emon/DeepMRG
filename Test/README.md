# Testing of DeepMRG
We've included the datasets (TEST.fasta, BacMet_Predicted_MRG_DB_partition_2.fasta, Fold1_Validation.fasta, Fold2_Validation.fasta, Fold3_Validation.fasta, Fold4_Validation.fasta, and Fold5_Validation.fasta) used for DeepMRG evaluation. To test DeepMRG and reproduce the results with these datasets, please follow the instructions below.<br><br>
At first, go inside the Test directory
<pre>
cd DeepMRG/Test
</pre>
# Evaluation on the test dataset (TEST.fasta)
Execute the following command:
<pre>
nextflow run test_pipeline.nf --test_set TEST
</pre>
As output, a txt file named <b>TEST_classification_report.txt</b> will be generated inside the Test directory which contains the classification results of DeepMRG on TEST.fasta.
# Evaluation on BacMet Predicted MRG DB partition 2 (BacMet_Predicted_MRG_DB_partition_2.fasta)
Execute the following command:
<pre>
nextflow run test_pipeline.nf --test_set LOW
</pre>
As output, a txt file named <b>LOW_classification_report.txt</b> will be generated inside the Test directory which contains the classification results of DeepMRG on BacMet_Predicted_MRG_DB_partition_2.fasta.
# Evaluation through 5-fold cross-validation
To evaluate DeepMRG on the 1st validation fold, execute the following command:
<pre>
nextflow run test_pipeline.nf --test_set fold1
</pre>
As output, a txt file named <b>fold1_classification_report.txt</b> will be generated inside the Test directory which contains the classification results of DeepMRG on Fold1_Validation.fasta.

To evaluate DeepMRG on the 2nd validation fold, execute the following command:
<pre>
nextflow run test_pipeline.nf --test_set fold2
</pre>
As output, a txt file named <b>fold2_classification_report.txt</b> will be generated inside the Test directory which contains the classification results of DeepMRG on Fold2_Validation.fasta.

To evaluate DeepMRG on the 3rd validation fold, execute the following command:
<pre>
nextflow run test_pipeline.nf --test_set fold3
</pre>
As output, a txt file named <b>fold3_classification_report.txt</b> will be generated inside the Test directory which contains the classification results of DeepMRG on Fold3_Validation.fasta.

To evaluate DeepMRG on the 4th validation fold, execute the following command:
<pre>
nextflow run test_pipeline.nf --test_set fold4
</pre>
As output, a txt file named <b>fold4_classification_report.txt</b> will be generated inside the Test directory which contains the classification results of DeepMRG on Fold4_Validation.fasta.

To evaluate DeepMRG on the 5th validation fold, execute the following command:
<pre>
nextflow run test_pipeline.nf --test_set fold5
</pre>
As output, a txt file named <b>fold5_classification_report.txt</b> will be generated inside the Test directory which contains the classification results of DeepMRG on Fold5_Validation.fasta.
