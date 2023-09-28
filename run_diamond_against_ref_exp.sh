#!/bin/bash

#--id                   minimum identity% to report an alignment
#--query-cover          minimum query cover% to report an alignment
#--subject-cover        minimum subject cover% to report an alignment
#-f 6 qtitle stitle pident bitscore evalue qlen slen qstart qend sstart send

# $1 = query_file;
ref_db=$2/ref_exp

chmod +x $2/diamond
$2/diamond blastp -q "$1" -d ${ref_db} -o diamond_with_ref_exp.tsv --id 20 --evalue 1e-7 -k 485 --query-cover 60 --subject-cover 60 -f 6 qtitle stitle pident bitscore evalue qlen slen qstart qend sstart send --very-sensitive --quiet
