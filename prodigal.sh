#!/bin/bash

fname_with_ext=$(basename "$1")
only_filename="${fname_with_ext%.*}"
$2/prodigal -a ${only_filename}.faa -p meta -i $1 -q
