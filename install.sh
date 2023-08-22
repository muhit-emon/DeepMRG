#!/bin/bash

# diamond installation
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
chmod +x diamond
rm diamond-linux64.tar.gz

# conda environment installation
#conda env create -f environment.yml
