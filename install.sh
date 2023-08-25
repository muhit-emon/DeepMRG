#!/bin/bash

# diamond installation
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
chmod +x diamond
rm diamond-linux64.tar.gz

# prodigal installation
wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux
mv prodigal.linux prodigal
chmod +x prodigal

# conda environment installation
conda env create -f environment.yml
