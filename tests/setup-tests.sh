#!/usr/bin/env bash

# Get and install anaconda for custom Python installation
echo -n "Downloading conda... "
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh 2> /dev/null > /dev/null
echo "Done!"
bash Miniconda-latest-Linux-x86_64.sh -b -p $HOME/miniconda2 -f
conda config --add channels r
conda config --add channels bioconda
conda config --add channels dakl

conda config --set always_yes yes --set changeps1 no
# - conda update -q conda
# Useful for debugging any issues with conda
conda info -a
conda install --file conda-list-tests.txt


TEST_FILE=$1
if [[ ${TEST_FILE} == integration* ]];
then

    # ref
    mkdir -p /tmp/genome-zipped/
    mkdir /tmp/genome
    wget --no-clobber -O /tmp/genome-zipped/test-genome.tar.gz https://export.uppmax.uu.se/b2010040/test-genome.tar.gz 2> /dev/null > /dev/null
    tar xvfz /tmp/genome-zipped/test-genome.tar.gz -C /tmp 2> /dev/null > /dev/null

    # test data
    wget --no-clobber -O /tmp/test-libraries.tar.gz https://export.uppmax.uu.se/b2010040/test-libraries.tar.gz 2> /dev/null > /dev/null
    tar xvfz /tmp/test-libraries.tar.gz -C /tmp 2> /dev/null > /dev/null

    # latex
    time sudo apt-get -y install texlive-latex-base 2> /dev/null > /dev/null
    time sudo apt-get -y install latex-xcolor texlive-fonts-recommended texlive-latex-extra 2> /dev/null > /dev/null

    # binaries for pipeline/s
    conda install --file conda-list.txt
fi