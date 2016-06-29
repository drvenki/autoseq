#!/usr/bin/env bash

TEST_FILE=$1
if [[ ${TEST_FILE} == integration* ]];
then

    # scripts
    wget -O /tmp/autoseq-scripts.zip https://github.com/dakl/autoseq-scripts/archive/master.zip
    unzip -o /tmp/autoseq-scripts.zip -d /tmp

    # ref
    mkdir -p /tmp/genome-zipped/
    mkdir /tmp/genome
    wget --no-clobber -O /tmp/genome-zipped/test-genome.tar.gz https://export.uppmax.uu.se/b2010040/test-genome.tar.gz
    tar xvfz /tmp/genome-zipped/test-genome.tar.gz -C /tmp

    # test data
    wget --no-clobber -O /tmp/test-libraries.tar.gz https://export.uppmax.uu.se/b2010040/test-libraries.tar.gz
    tar xvfz /tmp/test-libraries.tar.gz -C /tmp

      # latex
    time sudo apt-get -y install texlive-latex-base
    time sudo apt-get -y install latex-xcolor texlive-fonts-recommended texlive-latex-extra

    # test genome and data
    bash tests/setup-tests.sh

    # Get and install anaconda for custom Python installation
    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    bash Miniconda-latest-Linux-x86_64.sh -b -p $HOME/miniconda2 -f
    export PATH=$HOME/miniconda2/bin/:/tmp/autoseq-scripts-master/:$PATH
    conda config --add channels r
    conda config --add channels bioconda

    conda config --set always_yes yes --set changeps1 no
    # - conda update -q conda
    # Useful for debugging any issues with conda
    conda info -a
    conda install --file conda-list.txt

fi