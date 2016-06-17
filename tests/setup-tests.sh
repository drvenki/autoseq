#!/usr/bin/env bash

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
