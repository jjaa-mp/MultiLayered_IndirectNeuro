#!/bin/bash

## Create dir in desired location

mkdir -p /home/rstudio/jm_rstudio/data_indNeuro
export DATA_PATH=/home/rstudio/jm_rstudio/data_indNeuro


## Download Trevino 2021 RNA data. metadata, and supplementary data

wget https://atrev.s3.amazonaws.com/brainchromatin/rna_counts.tsv.gz -P ${DATA_PATH}
wget https://atrev.s3.amazonaws.com/brainchromatin/rna_cell_metadata.txt -P ${DATA_PATH}
wget https://www.cell.com/cms/10.1016/j.cell.2021.07.039/attachment/94887a22-4410-42a5-ab28-8c4fc093e99a/mmc1.xlsx -O ${DATA_PATH}/ST_trevino21_rna.xlsx

## Polioudakis 2019 data downloaded from http://solo.bmap.ucla.edu/shiny/webapp/

