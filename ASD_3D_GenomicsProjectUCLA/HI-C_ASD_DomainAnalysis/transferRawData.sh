#!/bin/bash

DIR=/u/project/geschwind/jennybea/HiC_Pro/
FILES_HiCPro=$DIR/ASD_800Mreads_HiCPro_BlankParameters

cd $FILES_HiCPro
mkdir -p rawdata_matrix_10kb iced_matrix_100kb iced_matrix_40kb iced_matrix_10kb
mkdir -p iced_matrix_2kb iced_matrix_3kb iced_matrix_4kb iced_matrix_5kb
mkdir -p rawdata_matrix_2kb rawdata_matrix_3kb rawdata_matrix_4kb rawdata_matrix_5kb

for i in B*
# for i in B5718
do
  ln -s $FILES_HiCPro/${i}/hic_results/matrix/rawdata/raw/10000/rawdata_10000.matrix rawdata_matrix_10kb/
  mv rawdata_matrix_10kb/rawdata_10000.matrix rawdata_matrix_10kb/${i}_rawdata_10000.matrix
done

## Transfer iced matrix
for i in B*
do
  ln -s $FILES_HiCPro/${i}/hic_results/matrix/rawdata/iced/100000/rawdata_100000_iced.matrix iced_matrix_100kb/
  mv iced_matrix_100kb/rawdata_100000_iced.matrix iced_matrix_100kb/${i}_rawdata_100000_iced.matrix
done
ln -s $FILES_HiCPro/${i}/hic_results/matrix/rawdata/raw/100000/rawdata_100000_abs.bed iced_matrix_100kb/
