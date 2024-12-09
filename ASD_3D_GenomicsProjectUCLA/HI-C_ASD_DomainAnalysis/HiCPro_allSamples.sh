#!/bin/bash
qrsh -l h_data=2G,h_rt=1:00:00

HOME=/u/project/geschwind/jennybea/HiC_Pro/
FILES=$HOME/ASD_800Mreads
mkdir -p $FILES/debug
FILES_HiCPro=$HOME/ASD_800Mreads_HiCPro_BlankParameters # Note the change here in contig-hicpro.txt
mkdir -p $FILES_HiCPro
# Copy config-hicpro_B5302.txt to $FILES

sample_array=("B5144" "B5666" "B4721" "B4787" "B5309" "B5242A" "B7769" "B5000" "B5352" "B7109" "B4498" "B6200" "B5302" "B5342B")
sample_array=("B5718")
cd $FILES
for i in B*
for i in ${sample_array[@]}
do
  FASTQ=$FILES/$i
  # Do not creat the following JY,JH,JG output directories in advance.
  OUTPUT=$FILES_HiCPro/$i
  CONFIG=$FILES/config-hicpro_allsamples.txt #config-hicpro_$i.txt
  qsub -cwd -o $FILES/debug -e $FILES/debug -V -S /bin/bash -N HiCpro_${i} -l h_data=20G,h_rt=120:00:00,highp -pe shared 8 $HOME/scripts/2_HiCpro_2ndround.sh $FASTQ $OUTPUT $CONFIG #Build the raw matrix:2_HiCpro_2ndround.sh, Include ice: 3_2_2_HiCpro_Dup15q2a3.sh
done

