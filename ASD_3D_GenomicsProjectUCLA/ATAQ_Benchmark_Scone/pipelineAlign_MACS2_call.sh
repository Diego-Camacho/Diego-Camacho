#!/bin/bash

DIR=/u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/
FASTQ=$DIR/fastq_raw
NGMERGE_OUTPUT=$DIR/NGmerge_output
SCRIPT=$DIR/scripts
Debug=$DIR/debug

cd $FASTQ

for file in B*_R1.fastq.gz;
do
  name=${file%_R1.fastq.gz}
  qsub -hold_jid NGmerge_$name -cwd -o $Debug -e $Debug -N p_align_$name \
  -l h_data=10G,h_rt=10:00:00,highp -pe shared 8 \
  $SCRIPT/2_2_2_p_alignment.sh $folder $name
  qsub -hold_jid p_align_$name -cwd -o $Debug -e $Debug -N p_delchrM_$name \
  -l h_data=3G,h_rt=4:00:00,highp \
  $SCRIPT/2_2_3_p_delchrM.sh $folder $name
  qsub -hold_jid p_align_$name,p_delchrM_$name \
  -cwd -o $Debug -e $Debug -N p_markdup_$name \
  -l h_data=10G,h_rt=3:00:00,highp \
  $SCRIPT/2_2_4_p_picard_markdup.sh $folder $name
  qsub -hold_jid p_align_$name,p_delchrM_$name,p_markdup_$name \
  -cwd -o $Debug -e $Debug -N p_index_$name \
  -l h_rt=1:00:00,highp \
  $SCRIPT/2_2_5_p_index_bam.sh $folder $name
  qsub -hold_jid p_align_$name,p_delchrM_$name,p_markdup_$name,p_index_$name \
  -cwd -o $Debug -e $Debug -N p_unique_$name \
  -l h_rt=3:00:00,highp \
  $SCRIPT/2_2_6_p_remove_NonUniqueAligned.sh $folder $name
  qsub -hold_jid p_align_$name,p_delchrM_$name,p_markdup_$name,p_index_$name,p_unique_$name \
  -cwd -o $Debug -e $Debug -N p_sam2bed_$name \
  -l h_data=25G,h_rt=10:00:00,highp \
  $SCRIPT/2_2_7_p_samtobed.sh $folder $name
  qsub -hold_jid p_align_$name,p_delchrM_$name,p_markdup_$name,p_index_$name,p_unique_$name,p_sam2bed_$name \
  -cwd -o $Debug -e $Debug -N macs2_$name \
  -l h_data=5G,h_rt=4:00:00,highp \
  $SCRIPT/2_2_8_p_macs2.sh $folder $name
done

