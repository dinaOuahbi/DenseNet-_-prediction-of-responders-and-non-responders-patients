#!/bin/bash
#set NUMEXPR_MAX_THREADS=68
#eval slide=$1
python 2_DenseNet_v6.py


#for slide in `ls /work/shared/ptbc/CNN_Pancreas_V2/Donnees/project_kernel03_scan21000500/scan_tiles/`
#do
#  echo $slide
#  qsub -q batch script/rein.sh $slide
#done

# while IFS=read line
#do
#qsub -q batch ../script/s7.sh $line 
#done < cnn1
