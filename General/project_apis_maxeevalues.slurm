#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20		                ##number of requested CPUs
#SBATCH --ntasks=8
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name test_maxerror    ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder
#sbatch --mem=10G
#SBATCH --error=%x.e%j
#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /scratch/bell/jmnixon/project_files/apis_caste_sequence_data/QC

pwd

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2022.8"
module load Qiime/2-2022.8

echo "Data Import"
#Make the sequences into a qiime recognizable format

#qiime tools import \
  #--type 'SampleData[PairedEndSequencesWithQuality]' \
  #--input-path apis-manifest1.tsv \
  #--output-path apis-paired-end-seqs.qza \
  #--input-format PairedEndFastqManifestPhred33V2
#echo "import complete"

##De-noise. Error corrections.2

#echo "dada2"
#date +"%d %B %Y %H:%M:%S"

#qiime demux emp-paired \
#  --i-seqs apis-paired-end-sequences.qza \
#  --m-barcodes-file apis-metadata.tsv \
#  --m-barcodes-column Sample-ID \
#  --p-rev-comp-mapping-barcodes \
#  --o-per-sample-sequences apis-demux.qza \
#  --o-error-correction-details apis-demux-details.qza

#qiime demux summarize \
#  --i-data apis-paired-end-seqs.qza \
#  --o-visualization apis-demux.qzv

#qiime dada2 denoise-paired \
#  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
#  --p-trim-left-f 25 \
#  --p-trim-left-r 25 \
#  --p-trunc-len-f 284 \
#  --p-trunc-len-r 198 \
#  --p-n-threads 2 \
#  --p-max-ee-f 2 \
#  --p-max-ee-r 2 \
#   --verbose \
#  --o-table Max_Error/apis2-table.qza \
#  --o-representative-sequences Max_Error/apis2-rep-seqs.qza \
#  --o-denoising-stats Max_Error/apis2-stats-dada2.qza \
 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 25 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 198 \
  --p-n-threads 2 \
  --p-max-ee-f 3 \
  --p-max-ee-r 3 \
  --verbose \
  --o-table Max_Error/apis3-table.qza \
  --o-representative-sequences Max_Error/apis3-rep-seqs.qza \
  --o-denoising-stats Max_Error/apis3-stats-dada2.qza 
 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 25 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 198 \
  --p-n-threads 2 \
  --p-max-ee-f 4 \
  --p-max-ee-r 4 \
  --verbose \
  --o-table Max_Error/apis4-table.qza \
  --o-representative-sequences Max_Error/apis4-rep-seqs.qza \
  --o-denoising-stats Max_Error/apis4-stats-dada2.qza 
 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 25 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 198 \
  --p-n-threads 2 \
  --p-max-ee-f 5 \
  --p-max-ee-r 5 \
  --verbose \
  --o-table Max_Error/apis5-table.qza \
  --o-representative-sequences Max_Error/apis5-rep-seqs.qza \
  --o-denoising-stats Max_Error/apis5-stats-dada2.qza 
 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 25 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 198 \
  --p-n-threads 2 \
  --p-max-ee-f 6 \
  --p-max-ee-r 6 \
  --verbose \
  --o-table Max_Error/apis6-table.qza \
  --o-representative-sequences Max_Error/apis6-rep-seqs.qza \
  --o-denoising-stats Max_Error/apis6-stats-dada2.qza 
 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 25 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 198 \
  --p-n-threads 2 \
  --p-max-ee-f 7 \
  --p-max-ee-r 7 \
  --verbose \
  --o-table Max_Error/apis7-table.qza \
  --o-representative-sequences Max_Error/apis7-rep-seqs.qza \
  --o-denoising-stats Max_Error/apis7-stats-dada2.qza 
 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 25 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 198 \
  --p-n-threads 2 \
  --p-max-ee-f 8 \
  --p-max-ee-r 8 \
  --verbose \
  --o-table Max_Error/apis8-table.qza \
  --o-representative-sequences Max_Error/apis8-rep-seqs.qza \
  --o-denoising-stats Max_Error/apis8-stats-dada2.qza 
 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 25 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 198 \
  --p-n-threads 2 \
  --p-max-ee-f 9 \
  --p-max-ee-r 9 \
  --verbose \
  --o-table Max_Error/apis9-table.qza \
  --o-representative-sequences Max_Error/apis9-rep-seqs.qza \
  --o-denoising-stats Max_Error/apis9-stats-dada2.qza 
 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs apis-paired-end-seqs.qza \
  --p-trim-left-f 25 \
  --p-trim-left-r 25 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 198 \
  --p-n-threads 2 \
  --p-max-ee-f 10 \
  --p-max-ee-r 10 \
  --verbose \
  --o-table Max_Error/apis10-table.qza \
  --o-representative-sequences Max_Error/apis10-rep-seqs.qza \
  --o-denoising-stats Max_Error/apis10-stats-dada2.qza 
 

echo "finish denoise"
date +"%d %B %Y %H:%M:%S"