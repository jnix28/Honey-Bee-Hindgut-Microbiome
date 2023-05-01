#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20		                ##number of requested CPUs
#SBATCH --ntasks=8
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name filter4    ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder
#sbatch --mem=10G
#SBATCH --error=%x.e%j
#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /scratch/bell/jmnixon/project_files/apis_caste_sequence_data/QC/filter4

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
  #--output-path apis4-paired-end-seqs.qza \
  #--input-format PairedEndFastqManifestPhred33V2
#echo "import complete"

##De-noise. Error corrections.

echo "dada2"
date +"%d %B %Y %H:%M:%S"

#qiime dada2 denoise-paired \
## --i-demultiplexed-seqs apis4-paired-end-seqs.qza \
## --p-trim-left-f 25 \
## --p-trim-left-r 25 \
## --p-trunc-len-f 284 \
## --p-trunc-len-r 198 \
## --p-n-threads 2 \
## --p-max-ee-f 4 \
## --p-max-ee-r 4 \
## --verbose \
## --o-table Max_Error/apis4-table.qza \
## --o-representative-sequences Max_Error/apis4-rep-seqs.qza \
## --o-denoising-stats Max_Error/apis4-stats-dada2.qza 
  

#echo "finish denoise"
#date +"%d %B %Y %H:%M:%S"

#echo "Classify and get .qzv"
#qiime feature-classifier classify-sklearn \
## --i-classifier silva-classifier-v1v3.qza \
## --i-reads apis4-rep-seqs.qza \
## --o-classification apis4-taxonomy.qza

#qiime metadata tabulate \
## --m-input-file apis4-taxonomy.qza \
## --o-visualization apis4-taxonomy.qzv

#echo "Remove mitochondria and chloroplast"

#qiime taxa filter-table \
## --i-table apis4-table.qza \
## --i-taxonomy apis4-taxonomy.qza \
## --p-exclude mitochondria,chloroplast \
## --o-filtered-table apis4-table-no-mitochondria-no-chloroplast.qza

#Rename files so that the filtered one is table.qza

#mv apis4-table.qza apis4-table-with-mitochondria-and-chloroplast.qza
#mv apis4-table-no-mitochondria-no-chloroplast.qza apis4-table.qza

#echo "finish classifying data"
#Convert .qza to .qzv format

#qiime metadata tabulate \
## --m-input-file apis4-stats-dada2.qza \
## --o-visualization apis4-stats-dada2.qzv

#Create a .qzv from the output of denoise

#qiime feature-table summarize \
# --i-table apis4-table.qza \
# --o-visualization apis4-table.qzv \
# --m-sample-metadata-file apis4-metadata.tsv

##qiime feature-table tabulate-seqs \
## --i-data apis4-rep-seqs.qza \
## --o-visualization apis4-rep-seqs.qzv

echo "Create a phylogenetic tree. (Phylogenetic method)"
#qiime phylogeny align-to-tree-mafft-fasttree \
# --i-sequences apis4-rep-seqs.qza \
# --o-alignment apis4-aligned-rep-seqs.qza \
# --o-masked-alignment apis4-masked-aligned-rep-seqs.qza \
# --o-tree apis4-unrooted-tree.qza \
# --o-rooted-tree apis4-rooted-tree.qza

echo "Rarefaction"
#qiime diversity alpha-rarefaction \
# --i-table apis4-table.qza \
# --i-phylogeny apis4-rooted-tree.qza \
# --p-max-depth 1275 \
# --m-metadata-file apis4-metadata.tsv \
# --o-visualization apis4-alpha-rarefaction.qzv

echo "Make a stacked barplot"
#qiime taxa barplot \
# --i-table apis4-table.qza \
# --i-taxonomy apis4-taxonomy.qza \
# --m-metadata-file apis4-metadata.tsv \
# --o-visualization apis4-taxa-bar-plots.qzv

#Calculate diversity metrics. You will need to designate a sampling depth. You might want to consult the rarifaction output prior to this step.

#qiime diversity core-metrics-phylogenetic \
# --i-phylogeny apis4-rooted-tree.qza \
# --i-table apis4-table.qza \
# --p-sampling-depth 130 \
# --m-metadata-file apis4-metadata.tsv \
# --output-dir apis4-core-metrics-results

#qiime diversity alpha \
# --i-table apis4-table.qza \
# --p-metric chao1 \
# --o-alpha-diversity apis4-core-metrics-results/chao1_vector.qza

#qiime diversity alpha \
# --i-table apis4-table.qza \
# --p-metric simpson \
# --o-alpha-diversity apis4-core-metrics-results/simpson_vector.qza

echo "Diversity testing complete"
#Statistical testing

#qiime diversity alpha-group-significance \
# --i-alpha-diversity apis4-core-metrics-results/faith_pd_vector.qza \
# --m-metadata-file apis4-metadata.tsv \
# --o-visualization apis4-core-metrics-results/faith-pd-group-significance.qzv

#qiime diversity alpha-group-significance \
# --i-alpha-diversity apis4-core-metrics-results/evenness_vector.qza \
# --m-metadata-file apis4-metadata.tsv \
# --o-visualization apis4-core-metrics-results/evenness-group-significance.qzv

#qiime diversity alpha-group-significance \
# --i-alpha-diversity apis4-core-metrics-results/shannon_vector.qza \
# --m-metadata-file apis4-metadata.tsv \
# --o-visualization apis4-core-metrics-results/shannon-group-significance.qzv

#qiime diversity alpha-group-significance \
# --i-alpha-diversity apis4-core-metrics-results/observed_features_vector.qza \
# --m-metadata-file apis4-metadata.tsv \
# --o-visualization apis4-core-metrics-results/observed_features-group-significance.qzv

#qiime diversity alpha-group-significance \
# --i-alpha-diversity apis4-core-metrics-results/chao1_vector.qza \
# --m-metadata-file apis4-metadata.tsv \
# --o-visualization apis4-core-metrics-results/chao1-group-significance.qzv

#qiime diversity alpha-group-significance \
# --i-alpha-diversity apis4-core-metrics-results/simpson_vector.qza \
# --m-metadata-file apis4-metadata.tsv \
# --o-visualization apis4-core-metrics-results/simpson-group-significance.qzv

##You will have to enter a categorical column name here from your metadata file.
#qiime diversity beta-group-significance \
# --i-distance-matrix apis4-core-metrics-results/weighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis4-metadata.tsv \
# --m-metadata-column Colony \
# --o-visualization apis4-core-metrics-results/weighted-unifrac-colony-significance.qzv \
# --p-pairwise

#qiime diversity beta-group-significance \
# --i-distance-matrix apis4-core-metrics-results/weighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis4-metadata.tsv \
# --m-metadata-column Colony \
# --p-method permdisp \
# --o-visualization apis4-core-metrics-results/weighted-unifrac-colony-permdisp.qzv \
# --p-pairwise

#qiime diversity beta-group-significance \
# --i-distance-matrix apis4-core-metrics-results/weighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis4-metadata.tsv \
# --m-metadata-column Caste \
# --o-visualization apis4-core-metrics-results/weighted-unifrac-caste-significance.qzv \
# --p-pairwise

#qiime diversity beta-group-significance \
# --i-distance-matrix apis4-core-metrics-results/weighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis4-metadata.tsv \
# --m-metadata-column Caste \
# --p-method permdisp \
# --o-visualization apis4-core-metrics-results/weighted-unifrac-caste-permdisp.qzv \
# --p-pairwise

#qiime diversity beta-group-significance \
# --i-distance-matrix apis4-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis4-metadata.tsv \
# --m-metadata-column Colony \
# --o-visualization apis4-core-metrics-results/unweighted-unifrac-colony-significance.qzv \
# --p-pairwise

#qiime diversity beta-group-significance \
# --i-distance-matrix apis4-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis4-metadata.tsv \
# --m-metadata-column Colony \
# --p-method permdisp \
# --o-visualization apis4-core-metrics-results/unweighted-unifrac-colony-disp.qzv \
# --p-pairwise
 
#qiime diversity beta-group-significance \
# --i-distance-matrix apis4-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis4-metadata.tsv \
# --m-metadata-column Caste \
# --o-visualization apis4-core-metrics-results/unweighted-unifrac-caste-significance.qzv \
# --p-pairwise

#qiime diversity beta-group-significance \
# --i-distance-matrix apis4-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis4-metadata.tsv \
# --m-metadata-column Caste \
# --p-method permdisp \
# --o-visualization apis4-core-metrics-results/unweighted-unifrac-caste-disp.qzv \
# --p-pairwise

#ANCOM  
qiime composition add-pseudocount \
  --i-table apis4-table.qza \
  --o-composition-table apis4-comp-table.qza
  
qiime composition ancom \
  --i-table apis4-comp-table.qza \
  --m-metadata-file apis4-metadata.tsv \
  --m-metadata-column Caste \
  --o-visualization apis4-ancom-caste.qzv


echo "Collapse levels and ancom again."

qiime taxa collapse \
  --i-table apis4-table.qza \
  --i-taxonomy apis4-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table apis4-table-l6.qza

qiime composition add-pseudocount \
  --i-table apis4-table-l6.qza \
  --o-composition-table apis4-comp-table-l6.qza

qiime composition ancom \
  --i-table apis4-comp-table-l6.qza \
  --m-metadata-file apis4-metadata.tsv \
  --m-metadata-column Caste \
  --o-visualization apis4-l6-ancom-caste.qzv

echo "Stats completed"
echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch project_apis_filter_0.slurm