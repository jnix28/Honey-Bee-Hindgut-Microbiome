#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20		                ##number of requested CPUs
#SBATCH --ntasks=8
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name filter0    ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder
#sbatch --mem=10G
#SBATCH --error=%x.e%j
#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /scratch/bell/jmnixon/project_files/apis_caste_sequence_data/QC/filter0

pwd

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2022.8"
module load Qiime/2-2022.8

#echo "Data Import"
#Make the sequences into a qiime recognizable format

#qiime tools import \
  #--type 'SampleData[PairedEndSequencesWithQuality]' \
  #--input-path apis-manifest1.tsv \
  #--output-path apis-paired-end-seqs.qza \
  #--input-format PairedEndFastqManifestPhred33V2
#echo "import complete"

##De-noise. Error corrections.

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
# --i-demultiplexed-seqs apis0-paired-end-seqs.qza \
# --p-trim-left-f 25 \
## --p-trim-left-r 25 \
# --p-trunc-len-f 284 \
# --p-trunc-len-r 198 \
# --p-n-threads 2 \
# --o-table apis0-table.qza \
# --o-representative-sequences apis0-rep-seqs.qza \
# --o-denoising-stats apis0-stats-dada2.qza


#echo "finish denoise"
#date +"%d %B %Y %H:%M:%S"

#echo "Classify and get .qzv"
#qiime feature-classifier classify-sklearn \
# --i-classifier silva-classifier-v1v3.qza \
# --i-reads apis0-rep-seqs.qza \
# --o-classification apis0-taxonomy.qza

#qiime metadata tabulate \
# --m-input-file apis0-taxonomy.qza \
# --o-visualization apis0-taxonomy.qzv

#echo "Remove mitochondria and chloroplast"

#qiime taxa filter-table \
# --i-table apis0-table.qza \
# --i-taxonomy apis0-taxonomy.qza \
# --p-exclude mitochondria,chloroplast \
# --o-filtered-table apis0-table-no-mitochondria-no-chloroplast.qza

#Rename files so that the filtered one is table.qza

#mv apis0-table.qza apis0-table-with-mitochondria-and-chloroplast.qza
#mv apis0-table-no-mitochondria-no-chloroplast.qza table.qza

#echo "finish classifying data"
#Convert .qza to .qzv format

#qiime metadata tabulate \
# --m-input-file apis0-stats-dada2.qza \
# --o-visualization apis0-stats-dada2.qzv

#Create a .qzv from the output of denoise

#qiime feature-table summarize \
# --i-table apis0-table.qza \
# --o-visualization apis0-table.qzv \
# --m-sample-metadata-file apis0-metadata.tsv

##qiime feature-table tabulate-seqs \
## --i-data apis0-rep-seqs.qza \
## --o-visualization apis0-rep-seqs.qzv

#echo "Create a phylogenetic tree. (Phylogenetic method)"
#qiime phylogeny align-to-tree-mafft-fasttree \
# --i-sequences apis0-rep-seqs.qza \
# --o-alignment apis0-aligned-rep-seqs.qza \
# --o-masked-alignment apis0-masked-aligned-rep-seqs.qza \
# --o-tree apis0-unrooted-tree.qza \
# --o-rooted-tree apis0-rooted-tree.qza

#echo "Rarefaction"
#qiime diversity alpha-rarefaction \
# --i-table apis0-table.qza \
# --i-phylogeny apis0-rooted-tree.qza \
# --p-max-depth 47553 \
# --m-metadata-file apis0-metadata.tsv \
# --o-visualization apis0-alpha-rarefaction.qzv

#echo "Make a stacked barplot"
#qiime taxa barplot \
# --i-table apis0-table.qza \
# --i-taxonomy apis0-taxonomy.qza \
# --m-metadata-file apis0-metadata.tsv \
# --o-visualization apis0-taxa-bar-plots.qzv

#Calculate diversity metrics. You will need to designate a sampling depth. You might want to consult the rarifaction output prior to this step.

#qiime diversity core-metrics-phylogenetic \
# --i-phylogeny apis0-rooted-tree.qza \
# --i-table apis0-table.qza \
# --p-sampling-depth 86 \
# --m-metadata-file apis0-metadata.tsv \
# --output-dir apis0-core-metrics-results

#qiime diversity alpha \
# --i-table apis0-table.qza \
# --p-metric chao1 \
# --o-alpha-diversity apis0-core-metrics-results/apis0-chao1_vector.qza

#qiime diversity alpha \
# --i-table apis0-table.qza \
# --p-metric simpson \
# --o-alpha-diversity apis0-core-metrics-results/apis0-simpson_vector.qza

#echo "Diversity testing complete"
#Statistical testing

#qiime diversity alpha-correlation \
# --i-alpha-diversity apis0-core-metrics-results/faith_pd_vector.qza \
# --m-metadata-file apis0-metadata.tsv \
# --o-visualization apis0-core-metrics-results/faith-pd-group-significance.qzv

#qiime diversity alpha-correlation \
# --i-alpha-diversity apis0-core-metrics-results/evenness_vector.qza \
# --m-metadata-file apis0-metadata.tsv \
# --o-visualization apis0-core-metrics-results/evenness-group-significance.qzv

#qiime diversity alpha-correlation \
# --i-alpha-diversity apis0-core-metrics-results/shannon_vector.qza \
# --m-metadata-file apis0-metadata.tsv \
# --o-visualization apis0-core-metrics-results/shannon-group-significance.qzv

#qiime diversity alpha-correlation \
# --i-alpha-diversity apis0-core-metrics-results/observed_features_vector.qza \
# --m-metadata-file apis0-metadata.tsv \
# --o-visualization apis0-core-metrics-results/observed_features-group-significance.qzv

#qiime diversity alpha-correlation \
# --i-alpha-diversity apis0-core-metrics-results/chao1_vector.qza \
# --m-metadata-file apis0-metadata.tsv \
# --o-visualization apis0-core-metrics-results/chao1-group-significance.qzv

#qiime diversity alpha-correlation \
# --i-alpha-diversity apis0-core-metrics-results/simpson_vector.qza \
# --m-metadata-file apis0-metadata.tsv \
# --o-visualization apis0-core-metrics-results/simpson-group-significance.qzv

#qiime diversity beta-group-significance \
# --i-distance-matrix apis0-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis0-metadata.tsv \
# --m-metadata-column Caste \
# --o-visualization apis0-core-metrics-results/unweighted-unifrac-Caste-significance.qzv \
# --p-pairwise

#qiime diversity beta-group-significance \
# --i-distance-matrix apis0-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis0-metadata.tsv \
# --m-metadata-column Caste \
# --p-method permdisp \
# --o-visualization apis0-core-metrics-results/unweighted-unifrac-Caste-significance.qzv

##You will have to enter a categorical column name here from your metadata file.
#qiime diversity beta-group-significance \
# --i-distance-matrix apis0-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis0-metadata.tsv \
# --m-metadata-column Colony \
# --o-visualization apis0-core-metrics-results/unweighted-unifrac-colony-significance.qzv \
# --p-pairwise

#qiime diversity beta-group-significance \
# --i-distance-matrix apis0-core-metrics-results/unweighted_unifrac_distance_matrix.qza \
# --m-metadata-file apis0-metadata.tsv \
# --m-metadata-column Colony \
# --p-method permdisp \
# --o-visualization apis0-core-metrics-results/unweighted-unifrac-colony-significance.qzv 
 
#ANCOM  
#qiime composition add-pseudocount \
#  --i-table apis0-table.qza \
#  --o-composition-table apis0-comp-table.qza
  
#qiime composition ancom \
#  --i-table apis0-comp-table.qza \
#  --m-metadata-file apis0-metadata.tsv \
#  --m-metadata-column Caste \
#  --o-visualization apis0-ancom-caste.qzv


#echo "Collapse levels and ancom again."

#qiime taxa collapse \
#  --i-table apis0-table.qza \
#  --i-taxonomy apis0-taxonomy.qza \
#  --p-level 6 \
#  --o-collapsed-table apis0-table-l6.qza

#qiime composition add-pseudocount \
#  --i-table apis0-table-l6.qza \
#  --o-composition-table apis0-comp-table-l6.qza

#qiime composition ancom \
#  --i-table apis0-comp-table-l6.qza \
#  --m-metadata-file apis0-metadata.tsv \
#  --m-metadata-column Caste \
#  --o-visualization apis0-l6-ancom-caste.qzv

#############################################################################
#Lefse
#qiime taxa collapse \
#--i-table apis0-table.qza \
#--o-collapsed-table collapse.table.qza \
#--p-level 6 \
#--i-taxonomy apis0-taxonomy.qza

#qiime feature-table relative-frequency \
#--i-table collapse.table.qza \
#--o-relative-frequency-table collapse.frequency.table.qza \
#--output-dir collapse.frequency/

#qiime tools export \
#--input-path collapse.frequency.table.qza \
#--output-path collapse.frequency/

#biom convert \
#-i collapse.frequency/feature-table.biom \
#-o collapse.frequency.table.txt \
#--header-key “taxonomy” --to-tsv

#############################################################################
#Random Forest

#qiime sample-classifier classify-samples \
#  --i-table apis0-table.qza \
#  --m-metadata-file apis0-metadata.tsv \
#  --m-metadata-column Caste \
#  --p-optimize-feature-selection \
#  --p-parameter-tuning \
#  --p-estimator RandomForestClassifier \
#  --p-n-estimators 20 \
#  --p-random-state 123 \
#  --output-dir apis0-classifier

#qiime metadata tabulate \
#  --m-input-file apis0-classifier/predictions.qza \
#  --o-visualization apis0-classifier/predictions.qzv

#qiime metadata tabulate \
#  --m-input-file apis0-classifier/probabilities.qza \
#  --o-visualization apis0-classifier/probabilities.qzv

#qiime metadata tabulate \
#  --m-input-file apis0-classifier/feature_importance.qza \
#  --o-visualization apis0-classifier/feature_importance.qzv

#qiime feature-table filter-features \
##  --i-table apis0-table.qza \
#  --m-metadata-file apis0-classifier/feature_importance.qza \
#  --o-filtered-table apis0-classifier/important-feature-table.qza

#qiime sample-classifier heatmap \
#  --i-table apis0-table.qza \
#  --i-importance apis0-classifier/feature_importance.qza \
#  --m-sample-metadata-file apis0-metadata.tsv \
#  --m-sample-metadata-column Caste \
#  --p-group-samples \
#  --p-feature-count 30 \
# --o-filtered-table apis0-classifier/important-feature-table-top-30.qza \
#  --o-heatmap apis0-classifier/important-feature-heatmap.qzv

#qiime sample-classifier heatmap \
#  --i-table apis0-table.qza \
#  --i-importance apis0-classifier/feature_importance.qza  \
#  --m-sample-metadata-file apis0-metadata.tsv  \
#  --m-sample-metadata-column Caste \
#  --m-feature-metadata-file apis0-taxonomy.qza \
#  --m-feature-metadata-column Taxon \
#  --p-group-samples  \
#  --p-feature-count 30  \
#  --o-filtered-table apis0-classifier/important-feature-table-taxonomy-top-30.qza \
#  --o-heatmap apis0-classifier/important-feature-heatmap-taxonomy.qzv

#########################################################################################################
#Load Picrust2

module use /depot/microbiome/apps/modulefiles
module load conda-env/qiime2-2021.11-py3.8.5

qiime picrust2 full-pipeline \
   --i-table apis0-table.qza \
   --i-seq apis0-rep-seqs.qza \
   --output-dir q2-picrust2_output \
   --p-placement-tool sepp \
   --p-threads 1 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose


echo "Stats completed"
echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch project_apis_filter_0.slurm
