#!/usr/bin/env bash
KSFINDER2_HOME_DIR=$(pwd)

# Create data directories
mkdir data
cd $KSFINDER2_HOME_DIR/data
mkdir kg
mkdir embeddings

## Download kg data
cd $KSFINDER2_HOME_DIR/data/kg
wget https://zenodo.org/records/14713589/files/kg_data.zip
unzip kg_data.zip
rm kg_data.zip

## Download embeddings
cd $KSFINDER2_HOME_DIR/data/embeddings
wget https://zenodo.org/records/14713589/files/embeddings.zip
unzip embeddings.zip
rm embeddings.zip

## Download classification assessment datasets, kg kinases and substrate_motifs
cd $KSFINDER2_HOME_DIR/data
wget https://zenodo.org/records/14713589/files/classifier_datasets.zip
unzip classifier_datasets.zip
rm classifier_datasets.zip

wget https://zenodo.org/records/14713589/files/assessment_datasets.zip
unzip assessment_datasets.zip
rm assessment_datasets.zip

wget https://zenodo.org/records/14713589/files/kg_ks.zip
unzip kg_ks.zip
rm kg_ks.zip

## Download assessment1 models
cd $KSFINDER2_HOME_DIR/assessments/assessment1
wget https://zenodo.org/records/14713589/files/kge_models_assess1.zip
unzip kge_models_assess1.zip
rm kge_models_assess1.zip

## Download assessment2 models
cd $KSFINDER2_HOME_DIR/assessments/assessment2
wget https://zenodo.org/records/14713589/files/models_assess2.zip
unzip models_assess2.zip
rm models_assess2.zip

## Download assessment3 models
cd $KSFINDER2_HOME_DIR/assessments/assessment3
wget https://zenodo.org/records/14713589/files/models_assess3.zip
unzip models_assess3.zip
rm models_assess3.zip

## Download datasets for assessment4 (comparative analysis)
cd $KSFINDER2_HOME_DIR/assessments/assessment4
wget https://zenodo.org/records/14713589/files/other_model_predictions.zip
unzip other_model_predictions.zip
rm other_model_predictions.zip
mv SER_THR_atlas.csv phosformer-ST
mkdir ksfinder/td2
mkdir link_phinder/td2
mkdir phosformer-ST/td2
mkdir predkinkg/td2

# Create output directory
cd $KSFINDER2_HOME_DIR
mkdir output
mkdir model

cd $KSFINDER2_HOME_DIR/output
## Download KSFinder 2.0's predictions if needed
wget https://zenodo.org/records/14713589/files/ksf2_predictions.zip

cd $KSFINDER2_HOME_DIR/model
## Download KSFinder 2.0
wget https://zenodo.org/records/14713589/files/model_ksfinder2.zip
unzip model_ksfinder2.zip
rm model_ksfinder2.zip