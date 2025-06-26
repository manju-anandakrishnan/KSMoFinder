# #!/usr/bin/env bash
KSFINDER2_HOME_DIR=$(pwd)
DATA_BASE_URI=https://zenodo.org/records/15730847/files

# Create data directories
mkdir data
cd $KSFINDER2_HOME_DIR/data
mkdir kg
mkdir embeddings

## Download kg data
cd $KSFINDER2_HOME_DIR/data/kg
wget $DATA_BASE_URI/kg_data.zip
unzip kg_data.zip
rm kg_data.zip
rm -r kg_data
rm -r __MACOSX

## Download embeddings
cd $KSFINDER2_HOME_DIR/data/embeddings
wget $DATA_BASE_URI/embeddings.zip
mv embeddings/*.csv .
unzip embeddings.zip
rm embeddings.zip
rm embeddings
rm -r __MACOSX

# ## Download classification assessment datasets, kg kinases and substrate_motifs
cd $KSFINDER2_HOME_DIR/data
wget $DATA_BASE_URI/classifier_datasets.zip
unzip classifier_datasets.zip
rm classifier_datasets.zip

wget $DATA_BASE_URI/assessments_data.zip
unzip assessments_data.zip
rm assessments_data.zip
mv assessments_data/* .
rm -r assessments_data   

wget $DATA_BASE_URI/kg_ks.zip
unzip kg_ks.zip
rm kg_ks.zip

wget $DATA_BASE_URI/other_datafiles.zip
unzip other_datafiles.zip
rm other_datafiles.zip

rm -r __MACOSX

## Download assessment1 models
cd $KSFINDER2_HOME_DIR/assessments/assessment1
wget $DATA_BASE_URI/kge_models_assess1.zip
unzip kge_models_assess1.zip
rm kge_models_assess1.zip
mv kge_models_assess1/* .
rm -r kge_models_assess1
rm -r __MACOSX

## Download assessment2 models
cd $KSFINDER2_HOME_DIR/assessments/assessment2
wget $DATA_BASE_URI/models_assess2.zip
unzip models_assess2.zip
rm models_assess2.zip
rm -r __MACOSX

## Download assessment3 models
cd $KSFINDER2_HOME_DIR/assessments/assessment3
wget $DATA_BASE_URI/models_assess3.zip
unzip models_assess3.zip
rm models_assess3.zip
rm -r __MACOSX

## Download datasets for assessment4 (comparative analysis)
cd $KSFINDER2_HOME_DIR/assessments/assessment4
wget $DATA_BASE_URI/other_model_predictions.zip
unzip other_model_predictions.zip
rm other_model_predictions.zip
mv SER_THR_atlas.csv phosformer-ST
mkdir ksfinder/td2
mkdir link_phinder/td2
mkdir phosformer-ST/td2
mkdir predkinkg/td2
rm -r __MACOSX

# Create output directory
cd $KSFINDER2_HOME_DIR
mkdir output
mkdir model

cd $KSFINDER2_HOME_DIR/output
## Download KSFinder 2.0's predictions if needed
wget $DATA_BASE_URI/ksf2_predictions.zip

cd $KSFINDER2_HOME_DIR/model
## Download KSFinder 2.0
wget $DATA_BASE_URI/model_ksfinder.zip
unzip model_ksfinder.zip
mv model_ksfinder/ksfinder2.pt .
rm model_ksfinder.zip
rm -r model_ksfinder
rm -r __MACOSX