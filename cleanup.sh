#!/usr/bin/env bash
KSFINDER2_HOME_DIR=$(pwd)

rm -r data

cd $KSFINDER2_HOME_DIR/assessments/assessment1
rm *.pt

cd $KSFINDER2_HOME_DIR/assessments/assessment2
rm *.pt

cd $KSFINDER2_HOME_DIR/assessments/assessment3
rm *.pt

cd $KSFINDER2_HOME_DIR/assessments/assessment4
rm ksfinder/ksfinder_predictions.csv
rm -r ksfinder/td2
rm link_phinder/linkphinder_predictions.csv
rm -r link_phinder/td2
rm phosformer-ST/phosformer_predictions.csv
rm -r phosformer-ST/td2
rm phosformer-ST/SER_THR_atlas.csv
rm predkinkg/predkinkg_predictions.csv
rm -r predkinkg/td2
rm license
find . -type f -name "*.png" -exec rm -f {} \;

cd $KSFINDER2_HOME_DIR/output
rm ksf2_predictions.zip

cd $KSFINDER2_HOME_DIR/model
rm ksfinder2.pt

cd $KSFINDER2_HOME_DIR/output/batch
rm *.csv
cd ..
rm *.csv

cd $KSFINDER2_HOME_DIR/output/handler
rm *.csv

