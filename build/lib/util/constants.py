import os

DIR_ASSESSMENTS = 'assessments'
DIR_ASSESSMENT1 = os.path.join(DIR_ASSESSMENTS, 'assessment1')
DIR_ASSESSMENT2 = os.path.join(DIR_ASSESSMENTS, 'assessment2')
DIR_ASSESSMENT3 = os.path.join(DIR_ASSESSMENTS, 'assessment3')
DIR_ASSESSMENT4 = os.path.join(DIR_ASSESSMENTS, 'assessment4')

MODEL_TRANSE = os.path.join(DIR_ASSESSMENT1, 'model_transE.pt')
MODEL_DISTMULT = os.path.join(DIR_ASSESSMENT1, 'model_distmult.pt')
MODEL_COMPLEX = os.path.join(DIR_ASSESSMENT1, 'model_complex.pt')
MODEL_EXPRESSIVE = os.path.join(DIR_ASSESSMENT1, 'model_expressivE.pt')

MODEL_ESM2 = os.path.join(DIR_ASSESSMENT2, 'model_esm2.pt')
MODEL_ESM3 = os.path.join(DIR_ASSESSMENT2, 'model_esm3.pt')
MODEL_PROTT5 = os.path.join(DIR_ASSESSMENT2, 'model_protT5.pt')
MODEL_KSF2_ASSESS2 = os.path.join(DIR_ASSESSMENT2, 'model_ksf2.pt')
MODEL_RANDOM = os.path.join(DIR_ASSESSMENT2, 'model_random.pt')

MODEL_STRUCT_SEQ = os.path.join(DIR_ASSESSMENT3,'model_struct_seq.pt')
MODEL_STRUCT_SEQ_FUNC = os.path.join(DIR_ASSESSMENT3,'model_struct_seq_func.pt')
MODEL_STRUCT_FUNC = os.path.join(DIR_ASSESSMENT3,'model_struct_func.pt')
MODEL_SEQ_FUNC = os.path.join(DIR_ASSESSMENT3,'model_seq_func.pt')
MODEL_FUNC = os.path.join(DIR_ASSESSMENT3,'model_func.pt')
MODEL_STRUCT_MOTIF = os.path.join(DIR_ASSESSMENT3,'model_struct_motif.pt')

DIR_MODEL = 'model'
MODEL_KSF2 = os.path.join(DIR_MODEL,'ksfinder2.pt')
DK_MODEL_KSF2 = os.path.join(DIR_MODEL,'dk_ksfinder2.pt')

DIR_DATA = 'data'
DIR_TD1 = 'td1_ratio_1.1'
DIR_TD2 = 'td2_tr_dist'

DIR_DATA_CLASSIFIER = os.path.join(DIR_DATA,'classifier')
CSV_CLF_TRAIN_DATA = os.path.join(DIR_DATA_CLASSIFIER,'clf_train_data.csv')
CSV_CLF_TEST_D1 = os.path.join(DIR_DATA_CLASSIFIER,DIR_TD1,'clf_test_data.csv')
CSV_CLF_TEST_D2 = os.path.join(DIR_DATA_CLASSIFIER,DIR_TD2,'clf_test_data.csv')

DIR_DATA_CLASSIFIER_ASSESS1 = os.path.join(DIR_DATA,'assessment1')
CSV_CLF_TRAIN_DATA_ASSESS1 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS1,'clf_train_data.csv')
CSV_CLF_TEST_D1_ASSESS1 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS1,DIR_TD1,'clf_test_data.csv')
CSV_CLF_TEST_D2_ASSESS1 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS1,DIR_TD2,'clf_test_data.csv')

DIR_DATA_CLASSIFIER_ASSESS2 = os.path.join(DIR_DATA,'assessment2')
CSV_CLF_TRAIN_DATA_ASSESS2 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS2,'clf_train_data.csv')
CSV_CLF_TEST_D1_ASSESS2 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS2,DIR_TD1,'clf_test_data.csv')
CSV_CLF_TEST_D2_ASSESS2 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS2,DIR_TD2,'clf_test_data.csv')

DIR_DATA_CLASSIFIER_ASSESS2_2 = os.path.join(DIR_DATA,'assessment2_2')
CSV_CLF_TEST_D1_ASSESS2_2 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS2_2,DIR_TD1,'clf_test_data.csv')
CSV_CLF_TEST_D2_ASSESS2_2 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS2_2,DIR_TD2,'clf_test_data.csv')

DIR_DATA_CLASSIFIER_ASSESS3 = os.path.join(DIR_DATA,'assessment3')
CSV_CLF_TRAIN_DATA_ASSESS3 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS3,'clf_train_data.csv')
CSV_CLF_TEST_D1_ASSESS3 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS3,DIR_TD1,'clf_test_data.csv')
CSV_CLF_TEST_D2_ASSESS3 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS3,DIR_TD2,'clf_test_data.csv')

DIR_DATA_CLASSIFIER_ASSESS4 = os.path.join(DIR_DATA,'assessment4')
CSV_CLF_TEST_D1_ASSESS4 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS4,DIR_TD1,'clf_test_data.csv')
CSV_CLF_TEST_D2_ASSESS4 = os.path.join(DIR_DATA_CLASSIFIER_ASSESS4,DIR_TD2,'clf_test_data.csv')

DIR_CLASSIFIER = 'classifier'
KSF2_ASSESS_ROC_CURVES = os.path.join(DIR_CLASSIFIER,'ksf2','roc_curve.png')
KSF2_ASSESS_PR_CURVES = os.path.join(DIR_CLASSIFIER,'ksf2','pr_curve.png')

DIR_EMBEDDINGS = os.path.join(DIR_DATA,'embeddings')
CSV_TRANSE_EMB = os.path.join(DIR_EMBEDDINGS,'transE_emb.csv')
CSV_DISTMULT_EMB = os.path.join(DIR_EMBEDDINGS,'distmult_emb.csv')
CSV_COMPLEX_EMB = os.path.join(DIR_EMBEDDINGS,'complex_emb.csv')
CSV_EXPRESSIVE_EMB = os.path.join(DIR_EMBEDDINGS,'expressivE_emb.csv')
CSV_ESM2_EMB = os.path.join(DIR_EMBEDDINGS,'esm2_emb.csv')
CSV_ESM3_EMB = os.path.join(DIR_EMBEDDINGS,'esm3_emb.csv')
CSV_PHOSFORMER_EMB = os.path.join(DIR_EMBEDDINGS,'phosformer_emb.csv')
CSV_PROTT5_EMB = os.path.join(DIR_EMBEDDINGS,'protT5_emb.csv')
CSV_PROSTT5_EMB = os.path.join(DIR_EMBEDDINGS,'prostT5_emb.csv')
CSV_RANDOM_EMB = os.path.join(DIR_EMBEDDINGS,'random_emb.csv')

CSV_SUBSTRATES_MOTIF = os.path.join(DIR_DATA,'substrates_motif.csv')
CSV_KG_KINASES = os.path.join(DIR_DATA,'kinases.csv')

CSV_PHOSFORMER_PREDICTIONS = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','phosformer_predictions.csv')
PHOS_ST_KINASES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','sup_file_x3_kinases.txt')
CSV_SER_THR_ATLAS = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','SER_THR_atlas.csv')
KSF2_PHOS_ST_ASSESS1_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','subassess1_roc_curve.png')
KSF2_PHOS_ST_ASSESS1_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','subassess1_pr_curve.png')
KSF2_PHOS_ST_TD2_ASSESS1_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','td2','subassess1_roc_curve.png')
KSF2_PHOS_ST_TD2_ASSESS1_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','td2','subassess1_pr_curve.png')
KSF2_PHOS_ST_ASSESS2_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','subassess2_roc_curve.png')
KSF2_PHOS_ST_ASSESS2_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','subassess2_pr_curve.png')
KSF2_PHOS_ST_TD2_ASSESS2_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','td2','subassess2_roc_curve.png')
KSF2_PHOS_ST_TD2_ASSESS2_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'phosformer-ST','td2','subassess2_pr_curve.png')

CSV_LINKPHINDER_PREDICTIONS = os.path.join(DIR_ASSESSMENT4,'link_phinder','linkphinder_predictions.csv')
KSF2_LINKPHINDER_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'link_phinder','roc_curve.png')
KSF2_LINKPHINDER_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'link_phinder','pr_curve.png')
KSF2_LINKPHINDER_TD2_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'link_phinder','td2','roc_curve.png')
KSF2_LINKPHINDER_TD2_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'link_phinder','td2','pr_curve.png')

CSV_KSFINDER_PREDICTIONS = os.path.join(DIR_ASSESSMENT4,'ksfinder','ksfinder_predictions.csv')
KSF2_KSF1_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'ksfinder','roc_curve.png')
KSF2_KSF1_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'ksfinder','pr_curve.png')
KSF2_KSF1_TD2_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'ksfinder','td2','roc_curve.png')
KSF2_KSF1_TD2_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'ksfinder','td2','pr_curve.png')

CSV_PREDKINKG_PREDICTIONS = os.path.join(DIR_ASSESSMENT4,'predkinkg','predkinkg_predictions.csv')
KSF2_PREDKINKG_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'predkinkg','roc_curve.png')
KSF2_PREDKINKG_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'predkinkg','pr_curve.png')
KSF2_PREDKINKG_TD2_ROC_CURVES = os.path.join(DIR_ASSESSMENT4,'predkinkg','td2','roc_curve.png')
KSF2_PREDKINKG_TD2_PR_CURVES = os.path.join(DIR_ASSESSMENT4,'predkinkg','td2','pr_curve.png')

CSV_KSF2_KG = os.path.join(DIR_DATA,'kg','kg_data.csv')
CSV_KSF2_KG_TRAIN = os.path.join(DIR_DATA,'kg','kg_train_data.csv')
CSV_KSF2_KG_TRAIN_VAL = os.path.join(DIR_DATA,'kg','kg_train_val.csv')
CSV_KSF2_KG_VAL = os.path.join(DIR_DATA,'kg','kg_val_data.csv')

DIR_KGE = 'kge'
RESULT_DIR_TRANSE = os.path.join(DIR_KGE,'transE')
RESULT_DIR_DISTMULT = os.path.join(DIR_KGE,'distmult')
RESULT_DIR_COMPLEX = os.path.join(DIR_KGE,'complex')
RESULT_DIR_EXPRESSIVE = os.path.join(DIR_KGE,'expressivE')

DIR_KSF2_PREDICTIONS = 'output'
DIR_KSF2_PREDICTIONS_BATCH = 'output/batch'

