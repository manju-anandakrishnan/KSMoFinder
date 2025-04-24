echo '(Assessment1) Assessing models developed using different KGE algorithms'

echo 'Evaluating model with transE embeddings'
python assessments/assessment1/nn_classifier_transE.py

echo 'Evaluating model with distmult embeddings'
python assessments/assessment1/nn_classifier_distmult.py

echo 'Evaluating model with complex embeddings'
python assessments/assessment1/nn_classifier_complex.py

echo 'Evaluating model with expressivE embeddings'
python assessments/assessment1/nn_classifier_expressivE.py

# ## To retrain these models, include argument retrain=T. E.g. python assessments/assessment1/nn_classifier_transE.py --retrain=T
echo '-----------------Assessment1 Ended---------------------'

echo '(Assessment2) KSFinder 2.0 KGE vs models developed using embeddings from other models'

echo 'Evaluating model with embeddings from ESM2'
python assessments/assessment2/nn_classifier_esm2.py

echo 'Evaluating model with embeddings from ESM3'
python assessments/assessment2/nn_classifier_esm3.py

echo 'Evaluating model with embeddings from ProtT5'
python assessments/assessment2/nn_classifier_protT5.py

echo 'Evaluating model with embeddings from KSFinder 2.0 KGE (transE)'
python assessments/assessment2/nn_classifier_ksf2.py

echo 'Evaluating model with embeddings from random generation'
python assessments/assessment2/nn_classifier_random.py

echo 'Compare KSMoFinder-KGE with other embeddings after removing easy test scenarios (similar kinase-motif)'
echo 'Evaluating model with embeddings from ESM2'
python assessments/assessment2/nn_classifier_esm2.py --test_data_path ASSESS_2_2

echo 'Evaluating model with embeddings from ESM3'
python assessments/assessment2/nn_classifier_esm3.py --test_data_path ASSESS_2_2

echo 'Evaluating model with embeddings from ProtT5'
python assessments/assessment2/nn_classifier_protT5.py --test_data_path ASSESS_2_2

echo 'Evaluating model with embeddings from KSFinder 2.0 KGE (transE)'
python assessments/assessment2/nn_classifier_ksf2.py --test_data_path ASSESS_2_2

# ## To retrain these models, include argument retrain=T. E.g. python assessments/assessment2/nn_classifier_esm2.py --retrain=T
echo '-----------------Assessment2 Ended---------------------'

echo '(Assessment3) Influence of embeddings from models - Phosformer and ProstT5'

echo 'Evaluating model ProsT5 (structure) and Phosformer (domain sequence & motif)'
python assessments/assessment3/prostT5_phosformer_only/nn_classifier.py

echo 'Evaluating model with KSFinder 2.0 embeddings (function and motif) and ProsT5 (structure)'
python assessments/assessment3/w_prostT5/nn_classifier.py

echo 'Evaluating model with KSFinder 2.0 embeddings (function and motif) and Phosformer (domain sequence & motif)'
python assessments/assessment3/w_phosformer/nn_classifier.py

echo 'Evaluating model with KSFinder 2.0 embeddings (function and motif), ProstT5(structure) and Phosformer (domain sequence & motif)'
python assessments/assessment3/w_prostT5_phosformer/nn_classifier.py

echo 'Evaluating model with KSFinder 2.0 embeddings only (function and motif)'
python assessments/assessment3/ksf2_only/nn_classifier.py

echo 'Evaluating model with KSFinder 2.0 embeddings (function and motif), ProstT5 (structure)'
python assessments/assessment3/w_prostT5_ksf2_motif_only/nn_classifier.py

# ## To retrain these models, include argument retrain=T. E.g. python assessments/assessment3/prostT5_phosformer_only/nn_classifier.py --retrain=T
echo '-----------------Assessment3 Ended---------------------'

echo '(Assessment4) Comparative evaluation of KSFinder 2.0 with other kinase-substrate tools'

echo 'Comparing KSFinder 2.0 and Phosformer-ST (subassessment 1)'
python assessments/assessment4/phosformer-ST/compare_subassess1.py

echo 'Comparing KSFinder 2.0 and Phosformer-ST (subassessment 2)'
python assessments/assessment4/phosformer-ST/compare_subassess2.py

echo 'Comparing KSFinder 2.0 and LinkPhinder'
python assessments/assessment4/link_phinder/compare.py

echo 'Comparing KSFinder 2.0 and ksfinder'
python assessments/assessment4/ksfinder/compare.py

echo 'Comparing KSFinder 2.0 and predkinkg'
python assessments/assessment4/predkinkg/compare.py

echo 'Compare KSFinder 2.0 after removing easy test scenarios (similar kinase-motif)'
echo 'KSFinder 2.0 vs Phosformer-ST (subassessment 1)'
python assessments/assessment4/phosformer-ST/compare_subassess1.py --test_data_path ASSESS4

echo 'KSFinder 2.0 vs Phosformer-ST (subassessment 2)'
python assessments/assessment4/phosformer-ST/compare_subassess2.py --test_data_path ASSESS4

echo 'KSFinder 2.0 vs LinkPhinder'
python assessments/assessment4/link_phinder/compare.py --test_data_path ASSESS4

echo '-----------------Assessment4 Ended---------------------'

echo '(Assessment5) Evaluation at group, family and kinase level'
python assessments/assessment5/evaluate_by_gr_fam_k.py