import pandas as pd
import random
import argparse

from classifier.predict.generator import predict
from util.metrics import Curve, Score
from util import constants
from util import data_util

random.seed(13)

def load_data():
    train_df = pd.read_csv(constants.CSV_CLF_TRAIN_DATA,sep='|')
    train_df['motif'] = train_df['tail'].apply(lambda x:x[x.find('_')+1:])
    train_df['substrate'] = train_df['tail'].apply(lambda x:x[:x.find('_')])
    train_df['km'] = train_df.apply(lambda x:(x['head'],x['motif']),axis=1)
    tr_pos_data = train_df[train_df['label']==1].copy()
    tr_pos_data_km = tr_pos_data['km'].unique()

    test_df = pd.read_csv(constants.CSV_CLF_TEST_D2,sep='|')
    test_df['motif'] = test_df['tail'].apply(lambda x:x[x.find('_')+1:])
    test_df['substrate'] = test_df['tail'].apply(lambda x:x[:x.find('_')])
    test_df['km'] = test_df.apply(lambda x:(x['head'],x['motif']),axis=1)
    test_pos_data = test_df[test_df['label']==1].copy()
    test_pos_data_km = test_pos_data['km'].unique()

    tr_pos_data_km = set(tr_pos_data_km)
    tr_pos_data_km.update(set(test_pos_data_km))

    print(len(tr_pos_data_km))

    with open('data/phosphorylation_km_pairs.txt') as ip_f:
        for line in ip_f:
            k,m = line.strip().split(':')
            tr_pos_data_km.add((k,m))
            
    print(len(tr_pos_data_km))
    return tr_pos_data_km


if __name__ == '__main__':

    tr_pos_data_km = load_data()

    parser = argparse.ArgumentParser()
    parser.add_argument('--test_data_path', type=str,default=False)
    args = parser.parse_args()

    test_data_path = constants.CSV_CLF_TEST_D2

    if (args.test_data_path) and (args.test_data_path == 'ASSESS4'):
        test_data_path = constants.CSV_CLF_TEST_D2_ASSESS4
    
    test_data = pd.read_csv(test_data_path,sep='|')
    test_data['htl'] = test_data.apply(lambda x:x['head']+x['tail']+str(x['label']),axis=1)
    test_data_htl = test_data['htl'].unique()

    test_df = pd.read_csv(constants.CSV_PHOSFORMER_PREDICTIONS,sep='|')
    test_df['htl'] = test_df.apply(lambda x:x['head']+x['tail']+str(x['label']),axis=1)
    test_df = test_df[test_df['htl'].isin(test_data_htl)].copy()

    test_df.drop(columns=['htl'],axis=1,inplace=True)

    pos_test_df = test_df[test_df['label']==1].copy()
    neg_test_df = test_df[test_df['label']==0].copy()
    
    print("Filtering of 300 serine-threonine kinases used in Phosformer-ST's work")
    phos_comp_kinases = None
    with open(constants.PHOS_ST_KINASES) as ip_f:
        phos_comp_kinases =  [x.strip() for x in ip_f]
    if phos_comp_kinases:
        test_df = test_df[test_df['head'].isin(phos_comp_kinases)].copy()

    pos_test_df = test_df[test_df['label']==1].copy()
    neg_test_df = test_df[test_df['label']==0].copy()
    print('Total test sample count::',pos_test_df.shape, neg_test_df.shape)

    test_df = pd.concat([pos_test_df,neg_test_df])
    test_df = test_df[['head','tail','label','motif','substrate','motif_15mer','phosST_pred']].copy()

    na_df = test_df[test_df.isnull().any(axis=1)].copy()
    
    test_df.drop_duplicates(inplace=True)
    test_df.dropna(axis=0,how='any',inplace=True)

    test_df['km'] = test_df.apply(lambda x:(x['head'],x['motif']),axis=1)
    pos_data = test_df[test_df['label']==1].copy()
    pos_data_km = pos_data['km'].unique()

    neg_data = test_df[test_df['label']==0].copy()
    neg_km_indices = neg_data[(neg_data['km'].isin(pos_data_km))
                          | (neg_data['km'].isin(tr_pos_data_km))].index.to_list()
    neg_data = neg_data[~neg_data.index.isin(neg_km_indices)].copy()
    
    test_df = pd.concat([pos_data,neg_data])
    test_df.drop_duplicates(inplace=True)
    test_df.dropna(axis=0,how='any',inplace=True)
    test_df.reset_index(drop=True,inplace=True)    

    ksf_pred = predict(test_df)
    test_df['ksf_pred'] = ksf_pred
    test_df.dropna(axis=0,how='any',inplace=True)

    # Testing dataset 2
    test_df2 = test_df.copy()

    # Testing dataset 1
    test_df1 = data_util.normalize_data_cnt(test_df)
    print(f"Total:{test_df1.shape[0]}....{test_df1['label'].value_counts().to_dict()}")
    
    y_true = test_df1['label'].to_list()
    phosST_pred = test_df1['phosST_pred'].to_list()
    ksf_pred = test_df1['ksf_pred'].to_list()

    roc_curve = Curve.get_roc_curves([y_true, y_true],
                                     [phosST_pred,ksf_pred],
                                     ['blue','magenta'],)
    pr_curve = Curve.get_pr_curves([y_true, y_true],[phosST_pred,ksf_pred],
                                   ['blue','magenta'],)
    
    roc_score, _,_,_ = Score.get_roc_score(y_true,ksf_pred)
    pr_score, _,_,_ = Score.get_pr_score(y_true,ksf_pred)
    print(f'KSFinder 2.0:: ROC-AUC: {roc_score} | PR-AUC: {pr_score}')
    print(f'KSFinder 2.0: CI (95%):: {data_util.get_confidence_interval(0.95,y_true,ksf_pred)}')
    
    roc_score, _,_,_ = Score.get_roc_score(y_true,phosST_pred)
    pr_score, _,_,_ = Score.get_pr_score(y_true,phosST_pred)
    print(f'Phosformer-ST:: ROC-AUC: {roc_score} | PR-AUC: {pr_score}')
    print(f'Phosformer-ST: CI (95%):: {data_util.get_confidence_interval(0.95,y_true,phosST_pred)}')

    roc_curve.savefig(constants.KSF2_PHOS_ST_ASSESS1_ROC_CURVES)
    pr_curve.savefig(constants.KSF2_PHOS_ST_ASSESS1_PR_CURVES)

    # Testing dataset 2 evaluation
    print(f"Total:{test_df2.shape[0]}....{test_df2['label'].value_counts().to_dict()}")
    
    y_true = test_df2['label'].to_list()
    phosST_pred = test_df2['phosST_pred'].to_list()
    ksf_pred = test_df2['ksf_pred'].to_list()

    roc_curve = Curve.get_roc_curves([y_true, y_true],
                                     [phosST_pred,ksf_pred],
                                     ['blue','magenta'],)
    pr_curve = Curve.get_pr_curves([y_true, y_true],[phosST_pred,ksf_pred],
                                   ['blue','magenta'],)
    
    roc_score, _,_,_ = Score.get_roc_score(y_true,ksf_pred)
    pr_score, _,_,_ = Score.get_pr_score(y_true,ksf_pred)
    print(f'KSFinder 2.0:: ROC-AUC: {roc_score} | PR-AUC: {pr_score}')
    print(f'KSFinder 2.0: CI (95%):: {data_util.get_confidence_interval(0.95,y_true,ksf_pred)}')

    roc_score, _,_,_ = Score.get_roc_score(y_true,phosST_pred)
    pr_score, _,_,_ = Score.get_pr_score(y_true,phosST_pred)
    print(f'Phosformer-ST:: ROC-AUC: {roc_score} | PR-AUC: {pr_score}')
    print(f'Phosformer-ST: CI (95%):: {data_util.get_confidence_interval(0.95,y_true,phosST_pred)}')

    roc_curve.savefig(constants.KSF2_PHOS_ST_TD2_ASSESS1_ROC_CURVES)
    pr_curve.savefig(constants.KSF2_PHOS_ST_TD2_ASSESS1_PR_CURVES)