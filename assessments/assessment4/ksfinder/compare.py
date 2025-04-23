import pandas as pd
import random
import argparse

from classifier.predict.generator import predict
from util.metrics import Curve, Score
from util import constants
from util import data_util

random.seed(13)

if __name__ == '__main__':

    test_df = pd.read_csv(constants.CSV_CLF_TEST_D2,sep='|')
    test_df['motif'] = test_df['tail'].apply(lambda x:x[x.find('_')+1:])
    test_df['substrate'] = test_df['tail'].apply(lambda x:x[:x.find('_')])
    test_df['kinase'] = test_df['head'].copy()
    test_pos_data = test_df[test_df['label']==1].copy()
    neg_data = test_df[test_df['label']==0].copy()
    test_pos_data['ks'] = test_pos_data.apply(lambda x:(x['kinase'],x['substrate']),axis=1)
    pos_ks = test_pos_data['ks'].to_list()
    neg_data['ks'] = neg_data.apply(lambda x:(x['kinase'],x['substrate']),axis=1)

    train_df = pd.read_csv(constants.CSV_CLF_TRAIN_DATA,sep='|')
    train_df['motif'] = train_df['tail'].apply(lambda x:x[x.find('_')+1:])
    train_df['substrate'] = train_df['tail'].apply(lambda x:x[:x.find('_')])
    train_df['kinase'] = train_df['head'].copy()
    train_pos_data = train_df[train_df['label']==1].copy()
    train_pos_data['ks'] = train_pos_data.apply(lambda x:(x['kinase'],x['substrate']),axis=1)
    pos_ks.extend(train_pos_data['ks'].to_list())
    test_neg_data = neg_data[~neg_data['ks'].isin(pos_ks)].copy()

    test_df = pd.concat([test_pos_data,test_neg_data])

    ksf1_pred_df = pd.read_csv(constants.CSV_KSFINDER_PREDICTIONS,sep=',')
    ksf1_pred_df = ksf1_pred_df[['kinase','substrate','prob']]

    test_df = test_df[['head','tail','label','substrate']]
    test_df.dropna(how='any',axis=0,inplace=True)
    test_df.drop_duplicates(inplace=True)

    test_df = test_df.merge(ksf1_pred_df,how='left',left_on=['head','substrate'],right_on=['kinase','substrate'])
    test_df.dropna(axis=0,how='any',inplace=True)
    test_df = test_df[['head','tail','label','prob']].copy()
    test_df.drop_duplicates(inplace=True)

    test_df2 = test_df.copy()

    # Testing dataset 1
    test_df1 = data_util.normalize_data_cnt(test_df)
    print(test_df1['label'].value_counts().to_dict())

    y_true = test_df1['label'].to_list()
    ksf1_pred = test_df1['prob'].to_list()
    test_df1.dropna(axis=0,how='any',inplace=True)

    ksf2_pred = predict(test_df1)
    test_df1['ksf_pred'] = ksf2_pred
    test_df1.dropna(axis=0,how='any',inplace=True)

    roc_curve = Curve.get_roc_curves([y_true, y_true],
                                        [ksf1_pred,ksf2_pred],
                                        ['blue','magenta'],)
    pr_curve = Curve.get_pr_curves([y_true, y_true],[ksf1_pred,ksf2_pred],
                                    ['blue','magenta'],)
        
    roc_score, _,_,_ = Score.get_roc_score(y_true,ksf2_pred)
    pr_score, _,_,_ = Score.get_pr_score(y_true,ksf2_pred)
    print(f'KSFinder 2.0:: ROC-AUC: {roc_score} | PR-AUC: {pr_score} | CI (95%):: {data_util.get_confidence_interval(0.95,y_true,ksf2_pred)}')

    roc_score, _,_,_ = Score.get_roc_score(y_true,ksf1_pred)
    pr_score, _,_,_ = Score.get_pr_score(y_true,ksf1_pred)
    print(f'KSFinder:: ROC-AUC: {roc_score} | PR-AUC: {pr_score} | CI (95%):: {data_util.get_confidence_interval(0.95,y_true,ksf1_pred)}')

    roc_curve.savefig(constants.KSF2_KSF1_ROC_CURVES)
    pr_curve.savefig(constants.KSF2_KSF1_PR_CURVES)

    # Testing dataset 2
    print(test_df2['label'].value_counts().to_dict())

    y_true = test_df2['label'].to_list()
    ksf1_pred = test_df2['prob'].to_list()
    test_df2.dropna(axis=0,how='any',inplace=True)

    ksf2_pred = predict(test_df2)
    test_df2['ksf_pred'] = ksf2_pred
    test_df2.dropna(axis=0,how='any',inplace=True)

    roc_curve = Curve.get_roc_curves([y_true, y_true],
                                        [ksf1_pred,ksf2_pred],
                                        ['blue','magenta'],)
    pr_curve = Curve.get_pr_curves([y_true, y_true],[ksf1_pred,ksf2_pred],
                                    ['blue','magenta'],)
        
    roc_score, _,_,_ = Score.get_roc_score(y_true,ksf2_pred)
    pr_score, _,_,_ = Score.get_pr_score(y_true,ksf2_pred)
    print(f'KSFinder 2.0:: ROC-AUC: {roc_score} | PR-AUC: {pr_score} | CI (95%):: {data_util.get_confidence_interval(0.95,y_true,ksf2_pred)}')
        
    roc_score, _,_,_ = Score.get_roc_score(y_true,ksf1_pred)
    pr_score, _,_,_ = Score.get_pr_score(y_true,ksf1_pred)
    print(f'KSFinder:: ROC-AUC: {roc_score} | PR-AUC: {pr_score} | CI (95%):: {data_util.get_confidence_interval(0.95,y_true,ksf1_pred)}')

    roc_curve.savefig(constants.KSF2_KSF1_TD2_ROC_CURVES)
    pr_curve.savefig(constants.KSF2_KSF1_TD2_PR_CURVES)