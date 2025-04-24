import pandas as pd
import random

from classifier.predict.generator import predict
from util.metrics import Curve, Score
from util import constants
from util import data_util

random.seed(13)

def _by_kinase_groups(test_df):
    for key, test_data in data_util.group_by_kinase_group(test_df):
        test_df1 = test_data.copy()
        y_true = test_df1['label'].to_list()
        ksf_pred = test_df1['ksf_pred'].to_list()
        roc_score, _,_,_ = Score.get_roc_score(y_true,ksf_pred)
        pr_score, _,_,_ = Score.get_pr_score(y_true,ksf_pred)
        label_cnts = test_df1['label'].value_counts().to_dict()
        unique_kinases = list(test_df1['head'].unique())
        print(f"{key}|(Positive: {label_cnts[1]}; Negative: {label_cnts[0]})|{len(unique_kinases)}|{','.join(unique_kinases)}|{roc_score}|{pr_score}")
        # for k in unique_kinases:
        #     print(k)


def _by_kinase_families(test_df):
    for key, test_data in data_util.group_by_kinase_family(test_df):
        test_df1 = test_data.copy()
        label_cnts = test_df1['label'].value_counts().to_dict()
        if (label_cnts.get(1) is None) or (label_cnts.get(0) is None): 
            #print(f"{key[0]}|{key[1]}|Skipped")
            continue
        y_true = test_df1['label'].to_list()
        ksf_pred = test_df1['ksf_pred'].to_list()
        roc_score, _,_,_ = Score.get_roc_score(y_true,ksf_pred)
        pr_score, _,_,_ = Score.get_pr_score(y_true,ksf_pred)
        
        unique_kinases = list(test_df1['head'].unique())
        print(f"{key[0]}|{key[1]}|(Positive: {label_cnts[1]}; Negative: {label_cnts[0]})|{len(unique_kinases)}|{','.join(unique_kinases)}|{roc_score}|{pr_score}")


def _by_kinase(test_df):
    skipped_kinases = set()
    for group, family, kinase, test_data in data_util.group_by_kinase(test_df):
        test_df1 = test_data.copy()
        label_cnts = test_df1['label'].value_counts().to_dict()
        if (label_cnts.get(1) is None) or (label_cnts.get(0) is None): 
            #print(f"{group}|{family}|{kinase}|Skipped")
            skipped_kinases.add(kinase)
            continue
        y_true = test_df1['label'].to_list()
        ksf_pred = test_df1['ksf_pred'].to_list()
        roc_score, _,_,_ = Score.get_roc_score(y_true,ksf_pred)
        pr_score, _,_,_ = Score.get_pr_score(y_true,ksf_pred)
        
        unique_kinases = list(test_df1['head'].unique())
        print(f"{group}|{family}|(Positive: {label_cnts[1]}; Negative: {label_cnts[0]})|{','.join(unique_kinases)}|{roc_score}|{pr_score}")
    else:
        #print(",".join(skipped_kinases))
        pass


def _avg_by_kinase(test_df,min_positives=5):
    skipped_kinases = set()
    sum_roc = cnt_roc = sum_pr = cnt_pr = 0
    min_positive_kinases = list()
    
    for group, family, kinase, test_data in data_util.group_by_kinase(test_df):
        test_df1 = test_data.copy()
        label_cnts = test_df1['label'].value_counts().to_dict()
        if (label_cnts.get(1) is None) or (label_cnts.get(0) is None): 
            #print(f"{group}|{family}|{kinase}|Skipped")
            skipped_kinases.add(kinase)
            continue

        y_true = test_df1['label'].to_list()
        ksf_pred = test_df1['ksf_pred'].to_list()
        roc_score, _,_,_ = Score.get_roc_score(y_true,ksf_pred)
        pr_score, _,_,_ = Score.get_pr_score(y_true,ksf_pred)
        
        if label_cnts[1] < min_positives:
            sum_roc += roc_score
            cnt_roc += 1
            sum_pr += pr_score
            cnt_pr +=1
            min_positive_kinases.extend(list(test_df1['head'].unique()))
        else:
            unique_kinases = list(test_df1['head'].unique())
            print(f"{group}|{family}|(Positive: {label_cnts[1]}; Negative: {label_cnts[0]})|{','.join(unique_kinases)}|{roc_score}|{pr_score}")
    else:
        #print(",".join(skipped_kinases))
        pass
    avg_roc = sum_roc/cnt_roc
    avg_pr = sum_pr/cnt_pr
    print(f"||NA|{','.join(min_positive_kinases)}|{avg_roc}|{avg_pr}")


def _avg_by_family(test_df, min_positives = 5):
    min_positive_families = list()
    min_positive_kinases = list()
    sum_roc = cnt_roc = sum_pr = cnt_pr = 0
    for key, test_data in data_util.group_by_kinase_family(test_df):
        test_df1 = test_data.copy()
        label_cnts = test_df1['label'].value_counts().to_dict()
        if (label_cnts.get(1) is None) or (label_cnts.get(0) is None): 
            #print(f"{key[0]}|{key[1]}|Skipped")
            continue
        y_true = test_df1['label'].to_list()
        ksf_pred = test_df1['ksf_pred'].to_list()
        roc_score, _,_,_ = Score.get_roc_score(y_true,ksf_pred)
        pr_score, _,_,_ = Score.get_pr_score(y_true,ksf_pred)
        group = key[0]
        family = key[1]
        
        if label_cnts[1] < min_positives:
            sum_roc += roc_score
            cnt_roc += 1
            sum_pr += pr_score
            cnt_pr +=1
            min_positive_families.append(family)
            min_positive_kinases.extend(test_df1['head'].unique())
        else:
            unique_kinases = list(test_df1['head'].unique())
            print(f"{group}|{family}|(Positive: {label_cnts[1]}; Negative: {label_cnts[0]})|{len(unique_kinases)}|{','.join(unique_kinases)}|{roc_score}|{pr_score}")
    else:
        avg_roc = sum_roc/cnt_roc
        avg_pr = sum_pr/cnt_pr
        print(f"||NA|{','.join(min_positive_families)}|{','.join(min_positive_kinases)}|{avg_roc}|{avg_pr}")

if __name__ == '__main__':
    test_df = pd.read_csv(constants.CSV_CLF_TEST_D2,sep='|')
    test_df['motif'] = test_df['tail'].apply(lambda x:x[x.find('_')+1:])
    test_df['substrate'] = test_df['tail'].apply(lambda x:x[:x.find('_')])

    substrates_df = pd.read_csv(constants.CSV_SUBSTRATES_MOTIF,sep='|')
    test_df = test_df.merge(substrates_df,how='left',left_on=['substrate','motif'],right_on=['Seq_Substrate','Motif'])
    
    test_df = test_df[['head','tail','label','substrate','Site']]
    test_df = test_df[['head','tail','label']].copy()
    test_df['ksf_pred'] = predict(test_df)
    test_df.drop_duplicates(inplace=True)

    print('********* By kinase group **********')
    _by_kinase_groups(test_df)
    print('********* By kinase family **********')
    _by_kinase_families(test_df)
    print('********* By kinases **********')
    _by_kinase(test_df)
    print('********* Averaged by kinases **********')
    _avg_by_kinase(test_df)
    print('********* Averaged by kinase families **********')
    _avg_by_family(test_df)
