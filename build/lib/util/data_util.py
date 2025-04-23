import random
from util import constants
import pandas as pd
from sklearn.utils import resample
from itertools import product
from util.metrics import Score
import scipy.stats as st
import numpy as np

random.seed(13)

def get_confidence_interval(ci, y_true,y_pred, n_boostrap_cnt=1000):
    
    bootstrap_roc_auc = []
    bootstrap_pr_auc = []

    label_1_indices = [idx for idx,y in enumerate(y_true) if y == 1]
    label_0_indices = [idx for idx,y in enumerate(y_true) if y == 0]

    for _ in range(n_boostrap_cnt):
        sampled_indices = resample(label_1_indices, n_samples=len(label_1_indices), 
                                replace=True)
        sampled_indices.extend(resample(label_0_indices, n_samples=len(label_0_indices), 
                                replace=True))
        sampled_y_true = [y_true[i] for i in sampled_indices]
        sampled_y_pred = [y_pred[i] for i in sampled_indices]

        roc_score, _, _, _ = Score.get_roc_score(sampled_y_true,sampled_y_pred)
        pr_score, _, _, _ = Score.get_pr_score(sampled_y_true,sampled_y_pred)

        bootstrap_roc_auc.append(roc_score)
        bootstrap_pr_auc.append(pr_score)

    roc_mean = round(np.mean(bootstrap_roc_auc),3)
    pr_mean = round(np.mean(bootstrap_pr_auc),3)
    roc_ci_lower, roc_ci_upper = st.norm.interval(confidence=ci, loc=roc_mean, scale=st.sem(bootstrap_roc_auc))
    roc_moe = round((roc_ci_upper - roc_ci_lower),3)
    
    pr_ci_lower, pr_ci_upper = st.norm.interval(confidence=ci, loc=pr_mean, scale=st.sem(bootstrap_pr_auc))
    pr_moe = round((pr_ci_upper - pr_ci_lower),3)

    return f'{roc_mean}+/-{roc_moe}', f'{pr_mean}+/-{pr_moe}'


def normalize_data_cnt(test_df):
    label_cnts = test_df['label'].value_counts()
    pos_label_cnt, neg_label_cnt = label_cnts[1], label_cnts[0]
    if pos_label_cnt > neg_label_cnt:
        select_indices = random.sample(test_df[test_df['label']==1].index.to_list(),neg_label_cnt)
        select_indices.extend(test_df[test_df['label']==0].index.to_list())
        test_df = test_df[test_df.index.isin(select_indices)].copy()
    elif neg_label_cnt > pos_label_cnt:
        select_indices = random.sample(test_df[test_df['label']==0].index.to_list(),pos_label_cnt)
        select_indices.extend(test_df[test_df['label']==1].index.to_list())
        test_df = test_df[test_df.index.isin(select_indices)].copy()
    return test_df

def get_kg_substrate_motifs():
    substrate_motif = pd.read_csv(constants.CSV_SUBSTRATES_MOTIF,sep='|')
    substrate_motif['substrate_motif']=substrate_motif.apply(lambda x:x['Seq_Substrate']+'_'+x['Motif'],axis=1)
    substrate_motifs = substrate_motif['substrate_motif'].unique()
    return substrate_motifs

def process_prediction_output(df):
    df.rename({'head':'kinase','tail':'substrate_motif','ksf_pred':'ksf2_pred'},axis=1,inplace=True)
    df = df[['kinase','substrate_motif','ksf2_pred']].copy()
    df.drop_duplicates(inplace=True)
    return df

def get_kg_kinases():
    kinase_df = pd.read_csv(constants.CSV_KG_KINASES,sep='|')
    return kinase_df['Kinase'].to_list()


def read_kinase_files(file_path):
    contents = []
    with open(file_path) as ip_f:
        for line in ip_f:
            contents.append(line.strip())
    return contents

def _get_kinases_by_type():
    tyr_kinases = read_kinase_files('data/kinases_tyr.txt')
    st_kinases = read_kinase_files('data/kinases_ser_thr.txt')
    other_kinases = read_kinase_files('data/kinases_other.txt')
    return tyr_kinases, st_kinases, other_kinases

def _get_substrate_motifs_by_type(substrate_motifs):
    tyr_substrate_motifs = []
    st_substrate_motifs = []
    for sm in substrate_motifs:
        substrate, motif = sm[:sm.index('_')], sm[sm.index('_')+1:]
        if motif[4] == 'Y':
            tyr_substrate_motifs.append(sm)
        else:
            st_substrate_motifs.append(sm)
    return tyr_substrate_motifs, st_substrate_motifs, substrate_motifs

def _get_valid_ksm_combinations(kinases, substrate_motifs):
    
    tyr_kinases, st_kinases, other_kinases = _get_kinases_by_type()
    tyr_substrate_motifs, st_substrate_motifs, substrate_motifs = _get_substrate_motifs_by_type(substrate_motifs)
    for kinase in kinases:
        if kinase in tyr_kinases:
            if len(tyr_substrate_motifs) == 0: return
            yield zip(*product([kinase], tyr_substrate_motifs))
        elif kinase in st_kinases:
            if len(st_substrate_motifs) == 0: return
            yield zip(*product([kinase], st_substrate_motifs))
        else:
            yield zip(*product([kinase], substrate_motifs))


def group_by_kinase_group(data):
    kinase_group_df = pd.read_csv('data/kinase_group.csv',sep='|')
    k_groups = kinase_group_df.groupby(['group'])['kinase'].apply(list).to_dict()

    for key in k_groups:
         kinases = k_groups[key]
         group_df = data[data['head'].isin(kinases)].copy().reset_index(drop=True)
         yield key, group_df


def group_by_kinase_family(data):
    kinase_group_df = pd.read_csv('data/kinase_group.csv',sep='|')
    kinase_family_df = pd.read_csv('data/kinase_family.csv',sep='|')
    kinase_group_family_df = pd.merge(kinase_group_df,kinase_family_df,on='kinase',how='left')
    k_family_groups = kinase_group_family_df.groupby(['group','family'])['kinase'].apply(list).to_dict()

    for key in k_family_groups:
         kinases = k_family_groups[key]
         group_df = data[data['head'].isin(kinases)].copy().reset_index(drop=True)
         yield key, group_df


def group_by_kinase(data):
    kinase_group_df = pd.read_csv('data/kinase_group.csv',sep='|')
    kinase_family_df = pd.read_csv('data/kinase_family.csv',sep='|')
    kinase_group_family_df = pd.merge(kinase_group_df,kinase_family_df,on='kinase',how='left')
    k_family_groups = kinase_group_family_df.groupby(['group','family','kinase'])['kinase'].apply(list).to_dict()
    
    processed_kinases = set()
    for key in k_family_groups:
         kinases = k_family_groups[key]
         for kinase in kinases:
             if kinase in processed_kinases: 
                 continue
             else:
                processed_kinases.add(kinase)
                group_df = data[data['head']==kinase].copy().reset_index(drop=True)
                group, family, kinase = key[0], key[1], key[2]
                yield group, family, kinase, group_df
    
    k_groups = kinase_group_family_df.groupby(['group','kinase'])['kinase'].apply(list).to_dict()
    for key in k_groups:
         kinases = k_groups[key]
         for kinase in kinases:
             if kinase in processed_kinases: 
                 continue
             else:
                processed_kinases.add(kinase)
                group_df = data[data['head']==kinase].copy().reset_index(drop=True)
                group, family, kinase = key[0], '', key[1]
                yield group, family, kinase, group_df

    test_df = pd.read_csv(constants.CSV_CLF_TEST_D2,sep='|')
    for key in test_df['head'].unique():
        kinase = key
        if kinase in processed_kinases: continue
        processed_kinases.add(kinase)
        group_df = data[data['head']==kinase].copy().reset_index(drop=True)
        group, family, kinase = '', '', key
        yield group, family, kinase, group_df