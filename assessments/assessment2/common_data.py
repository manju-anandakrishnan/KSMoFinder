import pandas as pd
import random
from util import constants

random.seed(13)

if __name__ == '__main__':
    train_data_csv = constants.CSV_CLF_TRAIN_DATA
    test_data_csv = constants.CSV_CLF_TEST_D2
        
    protT5_df = pd.read_csv(constants.CSV_PROTT5_EMB,sep='|')
    protT5_df = protT5_df[['entity']]
    
    esm2_df = pd.read_csv(constants.CSV_ESM2_EMB,sep='|')
    esm2_df = esm2_df[['entity']]
    
    esm3_df = pd.read_csv(constants.CSV_ESM3_EMB,sep='|')
    esm3_df = esm3_df[['entity']]

    kge_df = pd.read_csv(constants.CSV_TRANSE_EMB,sep='|')
    kge_df = kge_df[['entity']]

    train_data_df = pd.read_csv(train_data_csv,sep='|')
    test_data_df = pd.read_csv(test_data_csv,sep='|')
    clf_data_df = pd.concat([train_data_df,test_data_df],axis=0)
    clf_data_df['substrate'] = clf_data_df['tail'].apply(lambda x:x[:x.find('_')])
    clf_data_df['motif'] = clf_data_df['tail'].apply(lambda x:x[x.find('_')+1:])

    clf_data_df = clf_data_df.merge(protT5_df,how='inner',left_on='tail',right_on='entity')
    clf_data_df = clf_data_df.merge(protT5_df,how='inner',left_on='substrate',right_on='entity')
    clf_data_df = clf_data_df.merge(protT5_df,how='inner',left_on='head',right_on='entity')
    clf_data_df.dropna(axis=0,how='any',inplace=True)
    clf_data_df = clf_data_df[['head','tail','label','substrate','motif']]
    print(f'Merged ProtT5::{clf_data_df.shape}')
    print(clf_data_df['label'].value_counts())

    clf_data_df = clf_data_df.merge(esm3_df,how='inner',left_on='substrate',right_on='entity')
    clf_data_df = clf_data_df.merge(esm3_df,how='inner',left_on='motif',right_on='entity')
    clf_data_df = clf_data_df.merge(esm3_df,how='inner',left_on='head',right_on='entity')
    clf_data_df.dropna(axis=0,how='any',inplace=True)
    clf_data_df = clf_data_df[['head','tail','label','substrate','motif']]
    print(f'Merged ESM3::{clf_data_df.shape}')
    print(clf_data_df['label'].value_counts())

    clf_data_df = clf_data_df.merge(esm2_df,how='inner',left_on='substrate',right_on='entity')
    clf_data_df = clf_data_df.merge(esm2_df,how='inner',left_on='motif',right_on='entity')
    clf_data_df = clf_data_df.merge(esm2_df,how='inner',left_on='head',right_on='entity')
    clf_data_df.dropna(axis=0,how='any',inplace=True)
    clf_data_df = clf_data_df[['head','tail','label','substrate','motif']]
    print(f'Merged ESM2::{clf_data_df.shape}')
    print(clf_data_df['label'].value_counts())

    clf_data_df = clf_data_df.merge(kge_df,how='inner',left_on='substrate',right_on='entity')
    clf_data_df = clf_data_df.merge(kge_df,how='inner',left_on='motif',right_on='entity')
    clf_data_df = clf_data_df.merge(kge_df,how='inner',left_on='head',right_on='entity')
    clf_data_df.dropna(axis=0,how='any',inplace=True)
    clf_data_df = clf_data_df[['head','tail','label']]
    print(f'Merged KGE::{clf_data_df.shape}')

    clf_data_df = clf_data_df.sample(random_state=13, frac=1).reset_index(drop=True)
    pos_data_df = clf_data_df[clf_data_df['label']==1].copy().reset_index(drop=True)
    neg_data_df = clf_data_df[clf_data_df['label']==0].copy().reset_index(drop=True)

    pos_test_indices = list()
    neg_test_indices_td1 = list()
    neg_test_indices_td2 = list()
    pos_kinase_cnt = pos_data_df['head'].value_counts().to_dict()
    for kinase, cnt in pos_kinase_cnt.items():
        p_test_size = int(cnt*0.2)
        if p_test_size > 0:
            kinase_indices = pos_data_df[pos_data_df['head'] == kinase].index.to_list()
            pos_test_indices.extend(random.sample(kinase_indices,p_test_size))
        
        n_kinase_indices = neg_data_df[neg_data_df['head'] == kinase].index.to_list()        
        test_n_kinase_indices = []
        # TD1 (1:1 ratio)
        n_test_size = p_test_size
        if n_test_size < len(n_kinase_indices):
            test_n_kinase_indices = random.sample(n_kinase_indices,n_test_size)
            neg_test_indices_td1.extend(test_n_kinase_indices)
        
        # TD2 (similar dist as training)
        neg_test_indices_td2.extend(test_n_kinase_indices)
        n_test_size = int(len(n_kinase_indices)*0.2)
        diff = n_test_size-len(test_n_kinase_indices)
        if (diff > 0) & (diff < len(n_kinase_indices)):
            neg_test_indices_td2.extend(random.sample(n_kinase_indices,diff))

    pos_train_data_df = pos_data_df[~pos_data_df.index.isin(pos_test_indices)]
    pos_test_data_df = pos_data_df[pos_data_df.index.isin(pos_test_indices)]

    neg_train_data_df = neg_data_df[~neg_data_df.index.isin(neg_test_indices_td2)]
    neg_test_data_df2 = neg_data_df[neg_data_df.index.isin(neg_test_indices_td2)].copy().reset_index(drop=True)
    neg_test_data_df1 = neg_data_df[neg_data_df.index.isin(neg_test_indices_td1)].copy().reset_index(drop=True)
    
    test_data_df2 = pd.concat([pos_test_data_df,neg_test_data_df2])
    test_data_df1 = pd.concat([pos_test_data_df,neg_test_data_df1])
    train_data_df = pd.concat([pos_train_data_df,neg_train_data_df])
    
    print(train_data_df['label'].value_counts().to_dict())
    print(test_data_df1['label'].value_counts().to_dict())
    print(test_data_df2['label'].value_counts().to_dict())

    train_data_df.to_csv(constants.CSV_CLF_TRAIN_DATA_ASSESS2, sep='|',index=False)
    test_data_df1.to_csv(constants.CSV_CLF_TEST_D1_ASSESS2, sep='|',index=False)
    test_data_df2.to_csv(constants.CSV_CLF_TEST_D2_ASSESS2, sep='|',index=False)