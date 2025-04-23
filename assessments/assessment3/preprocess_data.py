import os
import pandas as pd

from util import constants

def load_motif_15mer(subs_motif_df):
    kg2_substrates = pd.read_csv(constants.CSV_SUBSTRATES_MOTIF,sep='|')
    kg2_substrates_unique = kg2_substrates.drop_duplicates(subset=['Seq_Substrate', 'Motif'])
    subs_motif_df = subs_motif_df.merge(kg2_substrates_unique,how='left',left_on=['substrate1','motif1'],right_on=['Seq_Substrate','Motif'])
    motif_15mer_list = subs_motif_df['motif_15mer'].to_list()
    return motif_15mer_list

def preprocess(df):
    df.drop_duplicates(inplace=True)
    df['kinase1'] = df['head'].copy()
    df['substrate1'] = df['tail'].apply(lambda x:x[:x.find('_')])
    df['motif1'] = df['tail'].apply(lambda x:x[x.find('_')+1:])
    df['kinase_struct'] = df['head'].copy()
    df['substrate_struct'] = df['substrate1'].copy()   
    df['kinase_domain'] = df['head'].copy()
    df['motif_11mer'] = load_motif_15mer(df)     
    return df

def _load_ksf_embeddings():
    emb_df = pd.read_csv(constants.CSV_TRANSE_EMB,sep='|')
    emb_df.set_index('entity', inplace=True)
    entity_embedding_dict = emb_df.apply(lambda row: row.values, axis=1).to_dict()
    return entity_embedding_dict

def _load_prostT5_embeddings():
    emb_df = pd.read_csv(constants.CSV_PROSTT5_EMB,sep='|')
    emb_df.set_index('Entity', inplace=True)
    entity_embedding_dict = emb_df.apply(lambda row: row.values, axis=1).to_dict()
    return entity_embedding_dict
    
def _load_phosformer_embeddings():
    emb_df = pd.read_csv(constants.CSV_PHOSFORMER_EMB,sep='|')
    emb_df.set_index('entity', inplace=True)
    entity_embedding_dict = emb_df.apply(lambda row: row.values, axis=1).to_dict()
    return entity_embedding_dict

if __name__ == '__main__':
    
    prostT5_entities = _load_prostT5_embeddings().keys()
    phosformer_entities = _load_phosformer_embeddings().keys()
    ksf_entities = _load_ksf_embeddings().keys()

    train_df = pd.read_csv(constants.CSV_CLF_TRAIN_DATA,sep='|')
    train_df = preprocess(train_df)
    train_df = train_df[(train_df['kinase1'].isin(ksf_entities))
                        & (train_df['substrate1'].isin(ksf_entities))
                        & (train_df['motif1'].isin(ksf_entities))
                        & (train_df['kinase_struct'].isin(prostT5_entities))
                        & (train_df['substrate_struct'].isin(prostT5_entities))
                        & (train_df['kinase_domain'].isin(phosformer_entities))
                        & (train_df['motif_11mer'].isin(phosformer_entities))
                        ].copy()
    train_df = train_df[['head','tail','label']]

    test_df2 = pd.read_csv(constants.CSV_CLF_TEST_D2,sep='|')    
    test_df2 = preprocess(test_df2)
    test_df2 = test_df2[(test_df2['kinase1'].isin(ksf_entities))
                        & (test_df2['substrate1'].isin(ksf_entities))
                        & (test_df2['motif1'].isin(ksf_entities))
                        & (test_df2['kinase_struct'].isin(prostT5_entities))
                        & (test_df2['substrate_struct'].isin(prostT5_entities))
                        & (test_df2['kinase_domain'].isin(phosformer_entities))
                        & (test_df2['motif_11mer'].isin(phosformer_entities))
                        ].copy()
    
    test_df2 = test_df2[['head','tail','label']]

    test_df1 = pd.read_csv(constants.CSV_CLF_TEST_D1,sep='|')    
    test_df1 = preprocess(test_df1)
    test_df1 = test_df1[(test_df1['kinase1'].isin(ksf_entities))
                        & (test_df1['substrate1'].isin(ksf_entities))
                        & (test_df1['motif1'].isin(ksf_entities))
                        & (test_df1['kinase_struct'].isin(prostT5_entities))
                        & (test_df1['substrate_struct'].isin(prostT5_entities))
                        & (test_df1['kinase_domain'].isin(phosformer_entities))
                        & (test_df1['motif_11mer'].isin(phosformer_entities))
                        ].copy()
    
    test_df1 = test_df1[['head','tail','label']]
    
    print(f"Training data shape:{train_df['label'].value_counts().to_dict()} \
          | Testing data 2 shape:{test_df2['label'].value_counts().to_dict()} \
          | Testing data 1 shape:{test_df1['label'].value_counts().to_dict()}")
    
    train_df.to_csv(constants.CSV_CLF_TRAIN_DATA_ASSESS3,index=False,sep='|')
    test_df1.to_csv(constants.CSV_CLF_TEST_D1_ASSESS3,index=False,sep='|')
    test_df2.to_csv(constants.CSV_CLF_TEST_D2_ASSESS3,index=False,sep='|')
