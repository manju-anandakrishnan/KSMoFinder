import pandas as pd
import numpy as np

import torch

from sklearn.preprocessing import MinMaxScaler
from util import constants

torch.manual_seed(13)
np.random.seed(13)

class DataLoader:

    def __init__(self):
        pass

    def _load_motif_11mer(self,subs_motif_df):
        kg2_substrates = pd.read_csv(constants.CSV_SUBSTRATES_MOTIF,sep='|')
        kg2_substrates_unique = kg2_substrates.drop_duplicates(subset=['Seq_Substrate', 'Motif'])
        subs_motif_df = subs_motif_df.merge(kg2_substrates_unique,how='left',left_on=['substrate1','motif1'],right_on=['Seq_Substrate','Motif'])
        motif_11mer_list = subs_motif_df['motif_15mer'].to_list()
        return motif_11mer_list
    
    def _preprocess(self,df):
        df.drop_duplicates(inplace=True)
        df['kinase1'] = df['head'].copy()
        df['substrate1'] = df['tail'].apply(lambda x:x[:x.find('_')])
        df['motif1'] = df['tail'].apply(lambda x:x[x.find('_')+1:])
        df['kinase_substrate'] = df.apply(lambda x:(x['kinase1'],x['substrate1']),axis=1)
        df['kinase_motif'] = df.apply(lambda x:(x['kinase1'],x['motif1']),axis=1)
        df['kinase_struct'] = df['head'].copy()
        df['substrate_struct'] = df['substrate1'].copy()   
        df['kinase_domain'] = df['head'].copy()
        df['motif_11mer'] = self._load_motif_11mer(df)     
        return df

    def get_training_data(self,train_data):
        df = pd.read_csv(train_data,sep='|')
        return self._preprocess(df)
    
    def get_testing_data(self,test_data):
        df = pd.read_csv(test_data,sep='|')
        return self._preprocess(df)

class KSMEmbeddings:

    def __init__(self):
        self.entity_embedding_dict = self._load_embeddings()
        self.prostT5_embedding_dict = self._load_prostT5_embeddings()
        self.phosformer_embedding_dict = self._load_phosformer_embeddings()
        self.features = ['kinase1','substrate1','motif1','kinase_struct','substrate_struct','kinase_domain','motif_11mer','label']

    def _scale_data(self, df):
        scaler = MinMaxScaler(feature_range=(-0.4,0.4))
        curr_index = df.index.to_list()
        scaled_df = pd.DataFrame(scaler.fit_transform(df),index=curr_index)
        return scaled_df

    def _load_embeddings(self):
        emb_df = pd.read_csv(constants.CSV_TRANSE_EMB,sep='|')
        emb_df.set_index('entity', inplace=True)
        emb_df = self._scale_data(emb_df)
        entity_embedding_dict = emb_df.apply(lambda row: row.values, axis=1).to_dict()
        return entity_embedding_dict

    def _load_prostT5_embeddings(self):
        emb_df = pd.read_csv(constants.CSV_PROSTT5_EMB,sep='|')
        emb_df.set_index('Entity', inplace=True)
        emb_df = self._scale_data(emb_df)
        entity_embedding_dict = emb_df.apply(lambda row: row.values, axis=1).to_dict()
        return entity_embedding_dict
    
    def _load_phosformer_embeddings(self):
        emb_df = pd.read_csv(constants.CSV_PHOSFORMER_EMB,sep='|')
        emb_df.set_index('entity', inplace=True)
        emb_df = self._scale_data(emb_df)
        entity_embedding_dict = emb_df.apply(lambda row: row.values, axis=1).to_dict()
        return entity_embedding_dict

    def _map_embeddings(self,df):
        for df_col in df.columns:
            if df_col == 'label': continue
            if df_col in ['kinase_struct','substrate_struct']: 
                df[df_col] = df[df_col].apply(lambda x:self.prostT5_embedding_dict.get(x))
            elif df_col in ['kinase_domain','motif_11mer']: 
                df[df_col] = df[df_col].apply(lambda x:self.phosformer_embedding_dict.get(x))
            else:
                df[df_col] = df[df_col].apply(lambda x:self.entity_embedding_dict.get(x))
        df.dropna(axis=0,how='any',inplace=True)  
        return tuple((pd.DataFrame(df[x].to_list())).to_numpy() for x in self.features)

    def get_training_data(self,raw_training_data):
        self._training_data = raw_training_data
        training_data = self._training_data[self.features].copy()
        emb_tup = self._map_embeddings(training_data)
        label = emb_tup[-1]
        print('Training data count::',label.shape[0],
              '| Positive data count::',label[label[:,-1]==1].shape[0],
              '| Negative data count::',label[label[:,-1]==0].shape[0])
        return emb_tup
    
    def get_testing_data(self,raw_testing_data):
        self._testing_data = raw_testing_data
        testing_data = self._testing_data[self.features].copy()
        emb_tup = self._map_embeddings(testing_data)        
        label = emb_tup[-1]
        print('Testing data count::',label.shape[0],
            '| Positive data count::',label[label[:,-1]==1].shape[0],
            '| Negative data count::',label[label[:,-1]==0].shape[0])
        return emb_tup
