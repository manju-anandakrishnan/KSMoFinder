import pandas as pd
import numpy as np
import os
import random
from util import constants

random.seed(13)   

class KGTriples:

    _instance = None
    
    VAL_RELATIONS = ['k_specific_motif', 'participating_pathway', 'part_of_complex', \
                                'bio_process', 'mol_func', 'cellular_comp', 'expressed_in', 'has_domain',\
                                    'belongs_to_family','homologous_superfamily',\
                                        'residue_1','residue_2','residue_3',
                                        'residue_4','residue_5','residue_6',
                                        'residue_7','residue_8','residue_9']

    def __init__(self):
        self._load_kg_data()

    def _load_kg_data(self):
        kg_df = pd.read_csv(constants.CSV_KSF2_KG,sep='|')
        phos_rel_kg = kg_df[kg_df['rel']=='phosphorylates(site)'].copy()
        phos_rel_kg['kinase_motif'] = phos_rel_kg.apply(lambda x:(x['head'],x['tail'][x['tail'].find('_')+1:]),axis=1)
        phos_rel_km = phos_rel_kg['kinase_motif'].unique()
        kg_df = kg_df[~kg_df['rel'].isin(['phosphorylates(site)','interacts_with','psite_motif','psite_of_substrate'])].copy()
        print(kg_df.shape)
        kg_df = kg_df[~kg_df.apply(
            lambda x:(x['head'],x['tail']) 
            if x['rel']=='k_specific_motif'
            else None,axis=1).isin(phos_rel_km)].copy()
        #shuffle the data
        self.kg_df = kg_df.sample(frac=1, random_state=13).reset_index(drop=True)

    def get_instance():
        if not KGTriples._instance:
            instance = KGTriples()
            instance._get_entity_cnt()
            instance._prepare_val_data()
            instance._prepare_train_data()
            instance._verify_duplicates()
            KGTriples._instance = instance
        return KGTriples._instance    
    
    def _get_entity_cnt(self):
        kg_np = self.kg_df.to_numpy()
        entities, entity_cnt = np.unique(np.concatenate([kg_np[:, 0], 
                                    kg_np[:, 2]]), return_counts=True)
        self.entity_cnt_dict = dict(zip(entities, entity_cnt))      

    def _prepare_val_data(self):
        self.idx_val = list()
        for relation in KGTriples.VAL_RELATIONS:
            rel_indices = self.kg_df[(self.kg_df['rel'] == relation)].index.to_list()
            random.shuffle(rel_indices)
            rel_val_cnt = int(len(rel_indices)*0.0025)
            idx_rel_val = list()
            for data_index in rel_indices:
                head, tail = self.kg_df.iloc[data_index,0], self.kg_df.iloc[data_index,2]                
                
                if len(idx_rel_val) == rel_val_cnt: break

                if ((self.entity_cnt_dict[head] > 3) & (self.entity_cnt_dict[tail] > 3)):
                    idx_rel_val.append(data_index)
                    self.entity_cnt_dict[head] = self.entity_cnt_dict[head] - 1
                    self.entity_cnt_dict[tail] = self.entity_cnt_dict[tail] - 1
            self.idx_val.extend(idx_rel_val)
            print('Relation:',relation,'Relation triple count:',len(rel_indices),'size of val:',len(idx_rel_val))
        print('Val length:',len(self.idx_val))
        self.val_triples_df = self.kg_df[self.kg_df.index.isin(self.idx_val)].copy()

    def _prepare_train_data(self):
        self.train_triples_df = self.kg_df[~self.kg_df.index.isin(self.idx_val)].copy()
        
    def _verify_duplicates(self):
        train_indices = self.train_triples_df.index.to_list()
        val_indices = self.val_triples_df.index.to_list()
        overlap_indices = set(train_indices).intersection(set(val_indices))
        if len(overlap_indices) > 0:
            print('Warning::: There are {n} duplicates between validation and test'.format(n=len(overlap_indices)))
        print('Number of triples for training: {n}'.format(n=len(train_indices)))
        print('Number of triples for validation: {n}'.format(n=len(val_indices)))

    def get_train_triples(self):
        return self.train_triples_df.reset_index(drop=True)
    
    def get_val_triples(self):
        return self.val_triples_df.reset_index(drop=True)
    
    def get_train_val_triples(self):
        train_val_triples_df = pd.concat([self.get_train_triples(),self.get_val_triples()])
        print('Train_Validation triple shape::',train_val_triples_df.shape)
        train_val_triples_df.drop_duplicates(inplace=True)
        print('Train_Validation triple shape::',train_val_triples_df.shape)
        return train_val_triples_df

if __name__ == '__main__':
    KGTriples.get_instance().get_train_triples().to_csv(constants.CSV_KSF2_KG_TRAIN,index=False,sep='|')
    KGTriples.get_instance().get_val_triples().to_csv(constants.CSV_KSF2_KG_VAL,index=False,sep='|')
    KGTriples.get_instance().get_train_val_triples().to_csv(constants.CSV_KSF2_KG_TRAIN_VAL,index=False,sep='|')







