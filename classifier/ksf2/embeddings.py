import pandas as pd
import numpy as np

import traceback

from util import constants

np.random.seed(13)

class DataLoader:

    """
    Class for loading data for classification models

    Parameters
    ----------
    None
    """
    def __init__(self):
        pass

    """
        Preprocesses data format to retrieve embeddings

        Parameters
        ----------
        df: pd.DataFrame

        Returns
        ----------
        dataframe formatted with kinase, substrate, motif, kinase_substrate, kinase_motif
    """
    def _preprocess(self,df):
        df['kinase1'] = df['head'].copy()
        df['substrate1'] = df['tail'].apply(lambda x:x[:x.find('_')])
        df['motif1'] = df['tail'].apply(lambda x:x[x.find('_')+1:])
        df['kinase_substrate'] = df.apply(lambda x:(x['kinase1'],x['substrate1']),axis=1)
        df['kinase_motif'] = df.apply(lambda x:(x['kinase1'],x['motif1']),axis=1)
        return df

    """
        Preprocess dataframe and shuffles the data

        Parameters
        ----------
        train_data: pd.DataFrame

        Returns
        ----------
        preprocessed dataframe after shuffling
    """
    def get_training_data(self,train_data):
        df = pd.read_csv(train_data,sep='|')
        df = self._preprocess(df)
        df = df.sample(random_state=13,frac=1).reset_index(drop=True)    
        return df
    
    """
        Preprocess dataframe 

        Parameters
        ----------
        test_data: pd.DataFrame

        Returns
        ----------
        preprocessed dataframe 
    """
    def get_testing_data(self,test_data):
        df = pd.read_csv(test_data,sep='|')
        df = self._preprocess(df)
        df = df.sample(random_state=13,frac=1).reset_index(drop=True)    
        return df

class KSMEmbeddings:

    """
    Class for mapping embeddings of entities

    Parameters
    ----------
    embedding_csv: pd.DataFrame
    
    """
    def __init__(self, embedding_csv):
        self.embedding_csv = embedding_csv
        self.entity_embedding_dict = self._load_embeddings()
        self.features = ['kinase1','substrate1','motif1','label']

    """
    Creates a dictionary of all entity:embedding

    Parameters
    ----------
    None

    Returns
    -------
    the dictionary of entity:embedding
    
    """
    def _load_embeddings(self):
        emb_df = pd.read_csv(self.embedding_csv,sep='|')
        emb_df.set_index('entity', inplace=True)
        entity_embedding_dict = emb_df.apply(lambda row: row.values, axis=1).to_dict()
        return entity_embedding_dict

    """
    Maps the data in the input dataframe to the embeddings

    Parameters
    ----------
    df: pd.DataFrame
    Input dataframe for which embeddings should be mapped

    Returns
    ----------
    embeddings of kinases, substrates, motifs as numpy array.
    
    """
    def _map_embeddings(self,df,for_predict=False):
        for df_col in df.columns:
            if df_col == 'label': continue
            df[df_col] = df[df_col].apply(lambda x:self.entity_embedding_dict.get(x))
        if for_predict:
            na_indices = df[df.isna().any(axis=1)].index.to_list()
            df.dropna(axis=0,how='any',inplace=True)
            emb_tup = tuple(np.concatenate([x for x in self.expand_np_object(df[feature].to_list())]) for feature in self.features)
            return emb_tup, na_indices
        else:
            df.dropna(axis=0,how='any',inplace=True)
            df_kinase = pd.DataFrame(df['kinase1'].to_list())
            df_substrate = pd.DataFrame(df['substrate1'].to_list())
            df_motif = pd.DataFrame(df['motif1'].to_list())
            df_label = pd.DataFrame(df['label'].to_list())            
            return df_kinase.to_numpy(), df_substrate.to_numpy(), df_motif.to_numpy(), df_label.to_numpy()

    def expand_np_object(self,x_series_list, subset_size=10000):
        for subset_i in range(0,len(x_series_list),subset_size):
            subset_x_list = x_series_list[subset_i:subset_i+subset_size]
            subset_x_array = np.array(subset_x_list, dtype=np.float32)
            yield subset_x_array

    """
    Maps the input training data to their embeddings

    Parameters
    ----------
    raw_training_data: pd.DataFrame
    training data dataframe is of format - kinase|substrate|motif|label

    Returns
    -------
    tuple of mapped embeddings of kinase,substrate and motif and label
    
    """
    def get_training_data(self,raw_training_data):
        self._training_data = raw_training_data
        training_data = self._training_data[['kinase1','substrate1','motif1','label']].copy()
        training_data_k, training_data_s, training_data_m, training_data_l = self._map_embeddings(training_data)
        print('Training data count::',training_data_l.shape[0],
              '| Positive data count::',training_data_l[training_data_l[:,-1]==1].shape[0],
              '| Negative data count::',training_data_l[training_data_l[:,-1]==0].shape[0])
        return (training_data_k, training_data_s, training_data_m, training_data_l)
    
    """
    Maps the input testing data to their embeddings

    Parameters
    ----------
    raw_testing_data: pd.DataFrame
    testing data dataframe is of format - kinase|substrate|motif|label

    Returns
    -------
    tuple of mapped embeddings of kinase,substrate and motif and label
    
    """
    def get_testing_data(self,raw_testing_data):
        self._testing_data = raw_testing_data
        testing_data = self._testing_data[['kinase1','substrate1','motif1','label']].copy()
        testing_data_k, testing_data_s, testing_data_m, testing_data_l = self._map_embeddings(testing_data)
        print('Testing data count::',testing_data_l.shape[0],
            '| Positive data count::',testing_data_l[testing_data_l[:,-1]==1].shape[0],
            '| Negative data count::',testing_data_l[testing_data_l[:,-1]==0].shape[0])
        return (testing_data_k, testing_data_s, testing_data_m, testing_data_l)
    
    """
    Maps the input data to their embeddings. This method is used for predictions after model's training and testing.

    Parameters
    ----------
    kinase: str
        kinase identifier (UniProtID)
    susbtrate: str
        susbtrate identifier (UniProtID)
    motif: str
        -/+4 mer motif

    Returns
    -------
    tuple of mapped embeddings of kinase,substrate and motif 
    
    """
    def get_data_embedding(self,kinase,substrate, motif):
        try:
            k_emb = self.entity_embedding_dict.get(kinase).reshape(1,-1)
            s_emb = self.entity_embedding_dict.get(substrate).reshape(1,-1)
            m_emb = self.entity_embedding_dict.get(motif).reshape(1,-1)
            return (k_emb, s_emb, m_emb)
        except:
            print(kinase, substrate, motif)
            print(traceback.print_exc())
            return None, None, None
    
    """
    Retrieves embeddings for the data in the input dataframe

    Parameters
    ----------
    pred_df: pd.DataFrame
        dataframe of kianse, subtrate and motif
        
    Returns
    -------
    tuple of mapped embeddings and position of null embeddings (cannot retrieve embeddings)
    
    """
    def get_embeddings(self,pred_df):
        pred_df = pred_df[self.features].copy()
        emb_tup, na_indices = self._map_embeddings(pred_df,for_predict=True)
        return emb_tup, na_indices
