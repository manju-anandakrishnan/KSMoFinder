import pandas as pd
import numpy as np
import os
import torch
import matplotlib.pyplot as plt
from util import constants

from pykeen.triples import TriplesFactory
from pykeen.models import TransE, DistMult, ComplEx
from pykeen.pipeline import pipeline
from kge.custom_models import ExpressivE

np.random.seed(13)
random_seed = 13

class DataLoader:

    def __init__(self):
        self.training_df = pd.read_csv(constants.CSV_KSF2_KG_TRAIN,sep='|')
        self.validation_df = pd.read_csv(constants.CSV_KSF2_KG_VAL,sep='|')
        self._load_entities()
        self._index_entities()
        self._load_triple_factory()

    def _load_entities(self):
        node_entities = set(self.training_df['head'].unique())
        node_entities.update(self.training_df['tail'].unique())
        node_entities = list(node_entities)
        self.node_entities = sorted(node_entities)
        self.rel_entities = sorted(list(self.training_df['rel'].unique()))

    def _index_entities(self):
        self.node_entity_index = {value:index for index,value in enumerate(self.node_entities)}
        self.rel_entity_index = {value:index for index,value in enumerate(self.rel_entities)}

    def _index_dataframe(self, df):
        df['head'] = df['head'].apply(lambda x:self.node_entity_index.get(x))
        df['tail'] = df['tail'].apply(lambda x:self.node_entity_index.get(x))
        df['rel'] = df['rel'].apply(lambda x:self.rel_entity_index.get(x))
        return df
    
    def get_node_indices(self):
        node_indices_df = pd.DataFrame.from_dict(self.node_entity_index, orient='index', columns=['emb_index'])
        node_indices_df.index.name = 'entity'
        node_indices_df.reset_index(inplace=True)
        return node_indices_df
    
    def _create_triple_factory(self,df):
        data_np = self._index_dataframe(df).to_numpy()
        np.random.shuffle(data_np)
        return TriplesFactory(torch.tensor(data_np),
                                       self.node_entity_index,
                                       self.rel_entity_index,
                                       create_inverse_triples=False)

    def _load_triple_factory(self):
        df = pd.concat([self.training_df,self.validation_df])
        self.training_triples = self._create_triple_factory(self.training_df)
        self.validation_triples = self._create_triple_factory(self.validation_df)
        self.all_triples = self._create_triple_factory(df)

    def get_training_triples(self):
        return self.training_triples
    
    def get_validation_triples(self):
        return self.validation_triples

    def get_triples(self):
        return self.all_triples

class KGE:

    algorithm = {'transE':TransE,
                    'expressivE':ExpressivE,
                    'distmult':DistMult,
                    'complex':ComplEx}
    
    def __init__(self, algorithm_nm, result_dir):
        self.kge_model = KGE.algorithm[algorithm_nm]
        self.result_dir = result_dir
        if not os.path.exists(self.result_dir):
            os.mkdir(self.result_dir)            

    def train_kge(self, training_triples, validation_triples):
        
        params = {
        'optimizer': 'Adam',
        'loss': 'pairwisehinge',
        'optimizer_kwargs': dict(lr=1e-03),
        'training_kwargs': dict(num_epochs=1000,
                                checkpoint_frequency=10,
                                checkpoint_name='checkpoint.pt',
                                checkpoint_on_failure=True,
                                batch_size=25000,
                                ),
        'evaluation_kwargs': dict(additional_filter_triples=None),
        'stopper_kwargs': dict(metric='mrr',frequency=20, 
                               patience=3,
                               relative_delta=0.01),
        'model_kwargs': dict(embedding_dim=100, 
                             random_seed=random_seed),
        'negative_sampler_kwargs': dict(filtered=True,
                                        num_negs_per_pos=1),
        }
        params.get('training_kwargs')['checkpoint_directory'] = self.result_dir

        params.get('evaluation_kwargs')['additional_filter_triples'] = [training_triples.mapped_triples,
                                                                             validation_triples.mapped_triples]
        result = pipeline(
            model = self.kge_model,
            model_kwargs = params.get('model_kwargs'),
            training=training_triples,
            validation = validation_triples,
            testing = validation_triples,
            optimizer = params.get('optimizer'),
            optimizer_kwargs = params.get('optimizer_kwargs'),
            loss = params.get('loss'),
            training_kwargs = params.get('training_kwargs'),
            evaluation_kwargs = params.get('evaluation_kwargs'),
            negative_sampler_kwargs= params.get('negative_sampler_kwargs'),
            stopper ='early',
            stopper_kwargs=params.get('stopper_kwargs'),
            use_testing_data=True,
            device='cuda:1',
            random_seed = params.get('model_kwargs').get('random_seed'),
            use_tqdm=True,
            evaluation_fallback = True,
        )
        self.save_results(result,params)
        return result, params
    
    def save_results(self,result,params):
        result.plot_losses()
        plt.savefig(os.path.join(self.result_dir,'training_loss.png'))
        plt.clf()
        result.plot_early_stopping()
        plt.savefig(os.path.join(self.result_dir,'validation_results.png'))
        plt.clf()
        result.save_to_directory(self.result_dir)
        with open(os.path.join(self.result_dir,'parameters.log'),'w') as op_f:
            print(params,file=op_f)

    def retrain_best_model(self,training_triples,validation_triples,node_indices_df,result,params):

        params['training_kwargs']['num_epochs'] = result._get_results()['stopper']['best_epoch']
        params.get('training_kwargs')['checkpoint_name'] = 'best_model_checkpoint.pt'
        best_result = pipeline(
            model = self.kge_model,
            model_kwargs = params.get('model_kwargs'),
            training = training_triples,
            testing = validation_triples,
            optimizer = params.get('optimizer'),
            optimizer_kwargs = params.get('optimizer_kwargs'),
            loss = params.get('loss'),
            training_kwargs = params.get('training_kwargs'),
            device='cuda:1',
            random_seed = params.get('model_kwargs').get('random_seed'),
            use_tqdm=True,
        )
        self.save_embeddings(node_indices_df,best_result)

    def save_embeddings(self, node_indices_df,best_result):
        entity_embeddings = best_result.model.entity_representations[0](indices=None)
        node_indices_df['embedding'] = entity_embeddings.tolist()
        node_indices_df.drop(columns=['emb_index'],axis=1,inplace=True)
        node_indices_df.to_csv(os.path.join(self.result_dir,'entity_emb.csv'),index=False,sep='|')

if __name__ == '__main__':
    kge_algorithms = ['transE','distmult','complex','expressivE']
    kge_result_dirs = {'transE':constants.RESULT_DIR_TRANSE,
                       'distmult':constants.RESULT_DIR_DISTMULT,
                       'complex':constants.RESULT_DIR_COMPLEX,
                       'expressivE':constants.RESULT_DIR_EXPRESSIVE}
    for kge_alg in kge_algorithms:
        data_loader = DataLoader()
        training_triples, validation_triples = data_loader.get_training_triples(), data_loader.get_validation_triples()
        node_indices_df = data_loader.get_node_indices()
        kge_instance = KGE(kge_alg,kge_result_dirs.get(kge_alg))
        result, params = kge_instance.train_kge(training_triples,
                            validation_triples)
        print(f'****************{kge_alg} :: Determined optimal epoch*******************')
        all_triples = data_loader.get_triples()    
        kge_instance.retrain_best_model(all_triples,validation_triples,node_indices_df,result,params)
        print(f'****************{kge_alg} :: Retrained the best {kge_alg} model with all triples*******************')
        del kge_instance, data_loader
        del all_triples,training_triples,validation_triples,node_indices_df,result,params


