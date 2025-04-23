import pandas as pd
import numpy as np
import os
from util import constants

if __name__ == '__main__':
    graph_alg_folder_csv = {constants.RESULT_DIR_TRANSE:constants.CSV_TRANSE_EMB,
                                  constants.RESULT_DIR_DISTMULT:constants.CSV_DISTMULT_EMB,
                                  constants.RESULT_DIR_COMPLEX:constants.CSV_COMPLEX_EMB,
                                  constants.RESULT_DIR_EXPRESSIVE:constants.CSV_EXPRESSIVE_EMB}
    for alg_folder in graph_alg_folder_csv.keys():
        ip_file = os.path.join(alg_folder,'entity_emb.csv')
        emb_df = pd.read_csv(ip_file,sep='|')
        entities = emb_df['entity'].to_list()
        embeddings = emb_df['embedding'].to_list()
        print(emb_df.head())
        print(emb_df.shape)
        emb_np_list = []
        for embedding in embeddings:
            emb_np = np.fromstring(embedding[1:-1], dtype=float, sep=',')
            emb_np_list.append(emb_np)

        df_ent_emb = pd.DataFrame(emb_np_list)
        df_ent_emb['entity'] = entities

        columns = ['entity']
        columns.extend(df_ent_emb.columns[:-1])
        df_ent_emb = df_ent_emb[columns]
        print(df_ent_emb.head())
        print(df_ent_emb.shape)
        df_ent_emb.to_csv(graph_alg_folder_csv[alg_folder],index=False,sep='|')