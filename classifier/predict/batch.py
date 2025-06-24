import pandas as pd
import os
from datetime import datetime

from util import data_util,constants
from classifier.predict.generator import predict
from itertools import product


def batch_process_pred(df_l_kinases, df_l_sms, batch_size=500000):
    for i in range(0,len(df_l_kinases),batch_size):
        batch_pred_df = pd.DataFrame()
        batch_pred_df['head']= df_l_kinases[i:i+batch_size]
        batch_pred_df['tail']= df_l_sms[i:i+batch_size]
        print(f'Initial size::{batch_pred_df.shape[0]}',end='|')
        batch_pred_df['ksf_pred'] = predict(batch_pred_df,include_label=True)
        yield i, batch_pred_df

output_dir = constants.DIR_KSF2_PREDICTIONS_BATCH

processed_cnt = 0

kinases = data_util.get_kg_kinases()
substrate_motifs = data_util.get_kg_substrate_motifs()


print(f'Total kinase count {len(kinases)}')
processed_kinase_cnt = 0
for df_l_kinases, df_l_substrate_motifs in data_util.get_valid_ksm_combinations(kinases,substrate_motifs):
    processed_kinase_cnt += 1
    for i, batch_pred_df in batch_process_pred(df_l_kinases, df_l_substrate_motifs):
        batch_pred_df.dropna(axis=0,how='any',inplace=True)
        batch_pred_df['ksf_pred'] = batch_pred_df['ksf_pred'].apply(lambda x:round(x,3))
        file_name = f'pred_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv'
        data_util.process_prediction_output(batch_pred_df).to_csv(os.path.join(output_dir,file_name),sep='|',index=False)
    print(f'Yet to process kinase count {len(kinases)-processed_kinase_cnt}')
    
print('Completed processing all')