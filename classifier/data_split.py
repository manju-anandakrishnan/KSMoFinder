from util import constants
import pandas as pd
import random

random.seed(13)

if __name__ == '__main__':

    df_train_data_assess1 = pd.read_csv(constants.CSV_CLF_TRAIN_DATA_ASSESS1,sep='|')
    df_test_data_assess1 = pd.read_csv(constants.CSV_CLF_TEST_DATA_ASSESS1,sep='|')

    df_train_pos = df_train_data_assess1[df_train_data_assess1['label']==1].copy()
    df_test_pos = df_test_data_assess1[df_test_data_assess1['label']==1].copy()

    df_train_neg = df_train_data_assess1[df_train_data_assess1['label']==0].copy()
    df_test_neg = df_test_data_assess1[df_test_data_assess1['label']==0].copy()

    train_kinases = df_train_neg['head'].to_list()
    test_kinases = df_test_neg['head'].to_list()

    unique_train_kinases = set(train_kinases)

    move_to_test_indices = []
    for kinase in unique_train_kinases:
        kinase_sample_cnt = train_kinases.count(kinase)+test_kinases.count(kinase)
        new_train_cnt = int(0.8*kinase_sample_cnt)
        new_test_cnt = kinase_sample_cnt-new_train_cnt
        if new_train_cnt < train_kinases.count(kinase):
            move_cnt = train_kinases.count(kinase)-new_train_cnt
            kinase_indices = df_train_neg[df_train_neg['head']==kinase].index.to_list()
            move_to_test_indices.extend(random.sample(kinase_indices,move_cnt))
    
    new_df_train_neg = df_train_neg[~df_train_neg.index.isin(move_to_test_indices)].copy().reset_index(drop=True)
    new_df_test_neg = df_train_neg[df_train_neg.index.isin(move_to_test_indices)].copy().reset_index(drop=True)
    new_df_test_neg = pd.concat([new_df_test_neg,df_test_neg])

    df_train = pd.concat([df_train_pos,new_df_train_neg])
    df_test = pd.concat([df_test_pos,new_df_test_neg])
    df_train.to_csv(constants.CSV_CLF_TRAIN_DATA,sep='|',index=False)    
    df_test.to_csv(constants.CSV_CLF_TEST_DATA_SUBASSESS1,sep='|',index=False)

    df_test = pd.concat([df_test_pos,df_test_neg])
    df_test.to_csv(constants.CSV_CLF_TEST_DATA_SUBASSESS2,sep='|',index=False)
    
    




