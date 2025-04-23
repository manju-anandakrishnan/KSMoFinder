import torch
import numpy as np
from classifier.ksf2.embeddings import KSMEmbeddings as ksm
from util import constants

ksf2 = torch.jit.load(constants.MODEL_KSF2)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

def _preprocess(df):
    df['substrate1'] = df['tail'].apply(lambda x:x[:x.index('_')])
    df['motif1'] = df['tail'].apply(lambda x:x[x.index('_')+1:])
    df['kinase1'] = df['head'].copy()
    return df

def _prepare_data(emb_tup, idx_row, subset_size):
    row_start_idx = idx_row; row_end_idx = idx_row+subset_size
    k_emb_test = torch.tensor(emb_tup[0][row_start_idx:row_end_idx],dtype=torch.float32)
    s_emb_test = torch.tensor(emb_tup[1][row_start_idx:row_end_idx],dtype=torch.float32)
    m_emb_test = torch.tensor(emb_tup[2][row_start_idx:row_end_idx],dtype=torch.float32)
    
    data_args = {'k_emb':k_emb_test, 's_emb':s_emb_test, 'm_emb':m_emb_test,}
    return data_args
    
def _predict(model,data):
    model.to(device)
    model.eval()    
    data = {k:v.to(device) for k,v in data.items()}
    y_pred = None
    # Forward pass to get predictions
    with torch.no_grad():
        outputs = model(data['k_emb'],
                        data['k_emb'],
                        data['s_emb'],
                        data['m_emb'],
                        )
        outputs = outputs.to('cpu')
        y_pred = outputs.numpy().flatten() if isinstance(outputs, torch.Tensor) else outputs.flatten()
    return y_pred

def get_predictions(emb_tup,data_size,subset_size=100000):
    for idx_row in range(0,data_size,subset_size):
        data = _prepare_data(emb_tup, idx_row, subset_size)
        subset_y_pred = _predict(ksf2,data)
        yield subset_y_pred


def predict(df,include_label=False):
    if include_label:
        df['label'] = 1
    df = _preprocess(df)
    data_size = df.shape[0]
    emb_tup, na_indices = ksm(constants.CSV_TRANSE_EMB).get_embeddings(df)

    y_pred = np.array([])
    for subset_y_pred in get_predictions(emb_tup,data_size):
        y_pred = np.concatenate((y_pred,subset_y_pred))

    y_pred_idx = 0
    df_pred = [None]*df.shape[0]
    for idx in range(len(df_pred)):
        if idx in na_indices: continue
        else:
            df_pred[idx] = y_pred[y_pred_idx]
            y_pred_idx += 1
    return df_pred