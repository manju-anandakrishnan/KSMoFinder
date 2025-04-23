import torch

from sklearn.metrics import roc_auc_score, average_precision_score
from util import data_util

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

"""
Evaluates the model with the input testing_data

Parameters
------------
testing_data: pd.DataFrame

Returns
------------
ROC_SCORE and PR_SCORE based on input testing_data predictions

"""
def evaluate_model(model,testing_data):   

    k_emb_test = torch.tensor(testing_data[0],dtype=torch.float32).to(device)
    m_emb_test = torch.tensor(testing_data[2],dtype=torch.float32).to(device)
    k_st_emb_test = torch.tensor(testing_data[5],dtype=torch.float32).to(device)
    s_st_emb_test = torch.tensor(testing_data[6],dtype=torch.float32).to(device)
    y_test = torch.tensor(testing_data[7],dtype=torch.float32).view(-1,1)
    model.eval()
    with torch.no_grad():        
        outputs = model(k_emb_test,m_emb_test,k_st_emb_test,s_st_emb_test)
        outputs = outputs.to('cpu')
        y_pred = outputs.numpy().flatten() if isinstance(outputs, torch.Tensor) else outputs.flatten()
        roc_score = round(roc_auc_score(y_test, y_pred),3)        
        pr_score = round(average_precision_score(y_test, y_pred),3)
        ci = data_util.get_confidence_interval(0.95,y_test,y_pred)
        return roc_score, pr_score, ci