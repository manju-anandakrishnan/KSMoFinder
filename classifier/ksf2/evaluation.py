import torch

from sklearn.metrics import roc_auc_score, average_precision_score

device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

def evaluate_model(model,testing_data):        
    k_emb_test = torch.tensor(testing_data[0],dtype=torch.float32).to(device)
    s_emb_test = torch.tensor(testing_data[1],dtype=torch.float32).to(device)
    m_emb_test = torch.tensor(testing_data[2],dtype=torch.float32).to(device)
    y_test = torch.tensor(testing_data[3],dtype=torch.float32).view(-1,1)
    model.eval()
    with torch.no_grad():
        outputs = model(k_emb_test,k_emb_test,s_emb_test, m_emb_test)
        outputs = outputs.to('cpu')
        y_pred = outputs.numpy().flatten() if isinstance(outputs, torch.Tensor) else outputs.flatten()
        roc_score = round(roc_auc_score(y_test, y_pred),3)
        pr_score = round(average_precision_score(y_test, y_pred),3)
        return roc_score, pr_score