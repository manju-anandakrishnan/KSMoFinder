import argparse

import torch
import torch.nn as nn
import torch.optim as optim

import numpy as np

from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold

from custom_nn import BilinearDNNModel
from embeddings import DataLoader, KSMEmbeddings
from evaluation import evaluate_model

from util import data_util
from util import constants
torch.manual_seed(13)

device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

class KSFinder2:

    '''
    Class for training KSFinder2
    
    Parameters
    ----------
    training_data:tuple of numpy arrays (kinase,substrate,motif)
    
    '''
    def __init__(self,training_data):
        
        k_emb_train = training_data[0]
        s_emb_train = training_data[1]
        m_emb_train = training_data[2]
        y_train = training_data[3]

        self.k_emb_train = torch.tensor(k_emb_train,dtype=torch.float32)        
        self.s_emb_train = torch.tensor(s_emb_train,dtype=torch.float32)
        self.m_emb_train = torch.tensor(m_emb_train,dtype=torch.float32)
        self.y_train = torch.tensor(y_train,dtype=torch.float32).view(-1,1)

    """
    Model training
    
    Returns
    ----------
    Best model after k-fold cross-validation    
    """
    def train_model(self):

        input_size = self.k_emb_train.shape[1]  
        output_size = 1 
        num_epochs = 5000

        num_folds = 5
        skf = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=13)
        
        # Initialize a list to store the ROC AUC scores for each fold
        results = []

        param_combinations = [     
            {'hidden_sizes': [40,50], 'learning_rate':0.001},
            #{'hidden_sizes': [50,], 'learning_rate':0.0001},
            #{'hidden_sizes': [100,110], 'learning_rate':0.0001}                
        ]

        for param_comb in param_combinations:
            model_fold_scores = [None] * num_folds
            model_fold_epochs = [None] * num_folds
            for fold, (train_index, test_index) in enumerate(skf.split(self.k_emb_train,self.y_train)):
                
                # Split the data into training and test sets for this fold
                X_fold_train_kemb, X_fold_train_semb, X_fold_train_memb = self.k_emb_train[train_index], self.s_emb_train[train_index], self.m_emb_train[train_index]
                X_fold_test_kemb, X_fold_test_semb, X_fold_test_memb = self.k_emb_train[test_index], self.s_emb_train[test_index], self.m_emb_train[test_index]
                y_fold_train, y_fold_test = self.y_train[train_index], self.y_train[test_index]

                try:
                    # Create the model
                    model = BilinearDNNModel(input_size, input_size, input_size, param_comb.get('hidden_sizes'), output_size)
                    model = nn.DataParallel(model,device_ids=[1])
                    model = model.to(device)
                    
                    # Define loss function and optimizer
                    criterion = nn.BCELoss()  # Binary Cross-Entropy Loss
                    optimizer = optim.Adam(model.parameters(), lr=param_comb.get('learning_rate'))

                    # Initialize variables for early stopping
                    best_score = 0
                    patience_count = 0

                    X_fold_train_kemb = X_fold_train_kemb.to(device)
                    X_fold_train_semb = X_fold_train_semb.to(device)
                    X_fold_train_memb = X_fold_train_memb.to(device)
                    y_fold_train_gpu = y_fold_train.to(device)
                    

                    X_fold_test_kemb = X_fold_test_kemb.to(device)
                    X_fold_test_semb = X_fold_test_semb.to(device)
                    X_fold_test_memb = X_fold_test_memb.to(device)
                    # Training loop
                    for epoch in range(num_epochs):
                        # Forward pass
                        tr_outputs = model(X_fold_train_kemb, X_fold_train_kemb, X_fold_train_semb, X_fold_train_memb)
                        tr_loss = criterion(tr_outputs, y_fold_train_gpu)

                        # Backward pass and optimization
                        optimizer.zero_grad()
                        tr_loss.backward()
                        optimizer.step()
                        
                        # Forward pass to get predictions on the test data
                        with torch.no_grad():
                            model.eval()
                            val_outputs = model(X_fold_test_kemb, X_fold_test_kemb, X_fold_test_semb, X_fold_test_memb)
                            val_outputs = val_outputs.to('cpu')
                            epoch_val_loss = criterion(val_outputs, y_fold_test)
                            epoch_val_loss = round(epoch_val_loss.item(),4)
                            # Convert outputs to numpy array and flatten if necessary
                            y_pred = val_outputs.numpy().flatten() if isinstance(val_outputs, torch.Tensor) else val_outputs.flatten()
                            # Compute ROC AUC score for this fold
                            roc_score = round(roc_auc_score(y_fold_test, y_pred),6)
                            pr_score = round(average_precision_score(y_fold_test, y_pred),6)
                        
                        if pr_score > best_score:
                            best_score = pr_score
                            model_fold_scores[fold]=best_score
                            model_fold_epochs[fold]=epoch
                            patience_count = 0
                        else:
                            patience_count += 1
                            if patience_count == 3:
                                model_fold_epochs[fold]=epoch-3
                                print(f"\tEarly stopping triggered at {epoch} epoch. ROC-AUC:{roc_score};PR-AUC:{pr_score}")
                                break
                except ValueError as ve:
                    print(ve.with_traceback())

            max_score = np.max(model_fold_scores)
            param_comb['epoch'] = model_fold_epochs[model_fold_scores.index(max_score)]
            results.append({'Fold':fold, 'params': param_comb, 'max_score':max_score})
        
        print(results)
        best_model = max(results, key=lambda x: x['max_score'])
        best_model_params = best_model.get('params')
        return self.train_best_model(best_model_params)

    ## Train the best model with the validation data as well
    def train_best_model(self,best_model_params):
        input_size = self.k_emb_train.shape[1]
        model = BilinearDNNModel(input_size,input_size, input_size, best_model_params.get('hidden_sizes'), self.y_train.shape[1]).to(device)
        criterion = nn.BCELoss()  # Binary Cross-Entropy Loss
        optimizer = optim.Adam(model.parameters(), lr=best_model_params.get('learning_rate'))
        k_emb_train = self.k_emb_train.to(device)
        s_emb_train = self.s_emb_train.to(device)
        m_emb_train = self.m_emb_train.to(device)
        y_train = self.y_train.to(device)
        for _ in range(best_model_params.get('epoch')):
            outputs = model(k_emb_train, k_emb_train, s_emb_train, m_emb_train)
            loss = criterion(outputs, y_train)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        return model

if __name__ == '__main__':

    # Loads values from input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--retrain',default=False,type=bool)
    args = parser.parse_args()

    # Sets the path of input embedding and model
    embedding_csv, model_path = constants.CSV_EXPRESSIVE_EMB, constants.MODEL_EXPRESSIVE

    # Instantiates classes for loading data and retrieving embeddings
    data_loader = DataLoader()
    ksm_embeddings = KSMEmbeddings(embedding_csv) 

    # Retraining model from scratch
    if args.retrain:        
        raw_training_data = data_loader.get_training_data(constants.CSV_CLF_TRAIN_DATA_ASSESS1)
        ksfinder = KSFinder2(ksm_embeddings.get_training_data(raw_training_data))
        best_model = ksfinder.train_model()
        scripted_model = torch.jit.script(best_model)
        scripted_model.save(model_path)
    
    # Evaluates model with testing dataset 1
    print('****Testing dataset 1 ****')
    raw_testing_data = data_loader.get_testing_data(constants.CSV_CLF_TEST_D1_ASSESS1)
    testing_data = ksm_embeddings.get_testing_data(raw_testing_data)
    model = torch.jit.load(model_path)
    
    roc_score, pr_score, ci = evaluate_model(model,testing_data)
    print("ROC Score:", roc_score, "PR Score:",pr_score, "CI (95%):",{ci})
    
    # Evaluates model with testing dataset 2
    print('****Testing dataset 2 ****')
    raw_testing_data = data_loader.get_testing_data(constants.CSV_CLF_TEST_D2_ASSESS1)
    testing_data = ksm_embeddings.get_testing_data(raw_testing_data)
    
    roc_score, pr_score, ci = evaluate_model(model,testing_data)
    print("ROC Score:", roc_score, "PR Score:",pr_score, "CI (95%):",{ci})