import numpy as np
import argparse

import torch
import torch.nn as nn
import torch.optim as optim

from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold

from util import constants
from embeddings import DataLoader, KSMEmbeddings
from evaluation import evaluate_model
torch.manual_seed(13)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class BilinearFeatureModule(nn.Module):
    def __init__(self, ip_dim1, ip_dim2, ip_dim3, op_dim):
        super(BilinearFeatureModule, self).__init__()
        self.bilinear1 = nn.Bilinear(ip_dim1, ip_dim1, op_dim)
        self.bilinear2 = nn.Bilinear(ip_dim1, ip_dim1, op_dim)
        self.bilinear3 = nn.Bilinear(ip_dim2, ip_dim2, op_dim)
        self.bilinear4 = nn.Bilinear(ip_dim3, ip_dim3, op_dim)
        self.bilinear5 = nn.Bilinear(op_dim, op_dim, op_dim)
        self.bilinear6 = nn.Bilinear(op_dim, op_dim, op_dim)
        self.bilinear7 = nn.Bilinear(op_dim, op_dim, op_dim)
        
    
    def forward(self, k_emb, s_emb, m_emb, kd_emb, s11_emb, k_st, s_st):
        out1 = self.bilinear1(k_emb, s_emb)
        out2 = self.bilinear2(k_emb, m_emb)
        out3 = self.bilinear3(kd_emb, s11_emb)
        out4 = self.bilinear4(k_st,s_st)
        combined1 = self.bilinear5(out1,out2)
        combined2 = self.bilinear6(combined1,out3)
        combined3 = self.bilinear7(combined2,out4)
        return combined3
    
class MultiLayerDNN(nn.Module):
    def __init__(self, hidden_sizes, output_size):
        super(MultiLayerDNN, self).__init__()
        layers = []
        in_features = 100
        for hidden_size in hidden_sizes:
            layers.append(nn.Linear(in_features, hidden_size))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(p=0.5))
            in_features = hidden_size
        layers.append(nn.Linear(in_features, output_size))
        self.network = nn.Sequential(*layers)
    
    def forward(self, x):
        x = self.network(x)
        return torch.sigmoid(x)


class BilinearDNNModel(nn.Module):
    def __init__(self, ip_dim1, ip_dim2, ip_dim3, op_dim, hidden_sizes, output_size):
        super(BilinearDNNModel, self).__init__()
        self.bilinear_module = BilinearFeatureModule(ip_dim1, ip_dim2, ip_dim3, op_dim)
        self.dnn = MultiLayerDNN(hidden_sizes, output_size)
    
    def forward(self, k_emb, s_emb, m_emb, kd_emb, s11_emb, k_st, s_st):
        combined_features = self.bilinear_module(k_emb, s_emb, m_emb, kd_emb, s11_emb, k_st, s_st)
        output = self.dnn(combined_features)
        return output

class KSFinder2:

    def __init__(self,training_data):        
        k_emb_train = torch.tensor(training_data[0],dtype=torch.float32)        
        s_emb_train = torch.tensor(training_data[1],dtype=torch.float32)
        m_emb_train = torch.tensor(training_data[2],dtype=torch.float32)
        kd_emb_train = torch.tensor(training_data[3],dtype=torch.float32)
        s11_emb_train = torch.tensor(training_data[4],dtype=torch.float32)
        k_st_emb_train = torch.tensor(training_data[5],dtype=torch.float32)
        s_st_emb_train = torch.tensor(training_data[6],dtype=torch.float32)
        y_train = torch.tensor(training_data[7],dtype=torch.float32).view(-1,1)

        self.train_data_args = {'k_emb':k_emb_train, 's_emb':s_emb_train, 'm_emb':m_emb_train,
                     'kd_emb':kd_emb_train, 's11_emb':s11_emb_train, 
                     'k_st_emb':k_st_emb_train, 's_st_emb':s_st_emb_train,
                     'y': y_train}
        
        self.dim_args = {'input_dim1':k_emb_train.shape[1],
                     'input_dim2':kd_emb_train.shape[1],
                     'input_dim3':k_st_emb_train.shape[1],
                     'output_dim':100}
        
        self.dnn_output_size = 1

    def _split_fold(self,args,index):
        fold_args = {k:v[index] for k,v in args.items()}
        return fold_args, fold_args['y']

    def train_model(self):

        input_dim1 = self.dim_args.get('input_dim1')
        input_dim2 = self.dim_args.get('input_dim2')
        input_dim3 = self.dim_args.get('input_dim3')
        output_dim = self.dim_args.get('output_dim')
        num_epochs = 5000

        num_folds = 5
        print('Cross validation fold:',num_folds)
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
            k_emb_train, y_train = self.train_data_args['k_emb'], self.train_data_args['y']

            for fold, (train_index, test_index) in enumerate(skf.split(k_emb_train,y_train)):
                
                # Split the data into training and test sets for this fold
                X_fold_train, y_fold_train = self._split_fold(self.train_data_args,train_index)
                X_fold_test, y_fold_test = self._split_fold(self.train_data_args,test_index)

                try:
                    # Create the model
                    model = BilinearDNNModel(input_dim1, input_dim2, input_dim3, output_dim, param_comb.get('hidden_sizes'), self.dnn_output_size)
                    model = nn.DataParallel(model)
                    model = model.to(device)
                    
                    # Define loss function and optimizer
                    criterion = nn.BCELoss()  # Binary Cross-Entropy Loss
                    optimizer = optim.Adam(model.parameters(), lr=param_comb.get('learning_rate'))

                    # Initialize variables for early stopping
                    best_score = 0
                    patience_count = 0

                    X_fold_train_cuda = {k:v.to(device) for k,v in X_fold_train.items()}
                    X_fold_test_cuda = {k:v.to(device) for k,v in X_fold_test.items()}
                    y_fold_train_cuda = y_fold_train.to(device)

                    # Training loop
                    for epoch in range(num_epochs):
                        # Forward pass
                        tr_outputs = model(X_fold_train_cuda['k_emb'],
                                                X_fold_train_cuda['s_emb'],
                                                X_fold_train_cuda['m_emb'],
                                                X_fold_train_cuda['kd_emb'],
                                                X_fold_train_cuda['s11_emb'],
                                                X_fold_train_cuda['k_st_emb'],
                                                X_fold_train_cuda['s_st_emb'])
                        tr_loss = criterion(tr_outputs, y_fold_train_cuda)

                        # Backward pass and optimization
                        optimizer.zero_grad()
                        tr_loss.backward()
                        optimizer.step()
                        
                        # Forward pass to get predictions on the test data
                        with torch.no_grad():
                            model.eval()
                            val_outputs = model(X_fold_test_cuda['k_emb'],
                                                X_fold_test_cuda['s_emb'],
                                                X_fold_test_cuda['m_emb'],
                                                X_fold_test_cuda['kd_emb'],
                                                X_fold_test_cuda['s11_emb'],
                                                X_fold_test_cuda['k_st_emb'],
                                                X_fold_test_cuda['s_st_emb'])
                            val_outputs = val_outputs.to('cpu')
                            epoch_val_loss = criterion(val_outputs, y_fold_test)
                            epoch_val_loss = round(epoch_val_loss.item(),4)
                            # Convert outputs to numpy array and flatten if necessary
                            y_pred = val_outputs.numpy().flatten() if isinstance(val_outputs, torch.Tensor) else val_outputs.flatten()
                            # Compute ROC AUC score for this fold
                            roc_score = round(roc_auc_score(y_fold_test, y_pred),6)
                            pr_score = round(average_precision_score(y_fold_test, y_pred),6)
                        
                        if (epoch+1)%4 == 0: 
                            # Check for early stopping
                            if pr_score > best_score:
                                best_score = pr_score
                                model_fold_scores[fold]=best_score
                                model_fold_epochs[fold]=epoch
                                patience_count = 0
                            else:
                                patience_count += 1
                                if patience_count >= 3:
                                    model_fold_epochs[fold]=epoch-3*4
                                    print("\tEarly stopping triggered.",epoch)
                                    break  
                except ValueError as ve:
                    print(ve.with_traceback())
                print(model_fold_scores)
                
            model_fold_scores = model_fold_scores[:fold+1]
            max_score = np.max(model_fold_scores)
            param_comb['epoch'] = model_fold_epochs[model_fold_scores.index(max_score)]
            results.append({'params': param_comb, 
                                'max_score':max_score})
        print(results)
        best_model = max(results, key=lambda x: x['max_score'])
        best_model_params = best_model.get('params')
        return self.train_best_model(best_model_params)

    ## Train the best model with the validation data as well
    def train_best_model(self,best_model_params):
        input_dim1 = self.dim_args.get('input_dim1')
        input_dim2 = self.dim_args.get('input_dim2')
        input_dim3 = self.dim_args.get('input_dim3')
        output_dim = self.dim_args.get('output_dim')
        model = BilinearDNNModel(input_dim1, input_dim2, input_dim3,output_dim, \
                                                best_model_params.get('hidden_sizes'), \
                                                self.dnn_output_size).to(device)
        criterion = nn.BCELoss()  # Binary Cross-Entropy Loss
        optimizer = optim.Adam(model.parameters(), lr=best_model_params.get('learning_rate'))
        train_data_args_cuda = {k:v.to(device) for k,v in self.train_data_args.items()}
        for _ in range(best_model_params.get('epoch')):
            outputs = model(train_data_args_cuda['k_emb'],
                            train_data_args_cuda['s_emb'],
                            train_data_args_cuda['m_emb'],
                            train_data_args_cuda['kd_emb'],
                            train_data_args_cuda['s11_emb'],
                            train_data_args_cuda['k_st_emb'],
                            train_data_args_cuda['s_st_emb'])
            loss = criterion(outputs, train_data_args_cuda['y'])
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        return model

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--retrain', type=bool,default=False)
    args = parser.parse_args()

    model_path = constants.MODEL_STRUCT_SEQ_FUNC
    data_loader = DataLoader()
    ksm_embeddings = KSMEmbeddings()

    if args.retrain:        
        raw_training_data = data_loader.get_training_data(constants.CSV_CLF_TRAIN_DATA_ASSESS3)
        ksfinder = KSFinder2(ksm_embeddings.get_training_data(raw_training_data))
        best_model = ksfinder.train_model()
        scripted_model = torch.jit.script(best_model)
        scripted_model.save(model_path)

    model = torch.jit.load(model_path)

    raw_testing_data = data_loader.get_testing_data(constants.CSV_CLF_TEST_D1_ASSESS3)
    testing_data = ksm_embeddings.get_testing_data(raw_testing_data)
    
    roc_score, pr_score, ci = evaluate_model(model,testing_data)
    print("ROC Score:", roc_score, "PR Score:",pr_score, "CI (95%):",{ci})

    raw_testing_data = data_loader.get_testing_data(constants.CSV_CLF_TEST_D2_ASSESS3)
    testing_data = ksm_embeddings.get_testing_data(raw_testing_data)
    
    roc_score, pr_score, ci = evaluate_model(model,testing_data)
    print("ROC Score:", roc_score, "PR Score:",pr_score, "CI (95%):",{ci})