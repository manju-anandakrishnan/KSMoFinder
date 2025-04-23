import torch
import torch.nn as nn

class BilinearFeatureModule(nn.Module):
    """
    Class for bilinear transformation.

    Parameters
    ----------
    k_emb_dim: int
        Dimension size of kinase embedding
    s_emb_dim: int
        Dimension size of susbtrate embedding
    m_emb_dim: int
        Dimension size of motif embedding
    out_dim
        Dimension size of output embedding
    """
    def __init__(self, k_emb_dim, s_emb_dim, m_emb_dim, out_dim):
        super(BilinearFeatureModule, self).__init__()
        self.bilinear1 = nn.Bilinear(k_emb_dim, m_emb_dim, out_dim)
        self.bilinear2 = nn.Bilinear(k_emb_dim, s_emb_dim, out_dim)
        self.bilinear3 = nn.Bilinear(out_dim, out_dim, out_dim)
    
    def forward(self, k1_emb, k2_emb, s_emb, m_emb):
        """
        Combines the embeddings via bilinear transformation

        Parameters
        ----------
        k1_emb: torch.Tensor
            kinase embedding
        k2_emb: torch.Tensor
            kinase embedding
        s_emb: torch.Tensor
            substrate embedding
        m_emb
            motif embedding

        Returns
        ----------
        Bilinearly transformed vector of kinase, substrate and motif
        """
        combined1 = self.bilinear1(k1_emb, s_emb)
        combined2 = self.bilinear2(k2_emb, m_emb)
        combined = self.bilinear3(combined1, combined2)
        return combined
    
class MultiLayerDNN(nn.Module):
    """
    Class for Deep Neural Network Architecture.

    Parameters
    ----------
    input_size: int
        size of input
    hidden_sizes: list
        size of hidden networks as list
    output_size: int
        size of output
    """
    def __init__(self, input_size, hidden_sizes, output_size):
        super(MultiLayerDNN, self).__init__()
        layers = []
        in_features = input_size
        for hidden_size in hidden_sizes:
            layers.append(nn.Linear(in_features, hidden_size))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(p=0.5))
            in_features = hidden_size
        layers.append(nn.Linear(in_features, output_size))
        self.network = nn.Sequential(*layers)
    
    """
        Forwards the input data through neural network layers

        Parameters
        ----------
        x: torch.Tensor
            input embedding (combined vector)

        Returns
        ----------
        outcome transformed using sigmoid activation (range 0 to 1)
        """
    def forward(self, x):
        x = self.network(x)
        return torch.sigmoid(x)


class BilinearDNNModel(nn.Module):
    """
    Class for Bilinear Deep Neural Network.

    Parameters
    ----------
    k_emb_dim: int
        size of kinase embedding
    s_emb_dim: int
        size of substrate embedding
    m_emb_dim: int
        size of motif embedding
    dnn_hidden_sizes: list
        size of hidden layers for DNN
    dnn_output_size: int
        size of output from DNN
    bilinear_output_size: int, default=100
        size of output after bilinear transformation
    """
    def __init__(self, k_emb_dim, s_emb_dim, m_emb_dim, dnn_hidden_sizes, dnn_output_size, bilinear_output_size=100):
        super(BilinearDNNModel, self).__init__()
        self.bilinear_module = BilinearFeatureModule(k_emb_dim, s_emb_dim, m_emb_dim, bilinear_output_size)
        self.dnn = MultiLayerDNN(bilinear_output_size, dnn_hidden_sizes, dnn_output_size)
    
    """
        Forwards the input data for bilinear transformation and then through neural network layers for prediction

        Parameters
        ----------
        k1_emb: torch.Tensor
            kinase embedding
        k2_emb: torch.Tensor
            kinase embedding
        s_emb: torch.Tensor
            substrate embedding
        m_emb
            motif embedding

        Returns
        ----------
        outcome transformed using sigmoid activation (range 0 to 1)
    """
    def forward(self, k1_emb, k2_emb, s_emb, m_emb):
        combined_features = self.bilinear_module(k1_emb, k2_emb, s_emb, m_emb)
        output = self.dnn(combined_features)
        return output