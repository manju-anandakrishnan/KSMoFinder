from pykeen.models import ERModel
from pykeen.nn import EmbeddingSpecification
from pykeen.nn.modules import Interaction
from pykeen.utils import broadcast_cat
import torch
import torch.nn as nn

'''
This class is implemented based on ExpressivE (paper and github repository) cited below. 
The code here must be used solely to apply ExpressivE algorithm for KSFinder 2.0 model training.
License restrictions imposed by authors of ExpressivE applies to copy or use this class.
'''

'''
ExpressivE custom model code is copied from https://github.com/AleksVap/ExpressivE.git
@inproceedings{
pavlovic2023expressive,
title={ExpressivE: A Spatio-Functional Embedding For Knowledge Graph Completion},
author={Aleksandar Pavlovi{\'c} and Emanuel Sallinger},
booktitle={International Conference on Learning Representations},
year={2023},
url={https://openreview.net/forum?id=xkev3_np08z}
}

MIT License

Copyright (c) 2023 Aleksandar Pavlović

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
class Utils:
    def detach_embeddings(es):
        detached_embeddings = []
        for i in es:
            detached_embeddings.append(i.detach().to('cpu'))

        return detached_embeddings

    def preprocess_entities(es, tanh_map=True):
        prep_entites = []

        if tanh_map:
            for e in es:
                prep_entites.append(torch.tanh(e))
            return prep_entites
        else:
            return es


    def positive_sign(x):
        s = torch.sign(x)
        s[s == 0] = 1
        return s


    def preprocess_relations(r, tanh_map=True, min_denom=0.5):
        d_h, d_t, c_h, c_t, s_h, s_t = r.tensor_split(6, dim=-1)

        ReLU = nn.ReLU()

        d_h = torch.abs(d_h)
        d_t = torch.abs(d_t)

        if tanh_map:
            d_h = torch.tanh(d_h)
            d_t = torch.tanh(d_t)
            c_h = torch.tanh(c_h)
            c_t = torch.tanh(c_t)
        else:
            d_h = ReLU(d_h)
            d_t = ReLU(d_t)

        # Set s_t to a value unequal to 0!
        s_t = s_t + Utils.positive_sign(s_t) * 1e-4

        # We have to clone s_h to s_h_c due to inplace operations
        s_h_c = s_h.clone()

        diag_denominator = 1 - s_h.mul(s_t)
        slope_update_mask = torch.abs(diag_denominator) < min_denom

        adjusted_min_denom = diag_denominator - Utils.positive_sign(diag_denominator) * min_denom
        s_h_c[slope_update_mask] = s_h[slope_update_mask] + adjusted_min_denom[slope_update_mask] / s_t[slope_update_mask]

        return d_h, d_t, c_h, c_t, s_h_c, s_t

class ExpressivE(ERModel):
    def __init__(self,embedding_dim: int = 50, **kwargs,) -> None:

        super().__init__(
                interaction=ExpressivEScoringFn(p=2),
                entity_representations=EmbeddingSpecification(
                    embedding_dim=embedding_dim,
                ),
                relation_representations=EmbeddingSpecification(
                    embedding_dim=6 * embedding_dim,
                ),  # d_h, d_t, c_h, c_t, s_h, s_t
                **kwargs,
            )

class ExpressivEScoringFn(Interaction):
    def __init__(self, p: int, tanh_map: bool = True, min_denom: float = 0.5):
        super().__init__()
        self.p = p  # Norm that shall be used
        self.tanh_map = tanh_map
        self.min_denom = min_denom

    def get_score(self, d_h, d_t, c_h, c_t, s_h, s_t, h, t):
        # Calculate the score of the triple

        d = torch.concat([d_h, d_t], dim=-1)  # distance
        c = torch.concat([c_h, c_t], dim=-1)  # centers
        s = torch.concat([s_t, s_h], dim=-1)  # slopes

        ht = broadcast_cat([h, t], dim=-1)
        th = broadcast_cat([t, h], dim=-1)

        contextualized_pos = torch.abs(ht - c - torch.mul(s, th))

        is_entity_pair_within_para = torch.le(contextualized_pos, d).all(dim=-1)

        w = 2 * d + 1

        k = torch.mul(0.5 * (w - 1), (w - 1 / w))
        dist = torch.mul(contextualized_pos, w) - k

        dist[is_entity_pair_within_para] = torch.div(contextualized_pos, w)[is_entity_pair_within_para]

        return -dist.norm(p=self.p, dim=-1)

    def forward(self, h, r, t):
        d_h, d_t, c_h, c_t, s_h, s_t = Utils.preprocess_relations(r, tanh_map=self.tanh_map,
                                                            min_denom=self.min_denom)
        h, t = Utils.preprocess_entities([h, t], tanh_map=self.tanh_map)
        score = self.get_score(d_h, d_t, c_h, c_t, s_h, s_t, h, t)
        return score
'''
End of code for ExpressivE
'''