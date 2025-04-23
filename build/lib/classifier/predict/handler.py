import argparse
import pandas as pd
import os 
from datetime import datetime
from itertools import product

from util import constants, data_util
from classifier.predict.generator import predict

def batch_process_pred(df_l_kinases, df_l_sms, batch_size=50000):
    for i in range(0,len(df_l_kinases),batch_size):
        batch_pred_df = pd.DataFrame()
        batch_pred_df['head']= df_l_kinases[i:i+batch_size]
        batch_pred_df['tail']= df_l_sms[i:i+batch_size]
        batch_pred_df['ksf_pred'] = predict(batch_pred_df,include_label=True)
        yield i, batch_pred_df


def _is_valid_input(kinase,substrate,motif):
    valid_k = (kinase == '') or (kinase in data_util.get_kg_kinases())
    valid_s = (substrate == '') or (substrate in [x[:x.index('_')] for x in data_util.get_kg_substrate_motifs()])
    valid_m = (motif == '') or (motif in [x[x.index('_')+1:] for x in data_util.get_kg_substrate_motifs()])
    return valid_k, valid_s, valid_m

# def read_kinase_files(file_path):
#     contents = []
#     with open(file_path) as ip_f:
#         for line in ip_f:
#             contents.append(line.strip())
#     return contents

# def _get_kinases_by_type():
#     tyr_kinases = read_kinase_files('data/kinases_tyr.txt')
#     st_kinases = read_kinase_files('data/kinases_ser_thr.txt')
#     other_kinases = read_kinase_files('data/kinases_other.txt')
#     return tyr_kinases, st_kinases, other_kinases

# def _get_substrate_motifs_by_type(substrate_motifs):
#     tyr_substrate_motifs = []
#     st_substrate_motifs = []
#     for sm in substrate_motifs:
#         substrate, motif = sm[:sm.index('_')], sm[sm.index('_')+1:]
#         if motif[4] == 'Y':
#             tyr_substrate_motifs.append(sm)
#         else:
#             st_substrate_motifs.append(sm)
#     return tyr_substrate_motifs, st_substrate_motifs, substrate_motifs

# def _get_valid_ksm_combinations(kinases, substrate_motifs):
    
#     tyr_kinases, st_kinases, other_kinases = _get_kinases_by_type()
#     tyr_substrate_motifs, st_substrate_motifs, substrate_motifs = _get_substrate_motifs_by_type(substrate_motifs)
#     for kinase in kinases:
#         if kinase in tyr_kinases:
#             if len(tyr_substrate_motifs) == 0: return
#             yield zip(*product([kinase], tyr_substrate_motifs))
#         elif kinase in st_kinases:
#             if len(st_substrate_motifs) == 0: return
#             yield zip(*product([kinase], st_substrate_motifs))
#         else:
#             yield zip(*product([kinase], substrate_motifs))

def _process_data(kinase=None,substrate=None,motif=None):
    valid_k, valid_s, valid_m = _is_valid_input(kinase, substrate, motif)
    if valid_k & valid_s & valid_m:
        kinases = [kinase] if kinase else data_util.get_kg_kinases()
        substrate_motifs = [substrate+'_'+motif] if substrate and motif else []
        if (not substrate) and motif:        
            substrate_motifs = [x for x in data_util.get_kg_substrate_motifs() if x[x.index('_')+1:] == motif]
        elif (not motif) and substrate:
            substrate_motifs = [x for x in data_util.get_kg_substrate_motifs() if x[:x.index('_')] == substrate]
        elif (not substrate) and (not motif):
            substrate_motifs = [x for x in data_util.get_kg_substrate_motifs()]
        processed_any = False
        for df_l_kinases, df_l_substrate_motifs in data_util._get_valid_ksm_combinations(kinases,substrate_motifs):
            processed_any = True
            for i, batch_pred_df in batch_process_pred(df_l_kinases, df_l_substrate_motifs):
                batch_pred_df.dropna(axis=0,how='any',inplace=True)
                batch_pred_df['ksf_pred'] = batch_pred_df['ksf_pred'].apply(lambda x:round(x,3))
                yield i, data_util.process_prediction_output(batch_pred_df)
        else:
            if not processed_any:
                print('No predictions for this combination of input. Kinase is ST and motif is Y (or) motif is Y and kinase is ST')
    else:
        if not valid_k: print(f'Input kinase:{kinase} is not supported by KSFinder 2.0')
        if not valid_s: print(f'Input substrate:{substrate} is not supported by KSFinder 2.0')
        if not valid_m: print(f'Input motif:{motif} is not supported by KSFinder 2.0')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--k',default='',help='Kinase UniProt ID',required=False)
    parser.add_argument('--sp',default='',help='Substrate Protein UniProt ID',required=False)
    parser.add_argument('--sm',default='',help='Substrate Motif (-/+ 4 mer)',required=False)

    args = parser.parse_args()
    kinase = args.k
    substrate = args.sp
    motif = args.sm
    output_dir = constants.DIR_KSF2_PREDICTIONS

    if (kinase) or (substrate) or (motif):
        print('Predictions are available at:')
        for i, result_df in _process_data(kinase,substrate,motif):
            file_name = f'pred_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv'
            result_df.to_csv(os.path.join(output_dir,file_name),index=False,sep='|')
            print(os.path.join(output_dir,file_name))
        
    else:
        print('''Either kinase or substrate or motif information is required. 
              Pass argument --k, --sp, --sm for kinase, susbtrate protein and motif respectively.
               If you need all prediction, download from KSFinder 2.0 repository.''')
        

