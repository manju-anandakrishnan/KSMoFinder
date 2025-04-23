'''
This is an utility module maintaining classes for plotting curves and computing performance metrics
'''

from sklearn.metrics import average_precision_score,roc_auc_score, roc_curve, precision_recall_curve, recall_score
from sklearn.metrics import balanced_accuracy_score, average_precision_score, confusion_matrix
import matplotlib.pyplot as plt

class Curve:

    def get_roc_curve(y_labels,y_probs):
        roc_score, fpr, tpr, thresholds = Score.get_roc_scores(y_labels, y_probs)
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        plt.plot(fpr, tpr, label='ROC Curve (AUC = %0.3f)' % roc_score)
        plt.plot([0, 1], [0, 1], linestyle='--', color='orange', label='Random classifier')  
        plt.plot([0, 0, 1], [0, 1, 1], linestyle=':', color='green', label='Perfect classifier')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.legend(loc="lower right")
        return fig

    def get_pr_curve(y_labels,y_probs):
        pr_score, precision, recall, thresholds = Score.get_pr_scores(y_labels, y_probs)
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        plt.plot(recall, precision, label='PR Curve (Avg precision = %0.3f)' % pr_score)
        plt.plot([0, 0, 1], [0, 0, 0], linestyle='--', color='orange', label='Random classifier')  
        plt.plot([0, 1, 1], [1, 1, 0], linestyle=':', color='green', label='Perfect classifier')
        plt.xlabel('recall')
        plt.ylabel('precision')
        plt.legend(loc="lower left")
        return fig
    
    # Use it when plotting multiple curves. Inputs - y_labels, y_probs, curve_colors are lists of list
    def get_roc_curves(y_labels_list,y_probs_list,curve_colors,curve_labels=[]):
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        for i in range(0,len(y_labels_list)):
            roc_score, fpr, tpr, thresholds = Score.get_roc_score(
                y_labels_list[i], 
                y_probs_list[i])            
            if len(curve_labels) > 0:
                plt.plot(fpr, tpr, color=curve_colors[i],label=f'{curve_labels[i]}(AUC={round(roc_score,3)})')
            else:
                plt.plot(fpr, tpr, color=curve_colors[i])
        plt.plot([0, 1], [0, 1], linestyle='--', color='orange', label='Random classifier')  
        plt.plot([0, 0, 1], [0, 1, 1], linestyle=':', color='green', label='Perfect classifier')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.legend(loc="lower right")
        return fig
    
    def get_pr_curves(y_labels_list,y_probs_list,curve_colors,curve_labels=[]):
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        for i in range(0,len(y_labels_list)):
            pr_score, precision, recall, thresholds = Score.get_pr_score(y_labels_list[i], y_probs_list[i])        
            if len(curve_labels) > 0:
                plt.plot(recall, precision, color=curve_colors[i],label=f'{curve_labels[i]}(AUC={round(pr_score,3)})')
            else:
                plt.plot(recall, precision, color=curve_colors[i])
        plt.plot([0, 0, 1], [0, 0, 0], linestyle='--', color='orange', label='Random classifier')  
        plt.plot([0, 1, 1], [1, 1, 0], linestyle=':', color='green', label='Perfect classifier')
        plt.xlabel('recall')
        plt.ylabel('precision')
        plt.legend(loc=(0.3, 0.1))
        return fig

    
    def reset_plt():
        plt.clf()
    
class Score:

    def get_roc_score(y_labels,y_probs):
        roc_score = roc_auc_score(y_labels, y_probs,average='macro')
        roc_score = round(roc_score,3)
        fpr, tpr, thresholds = roc_curve(y_labels, y_probs)
        return roc_score, fpr, tpr, thresholds

    def get_pr_score(y_labels,y_probs):
        pr_score = average_precision_score(y_labels, y_probs,average='macro')
        pr_score = round(pr_score,3)
        precision, recall, thresholds = precision_recall_curve(y_labels, y_probs)
        return pr_score, precision, recall, thresholds
    
    def accuracy_score(y_labels,y_pred):
        return round(balanced_accuracy_score(y_labels, y_pred),3)

    def precision_score(y_labels,y_pred):
        return round(average_precision_score(y_labels, y_pred),3)
    
    def confusion_matrix(y_labels,y_pred):
        tn, fp, fn, tp = confusion_matrix(y_labels, y_pred).ravel()
        return tn,fp,fn,tp
    
    def recall_score(y_labels,y_pred):
        return round(recall_score(y_labels,y_pred),3)


