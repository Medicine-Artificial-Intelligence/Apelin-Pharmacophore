import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.append('../Ph4')
from Postprocessing import pharmacophore_postprocess
from Validation import pharmacophore_validation
from Autoresult import autoph4result
from sklearn.metrics import confusion_matrix, recall_score, precision_score,  roc_curve, auc, precision_recall_curve,average_precision_score
from sklearn.metrics import accuracy_score, make_scorer



def GH_score(y_true, y_pred):
    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)
    specificity = recall_score(y_true, y_pred, pos_label=0)
    GH = (0.75*precision + 0.25*recall)*specificity
    return GH


def EF1_score(actives_list, score_list):
    """ Calculates enrichment factor.
    Parameters:
    actives_list - binary array of active/decoy status.
    score_list - array of experimental scores.
    n_percent - a decimal percentage.
    """
    total_actives = len(actives_list[actives_list == 1])
    total_compounds = len(actives_list)
    # Sort scores, while keeping track of active/decoy status
    # NOTE: This will be inefficient for large arrays
    labeled_hits = sorted(zip(score_list, actives_list), reverse=True)
    # Get top n percent of hits
    num_top = int(total_compounds * 0.01)
    top_hits = labeled_hits[0:num_top]    
    num_actives_top = len([value for score, value in top_hits if value == 1])
    # Calculate enrichment factor
    EF1 = num_actives_top / (total_actives * 0.01)
    return EF1


def get_external(y_true, y_pred, score):
    f1 = f1_score(y_true, y_pred)
    auc = roc_auc_score(y_true, score)
    gh = GH_score(y_true, y_pred)
    ef1 = EF1_score(y_true, y_pred)
    return f1, auc, gh ,ef1  

def internal_validation(data, cv, active = 'Active', model = "rhha_52",  
                        predict = "rhha_52_predict", scores = "rhha_52_rescore"):
    table = pd.DataFrame()
    plt.figure(figsize = (14,10))
    for train_index, test_index in cv.split(data.drop([active], axis =1), data[active]):
        test = data.iloc[test_index,:]
        ph4 = Validation.pharmacophore_validation(data = test, active = active , predict = predict, 
                                   scores = scores, model = model, auc_thresh = 0.6)
        ph4.validation() 
        table = pd.concat([table, ph4.table], axis = 0)
    #table['F1'] = (table['Precision']*table['Sensitivity'])/(table['Precision']+table['Sensitivity'])
    #table['GH'] = (0.75*table['Precision']+0.25*table['Sensitivity'])*table['Specificity']
    plt.plot
    return table