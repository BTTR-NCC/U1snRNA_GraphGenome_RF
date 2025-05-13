#clone from U1analysis.random.forest_modeleval_list.U1est.v1.1.py
import argparse
import os
import random
import pandas as pd
import numpy as np
import pickle
######evalusation section
##Stratified cros validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_validate
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, confusion_matrix

##Parameters
k = 5 #number of subset
# arg definition
parser = argparse.ArgumentParser(
    prog = 'U1analysis.random.forest_modeleval_list.py',
    usage = 'Random forest model evaluation for the estimation of U1mutation from psi and TAG list.',
    description = 'description',
    epilog='end',
    add_help=True,
)
# set arguments
parser.add_argument('listF', type=str, help="Path to the input list file; TAB-delemited information experimentTAG<TAB>psitable.")
parser.add_argument('output', type=str, help="Abs path to output")
parser.add_argument('seed', type=int, help="Seed number")
parser.add_argument('-t','--tree', type=int, help="Number of decision tree, intiger, default='100'")
# parse arguments
args = parser.parse_args()
trainlistF = args.listF
seed = args.seed
outtag = args.output
Ntree = args.tree

##body
with open(trainlistF, 'r') as file:
    for line in file:  ##set colors for each line
        parts = line.strip().split('\t')  # separate row with tab
        if len(parts) == 3:  # check 3 contents
            TAG, psitrainF, psitestF = parts[0], parts[1], parts[2]
            #read training datasets
            psiTR = pd.read_csv(psitrainF, header = 0, sep = "\t")
            ##pick up subtype (explanatory variable) from training set
            labels = psiTR['genotype']
            #pick up feature value (psi)
            features = psiTR.drop(['sample', 'genotype'], axis=1)

            ##read test dataset
            psiTE = pd.read_csv(psitestF, header = 0, sep = "\t")
            #pick up feature values of test data
            test_features = psiTE.drop(['sample', 'genotype'], axis=1)
            test_labels = psiTE['genotype']

            #select features commonly included in test/train feature
            common_features = list(set(test_features.columns).intersection(features.columns))
            train_features_common = features[common_features]
            test_features_common = test_features[common_features]

            ##Cross Validation for Training dataset
            X, y = train_features_common, labels

            ######training section
            #set randomforest classifier and proceed training
            rfc = RandomForestClassifier(random_state = seed, n_estimators = Ntree, oob_score=True)

            X = X.values  # Convert X to NumPy array
            y = y.values  # Convert y to NumPy array

            #loop by prediction threshold
            for thresh in np.arange(0.1, 1.0, 0.1):
                thresh = round(thresh,1)
                #stratified k fold split
                skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
                splitsetting = []
                accuracy_scores = []
                f1_Mut = []
                f1_WT = []
                precision_Mut = []
                precision_WT = []
                recall_Mut = []
                recall_WT = []
                youden_indices = []
                i = 1
                #loop prediction and evaluation
                #Interpolate for average FPR / TPR
                avg_fpr = np.linspace(0, 1, 100)
                avg_tpr = np.zeros_like(avg_fpr)
                for train_index, test_index in skf.split(X, y):
                    #train
                    X_train, X_test = X[train_index], X[test_index]
                    y_train, y_test = y[train_index], y[test_index]
                    rfc.fit(X_train, y_train)
                    setting = f"{TAG}\t{i}\t{sum(y_train == 'WT')}\t{sum(y_train == 'Mut')}\n"
                    splitsetting.append(setting)
                    # Calculate AUC and store it in the list
                    labs = rfc.classes_
                    if labs[0] == 'Mut':
                        y_scores = rfc.predict_proba(X_test)[:, 0]  ###alphabed order
                    else:
                        y_scores = rfc.predict_proba(X_test)[:, 1]  ###alphabed order
                    #prediction by threshold
                    y_pred = np.where(y_scores >= thresh, 'Mut', 'WT')
                    #evaluation
                    accuracy_scores.append(accuracy_score(y_test, y_pred))
                    f1_Mut.append(f1_score(y_test, y_pred, pos_label='Mut'))
                    f1_WT.append(f1_score(y_test, y_pred, pos_label='WT'))
                    precision_Mut.append(precision_score(y_test, y_pred, pos_label='Mut'))
                    precision_WT.append(precision_score(y_test, y_pred, pos_label='WT'))
                    recall_Mut.append(recall_score(y_test, y_pred, pos_label='Mut'))
                    recall_WT.append(recall_score(y_test, y_pred, pos_label='WT'))
                    fpr, tpr, thresholds = roc_curve(y_test, y_scores, pos_label='Mut')
                    roc_auc = auc(fpr, tpr)
                    #Interpolate
                    avg_tpr += np.interp(avg_fpr, fpr, tpr)
                    avg_tpr[0] = 0.0
                    #confusion matrix
                    tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=['Mut', 'WT']).ravel()
                    sensitivity = tp / (tp + fn)
                    specificity = tn / (tn + fp)
                    youden = sensitivity + specificity - 1
                    youden_indices.append(youden)
                    i += 1

                #calc average score
                avg_accuracy = sum(accuracy_scores) / len(accuracy_scores)
                avg_f1_Mut = sum(f1_Mut) / len(f1_Mut)
                avg_f1_WT = sum(f1_WT) / len(f1_WT)
                avg_precision_Mut = sum(precision_Mut) / len(precision_Mut)
                avg_precision_WT = sum(precision_WT) / len(precision_WT)
                avg_recall_Mut = sum(recall_Mut) / len(recall_Mut)
                avg_recall_WT = sum(recall_WT) / len(recall_WT)
                avg_youden = sum(youden_indices) / len(youden_indices)

                avg_tpr /= k  # average the TPRs
                avg_auc = auc(avg_fpr, avg_tpr)

                ###output Cross Validation
                #1 cross validation setting
                setoutput = outtag + '_CV_modelsummary.tsv'
                with open(setoutput, 'a') as out:
                    out.writelines(splitsetting)
                #2 result
                resultout = outtag + '_CV_evaluationsummary.tsv'
                with open(resultout, 'a') as out2:
                    out2.write(f"{TAG}\t{str(len(common_features))}\t{thresh}\t{avg_youden}\t{avg_accuracy}\t{avg_precision_Mut}\t{avg_recall_Mut}\t{avg_f1_Mut}\t{avg_precision_WT}\t{avg_recall_WT}\t{avg_f1_WT}\t{', '.join(map(str, avg_fpr))}\t{', '.join(map(str, avg_tpr))}\t{avg_auc}\n")

                #Prediction of test data
                #set randomforest classifier and proceed training
                rfc.fit(train_features_common, labels) ##train model
                predictions = rfc.predict(test_features_common)
                #predicted_probabilities = rfc.predict_proba(test_features_common)
                # Calculate AUC and store it in the list
                labs = rfc.classes_
                if labs[0] == 'Mut':
                    predicted_probabilities = rfc.predict_proba(test_features_common)[:, 0]  ###alphabed order
                else:
                    predicted_probabilities = rfc.predict_proba(test_features_common)[:, 1]  ###alphabed order
                #evaluation for test data
                y_test = test_labels
                y_pred = np.where(predicted_probabilities >= thresh, 'Mut', 'WT')
                accuracyT = accuracy_score(y_test, y_pred)
                f1_MutT = f1_score(y_test, y_pred, pos_label='Mut')
                f1_WTT = f1_score(y_test, y_pred, pos_label='WT')
                precision_MutT = precision_score(y_test, y_pred, pos_label='Mut')
                precision_WTT = precision_score(y_test, y_pred, pos_label='WT')
                recall_MutT = recall_score(y_test, y_pred, pos_label='Mut')
                recall_WTT = recall_score(y_test, y_pred, pos_label='WT')
                fprT, tprT, thresholdsT = roc_curve(y_test, predicted_probabilities, pos_label='Mut')
                roc_aucT = auc(fprT, tprT)
                #confusion matrix
                tnT, fpT, fnT, tpT = confusion_matrix(y_test, y_pred, labels=['Mut', 'WT']).ravel()
                sensitivity = tpT / (tpT + fnT)
                specificity = tnT / (tnT + fpT)
                youdenT = sensitivity + specificity - 1

                #2 result.test.prediction
                resultouttest = outtag + '_TestCohort_evaluationsummary.tsv'
                with open(resultouttest, 'a') as out3:
                    out3.write(f"{TAG}\t{str(len(common_features))}\t{thresh}\t{youdenT}\t{accuracyT}\t{precision_MutT}\t{recall_MutT}\t{f1_MutT}\t{precision_WTT}\t{recall_WTT}\t{f1_WTT}\t{', '.join(map(str, fprT))}\t{', '.join(map(str, tprT))}\t{roc_aucT}\n")

                #3 Save RFC model
                modelout = outtag + '_thresh' + str(thresh) + '_RFCmodel.pkl'
                with open(modelout, 'wb') as f:
                    pickle.dump(rfc,f)
        else:
            raise ValueError("Error: The number of contents in a row is not three.")
