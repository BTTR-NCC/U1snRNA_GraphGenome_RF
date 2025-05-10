#Based on U1analysis.random.forest.v2.py for github
##source /home/ha6434/virturalenv3//scikitlearn/bin/activate
import argparse
import os
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
from sklearn.metrics import accuracy_score

##Parameters
seed = 1234
# arg definition
parser = argparse.ArgumentParser(
    prog = 'U1analysis.random.forest.py',
    usage = 'Random forest estimation of U1mutation.',
    description = 'description',
    epilog='end',
    add_help=True,
)
# set arguments
parser.add_argument('psitableTR', type=str, help="Path to the input psi table for training.")
parser.add_argument('psitableTE', type=str, help="Path to the input psi table for test.")
parser.add_argument('output', type=str, help="Abs path to output")
parser.add_argument('-t','--tree', type=int, help="Number of decision tree, intiger, default='100'")
# parse arguments
args = parser.parse_args()
psitrainF = args.psitableTR   ##psitrainF = "/home/ha6434/command/leafcutter/output/20220604_MB-SHH/20220604_MB-SHH_logef1.5_FDR0.1_psitable_TRth10ASth6_400x221_corrected_RFinput.tsv"
psitestF = args.psitableTE  ##psitestF = "/home/ha6434/command/leafcutter/output/20230426_JCCGMB-SHH_20220604_MB-SHH_pooled/20230426_JCCGMB-SHH_20220604_MB-SHH_pooled_perind.counts_20220604_MB-SHH_perind.counts_logef1.5_FDR0.1.junctions.gz_psitable_TRth10ASth6_400x48_corrected.txt"
outputF = args.output ##outputF = "/home/ha6434/command/leafcutter/output/20230426_JCCGMB-SHH_20220604_MB-SHH_pooled/20230426_JCCGMB-SHH_20220604_MB-SHH_pooled_perind.counts_20220604_MB-SHH_perind.counts_logef1.5_FDR0.1.junctions.gz_psitable_TRth10ASth6_400x48_corrected_RFestimated_subtype.tsv"
Ntree = args.tree   ##Ntree = 100

##body
##read training datasets
psiTR = pd.read_csv(psitrainF, header = 0, sep = "\t")
##pick up subtype (explanatory variable) from training set
labels = psiTR['genotype']
#pick up feature value (psi)
features = psiTR.drop(['sample', 'genotype'], axis=1)

##read test dataset
psiTE = pd.read_csv(psitestF, index_col = 0, header = 0, sep = "\t")
psiTE = psiTE.T
#pick up feature values of test data
test_features = pd.DataFrame(psiTE)

#select features commonly included in test/train feature
common_features = list(set(test_features.columns).intersection(features.columns))
train_features_common = features[common_features]
test_features_common = test_features[common_features]
#set randomforest classifier and proceed training
rfc = RandomForestClassifier(random_state = seed, n_estimators = Ntree, oob_score=True)
rfc.fit(train_features_common, labels) ##train model

#model evaluation
rfc_cv = cross_validate(rfc, train_features_common, labels, cv=3)  ##by cross validate
oob_score = rfc.oob_score_ #by Out-of-Bag

#prediction
predictions = rfc.predict(test_features_common)
predicted_probabilities = rfc.predict_proba(test_features_common)

with open(outputF, 'w') as out:
    out.write("sample\tprediction\tprobability_Mut\tprobability_WT\n")
    for sample, prediction, probabilities in zip(psiTE.index.tolist(), predictions, predicted_probabilities):
        out.write(f"{sample}\t{prediction}\t{probabilities[0]}\t{probabilities[1]}\n")
