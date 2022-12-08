## Compute prediction of recurrence using cell_proportion + clinical features (consider all available clinical data)
## Test various ML models
import pandas as pd
import scipy.stats as stat
import numpy as np
from collections import defaultdict
import time, os, random
from sklearn import metrics
from sklearn.model_selection import train_test_split, KFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.ensemble import VotingClassifier, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import xgboost
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from sklearn.feature_selection import SelectKBest, chi2
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings(action='ignore', category=DeprecationWarning)
warnings.filterwarnings(action='ignore', category=FutureWarning)
start = time.ctime()

c_dir = os.getcwd()
os.chdir('../utilities')
execfile('parse_patient_data.py', globals())
execfile('pathway_utilities.py', globals())
execfile('ML_utilities.py', globals())
os.chdir(c_dir)



## -------------------------------------------------------------
## Import immune cell proportions
def import_immune_cell_proportion(reference_set):
	if reference_set == 'LM22':
		fdir = '%s/CIBERSORT_LM22.txt' % fo_dir
	elif reference_set == 'FANTOM_TM_signature':
		fdir = '%s/FANTOM_string_network_LM22.txt' % fo_dir
	elif reference_set == 'FANTOM_TM_signature_no_network':
		fdir = '%s/FANTOM_LM22.txt' % fo_dir
	elif reference_set == 'methylCIBERSORT':
		fdir = '%s/MethylCIBERSORT_result.txt' % fo_dir
	
	df = pd.read_csv(fdir, sep='\t')
	tmp_cells = df.columns; cells = []
	for c in tmp_cells:
		if (not 'sample' in str(c).lower()) and (not 'P-value' in c) and (not 'Pearson' in c) and (not 'RMSE' in c) and (not 'Recurrent' in c):
			cells.append(c)
	return df, cells


## Import clinicopathological data
def import_clinico():
	df = pd.read_csv('%s/clinical_SGI_cohort.txt' % fi_dir, sep='\t')
	tmp_cp = df.columns; cp = []
	for c in tmp_cp:
		if (not 'Bank.ID' in c) and (not 'sampleID' in c):
			cp.append(c)
	return df, cp


## -------------------------------------------------------------
## Initialize
fi_dir = '../data'
fo_dir = '../result/SGI_cohort'
references = ['LM22', 'methylCIBERSORT', 'FANTOM_TM_signature_no_network', 'FANTOM_TM_signature']
recur_dic = parse_updated_retrospective_cohort_recurrence() # { sample : recurrence status }
ML_algorithms = ['RandomForest', 'ExtraTrees', 'LogisticRegression', 'SVM']
training_set = 0.8
number_of_cv_iterations = 2
n_jobs = 10
test_type = 'all_TNM_only' # 'all_TNM_only', 'clinical_TNM_only', 'immune_only'




## Load data
samples = []
T_ip, cell_types = import_immune_cell_proportion('LM22')
M_ip, cell_types_m = import_immune_cell_proportion('methylCIBERSORT')
TM_ip, _ = import_immune_cell_proportion('FANTOM_TM_signature_no_network')
TM_NET_ip, _ = import_immune_cell_proportion('FANTOM_TM_signature')

sample_column = {'Input Sample':'sampleID', 'Sample':'sampleID'}
T_ip = T_ip.rename(columns=sample_column)
M_ip = M_ip.rename(columns=sample_column)
TM_ip = TM_ip.rename(columns=sample_column)
TM_NET_ip = TM_NET_ip.rename(columns=sample_column)

M_ip['sampleID'] = M_ip['sampleID'].astype(str) + '-TT-R'

cp, cp_types = import_clinico() # ClinicoPathological data (CP)
cp['sampleID'] = cp['sampleID'].astype(str) + '-TT-R'
cp2 = pd.read_csv('/home/junghokong/PROJECT/co_work/SGI_CRC_seq/200213_clinical_data_InvasionMetastasisOnly.txt', sep='\t')
cp2['sampleID'] = cp2['sampleID'].astype(str) + '-TT-R'

## samples
samples = list(set(T_ip['sampleID']) & set(M_ip['sampleID']) & set(cp['sampleID']))



# cross validation
output = defaultdict(list)
output2 = defaultdict(list)
random.seed(42)
for i in range(number_of_cv_iterations):
	iternum = i+1
	train_samples = random.sample(samples, np.int(len(samples)*training_set))
	test_samples = []
	for sample in samples:
		if not sample in train_samples:
			test_samples.append(sample)
	print('##---------------------------------\ntestType: %s / iternum: %s, %s'%(test_type, iternum, time.ctime()))



	for ref, df in zip(references, [T_ip, M_ip, TM_ip, TM_NET_ip]):
		# input feature dataFrame
		if test_type == 'all_TNM_only':
			clinical_list = ['TNM.stage']
			merged = pd.merge(df, pd.DataFrame(data=cp, columns=list(set(cp.columns) & set(np.append(['sampleID'], clinical_list)))), how='inner')
			merged = pd.merge(merged, pd.DataFrame(data=cp2, columns=list(set(cp2.columns) & set(np.append(['sampleID'], clinical_list)))), how='inner')
		if test_type == 'immune_only':
			merged = df
		if test_type == 'clinical_TNM_only':	
			clinical_list = ['TNM.stage']
			merged = pd.merge(pd.DataFrame(data=cp, columns=list(set(cp.columns) & set(np.append(['sampleID'], clinical_list)))), pd.DataFrame(data=cp2, columns=list(set(cp2.columns) & set(np.append(['sampleID'], clinical_list)))), how='inner')

		for col in ['P-value', 'Pearson Correlation', 'RMSE', 'Recurrent']:
			if col in merged.columns:
				merged = merged.drop(columns=col)
		
			
		
		## recurrence dataframe
		rdf = pd.DataFrame(data=cp, columns=['sampleID', 'Recur'])
		rdf['Recur'] = rdf['Recur'].map({'R':1, 'N':0})
		
	
		
		## train, test datasets
		X_train, X_test, y_train, y_test = [], [], [], []
		for sample in train_samples:
			X_train.append(merged.loc[merged['sampleID']==sample,:].drop(columns='sampleID').values[0])
			y_train.append(rdf.loc[rdf['sampleID']==sample,:]['Recur'].tolist()[0])
		for sample in test_samples:
			X_test.append(merged.loc[merged['sampleID']==sample,:].drop(columns='sampleID').values[0])
			y_test.append(rdf.loc[rdf['sampleID']==sample,:]['Recur'].tolist()[0])
		X_train, X_test, y_train, y_test = np.array(X_train), np.array(X_test), np.array(y_train), np.array(y_test)
		
		
		
		## Cross-Validation
		# voting estimator list
		voting_estimators = []
		voting_weights = [] # average accuracy

		# Train ML
		for ML in ML_algorithms:
			# ML parameters
			if ML == 'SVM':
				model = SVC()
				param_grid = {'kernel':['rbf'], 'gamma':[0.001, 0.01, 1, 2, 3, 4, 5, 10, 20, 50, 100], 'C':[0.001, 0.01, 0.1, 1, 2, 3, 4, 5, 10, 20, 50, 100]}
			
			if ML == 'RandomForest':
				model = RFC()
				param_grid = {'n_estimators':[500, 1000, 2000]}
			
			if ML == 'ExtraTrees':
				model = ExtraTreesClassifier()
				param_grid = {'n_estimators':[500, 1000, 2000]}

			if ML == 'LogisticRegression':
				model = LogisticRegression()
				#param_grid = {'penalty':['l1','l2'], 'C':[0.001, 0.1, 0.1, 1, 2, 3, 4, 5, 10, 20, 50, 100]}
				param_grid = {'penalty':['l2'], 'max_iter':[1e10], 'solver':['lbfgs'], 'C':np.arange(0.1, 1, 0.1), 'class_weight':['balanced'] }


			# make predictions using best hyperparamters from training-set cross-validation
			gcv = []
			gcv = GridSearchCV(model, param_grid=param_grid, cv=5, n_jobs=10).fit(X_train, y_train)
			pred_status = gcv.best_estimator_.predict(X_test)

			# Voting Classifier
			voting_estimators.append((ML, gcv.best_estimator_))
			voting_weights.append(gcv.best_estimator_.score(X_train,y_train))

			# AUC
			fpr, tpr, threshold = metrics.roc_curve(y_test, pred_status, pos_label=1)
			AUC = metrics.auc(fpr, tpr)
			print '\t'.join(map(str, [ref, ML, AUC]))
			accuracy = accuracy_score(y_test, pred_status)
			precision = precision_score(y_test, pred_status)
			recall = recall_score(y_test, pred_status)
			F1_score = f1_score(y_test, pred_status)
			for metric, score in zip(['AUC', 'accuracy', 'precision', 'recall', 'F1'], [AUC, accuracy, precision, recall, F1_score]):
				output['%s_%s'%(ML, metric)].append(score)
		
		# Voting Classifier AUC
		votingC = VotingClassifier(estimators=voting_estimators, voting='hard').fit(X_train, y_train)
		pred_status = votingC.predict(X_test)

		fpr, tpr, threshold = metrics.roc_curve(y_test, pred_status, pos_label=1)
		AUC = metrics.auc(fpr, tpr)
		print '\t'.join(map(str, [ref, 'Voting Classifier', AUC]))
		# TP, FP, TN, FN, accuracy, precision, recall, F1_score
		#TP, FP, TN, FN, accuracy, precision, recall, F1_score = performance_measurements(y_test, pred_status)
		accuracy = accuracy_score(y_test, pred_status)
		precision = precision_score(y_test, pred_status)
		recall = recall_score(y_test, pred_status)
		F1_score = f1_score(y_test, pred_status)
		for metric, score in zip(['AUC', 'accuracy', 'precision', 'recall', 'F1'], [AUC, accuracy, precision, recall, F1_score]):
			output['%s_%s'%('VotingClassifier', metric)].append(score)
			if metric == 'AUC':
				output2['%s_score'%ref].append(score)
				print '\t'.join(map(str, [ref, 'average VotingClassifier AUC: %s'%np.mean(output2['%s_score'%ref])]))
			if metric == 'accuracy':
				output2['%s_accuracy_score'%ref].append(score)
				print '\t'.join(map(str, [ref, 'average VotingClassifier accuracy: %s'%np.mean(output2['%s_accuracy_score'%ref])]))


# output dataframe
output_col = ['reference', 'iternum']
for ML in np.append(ML_algorithms, ['VotingClassifier']):
	for metric in ['AUC', 'accuracy', 'precision', 'recall', 'F1']:
		output_col.append('%s_%s'%(ML, metric))
output = pd.DataFrame(data=output, columns=output_col)


# print output
output.to_csv('%s/FINAL_cross_validation_performance_SGI_cohort_%s.txt'%(fo_dir, test_type), index=False, sep='\t')

# finish
end = time.ctime()
print 'start: %s\nend: %s'%(start,end)
