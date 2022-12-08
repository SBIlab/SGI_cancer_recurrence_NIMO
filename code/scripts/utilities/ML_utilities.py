## Utilities for machine-learning
import scipy.stats as stat
import numpy as np
from collections import defaultdict

def performance_measurements(y_actual, y_predicted):
	'''
	TP, FP, TN, FN, accuracy, precision, recall, F1_score
	'''
	TP, FP, TN, FN = 0, 0, 0, 0

	for real, pred in zip(y_actual, y_predicted):
		if real == pred == 1:
			TP += 1
		if real == pred == 0:
			TN += 1
		if (pred == 1) and (real != pred):
			FP += 1
		if (pred == 0) and (real != pred):
			FN += 1
	TP, FP, TN, FN = float(TP), float(FP), float(TN), float(FN)
	try:
		accuracy = (TP+TN)/(TP+FP+TN+FN)
	except:
		accuracy = 'na'
	try:
		precision = TP/(TP+FP)
	except:
		precision = 'na'
	try:
		recall = TP/(TP+FN)
	except:
		recall = 'na'
	try:
		F1_score = 2*(recall * precision)/(recall + precision)
	except:
		F1_score = 'na'
	return TP, FP, TN, FN, accuracy, precision, recall, F1_score
