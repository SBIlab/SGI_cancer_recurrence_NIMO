## 
import pandas as pd
import os, time
import matplotlib.pyplot as plt
import scipy.stats as stat


ML = 'VotingClassifier'
metrics = ['accuracy', 'precision', 'recall', 'F1', "AUC"]
fi_dir = "../../result/GSE107422"




## load data
idf = pd.read_csv('%s/FINAL_cross_validation_performance_GSE107422_immune_only.txt'%fi_dir, sep='\t')


## boxplot #1 
## using only immune cell proportions
for metric in ['accuracy', 'precision', 'recall', 'F1', "AUC"]:
	T = idf.loc[idf['reference']=='LM22',:]['%s_%s'%(ML, metric)].tolist()
	TM_NET = idf.loc[idf['reference']=='FANTOM_TM_signature',:]['%s_%s'%(ML, metric)].tolist()
	_, T_TMNET = stat.mannwhitneyu(T, TM_NET)
	
	plt.figure(figsize=(8,8))
	plt.boxplot([T, TM_NET])
	plt.title('%s\nT vs TM NET = %s'%(ML, T_TMNET))
	plt.ylabel(metric)
	plt.xticks([1,2], ['LM22', 'NIMO'])
	plt.tight_layout()
	#plt.show()
	plt.savefig(fi_dir + "/result_immune_GSE107422_" + metric + ".png")
	plt.savefig(fi_dir + "/result_immune_GSE107422_" + metric + ".eps", format="eps")








