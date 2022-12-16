## 
import pandas as pd
import os, time
import matplotlib.pyplot as plt
import scipy.stats as stat


ML = 'VotingClassifier'
metrics = ['accuracy', 'precision', 'recall', 'F1', "AUC"]
# fi_dir = '/data/user/junghokong/co_work/SGI_CRC_seq/result/4_clinicopathologic_analysis'
fi_dir = "../../result/SGI_cohort"



## load data
adf = pd.read_csv('%s/FINAL_cross_validation_performance_SGI_cohort_all_TNM_only.txt'%fi_dir, sep='\t')
cdf = pd.read_csv('%s/FINAL_cross_validation_performance_SGI_cohort_clinical_TNM_only.txt'%fi_dir, sep='\t')
idf = pd.read_csv('%s/FINAL_cross_validation_performance_SGI_cohort_immune_only.txt'%fi_dir, sep='\t')


## boxplot #1 
## using only immune cell proportions
for metric in ['accuracy', 'precision', 'recall', 'F1', 'AUC']:
	T = idf.loc[idf['reference']=='LM22',:]['%s_%s'%(ML, metric)].tolist()
	M = idf.loc[idf['reference']=='methylCIBERSORT',:]['%s_%s'%(ML, metric)].tolist()
	TM_NET = idf.loc[idf['reference']=='FANTOM_TM_signature',:]['%s_%s'%(ML, metric)].tolist()
	#_, T_TMNET = stat.ttest_rel(T, TM_NET)
	#_, M_TMNET = stat.ttest_rel(M, TM_NET)
	_, T_TMNET = stat.mannwhitneyu(T, TM_NET)
	_, M_TMNET = stat.mannwhitneyu(M, TM_NET)

	plt.figure(figsize=(8,8))
	plt.boxplot([T, M, TM_NET])
	plt.title('%s\nT vs TM NET = %s\nM vs TM NET = %s'%(ML, T_TMNET, M_TMNET))
	plt.ylabel(metric)
	plt.xticks([1,2,3], ['T', 'M', 'TM NET'])
	plt.tight_layout()
	#plt.show()
	plt.savefig(fi_dir + "/result_immune_" + metric + ".png")
	plt.savefig(fi_dir + "/result_immune_" + metric + ".eps", format="eps")





## boxplot #2
## using immune cell proportions and/or TNM stage
ref = 'FANTOM_TM_signature'
for metric in ['accuracy', 'precision', 'recall', 'F1', 'AUC']:
	a = adf.loc[adf['reference']==ref,:]['%s_%s'%(ML,metric)].tolist()
	c = cdf.loc[cdf['reference']==ref,:]['%s_%s'%(ML, metric)].tolist()
	i = idf.loc[idf['reference']==ref,:]['%s_%s'%(ML, metric)].tolist()
	
	#_, a_c = stat.ttest_rel(a,c)
	#_, i_c = stat.ttest_rel(i,c)
	#_, i_a = stat.ttest_rel(i,a)
	_, a_c = stat.mannwhitneyu(a,c)
	_, i_c = stat.mannwhitneyu(i,c)
	_, i_a = stat.mannwhitneyu(i,a)

	plt.figure(figsize=(8,8))
	plt.boxplot([c,i,a])
	plt.title('%s\nclinical vs immune = %s\nclinical vs all = %s\nimmune vs all = %s'%(ML,i_c, a_c, i_a))
	plt.ylabel(metric)
	plt.xticks([1,2,3], ['clinical', 'immune', 'all'])
	plt.tight_layout()
	#plt.show()
	plt.savefig(fi_dir + "/result_withClinical_" + metric + ".png")
	plt.savefig(fi_dir + "/result_withClinical_" + metric + ".eps", format="eps")

