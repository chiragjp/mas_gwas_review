#get prevalances for UK Biobank data

import os
import pandas as pd
from collections import Counter

#load uids to ignore 
toremove=[]
with open('ID_non-Europeans.tab') as f:
	for line in f:
		toremove.append(line.rstrip().split('\t')[0])

with open('bolt.in_plink_but_not_imputed.FID_IID.642.txt') as f:
	for line in f:
		toremove.append(line.rstrip().split(' ')[0])

files=os.listdir('.')
files=[x for x in files if 'icd_phenotype_encoding' in x]
output=[]
for x in files:
	d=pd.read_csv(x,sep='\t')
	d=d[~d.IID.isin(toremove)]
	if 'obesity' in x:
		pheno='Obesity'
		data=Counter(d.obesity)
		cases=data[1]
		controls=data[0]
		output.append([pheno,cases,controls])	
	else:
		d=pd.read_csv(x,sep='\t')
		pheno=d.columns[-1]
		data=Counter(d.iloc[:,-1])
		cases=data[1]
		controls=data[0]
		output.append([pheno,cases,controls])

output=pd.DataFrame(output)
output.columns=['Phenotype','Cases','Controls']
output.loc[:,'Prevalence']=output.Cases/(output.Cases+output.Controls)
output.to_csv('phenotype_prevalence_data.tsv',sep='\t')


