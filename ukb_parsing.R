#compute prevalences from UKB data


####to generate dataframe from ukb34521_phesant.tsv file
"""
head -1 ukb34521_phesant.tsv  | sed 's/,/\t/g' | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' | grep -En 'userId|20002_1111|20002_1158|20001_1020|20002_1462|20002_1509|20002_1439|20002_1065|20002_1463|20002_1462|20002_1262|20001_1044|20002_1464|20002_1289|20002_1222|20002_1223|20002_1463' | cut -f1 -d: > line_numbers_for_phenotypes

touch phenotypes_of_interest.tsv

while read p; do paste <(cut -f${p} ukb34521_phesant.tsv) phenotypes_of_interest.tsv > foo; mv foo phenotypes_of_interest.tsv; echo $p; done<line_numbers_for_phenotypes

"""

library(tidyverse)

data = as_tibble(read.csv('phenotypes_of_interest.tsv',sep='\t'))

data = data %>% mutate(IID=userId)

data = data %>% select(userId,IID,X20002_1111,X20002_1158,X20001_1020,X20002_1462,X20002_1509,X20002_1439,X20002_1065,X20002_1262,X20001_1044,X20002_1464,X20002_1289,X20002_1222,X20002_1223,X20002_1463)

tofilter=as_tibble(read.csv('/n/groups/patel/yixuan/lung_function/data_frames/ID_non-Europeans.tab',sep='\t'))

tofilter2=as_tibble(read.csv('/n/groups/patel/yixuan/lung_function/data_frames/bolt.in_plink_but_not_imputed.FID_IID.642.txt',sep=' ',header=FALSE,row.names=NULL,stringsAsFactors=FALSE))

tofilter=bind_rows(tofilter,tofilter2)

data = data %>% filter(!(userId %in% tofilter$FID))

colnames(data)=c('FID','IID',"Asthma","Cirrhosis","Colorectal_Cancer","Crohns_Disease","enteric_diarrheal_disease","HIV","Hypertension","Parkinsons_Disease","Prostate_Cancer","Rheumatoid Arthritis","Schizophrenia","Type_1_Diabetes","Type_2_Diabetes","Ulcerative Colitis")

counts=data %>% select(-FID,-IID) %>% map(function(x) table(x))
counts=as.data.frame(t(bind_rows(counts))) %>% rownames_to_column()
counts=as_tibble(counts)
colnames(counts) = c('Disease','Controls','Cases')
counts=counts %>% mutate(Total=rowSums(counts[,2:3]),Prevalence=Cases/Total)
ibd_case=counts %>% filter(Disease == "Crohn’s Disease" | Disease == 'Ulcerative Colitis') %>% select(Cases) %>% sum()
ibd_control=counts %>% filter(Disease == "Crohn’s Disease" | Disease == 'Ulcerative Colitis') %>% select(Controls) %>% sum()
ibd_prev=ibd_case/501665
counts=add_row(counts,Disease='Inflammatory Bowel Disease',Controls=ibd_control,Cases=ibd_case,Total=501665,Prevalence=ibd_prev)

write.csv(counts,'/n/groups/patel/tierney/table_paper/disease_prevalence_data.tsv')

phenotypes = counts %>% select(Disease)

write.csv(phenotypes,'/n/groups/patel/tierney/table_paper/phenotypes_for_iteration',row.names=FALSE)
