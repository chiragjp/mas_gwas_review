
library(metafor)
library(rmeta)
library(stringr)
library(ggplot2)
library(cowplot)
library(tidyverse)

theme_set(theme_cowplot())
###FUNCTIONS
meta_analyze<-function(df){
  vals=unique(df$Disease_2)
  output=list()
  for(x in vals){
    df1=df[df$Disease_2==x,]
    l=nrow(df1)
    if(nrow(df1)==1){
      eff=df1[1,3]
      se=df1[1,4]
      output[[x]]=c(eff,se,as.integer(l))
    }
    else{
      meta=rma(yi=df1$AUC,sei=df1$SE, measure="OR", method="REML")
      eff=meta$beta[1]
      se=meta$se[1]
      output[[x]]=c(eff,se,as.integer(l))
    }
  }
  output=t(as.data.frame(output))
  colnames(output)=c('AUC','SE','N')
  return(output)
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  return(df)
}

###
metagenomic=read.csv('metagenomic.csv')
genomic=read.csv('genomic.csv')
tokeep=intersect(metagenomic$Disease_2,genomic$Disease_2)
d=rbind(metagenomic,genomic)

#########generate methods plot
barplotdata_method=metagenomic %>% filter(Disease_2 %in% tokeep) %>% filter(Disease_2!='Ulcerative Colitis') %>% filter(Disease_2 !="Crohn's Disease") %>% group_by(Type,Method) %>% tally()
ggplot(barplotdata_method, aes(y=n, x=reorder(Method, -n))) + theme(legend.title = element_blank()) +geom_bar(position="dodge", stat="identity")+ylab('Number of models')+theme(axis.text.x = element_text(angle = 60, hjust = 1))+ggtitle('Methodological Choices by Study Type') + theme(axis.title.x=element_blank())+theme(legend.position="bottom")
ggsave('methods_plot.pdf',height=6,width=6,units='in')

#########generate disease count plot
barplotdata_disease=metagenomic  %>% filter(Disease_2 %in% tokeep) %>% filter(Disease_2!='Ulcerative Colitis') %>% filter(Disease_2 !="Crohn's Disease") %>% group_by(Type,Disease_2) %>% tally() 
ggplot(barplotdata_disease, aes(y=n, x=reorder(Disease_2, -n))) + theme(legend.title = element_blank()) +geom_bar(position="dodge", stat="identity")+ylab('Number of models')+theme(axis.text.x = element_text(angle = 60, hjust = 1))+ggtitle('Diseases by Study Type') + theme(axis.title.x=element_blank())+theme(legend.position="bottom")
ggsave('disease_plot.pdf',height=6,width=6,units='in')

barplotdata_disease=metagenomic %>% group_by(Type,Disease_2) %>% tally() 
write.csv(barplotdata_disease,'metagenomic_model_counts.csv')

barplotdata_disease=genomic  %>% filter(Disease_2 %in% tokeep)%>% group_by(Type,Disease_2) %>% tally() 
write.csv(barplotdata_disease,'genomic_model_counts.csv')

d1=genomic[,c('Disease','Disease_2','AUC','SE')]
d2=metagenomic[,c('Disease','Disease_2','AUC','SE')]

d1_temp=d1 %>% filter(Disease_2 %in% tokeep) %>% filter(Disease_2!='Ulcerative Colitis') %>% filter(Disease_2 !="Crohn's Disease")
d2_temp=d2 %>% filter(Disease_2 %in% tokeep) %>% filter(Disease_2!='Ulcerative Colitis') %>% filter(Disease_2 !="Crohn's Disease")

###funnelplots
pdf('funnel_plot_genome.pdf')
funnel(rma(yi=d1_temp$AUC,sei=d1_temp$SE, measure="OR", method="REML"))
dev.off()
pdf('funnel_plot_metagenome.pdf')
funnel(rma(yi=d2_temp$AUC,sei=d2_temp$SE, measure="OR", method="REML"))
dev.off()

###egger regression
regtest(rma(yi=d1_temp$AUC,sei=d1_temp$SE, measure="OR", method="REML"))
regtest(rma(yi=d2_temp$AUC,sei=d2_temp$SE, measure="OR", method="REML"))

###meta analysis
d1=data.frame(meta_analyze(d1))
d2=data.frame(meta_analyze(d2))

d1$Type='GENOMIC'
d2$Type='METAGENOMIC'

d3=merge(d1,d2,all=TRUE,by="row.names")
d3=d3[,-which(names(d3) %in% c('Type.x','Type.y'))]
colnames(d3)=c('Phenotype','AUC_GENOMIC','SE_GENOMIC','N_GENOMIC','AUC_METAGENOMIC','SE_METAGENOMIC','N_METAGENOMIC')
d3=d3[which(!is.na(d3$AUC_METAGENOMIC)),]

#########generate overall plot
#max_aucs=read.csv('max_auc_g_values.csv',stringsAsFactors = FALSE)
#max_aucs$Phenotype=gsub(' ','.',max_aucs$Phenotype)
#max_aucs$Phenotype=gsub("'",'.',max_aucs$Phenotype)
#max_aucs = max_aucs %>% select(Phenotype,AUC_g,SE_AUC_g)

#d4=merge(d3,max_aucs,by.x='Phenotype',by.y='Phenotype',all=TRUE)
#d4=d4[which(!is.na(d4$AUC_GENOMIC) | !is.na(d4$AUC_g)),]
d4=round_df(d3,3)

d4=d4 %>% filter(Phenotype %in% gsub("'",'.',gsub(' ','.',tokeep)))
write.csv(d4,'auc_comparison.csv')

#max_aucs$Type=rep('Max_AUC_Genomic',nrow(max_aucs))
#max_aucs_forplot=max_aucs[,c('AUC_g','Type')]
#colnames(max_aucs_forplot)[1]='AUC'
d_temp = d %>% filter(Disease_2 %in% tokeep) %>% filter(Disease_2!='Ulcerative Colitis') %>% filter(Disease_2 !="Crohn's Disease") %>% select(AUC,Type,Disease_2)
#d_temp=rbind(max_aucs_forplot,d_temp)
ggplot(d_temp, aes(y=AUC, x=Type)) +ylim(0.4,1)+theme(legend.title = element_blank()) +geom_jitter(width = .25,aes(color=Disease_2))+geom_boxplot(alpha=.5)+ylab('AUC')+ggtitle("Overall Model Performance Across Disease") + theme(axis.title.x=element_blank())
ggsave('performance_plot.pdf',height=6,width=9,units='in')

#get overall predictive power
d4=d4 %>% filter(Phenotype!='Ulcerative.Colitis') %>% filter(Phenotype !="Crohn's.Disease")
d5=d4[which(!is.na(d4$AUC_GENOMIC)),]
rma(yi=d5$AUC_GENOMIC,sei=d5$SE_GENOMIC, measure="OR", method="REML")
rma(yi=d5$AUC_METAGENOMIC,sei=d5$SE_METAGENOMIC, measure="OR", method="REML")
#d6=d4[which(!is.na(d4$AUC_g)),]
#rma(yi=d6$AUC_g,sei=d6$SE_AUC_g, measure="OR", method="REML")

########compute b2 h2 comparison

b2h2=read.csv('b2_h2.csv') %>% drop_na()
h2= b2h2 %>% filter(type=='h2')
b2= b2h2 %>% filter(type=='b2') %>% filter(phenotype %in% h2$phenotype)

rma(yi=h2$estimate,sei=h2$SE, measure="OR", method="REML")
rma(yi=b2$estimate,sei=b2$SE, measure="OR", method="REML")


#######compute AUC PRS comparison
prs=read.csv('auc_g_prs.csv') %>% drop_na()
rma(yi=prs$AUC_PRS,sei=prs$AUC_PRS_SE, measure="OR", method="REML")


#####create forest plot
d1=d1 %>% rownames_to_column() %>% filter(rowname %in% rownames(d2))
d2=d2 %>% rownames_to_column() %>% filter(rowname %in% d1$rowname)
fpd=rbind(d1,d2)
colnames(fpd)[1]='Phenotype'
fpd$Phenotype=gsub("\\.", " ", fpd$Phenotype)
fpd$Phenotype=gsub("Crohn s", "Crohn's", fpd$Phenotype)
fpd$Phenotype=gsub("Parkinson s", "Parkinson's", fpd$Phenotype)
getorder = fpd %>% filter(Type=='METAGENOMIC') %>% arrange((AUC))%>% select(Phenotype) %>% unlist() %>% unname()
fpd$Phenotype <- factor(fpd$Phenotype, levels = getorder)

ggplot(fpd, aes(y=AUC, x=Phenotype, colour=Type)) +geom_point(position = position_dodge(width = .5)) +geom_errorbar(aes(ymin=AUC-1.96*SE, ymax=AUC+1.96*SE),position = position_dodge(width = 0.5)) +coord_flip() +theme(axis.text=element_text(size=8))  +scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) + theme(axis.title = element_text()) + ylab('AUC') + xlab('Disease')

ggsave('forestplot_ma.pdf',height=6,width=9,units='in')








