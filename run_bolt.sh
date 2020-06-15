#!/bin/bash


./BOLT-LMM_v2.3.4/bolt --reml --bed=/n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--fam=/n/groups/patel/yixuan/bolt/ukb_bolt_lmm.fam \
--bim=/n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--remlNoRefine \
--geneticMapFile=/n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--phenoFile="${2}" \
--phenoCol="${1}" \
--remove=./ID_non-Europeans.tab \
--remove=./bolt.in_plink_but_not_imputed.FID_IID.642.txt \
--numThreads=10 \
--LDscoresFile=/n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--verboseStats 2>&1 | tee ./"${1}"medium.log
