
rm(list = ls())
options("warn"=0)  #print max 10 warnings on screen
library(limma)      # for package 'limma'
library(edgeR)
library(DESeq2)
library(vsn)


##############################################  Parameters to set


				
work_path ="D:/disco_H/bandi_grants/Regione_Lazio_FILAS_2016_Confessore/E_LIFE_submitted_13dic2016/tesi_laurea/Zainab_Nazari/PPMI_analysis/R_limma"

#min_raw_counts_per_row=24  #minimum row counts per row(gene) for filtering , select only >= this limit
min_average_LOG_cpm_per_row= +1.0 #minimum Log2 cpm per row(gene) for filtering , select only >= this limit
min_raw_counts_per_row=1

min_raw_counts_per_sample=5		
percent_with_min_raw_counts_per_sample=0.1

Log2FC_threshold=log(2.0,2)
Pvalue_threshold=0.05
FDR_threshold=0.05
label_FDR_FC_threshold="_FDR_5pc_2fc"
label_FDR_threshold="_FDR_5pc"
label_pvalue_FC_threshold="_Pval_5pc_2fc"
label_pvalue_threshold="_Pval_5pc"


name_file_input="ir3_rna_step1_no_normalization.txt"

name_batch_factor_file_input="factors_PD_CTR.txt"


##############################################
   
setwd(work_path)



mydata<-read.table(file=name_file_input, sep = "\t", quote = "\"",row.names=1,header=TRUE, fill=TRUE)  # re-read data into a dataframe with just numbers as real data
data_rows=nrow(mydata)

myfactors<-read.table(file=name_batch_factor_file_input, sep = "\t", quote = "\"",header=TRUE, fill=TRUE)  # re-read data into a dataframe 



mydata_original=mydata
myfactors_original=myfactors



fac = factor(myfactors$Diagnosis,levels=c("PD","CTR"))                  #  factor() discretizza i gruppi, creando dei fattori
sex_fac = factor(myfactors$Sex,levels=c("M","F"))
Clinical_center_fac = factor(myfactors$Clinical_center)
RIN_covariate=as.vector(myfactors$RIN)


data_rows=nrow(mydata)




				# filtering by minimum percentage of raw counts > threshold
mydata_gt_min_counts_boolean=mydata>min_raw_counts_per_sample
rows_min_raw_counts_gt_fixed_percent_boolean= (rowSums(mydata_gt_min_counts_boolean)/ncol(mydata))>percent_with_min_raw_counts_per_sample
mydata_filtered=mydata[which(rows_min_raw_counts_gt_fixed_percent_boolean),]
dim(mydata_filtered)





dge<- DGEList(counts = mydata_filtered, group = fac)  # compute Differential IDs  group2-group1


missing_genes = setdiff(rownames(mydata),rownames(mydata_filtered))  ## filtered out genes


###############################################
#  here compute Limma test
###############################################
# Limma utilizza i modelli lineari per analizzare gli esperimenti di microarray o RNAseq
# calcola dei coefficienti per ogni gene, sono necessari due dati in input per calcolare i coefficienti:
#   1) una matrice "design", le cui colonne corrispondono alla condizione sperimentale
#      e a come viene disegnato l'esperimento, utilizzata per calcolare i fit, cioè i coefficienti
#      ogni gene, grazie alla funzione lmFit(matrice,design)
#   2) una matrice "contrast"  

		# disegnare/creare la matrice factor, with sex correction
# design <- model.matrix(~0 + fac + sex_fac  + Clinical_center_fac + RIN_covariate)       #  genera una matrice 'design' ,                               le righe corrispondono ai parametri da stimare, le colonne alla condizione sperimentale
design <- model.matrix(~0 + fac )       #  genera una matrice 'design' ,                               le righe corrispondono ai parametri da stimare, le colonne alla condizione sperimentale

colnames(design)[1:2]=c("PD","CTR")

head(design)

dge<- DGEList(counts = mydata_filtered, group = fac)  # compute Differential IDs  group2-group1


dge <- calcNormFactors(dge,method="TMM")

matrix_norm_factors= dd=t(replicate(nrow(mydata_filtered),dge[[2]]$lib.size*dge[[2]]$norm.factors))

dge_data_libsize_corrected=round( dge[[1]]/matrix_norm_factors*mean(dge[[2]]$lib.size*dge[[2]]$norm.factors) )

dge = DGEList(counts = dge_data_libsize_corrected, group = fac)  # compute Differential IDs  group2-group1

dge <- calcNormFactors(dge,method="TMM")




logCPM <- cpm(dge, log=TRUE, prior.count=2)  #The prior count is used here to avoid log(0). The logCPM values can then be used in any standard limma pipeline, using the trend=TRUE argument when running eBayes or treat. 



				##filtering by minimum LOG_mean(cpm ) per row
mydata_mean_LOG_cpm=rowMeans(logCPM)  ## log2 value of mean per row
logCPM_filtered=logCPM[mydata_mean_LOG_cpm >= 1.0,]
dim(logCPM_filtered)



logCPM_filtered_batch_effect_removed = removeBatchEffect(logCPM_filtered, batch=Clinical_center_fac, batch2=sex_fac, 
														covariates=RIN_covariate, design=design)
														
colnames(logCPM_filtered_batch_effect_removed)=myfactors$ID


write.table(logCPM_filtered_batch_effect_removed, "mydata_Log_CPM_filtered_bact_sex_effect_removed_RIN_covariate.txt", sep="\t",row.names=TRUE, col.names=TRUE)

