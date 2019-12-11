library(tidyverse)
library(GSVA)
library(DBI)
library(odbc)
library(reshape)

rm(list=ls())

#Load and prepare gene expression matrix
myinf1 <- "/projects/verhaak-lab/SCGP-analysis/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"

expr <- read.delim(myinf1,row.names=1)
colnames(expr) <- gsub("\\.","-",colnames(expr))
expr <- as.matrix(expr)

#only use genes expressed in at least half of the cohort
sums <- apply(expr,1,function(x)sum(x>0))
#means <- apply(expr,1,mean)
expr <- expr[which(sums>(0.5*max(sums))),]

#Load Davoli gene sets
con <- DBI::dbConnect(odbc::odbc(), "scgp")

q <- "SELECT * FROM ref.immune_signatures WHERE signature_set = 'Davoli'"

imm_sigs <- dbGetQuery(con, q)

#Convert table to signature list
immune_cell <- unique(imm_sigs[,"signature_name"])
gene_list <- list()
for(i in 1:length(immune_cell))
{
	sub_imm_sig <- imm_sigs[which(imm_sigs[,"signature_name"] == immune_cell[i]),]
	sig_genes <- sub_imm_sig[,"gene_symbol"]
	gene_list[[i]] <- sig_genes
}
names(gene_list) <- immune_cell

res <- gsva(expr, gene_list, method="ssgsea")

cor(t(res),method="s")

#                       CD4.mature CD8.effector  NK.cells   B.cells      T.reg
# CD4.mature             1.0000000    0.8095238 0.5000000 0.8809524 0.33333333
# CD8.effector           0.8095238    1.0000000 0.8333333 0.9761905 0.52380952
# NK.cells               0.5000000    0.8333333 1.0000000 0.7380952 0.50000000
# B.cells                0.8809524    0.9761905 0.7380952 1.0000000 0.40476190
# T.reg                  0.3333333    0.5238095 0.5000000 0.4047619 1.00000000
# Dendritic              0.8571429    0.9523810 0.6904762 0.9761905 0.33333333
# CD8.effector.NK.cells  0.6666667    0.8095238 0.8333333 0.7142857 0.76190476
# Macrophages            0.8333333    0.8333333 0.6190476 0.9047619 0.07142857
# Macrophages.M2         0.9523810    0.8095238 0.5476190 0.9047619 0.16666667
# Macrophages.M1         0.7619048    0.7142857 0.5476190 0.7857143 0.07142857
#                       Dendritic CD8.effector.NK.cells Macrophages
# CD4.mature            0.8571429             0.6666667  0.83333333
# CD8.effector          0.9523810             0.8095238  0.83333333
# NK.cells              0.6904762             0.8333333  0.61904762
# B.cells               0.9761905             0.7142857  0.90476190
# T.reg                 0.3333333             0.7619048  0.07142857
# Dendritic             1.0000000             0.6666667  0.88095238
# CD8.effector.NK.cells 0.6666667             1.0000000  0.54761905
# Macrophages           0.8809524             0.5476190  1.00000000
# Macrophages.M2        0.8809524             0.5714286  0.92857143
# Macrophages.M1        0.7380952             0.5476190  0.95238095
#                       Macrophages.M2 Macrophages.M1
# CD4.mature                 0.9523810     0.76190476
# CD8.effector               0.8095238     0.71428571
# NK.cells                   0.5476190     0.54761905
# B.cells                    0.9047619     0.78571429
# T.reg                      0.1666667     0.07142857
# Dendritic                  0.8809524     0.73809524
# CD8.effector.NK.cells      0.5714286     0.54761905
# Macrophages                0.9285714     0.95238095
# Macrophages.M2             1.0000000     0.85714286
# Macrophages.M1             0.8571429     1.00000000

########################################################
#Correlate these scores with flow sorted scores
########################################################

myinf2 <- "/projects/verhaak-lab/scgp/data/10x/scgp-10x-immune-proportions.txt"
imm_pct <- read.delim(myinf2,stringsAsFactor=FALSE)
imm_pct[,"overall_prop"] <- imm_pct[,"immune_cell_proportion"] * imm_pct[,"cd45_parental"]

#Make imm_pct a matrix
cells <- unique(imm_pct[,"cell_type"])
subject <- unique(imm_pct[,"subject_id"])
imm_mat <- matrix(NA,nrow=length(subject),ncol=length(cells))
rownames(imm_mat) <- subject
colnames(imm_mat) <- cells
for(i in 1:nrow(imm_pct))
{
	mysubj <- imm_pct[i,"subject_id"]
	mycell <- imm_pct[i,"cell_type"]
	imm_mat[mysubj,mycell] <- imm_pct[i,"overall_prop"]
}

cd45_parental <- imm_pct[,c("subject_id","cd45_parental")]
cd45_parental <- cd45_parental[-which(duplicated(cd45_parental)),2]
imm_mat <- cbind(imm_mat,cd45_parental)

#Change aliquot barcodes to case_barcodes
full_res <- t(res)
rownames(full_res) <- sapply(strsplit(rownames(full_res),"-"),function(x)paste(x[1:3],collapse="-"))

comxx <- intersect(rownames(imm_mat),rownames(full_res))
full_res <- full_res[comxx,]
imm_mat <- imm_mat[comxx,]

full_res <- cbind(full_res,imm_mat)

########################################################
#Correlate with methylation purity
########################################################

myinf3 <- "/projects/verhaak-lab/scgp/data/epic/scgp-mnp-dkfz-predict.tsv"
meth_pur <- read.delim(myinf3)
meth_pur[,"sample_id"] <- paste("SCGP-",meth_pur[,"sample_id"],sep="")

purity_absolute <- meth_pur[,"purity_absolute"]
purity_estimate <- meth_pur[,"purity_estimate"]
names(purity_absolute) <- names(purity_estimate) <- meth_pur[,"sample_id"]

purity_absolute <- purity_absolute[rownames(full_res)]
purity_estimate <- purity_estimate[rownames(full_res)]

full_res <- cbind(full_res,purity_absolute,purity_estimate)

cor(full_res,method="s",use="complete")

########################################################
#Try calculating cell values * purity instead of * CD45%
########################################################

purity_absolute <- meth_pur[,"purity_absolute"]
purity_estimate <- meth_pur[,"purity_estimate"]
names(purity_absolute) <- names(purity_estimate) <- meth_pur[,"sample_id"]

imm_pct[,"subject_id"] <- sapply(strsplit(imm_pct[,"subject_id"],"-"),function(x)paste(x[1:3],collapse="-"))
long_absolute <- purity_absolute[imm_pct[,"subject_id"]]
long_estimate <- purity_estimate[imm_pct[,"subject_id"]]
imm_pct <- cbind(imm_pct,long_absolute,long_estimate)

imm_pct[,"purity_prop"] <- imm_pct[,"immune_cell_proportion"] * imm_pct[,"long_absolute"]

#Make imm_pct a matrix
cells <- unique(imm_pct[,"cell_type"])
subject <- unique(imm_pct[,"subject_id"])
imm_mat <- matrix(NA,nrow=length(subject),ncol=length(cells))
rownames(imm_mat) <- subject
colnames(imm_mat) <- cells
for(i in 1:nrow(imm_pct))
{
	mysubj <- imm_pct[i,"subject_id"]
	mycell <- imm_pct[i,"cell_type"]
	imm_mat[mysubj,mycell] <- imm_pct[i,"overall_prop"]
}
full_res <- t(res)
rownames(full_res) <- sapply(strsplit(rownames(full_res),"-"),function(x)paste(x[1:3],collapse="-"))

comxx <- intersect(rownames(imm_mat),rownames(full_res))
full_res <- full_res[comxx,]
imm_mat <- imm_mat[comxx,]

full_res <- cbind(full_res,imm_mat)

########################################################
#Upload to db
########################################################

res <- as.data.frame(res)
res[,"signature_name"] <- rownames(res)
res <- melt(res, id="signature_name")

res <- res[,c(2,1,3)]
colnames(res) <- c("aliquot_barcode","signature_name","enrichment_score")
head(res)

dbWriteTable(con, Id(schema="analysis", table="davoli_immune_score"), res, overwrite=TRUE)

