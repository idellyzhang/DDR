##########HEATMAP for TCGA-LUAD ################
library(edgeR)
LUAD <- read.csv("../data/TCGA-LUAD.csv", header = T,check.names = F,row.names = 1)
grp <- as.factor(rep(1, ncol(LUAD)))
y <- DGEList(LUAD,group=grp)
y <- calcNormFactors(y,method="TMM")

expr_norm <- cpm(y)

marker <- read.csv("../data/TCGA-LUAD-DDR_05.csv", header = T)
Biomarkers <- as.character(marker$row.names.cancer_cpm.[1:20])
biomarker_cpm <- expr_norm[Biomarkers,]
biomarker_cancer <- biomarker_cpm[,1:535]
biomarker_normal <- biomarker_cpm[,536:594]
biomarker_label <- cbind(biomarker_cancer, biomarker_normal)
row.names(biomarker_label) <- c("AGER", "CAV1", "CLDN18", "STX11", "PECAM1", "GPM6A", "EMP2", "FAM107A", "SPOCK2", "CLIC5", 
                                "TNNC1", "RTKN2", "FHL1", "FMO2", "FABP4", "ADRB2", "LAMP3", "CAVIN2", "ITLN2", "SLC6A4")

annot <- data.frame(condition=c(rep("LUAD", 535), rep("normal",59)))
rownames(annot) <- colnames(biomarker_label)
pheatmap(log2(as.matrix(biomarker_label)+0.1),cluster_rows = F,show_colnames=F,cluster_cols = F,annotation_col = annot,cellheight = 15)

######HEATMAP for GSE40419 Validation dataset ################reference_gene <- c("ZUFSP", "PPIL4", "BCLAF1", "TMED10", "RPS6")


################EXPRESSION#######################
LUAD5_test_clean <- read.csv("../data/GSE40419_LUAD.csv", 
                             header = T,check.names = F,row.names = 1)
cancer_test <- LUAD5_test_clean[, 1:87]
normal_test <- LUAD5_test_clean[, 88:164]

Biomarkers <- c("AGER", "CAV1", "CLDN18", "STX11", "PECAM1", "GPM6A", "EMP2", "FAM107A", "SPOCK2", "CLIC5", 
                "TNNC1", "RTKN2", "FHL1", "FMO2", "FABP4", "ADRB2", "LAMP3", "ITLN2", "SLC6A4") ###### MISSING "CAVIN2"

LUAD_clean <- cbind(cancer_test, normal_test)
marker_counts <- LUAD_clean[Biomarkers,]
cancer_marker_counts <- marker_counts[, 1:87]
normal_marker_counts <- marker_counts[, 88:164]
marker_counts_clean <- cbind(cancer_marker_counts ,normal_marker_counts)

annot <- data.frame(condition=c(rep("LUAD", 87), rep("normal",77)))
rownames(annot) <- colnames(marker_counts_clean)
pheatmap(log2(as.matrix(marker_counts_clean)+0.1),cluster_rows = F,show_colnames=F,cluster_cols = F,annotation_col = annot,cellheight = 15)
