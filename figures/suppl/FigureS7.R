MR <- read.csv("../data/GSE37418.csv", header = T,check.names = F,row.names = 1)


########HEATMAP of G3 vs G4 for GSE37418 #######################
G4 <- c(1,2,4,6,7,8,13,14,15,18,19,20,21,22,25,27,28,33,36, 38,39,40,45,46,49,50, 53,54,55,56,57,58,68, 69, 70, 71, 73, 75, 76)
G3 <- c(9, 10, 11, 12, 17, 23, 24,34, 35, 37, 47, 52, 59, 60, 62, 63)
G4_m <-MR[, G4]
G3_m <- MR[, G3]
#WNT <- c(3, 16, 26, 48, 61, 64, 65, 66)
#SHH <- c(5, 29, 30, 31, 32, 42, 43, 44, 51, 72, 74)
MR_clean <- cbind( G4_m, G3_m)
Biomarkers <- c("HLX", "TRIM58", "MFAP4", 
               "EOMES", "RBM24", "EN2")
biomarker_cpm <- MR_clean[Biomarkers,]
biomarker_G4 <- biomarker_cpm[,1:39]
biomarker_G3 <- biomarker_cpm[,40:55]
biomarker_label <- cbind( biomarker_G3, biomarker_G4)
annot <- data.frame(condition=c(rep("G3",16), rep("G4",39)))
rownames(annot) <- colnames(biomarker_label)
pheatmap(as.matrix(biomarker_label),cluster_rows = F,show_colnames=F,cluster_cols = F,annotation_col = annot,cellheight = 25)

########HEATMAP of G3 vs G4 for GSE21140 #######################
MR_test_clean <- read.csv("../data/GSE21140.csv", header = T,check.names = F,row.names = 1)
G4 <- c(5,10,12,13,15,18,22,24,25,26,27,30,32,33,41,44,49,55,57,59,61,63,65,71,73,78,80,86,87,92,93,94,98,99,103)
G3 <- c(2,14,17,20,35,38,39,40,42,43,51,54,58,60,66,67,81,82,83,84,85,90,91,95,96,97,102)
WNT <- c(8,19,34,36,37,62,64,74)
SHH <- c(1,3,4,6,7,9,11,16,21,23,28,29,31,45,46,47,48,50,52,53,56,68,69,70,72,75,76,77,79,88,89,100,101)
G4_m <-MR_test_clean[, G4]
G3_m <- MR_test_clean[, G3]
WNT_m <- MR_test_clean[, WNT]
SHH_m <- MR_test_clean[, SHH]
Biomarkers <- c("HLX", "TRIM58", "MFAP4", 
                "EOMES", "RBM24", "EN2")
MR_clean <- cbind(WNT_m, SHH_m, G3_m, G4_m)
marker_counts <- MR_clean[Biomarkers,]
G3_marker_counts <- marker_counts[, 42:68]
G4_marker_counts <- marker_counts[, 69:103]
marker_counts_clean <- cbind(G3_marker_counts,G4_marker_counts)
annot <- data.frame(condition=c(rep("G3",27), rep("G4",35)))
rownames(annot) <- colnames(marker_counts_clean)
pheatmap(log2(as.matrix(marker_counts_clean)),cluster_rows = F,show_colnames=F,cluster_cols = F,annotation_col = annot,cellheight = 20)

########HEATMAP of G3 vs G4 for GSE37382 #######################
MR_test2_clean <- read.csv("../data/GSE37382.csv", header = T,check.names = F,row.names = 1)
SHH_m <- MR_test2_clean[, 1:51]
G3_m <- MR_test2_clean[, 52:97]
G4_m <-MR_test2_clean[, 98:285]
Biomarkers <- c("HLX", "TRIM58", "MFAP4", 
                "EOMES", "RBM24", "EN2")

MR_clean <- cbind(SHH_m, G3_m, G4_m)
marker_counts <- MR_clean[Biomarkers,]
#SHH_marker_counts <- marker_counts[, 1:51]
G3_marker_counts <- marker_counts[, 52:97]
G4_marker_counts <- marker_counts[, 98:285]
marker_counts_clean <- cbind(G3_marker_counts,G4_marker_counts)

annot <- data.frame(condition=c(rep("G3",46), rep("G4",188)))
rownames(annot) <- colnames(marker_counts_clean)
pheatmap(log2(as.matrix(marker_counts_clean)),cluster_rows = F,show_colnames=F,cluster_cols = F,annotation_col = annot,cellheight = 40)
