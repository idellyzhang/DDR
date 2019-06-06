###HEATMAP for 
BRCA_raw <- read.csv("../data/TCGA-BRCA.csv", header = T,check.names = F,row.names = 1)
BRCA_clean <- BRCA_raw[,-grep("11A$|11B$",colnames(BRCA_raw))]
colnames(BRCA_clean) <- substring(colnames(BRCA_clean),1,12)

bc_er_pr_her_status <- read.table(file="../data/bc_er_pr_her_status.txt",header=F,sep=" ",stringsAsFactors = F)
TNindex <- which(bc_er_pr_her_status[,2] == "Negative" & bc_er_pr_her_status[,3] == "Negative" & bc_er_pr_her_status[,4] == "Negative")
TN.barcodes <- as.character(bc_er_pr_her_status[TNindex,1])

OTindex <- which(bc_er_pr_her_status[,2] == "Positive" | bc_er_pr_her_status[,3] == "Positive" | bc_er_pr_her_status[,4] == "Positive")
OT.barcodes <- as.character(bc_er_pr_her_status[OTindex,1])
TN.index <- match(as.character(TN.barcodes),colnames(BRCA_clean))
TN.index <- TN.index[!is.na(TN.index)]
TN <- BRCA_clean[,TN.index]
OT.index <- match(as.character(OT.barcodes),colnames(BRCA_clean))
OT.index <- OT.index[!is.na(OT.index)]
OT <- BRCA_clean[,OT.index]
complete <- cbind(TN, OT)
FET_sig <- read.csv("../data/BRCA_TCGA_top20_ED_sort.csv", header = T, check.names = F,row.names = 1)


#### Heatmap for Top20 base adjusted p-values ######################
Biomarkers <- as.character(row.names(FET_sig))[1:20]
dic <- read.csv("../data/ensg2symbol.csv", header = T,check.names = F,row.names = 1)
Biomarker_name <- as.character(dic[Biomarkers,])
marker_counts <- complete[Biomarkers,]
row.names(marker_counts) <- as.character(Biomarker_name)
annot <- data.frame(condition=c(rep("TNBC",115),rep("OTHER",858)))
rownames(annot) <- colnames(marker_counts)

pheatmap(log2(as.matrix(marker_counts)+1),cluster_rows = F,show_colnames=F,cluster_cols = F,annotation_col = annot,cellheight = 12)


####### Heatmap for Top10 up-regulated genes################
positive <- read.csv("../data/BRCA_TCGA_top10_up.csv", header = T, check.names = F,row.names = 1)
Biomarkers <- as.character(c(row.names(positive)))#), row.names(negative)))
dic <- read.csv("../data/ensg2symbol.csv", header = T,check.names = F,row.names = 1)
Biomarker_name <- as.character(dic[Biomarkers,])
marker_counts <- complete[Biomarkers,]
row.names(marker_counts) <- as.character(Biomarker_name)
annot <- data.frame(condition=c(rep("TNBC",115),rep("OTHER",858)))
rownames(annot) <- colnames(marker_counts)

pheatmap(log2(as.matrix(marker_counts)+1),cluster_rows = F,show_colnames=F,cluster_cols = F,annotation_col = annot,cellheight = 20)
