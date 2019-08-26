
BRCA_raw <- read.table(gzfile("data/allhtseqcounts.brca.csv.gz"), header = T,check.names = F,row.names = 1,sep=" ")
BRCA_clean <- BRCA_raw[,-grep("11A$|11B$",colnames(BRCA_raw))]
colnames(BRCA_clean) <- substring(colnames(BRCA_clean),1,12)

bc_er_pr_her_status <- read.table(file="data/bc_er_pr_her_status.txt",header=F,sep=" ",stringsAsFactors = F)
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

write.csv(complete, "processedTable.csv")
# complete is the gene expression table that has first n columns as TN samples and remaining samples as OT (other)



