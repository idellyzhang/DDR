#####EdgeR###############
library(edgeR)
setwd("~/Box/lzhang/githubddr/DDR")
args <- commandArgs(trailingOnly = T)
PYTHON.PATH = "/usr/bin/python"
if (length(args[1]) > 0) PYTHON.PATH = args[1]

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
grp <- as.factor(rep(1, ncol(complete)))
y <- DGEList(complete,group=grp)
y <- calcNormFactors(y,method="TMM")

BRCA_norm <- cpm(y)

cov <- apply(BRCA_norm,1, sd)/apply(BRCA_norm,1, mean)
mean <- apply(BRCA_norm,1, mean)
std <- apply(BRCA_norm,1, sd)
MFC <- apply(BRCA_norm,1, max)/apply(BRCA_norm,1, min)


out <- data.frame(cov, mean, std, MFC)
out <- out[!is.infinite(out$MFC),]

out <- cbind("genes"= rownames(out),out)

write.csv(out, "final_out.csv",row.names=F)

###############REFERENCE######################

reference_gene = system2(PYTHON.PATH,args="src/ref_sel_test.py",stdout=T)
print(paste0("the reference genes:",paste(reference_gene,collapse = ",")))
ref_cpm <- BRCA_norm[reference_gene,]
write.csv(ref_cpm, "ref_cpm_BRCA.csv")



###############GROUP############
TN_cpm <- BRCA_norm[, 1:ncol(TN)]
OT_cpm <- BRCA_norm[, (ncol(TN)+1):ncol(complete)]

#################Caterories################################
TN_ref <- TN_cpm[reference_gene, ] 
OT_ref <- OT_cpm[reference_gene, ] 

##############Counting distribution in each category#################
##################Computing Categorized No. #####################
TN_C0 <- c()
for (i in 1:nrow(TN_cpm)){
  test <- c()
  for (j in 1:ncol(TN_cpm)){
    
    if(as.numeric(TN_cpm[i,j]) < as.numeric(TN_ref[1,j])){ 
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  TN_C0[i] <- sum(test)
}

TN_C1 <- c()
for (i in 1:nrow(TN_cpm)){
  test <- c()
  for (j in 1:ncol(TN_cpm)){
    
    if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[1,j]) && as.numeric(TN_cpm[i,j])< as.numeric(TN_ref[2,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  TN_C1[i] <- sum(test)
}

TN_C2 <- c()
#out <- c()
for (i in 1:nrow(TN_cpm)){
  test <- c()
  for (j in 1:ncol(TN_cpm)){
    
    if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[2,j]) && as.numeric(TN_cpm[i,j])< as.numeric(TN_ref[3,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  TN_C2[i] <- sum(test)
}

TN_C3 <- c()
for (i in 1:nrow(TN_cpm)){
  test <- c()
  for (j in 1:ncol(TN_cpm)){
    
    if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[3,j]) && as.numeric(TN_cpm[i,j])< as.numeric(TN_ref[4,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  TN_C3[i] <- sum(test)
}

TN_C4 <- c()
for (i in 1:nrow(TN_cpm)){
  test <- c()
  for (j in 1:ncol(TN_cpm)){
    
    if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[4,j]) && as.numeric(TN_cpm[i,j])< as.numeric(TN_ref[5,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  TN_C4[i] <- sum(test)
}


TN_C5 <- c()
for (i in 1:nrow(TN_cpm)){
  test <- c()
  for (j in 1:ncol(TN_cpm)){
    
    if(as.numeric(TN_cpm[i,j]) >= as.numeric(TN_ref[5,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  TN_C5[i] <- sum(test)
}

OT_C0 <- c()
for (i in 1:nrow(OT_cpm)){
  test <- c()
  for (j in 1:ncol(OT_cpm)){
    
    if(as.numeric(OT_cpm[i,j]) < as.numeric(OT_ref[1,j])){ 
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  OT_C0[i] <- sum(test)
}

OT_C1 <- c()
for (i in 1:nrow(OT_cpm)){
  test <- c()
  for (j in 1:ncol(OT_cpm)){
    
    if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[1,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[2,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  OT_C1[i] <- sum(test)
}

OT_C2 <- c()
for (i in 1:nrow(OT_cpm)){
  test <- c()
  for (j in 1:ncol(OT_cpm)){
    
    if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[2,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[3,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  OT_C2[i] <- sum(test)
}

OT_C3 <- c()
for (i in 1:nrow(OT_cpm)){
  test <- c()
  for (j in 1:ncol(OT_cpm)){
    
    if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[3,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[4,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  OT_C3[i] <- sum(test)
}

OT_C4 <- c()
for (i in 1:nrow(OT_cpm)){
  test <- c()
  for (j in 1:ncol(OT_cpm)){
    
    if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[4,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[5,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  OT_C4[i] <- sum(test)
}


OT_C5 <- c()
for (i in 1:nrow(OT_cpm)){
  test <- c()
  for (j in 1:ncol(OT_cpm)){
    
    if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[5,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  OT_C5[i] <- sum(test)
}



cat_count_total <- data.frame(row.names(TN_cpm), TN_C0, TN_C1, TN_C2, TN_C3, TN_C4,TN_C5, OT_C0, OT_C1, OT_C2,
                              OT_C3, OT_C4,OT_C5)

overlap_p <- c()
for (i in 1: nrow(cat_count_total)){
  A <- as.numeric(cat_count_total[i,2:7])
  B <- as.numeric(cat_count_total[i,8:13])
  tab=as.table(rbind(A,B))
  row.names(tab)=c('TN','OT')
  c <- fisher.test(tab, workspace=2e+09,hybrid=TRUE)
  overlap_p[i] <- c$p.value
}

ED <- c()
for (i in 1: nrow(cat_count_total)){
  A <- as.numeric(cat_count_total[i,2:7])
  B <- as.numeric(cat_count_total[i,8:13])
  D <- sum(A * seq_along(A))/ncol(TN_cpm)
  E <- sum(B * seq_along(B))/ncol(OT_cpm)
  c <- (D-E)
  ED[i] <- c
}

overlap_test <- data.frame(row.names(TN_cpm),overlap_p, ED )
overlap_test_padjust <- p.adjust(overlap_p, method = "fdr", n = length(overlap_p))
overlap_test_fdr <- data.frame(row.names(TN_cpm),overlap_test_padjust, ED, overlap_p)
overlap_pvalue_fdr <- overlap_test_fdr[order(overlap_test_padjust), ]
keep <- (overlap_pvalue_fdr$overlap_test_padjust <= 0.1) 
write.csv(overlap_pvalue_fdr[keep, ], "overlap_test_fdr_1_BRCA_DEC.csv")

keep <- (overlap_pvalue_fdr$overlap_test_padjust <= 0.05) 
write.csv(overlap_pvalue_fdr[keep, ], "overlap_test_fdr_05_BRCA.csv")


