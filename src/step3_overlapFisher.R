args <- commandArgs(trailingOnly = T)
group1.sampleSize <- args[1]
# group1.sampleSize <- 115

load("ref.gene.bin")
norm.table <- read.csv("normalized_table.csv",row.names = 1)
print(dim(norm.table))

###############GROUP############
TN_cpm <- norm.table[, 1:group1.sampleSize]
OT_cpm <- norm.table[, (group1.sampleSize+1):ncol(norm.table)]

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
write.csv(overlap_pvalue_fdr[keep, ], "overlap_test_fdr_1_DEC.csv")

keep <- (overlap_pvalue_fdr$overlap_test_padjust <= 0.05) 
write.csv(overlap_pvalue_fdr[keep, ], "overlap_test_fdr_05.csv")


