#####EdgeR###############
library(edgeR)

# complete is the gene expression table that has first n columns as TN samples and remaining samples as OT (other)
complete <- read.csv("processedTable.csv",row.names=1)
print(dim(complete))
grp <- as.factor(rep(1, ncol(complete)))
y <- DGEList(complete,group=grp)
y <- calcNormFactors(y,method="TMM")

norm.table <- cpm(y)

cov <- apply(norm.table,1, sd)/apply(norm.table,1, mean)
mean <- apply(norm.table,1, mean)
std <- apply(norm.table,1, sd)
MFC <- apply(norm.table,1, max)/apply(norm.table,1, min)
print(dim(norm.table))

write.csv(norm.table,"normalized_table.csv")
out <- data.frame(cov, mean, std, MFC)
out <- out[!is.infinite(out$MFC),]

out <- cbind("genes"= rownames(out),out)

write.csv(out, "final_out.csv",row.names=F)

