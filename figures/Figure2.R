########Plot for LUAD#####################
######plot precision########

prec <- read.csv("../data/precision_recall_results_LUAD.csv")
prec$size <- as.character(prec$size)
library(ggplot2)
p<- ggplot(prec, aes(x=size, y=precision, group=methods, color=methods)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=precision-precision_std, ymax=precision+precision_std), width=.2,
                position=position_dodge(0.05))

# Finished line plot
p <- p+labs(title="Precision of the five methods in LUAD Data", x="Number of samples per group", y = "Precision")+theme_classic() 
print(p)
###RECALL######
p<- ggplot(prec, aes(x=size, y=recall, group=methods, color=methods)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=recall-recall_std, ymax=recall+recall_std), width=.2,
                position=position_dodge(0.05))

# Finished line plot
p <- p+labs(title="Recall of the five methods in LUAD Data", x="Number of samples per group", y = "Recall")+theme_classic() 
print(p)
####### BOXPLOT#######

methods <- c("DDR", "DDR", "DDR", "DDR", "DDR", "DDR", "DDR", "DDR", "DDR", "DDR", 
             "EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM",
             "EdgeR_EXACT","EdgeR_EXACT","EdgeR_EXACT","EdgeR_EXACT", "EdgeR_EXACT","EdgeR_EXACT","EdgeR_EXACT","EdgeR_EXACT", "EdgeR_EXACT","EdgeR_EXACT",
             "DESeq", "DESeq","DESeq","DESeq",  "DESeq", "DESeq","DESeq","DESeq",  "DESeq", "DESeq",
             "DESeq2", "DESeq2","DESeq2","DESeq2",  "DESeq2", "DESeq2", "DESeq2","DESeq2",  "DESeq2", "DESeq2")

# calculations for FDR is shown in simulation section 
FDR <- c(0, 0, 0, 0, 0, 0, 0,0,0,0,
         0.01638471, 0.0370321, 0.01971493,  0.009058212, 0.01838284, 0.01971493, 0.04988677, 0.01731717, 0.01225523, 0.009990675,
         0.01653386,  0.03472776, 0.01992032, 0.008698539, 0.01746348, 0.0189907, 0.06088977, 0.01739708,  0.012417, 0.009694555, 
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0.002139205, 0.01330481, 0.003443598, 0.002400083, 0.004539288, 0.003391422, 0.05248878, 0.004382761,  0.004382761,  0.0008348116)

FDR_data <- data.frame(methods, FDR)
p<-ggplot(FDR_data, aes(x=methods, y=FDR, fill=methods)) + geom_boxplot()
print(p)

########Plot for BRCA#####################

######plot precision########

prec <- read.csv("../data/precision_recall_results_BRCA.csv")
prec$size <- as.character(prec$size)
#library(ggplot2)
p<- ggplot(prec, aes(x=size, y=precision, group=methods, color=methods)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=precision-precision_std, ymax=precision+precision_std), width=.2,
                position=position_dodge(0.05))

# Finished line plot
p <- p+labs(title="Precision of the five methods in TNBrCa Data", x="Number of samples per group", y = "Precision")+theme_classic() 
print(p)
###RECALL######
p<- ggplot(prec, aes(x=size, y=recall, group=methods, color=methods)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=recall-recall_std, ymax=recall+recall_std), width=.2,
                position=position_dodge(0.05))
# Finished line plot
p <- p+labs(title="Recall of the five methods in TNBrCa Data", x="Number of samples per group", y = "Recall")+
  theme_classic() 
print(p)
####### BOXPLOT#######

methods <- c("DDR", "DDR", "DDR", "DDR", "DDR", "DDR", 
             "EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM","EdgeR_GLM",
             "EdgeR_EXACT","EdgeR_EXACT","EdgeR_EXACT","EdgeR_EXACT", "EdgeR_EXACT","EdgeR_EXACT",
             "DESeq", "DESeq","DESeq","DESeq",  "DESeq", "DESeq",
             "DESeq2", "DESeq2","DESeq2","DESeq2",  "DESeq2", "DESeq2")
FDR <- c(0, 0, 0, 0, 0, 0, 0.039748819,0.034116657,0.062860102,0.045963618,0.041237781,0.042338318,
         0.039147956,0.034158511,0.063967249,0.045288812,0.042410286,0.042410286,
         0,0,0,0,0,0, 0.006532977,0.004431063,0.025620633,0.025620633,0.010623189,0.010679998)

FDR_data <- data.frame(methods, FDR)
p<-ggplot(FDR_data, aes(x=methods, y=FDR, fill=methods)) +
  geom_boxplot()
print(p)

