## Benchmark with limma using GSE62872 #######

library(ggplot2)
######plot precision########
size <- c(100, 80, 60, 40, 100, 80, 60, 40) 
methods <- c("DDR", "DDR", "DDR", "DDR", "limma","limma","limma","limma")
precision <- c(0.8054627, 0.82580569, 0.8443282, 0.82973703, 0.94493036, 0.94791051, 0.9657599, 0.87581541)
precision_std <- c(0.151107089, 0.152365062, 0.128307112, 0.115292093, 0.032143796, 0.024386778, 0.021204772, 0.308804092)

prec <- data.frame(size, methods, precision, precision_std)

p<- ggplot(prec, aes(x=size, y=precision, group=methods, color=methods)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=precision-precision_std, ymax=precision+precision_std), width=.2,
                position=position_dodge(0.05)) + labs(title="Precision of prostate cancer (microarray)", x="Number of samples per group", y = "Precision")+
  theme(text = element_text(size=30), axis.title=element_text(size=30))+theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(p)
# Finished line plot
p+labs(title="Precision of prostate cancer (microarray)", x="Number of samples per group", y = "Precision")+
  theme_classic() 

###RECALL######
recall <- c(0.64973146, 0.44017186, 0.206874341, 0.061761548, 0.70196535, 0.50709527, 0.231145896, 0.079413718)
recall_std <- c(0.054179682, 0.106960263, 0.098097445, 0.036894422, 0.054284526, 0.107355115, 0.117396954, 0.055602681)
rec<- data.frame(size, methods, recall, recall_std)

p<- ggplot(rec, aes(x=size, y=recall, group=methods, color=methods)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=recall-recall_std, ymax=recall+recall_std), width=.2,
                position=position_dodge(0.05))+labs(title="Recall of the five methods in TNBrCa Data", x="Number of samples per group", y = "Recall") + theme(text = element_text(size=30), axis.title=element_text(size=30))+theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
print(p)
# Finished line plot
p+labs(title="Recall of the five methods in TNBrCa Data", x="Number of samples per group", y = "Recall")+
  theme_classic() 

####### BOXPLOT for FDR#######

methods <- c("DDR", "DDR", "DDR", "DDR", "DDR", "DDR", "DDR", "DDR", "DDR", "DDR", 
             "limma","limma","limma","limma","limma","limma", "limma", "limma", "limma", "limma"
)
FDR <- c(0.06766917, 0, 0, 0, 0, 0, 0, 0, 0,0,
         0,0,0,0,0,0,0,0,0,0)

FDR_data <- data.frame(methods, FDR)
p<-ggplot(FDR_data, aes(x=methods, y=FDR, fill=methods)) +
  geom_boxplot() + theme(text = element_text(size=30), axis.title=element_text(size=30))
print(p)



