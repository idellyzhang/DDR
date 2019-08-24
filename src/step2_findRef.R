###### The types of data: Microarray and RNA-Seq #######################
args <- commandArgs(trailingOnly = T)
data <- args[1]
###############REFERENCE######################
PYTHON.PATH = "/Users/anaconda3/bin/python3.6"

norm.table <- read.csv("normalized_table.csv",row.names = 1)
print(dim(norm.table))
if(data == "microarray"){
  reference_gene = system2(PYTHON.PATH,args="src/ref_sel_test_symbol.py",stdout=T)
} else if (data == "RNASeq"){
  reference_gene = system2(PYTHON.PATH,args="src/ref_sel_test.py",stdout=T)
}

print(paste0("the reference genes:",paste(reference_gene,collapse = ",")))
ref_cpm <- norm.table[reference_gene,]
write.csv(ref_cpm, "ref_cpm.csv")
save(reference_gene,file = "ref.gene.bin")


