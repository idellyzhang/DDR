

###############REFERENCE######################
PYTHON.PATH = "/usr/bin/python"

norm.table <- read.csv("normalized_table.csv",row.names = 1)
print(dim(norm.table))

reference_gene = system2(PYTHON.PATH,args="src/ref_sel_test.py",stdout=T)
print(paste0("the reference genes:",paste(reference_gene,collapse = ",")))
ref_cpm <- norm.table[reference_gene,]
write.csv(ref_cpm, "ref_cpm.csv")
save(reference_gene,file = "ref.gene.bin")
