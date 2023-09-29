library("GENIE3")
library("tidyverse")
set.seed(42)
demo <- hpv_pdx_fold_change[complete.cases(hpv_pdx_fold_change),]
demo <- demo[complete.cases(demo),]
demo<- demo%>%select(-geneSymbol)
row.names(demo) <- demo[,1]
demo <- demo[,-1]
demo_gene_matrix <- as.matrix(demo)
demo_gene_matrix <- demo[,1:5]

# run Genie
# genei failed to run with a complete gene matrix list
demo_gene_matrix <- demo_gene_matrix[sample(1:nrow(demo_gene_matrix),4000),]
demo_gene_matrix <- as.matrix(demo_gene_matrix)
weightMat <- GENIE3(demo_gene_matrix)
# get list of linked genes
linkList <- getLinkList(weightMat)
dim(linkList)
demoList <- linkList[sample(1:nrow(linkList),50),]

write.csv(demoList,"demolist.csv",quote=F,row.names = F)
