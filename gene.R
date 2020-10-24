library(pcalg)
setwd("C:/Users/Administrator/Desktop")
GenePAV.matrix <- read.delim("./GenePAV.matrix.txt", row.names=1)
suffStat <-list(C = cor(GenePAV.matrix), n = nrow(GenePAV.matrix))
score <- new("GaussL0penObsScore", suffStat$C)
ges.fit <- ges(score)
library("marray")
write.list(ges.fit[["essgraph"]][[".->.in.edges"]], filename = "./samples.txt")

GenePAV.matrix=t(GenePAV.matrix[1:35533,])
A <- colSums(GenePAV.matrix)
GenePAV.matrix <- GenePAV.matrix[,-which(A<20|A==453)]
ncol(GenePAV.matrix)
suffStat <-list(C = cor(GenePAV.matrix), n = nrow(GenePAV.matrix))
score <- new("GaussL0penObsScore", suffStat$C)
ges.fit <- ges(score)

write.list(ges.fit[["essgraph"]][[".->.in.edges"]], filename = "./gene_interactions.txt")

