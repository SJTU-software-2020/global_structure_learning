setwd("C:/Users/Administrator/Desktop")
rice_453_flag_leaf_angle <- read.table("./phe/rice_453_flag_leaf_angle.phe", row.names=1, quote="\"", comment.char="")
rice_453_grain_length <- read.table("./phe/rice_453_grain_length.phe", row.names=1, quote="\"", comment.char="")
rice_453_grain_weight <- read.table("./phe/rice_453_grain_weight.phe", row.names=1, quote="\"", comment.char="")
rice_453_grain_width <- read.table("./phe/rice_453_grain_width.phe", row.names=1, quote="\"", comment.char="")
rice_453_height <- read.table("./phe/rice_453_height.phe", row.names=1, quote="\"", comment.char="")
rice_453_leaf_angle <- read.table("./phe/rice_453_leaf_angle.phe", row.names=1, quote="\"", comment.char="")
rice_453_leaf_length <- read.table("./phe/rice_453_leaf_length.phe", row.names=1, quote="\"", comment.char="")
rice_453_leaf_width <- read.table("./phe/rice_453_leaf_width.phe", row.names=1, quote="\"", comment.char="")
ricedata <- data.frame(
  flag_leaf_angle = rice_453_flag_leaf_angle$V3,
  grain_length = rice_453_grain_length$V3,
  grain_weight = rice_453_grain_weight$V3,
  grain_width = rice_453_grain_width$V3,
  height = rice_453_height$V3,
  leaf_angle = rice_453_leaf_angle$V3,
  leaf_length = rice_453_leaf_length$V3,
  leaf_width = rice_453_leaf_width$V3
)

library(pcalg)
tempdata=NULL
tdata=read.delim("./GenePAV.matrix.txt", row.names=1)
tdata=t(tdata[1:35533,])
A <- colSums(tdata)
tdata <- tdata[,-which(A<20|A==453)]
ncol(tdata)

temp=data.frame(ricedata$flag_leaf_angle,tdata)


suffStat <-list(C = cor(temp), n = nrow(temp))
score <- new("GaussL0penObsScore", suffStat$C)
ges.fit <- ges(score)
# ges.fit$repr$.in.edges
trytry=unlist(ges.fit$essgraph$`.->.in.edges`[lapply(ges.fit$essgraph$`.->.in.edges`,length)>0])
as.numeric(trytry)
names(ges.fit$essgraph$`.->.in.edges`[trytry])
try_again=c(names(ges.fit$essgraph$`.->.in.edges`[lapply(ges.fit$essgraph$`.->.in.edges`,length)>0]),names(ges.fit$essgraph$`.->.in.edges`[trytry]))
try_again=try_again[!duplicated(try_again)]
tempdata=cbind(tempdata,try_again)
smalldata <-list(C = cor(try[try_again]), n = nrow(try))
# suffStat <-list(C = cor(temp), n = nrow(temp))
# fciPlus.gmL <- fciPlus(suffStat, indepTest=gaussCItest,
#                        alpha = 0.0001, p=108)
# plot(fciPlus.gmL)
# tempnames=c(colnames(fciPlus.gmL@amat)[colSums(fciPlus.gmL@amat)!=0], 
#             rownames(fciPlus.gmL@amat)[rowSums(fciPlus.gmL@amat)!=0])
# # rownames(fciPlus.gmL@amat)[rowSums(fciPlus.gmL@amat)!=0])
# tempnames <- as.numeric(tempnames[!duplicated(tempnames)])
# colnames(temp)[tempnames]
# smalldata <-list(C = cor(temp[tempnames]), n = nrow(temp[tempnames]))
fciPlus.gmL <- fciPlus(smalldata, indepTest=gaussCItest,
                       alpha = 0.0001, labels =names(temp[tempnames]))

tempnames=c(colnames(fciPlus.gmL@amat)[which(fciPlus.gmL@amat[1,]!=0)], 
            rownames(fciPlus.gmL@amat)[which(fciPlus.gmL@amat[,1]!=0)])
tempnames <- tempnames[!duplicated(tempnames)]
write.csv(tempnames,"./flag_leaf_angle.csv", row.names = F,na = "NA")
# for (i in 0:289) {
#   left = 50*i+2
#   right = 50*(i+1)+1
#   temp=data.frame(ricedata$flag_leaf_angle,tdata[,left:right])
#   temp=na.omit(temp)
#   PCt=Phase_II(temp[,1],temp[,-1],c())
#   PCname=colnames(temp)[c(2:ncol(temp))[PCt[[2]]]]
#   tempdata=rbind.fill(data.frame(tempdata),data.frame(t(PCname)))
# }
# temp=data.frame(ricedata$flag_leaf_angle,tdata[,14552:14581])
# temp=na.omit(temp)
# PCt=Phase_II(temp[,1],temp[,-1],c())
# PCname=colnames(temp)[c(2:ncol(temp))[PCt[[2]]]]
# 
# tempdata=rbind.fill(data.frame(tempdata),data.frame(t(PCname)))
# 

