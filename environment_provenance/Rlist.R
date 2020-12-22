args = commandArgs(trailingOnly=TRUE)

ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
write.table(ip, file=args[1], sep=' ', row.names=F, col.names=F, quote=F)