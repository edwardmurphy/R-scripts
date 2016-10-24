############################################################################
# Program	 	: TI.R
# Purpose	 	: Determine tolerance intervals for DS and DP for tests using 
#		   	release data
# Client	 	: aTYR Pharma
# Project	 	: Consulting
# Programmer 	: E.Murphy
# Ver/Date	 	: v0.1 02Jun2015
# Comments	 	: 
#
# Modifications	:
#############################################################################
library(xlsx)
library(tolerance)

setwd("C:\\Users\\emurphy\\Dropbox\\aTYR Pharma")

### read data ###
res.all = read.xlsx("DS_DP_results.xlsx",1)
res.all
res.nondevel = res.all[-which(res.all$QASTATUS=="Development"),]
res.nondevel
res.gmp = res.all[which(res.all$QASTATUS=="GMP"),]
res.gmp

### normal tolerance interval ###
cols = c("PH","OSMO","SEMAIN","SELOW","RPMAIN","RPTOTIMP","RPPRE2","CEXMAIN","CEXTOTIMP","CEXDIMER","ACTIVITY","CONC","PS20")

res1 = NULL
res2 = NULL
j = 0
for (i in cols){
	j=j+1
	res1[[j]] = normtol.int(na.omit(res.all[,i]),P=c(0.90,0.95,0.99),side=1)
	res2[[j]] = normtol.int(na.omit(res.all[,i]),P=c(0.90,0.95,0.99),side=2)
	out = data.frame(append(res1[[j]],res2[[j]][4:5]))
	out[is.na(out)] <- -1000
	colnames(out) <- c("alpha","prob","mean","onesidelow","onesidehigh","twosidelow","twosidehigh")
	write.xlsx(out,file="tolints25JUN2015.xlsx",sheetName=i,append=TRUE,row.names=F)
}
