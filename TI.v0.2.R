############################################################################
# Program	 	: TI.R
# Purpose	 	: Determine tolerance intervals for DS and DP for tests using 
#		   	  release data
# Client	 	: aTYR Pharma
# Project	 	: Consulting
# Programmer 	: E.Murphy
# Ver/Date	 	: v0.1 02Jun2015
#			  v0.2 28Jun2015
# Comments	 	: 
#
# Modifications	:
#     E.Murphy; 28Jun2015; 
#       - added wrapper function to handle multiple datasets with one call
#############################################################################

# Preliminaries
library(xlsx,lib.loc="C:\\Users\\emurphy\\Dropbox\\R\\R-3.2.0\\library")
library(tolerance,lib.loc="C:\\Users\\emurphy\\Dropbox\\R\\R-3.2.0\\library")

setwd("C:\\Users\\emurphy\\Dropbox\\aTYR Pharma")

############## Read data ##################
# read data
res.all = read.xlsx("DS_DP_results.xlsx",1)
###########################################

################## Subsets ###################
# subsets (DS/DP)
res.nondevel = res.all[-which(res.all$QASTATUS=="Development"),]
res.gmp = res.all[which(res.all$QASTATUS=="GMP"),]
# subsets (DS only)
res.ds.all = res.all[which(res.all$TYPE=="DS"),]
res.ds.nondevel = res.nondevel[which(res.nondevel$TYPE=="DS"),]
res.ds.gmp = res.gmp[which(res.gmp$TYPE=="DS"),]
##############################################

############### Datasets for Analysis ################
## NAs retained, DS/DP pooled
# all lots
res.all	
# Non-GMP and GMP (no Development)
res.nondevel
# GMP
res.gmp

## DS only, NAs retained
# all DS lots
res.ds.all
# Non-GMP and GMP DS lots
res.ds.nondevel
# GMP DS lots
res.ds.gmp

### all columns 
write.table(t(colnames(res.all)),file="clipboard",row.names=F,col.names=F,sep=",")
# "TYPE","LOT","QASTATUS","PH","OSMO","SEMAIN","SEHIGH","SELOW","RPMAIN","RPTOTIMP","RPPRE2",
# "CEXMAIN","CEXTOTIMP","CEXDIMER","HCP","RESIDDNA","LAL","BIOBURDEN","ACTIVITY","CONC",
# "PS20","KANA","IPTG","PEI"
###

# columns for analysis 
cols = c("PH","OSMO","SEMAIN","SELOW","RPMAIN",
	   "RPTOTIMP","RPPRE2","CEXMAIN","CEXTOTIMP",
	   "CEXDIMER","ACTIVITY","CONC","PS20")

# put all datasets in list
res.list = list(all=res.all[,cols],nondevel=res.nondevel[,cols],gmp=res.gmp[,cols],
		     dsall = res.ds.all[,cols],dsnondevel=res.ds.nondevel[,cols],dsgmp=res.ds.gmp[,cols])

######################################################

############### FUNCTION TO CALCULATE TOLERANCE INTERVALS AND OUTPUT TO XLSX FORMAT ###############
NormInt <- function(INDATA){
  # get one-sided and two-sided tolerance intervals
  oneside = rapply(na.omit(INDATA),normtol.int,how="list",P=c(0.90,0.95,0.99),side=1)
  twoside = rapply(na.omit(INDATA),normtol.int,how="list",P=c(0.90,0.95,0.99),side=2)
  # system date for filename
  d = format(Sys.Date(),format="%d%b%Y")

  # use sapply function to loop through each member of list, then through another call to sapply
  #  through each submember (column)

  OUTDATA <- sapply(names(oneside),function(x) {
    # looping through each dataset (which are names of the list "oneside")
    sapply(attributes(oneside[[x]])$names, function(y){
      # looping through each column in each dataset
      tmp = cbind(oneside[[x]][[y]],twoside[[x]][[y]][,4:5])
	colnames(tmp) <- c("alpha","prob","mean","onesidelow","onesidehigh","twosidelow","twosidehigh")
      write.xlsx(tmp,file=paste0(x,'_tolints_',d,'.xlsx'),sheetName=y,append=TRUE,row.names=F,showNA=F)
      }
    , simplify=FALSE,USE.NAMES=TRUE)
    }
  ,simplify=FALSE,USE.NAMES=TRUE)

  return(OUTDATA)
}
###################################################################################################

########## Call function ##########

tolints = NormInt(res.list)




