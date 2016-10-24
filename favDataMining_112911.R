library(randomForest)

data<-read.delim("C:\\Users\\Owner\\Documents\\Favrille\\Fav.txt",header=TRUE)
data<-data[,-c(2,11,24,25,26)]	#remove date index,C2,D2,D3,fill vol,

############################################################################
> colnames(data)
 [1] "LotNumber"          "IgClassHCFam"       "IgClassHCGene"     
 [4] "IgClassLCFam"       "IgClassLCGene"      "NumHavestBags"     
 [7] "Amp1ELISA"          "AvgProteinHarvBags" "C1A280"            
[10] "C4SECPurity"        "RT"                 "C6SDSRed"          
[13] "C7SDSNonRed"        "C10pH"              "HC1PropDetected"   
[16] "HC2PropDetected"    "HC3PropDetected"    "LC1PropDetected"   
[19] "LC2PropDetected"    "LC3PropDetected"    "D1BCA"             
[22] "D5ParticleSize"     "D6pH"               "D7Osmo"            
[25] "D9Endo"             "D10ResGlut"         "pI"                
[28] "MolWtHC"            "MolWtLC"            "GlycoSitesHC"      
[31] "GlycoSitesLC"       "LysineRes"          "ExtCoeff" 
############################################################################

rf.tune <- tuneRF(data[ , ], 
   baseball[,22], ntreeTry=500, stepFactor=2, 
    improve=0.05, trace=TRUE, plot=TRUE, dobest=FALSE)
bb.rf