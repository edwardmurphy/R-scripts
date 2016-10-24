library(affy)
library(limma)
library(gplots)
library(xtable)
library(BHC)

setwd("/home/ed/Dropbox/MATH365/Project")
load("ketoData.RData")

annota <- read.csv("RAE230A.csv", header=T, sep="\t")
paper.genes <- read.table("paperGenes.txt", header=T)


