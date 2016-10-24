#####################################################
# R Assignment
# Created: 2011-01-20
# Author: G. Welch
# email: gwelch@sandiego.gov
# Purpose: Assignment # 1 for CoSD R workgroup
#####################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
# does not make strings into factors - leaves them as strings

################# USER DEFINED SETUP START ###############################

#  Set your path and filename
setwd("H:\\R\\2010-12-21 Assignment")

#I didn't use this - Why bring it in twice???
#file <- "H:\\R\\2010-12-21 Assignment\\WaterQualityData.csv"

# input raw data

DATA <- read.csv("WaterQualityData.csv", header = TRUE)
#might be able to use flush=T to not look at the last two lines where SQL gives you the number of lines it pulled.

str(DATA)

# 1, 2, 3 Find all stations and dates which have exceedances for the Parmcode:
#1. TOTAL (>10000 mL/L), 
#2. FECAL (>400 mL/L),  
#3. ENTERO (>104 mL/L)
#and put these into a data frame with the columns station, sample date, depth, value

x.total <- subset(DATA, PARMCODE == "TOTAL" & DATA$VALUE > 10000, c("STATION", "SAMPLE_DATE", "DEPTH_METER", "PARMCODE", "VALUE"))
x.fecal <- subset(DATA, PARMCODE == "FECAL" & DATA$VALUE > 400, c("STATION", "SAMPLE_DATE", "DEPTH_METER", "PARMCODE", "VALUE"))
x.entero <- subset(DATA, PARMCODE == "ENTERO" & DATA$VALUE > 104, c("STATION", "SAMPLE_DATE", "DEPTH_METER", "PARMCODE", "VALUE"))
str(x.total)
str(x.fecal)
str(x.entero)

# 4. Put all of the dataframes into a list and save them into an R workspace

all.x.lst <- list(x.total, x.fecal, x.entero)
save.image("H:\\R\\2010-12-21 Assignment\\xceeds.rdata")
#to bring it back - load("H:\\R\\2010-12-21 Assignment\\xceeds.rdata")

# 5. Merge the total, fecal and entro exceedances into one dataframe

xceeds.df <- rbind(x.total, x.fecal, x.entero)
summary(xceeds.df)
str(xceeds.df)

# 6. Make new columns identifying the stations where the maximum value
#for total, fecal and entro exceedances exists

max.t <- max(xceeds.df$VALUE[xceeds.df$PARMCODE == "TOTAL" & !is.na(xceeds.df$VALUE)])
max.f <- max(xceeds.df$VALUE[xceeds.df$PARMCODE == "FECAL" & !is.na(xceeds.df$VALUE)])
max.e <- max(xceeds.df$VALUE[xceeds.df$PARMCODE == "ENTERO" & !is.na(xceeds.df$VALUE)])

xceeds.df$MaxVals <-	(xceeds.df$PARMCODE == "TOTAL" & !is.na(xceeds.df$VALUE) & xceeds.df$VALUE==max.t)|
				(xceeds.df$PARMCODE == "FECAL" & !is.na(xceeds.df$VALUE) & xceeds.df$VALUE==max.f)|
				(xceeds.df$PARMCODE == "ENTERO" & !is.na(xceeds.df$VALUE) & xceeds.df$VALUE==max.e)
summary(xceeds.df)
str(xceeds.df)


# 7. Write a .csv file that contains the total, fecal, and entro exceedances

write.csv(xceeds.df, "H:\\R\\2010-12-21 Assignment\\xceeds.csv", row.names = FALSE)

# 8. BONUS: Write a function that uses the original dataset (DATA) that checks for
#compliance of total, fecal and entero and prints one data frame with
#the station, date, depth, and value of total, fecal and entro exceedances
#Again, exceedance values are Total > 10000 mL/L, Fecal > 400 mL/L, and Entero > 104 mL/L

xceeds.fun <- function(data){

	xceeds.logic <- 	(data$PARMCODE == "TOTAL" & data$VALUE > 10000 & !is.na(data$VALUE))|
				(data$PARMCODE == "FECAL" & data$VALUE > 400 & !is.na(data$VALUE))|
				(data$PARMCODE == "ENTERO" & data$VALUE > 104 & !is.na(data$VALUE))

	xceeds.data <- data[xceeds.logic,]
				xceeds.data
	}

#To run it, do:

xceeds.fun.df <-xceeds.fun(DATA)

summary(xceeds.fun.df)
str(xceeds.fun.df)







