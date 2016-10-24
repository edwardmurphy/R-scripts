rm(list=ls(all=TRUE)) 
#set(stringsAsChar=T  not right look it up

################# USER DEFINED SETUP START ###############################

#  Set your path and filename
file <- "H:\\Jen\\R\\RUserGroup\\Practice1\\WaterQualityData.csv"

# input raw data

DATA<-read.csv(file, header=T)  #might be able to use flush=T to not look at the last two lines where SQL gives you the number of lines it pulled.

# 1. Find all stations and dates which have exceedances for the Parmcode TOTAL (>10000 mL/L) and put these into a data frame with the columns station, sample date, depth, value 

# 2. Find all stations and dates which have exceedances for Parmcode FECAL (>400 mL/L) and put these into a data frame with the columns station, sample date, depth, value 

# 3. Find all stations and dates which have exceedances for Parmcode ENTERO (>104 mL/L) and put these into a data frame with the columns station, sample date, depth, value 

# 4. Put all of the dataframes into a list and save them into an R workspace

# 5. Merge the total, fecal and entro exceedances into one dataframe

# 6. Make new columns identifying the stations where the maximum value for total, fecal and entro exceedances exists

# 7. Write a .csv file that contains the total, fecal, and entro exceedances

# 8. BONUS: Write a function that uses the original dataset (DATA) that checks for compliance of total, fecal and entero and prints one data frame with the station, date, depth, and value of total, fecal and entro exceedances
#    Again, exceedance values are Total > 10000 mL/L, Fecal > 400 mL/L, and Entero > 104 mL/L   


