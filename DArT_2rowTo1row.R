# DArT_2rowTo1row.R 
#
# by Steven Myers
# 07.vii.2016
# 
# this script is for use in the DArT GBS processing pipeline
# description: after running 'DArT_hamFilter.R' 
# 
# input: '_hamFiltered7.csv' output from 'DArT_hamFilter.R' and using DArT Report file save sheet 'ScoringData_SNP_1row' as comma delimited (csv) file with extension "_1row.csv"
# 
# 
# set the working directory
setwd("C:/path/to/data")
# increase the print limit (just in case)
options(max.print=99999999)  #increases the print capacity for demanding print requirements; required for next step

# load the 1row SNP data
data <- read.table("Report-<project_code>_1row.csv", sep=",", header=FALSE, row.names=NULL) 
# load the 2row SNP data
data2 <- read.table('Report-<project_code>_2row.csv', sep=',', header=F, row.names=NULL)
# load the ham filtered SNP data
filteredData <- read.table("Report-<project_code>_2Row_hamFiltered7.csv", sep=",", header=FALSE, row.names=NULL)

#assign colnames for ease of use
colnames(data) <- as.character(unlist(data[7,]))   
colnames(filteredData) <- as.character(unlist(data2[7,]))
rm(data2)

# remove the top 7 rows and save for later
header <- data[1:7,1:ncol(data)]
data <- data[8:nrow(data),1:ncol(data)]

# extract from data CloneIDs that match CloneIDs in filteredData
data <- data[ data$CloneID %in% filteredData$CloneID,]

#add header back in (although we'll probably just be removing this later anyway, lol)
data <- rbind(header,data)
rm(filteredData,header)

# export data in same format that DArT_filter.pl imports and exports.
write.table(data, "Report-<project_code>_1row_hamFiltered.csv",row.names=F, col.names=F, quote=F, sep=",")
