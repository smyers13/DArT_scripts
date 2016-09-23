# DArTtoR_v2.R
#
# written by Steven Myers 
# 4.vii.2016
#
# modified 7.vii.2016 to remove the need to use excel
#
# description: convert DArT report data (raw) to a format that will feed into my analysis R-scripts (e.g. 'RAD_analysis_pipeline_project.R' and 'prelimPCA.R')
#
# input: DArT submission file and output from 'DArT_2rowTo1row.R' 

# increase the print limit (just in case)
options(max.print=99999999)  

setwd("C:/path/to/DArT/submission/document")
sub <- read.table("Diversity_Arrays_Technology-<custom_descriptive_name>.csv", sep=",", header=FALSE, row.names=NULL) 

setwd("C:/path/to/data")
data <- read.table("Report-<project_code>_1row_hamFiltered.csv", sep=",", header=FALSE, row.names=NULL) 

# remove superfluous header rows from data
data <- data[7:nrow(data),]
# remove superfluous columns from data
data <- data[,c(1,23:ncol(data))]

# reformat cloneID to remove everything after and including first pipe "|", and rename variable as "locus"
data[,1] <- gsub("\\|.*","",data[,1])

# transpose data #ready to rename samples
data <- as.data.frame(t(data))

# add some column names
colnames(data)[1] <- "Genotype"
colnames(sub) <- as.character(unlist(sub[1,]))

# remove any samples from sub not in data (should be able to match sub$Genotype and data[,1])
sub <- sub[sub$Genotype %in% data$Genotype,]

# make the new variable in sub combining "plate_RowColumn_Genotype". (perhaps the cuttlefish script can aid with this)
RowCol <- as.data.frame(paste(sub$Row,sub$Column,sep=""))
sub$SampleID <- paste(sub$Plate,RowCol[,1],sub$Genotype,sep="_")
rm(RowCol)

# sort data and sub by Genotype so they match (first remove header row from data)
header <- as.data.frame(data[1,])
data <- data[2:nrow(data),1:ncol(data)]

data <- data[order(data$Genotype),]
sub <- sub[order(sub$Genotype),]
data$Genotype <- sub$SampleID
rm(sub)

#add header back in and transpose
data <- rbind(header,data)
data <- as.data.frame(t(data))
rm(header)


#### convert data ####

# make header column names and remove
colnames(data) <- as.character(unlist(data[1,]))
data <- data[2:nrow(data),1:ncol(data)]

#Make sure there are no columns with contents including NA. 
data[data=="-"]<-NA

rownames(data) = data[1:nrow(data),1]
#data = t(data[1:nrow(data),2:ncol(data)])   # don't traspose now because we want less columns than rows as it makes writing and reading the data a lot quicker
data <- data[1:nrow(data),2:ncol(data)]

write.csv(data, "data_Filtered1row.csv")       
data = read.table('data_Filtered1row.csv', sep=',', header=F, row.names=NULL, stringsAsFactors=FALSE) # Load the data
colnames(data) <- as.character(unlist(data[1,]))
data <- t(data[2:nrow(data),1:ncol(data)])
colnames(data) <- as.character(unlist(data[1,]))
data <- data[2:nrow(data),1:ncol(data)]


# --> Now it's ready (in the console) for use with my other pipelines, e.g. "RAD_analysis_pipeline_project.R" and "prelimPCA.R".