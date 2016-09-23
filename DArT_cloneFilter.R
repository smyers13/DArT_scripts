# by Steven Myers
# 06.vii.2016
# 
# Description: this script is for use in the DArT GBS processing pipeline.
# It should be the first filter run on DArT Report, before running the 'DArT_filter_v2.pl' filter
# it makes a list of loci to remove from the filtered DArT Rport based on a user set threshold for hamming distance (default <4) and call rate data, it then removes those loci.
# 
# input: using DArT Report file save sheet 'Report-<project_code>_2row.csv' as comma delimited (csv) file.
# 
# set the working directory
setwd("C:/path/to/data")
# increase the print limit (just in case)
options(max.print=99999999)  #increases the print capacity for demanding print requirements; required for next step

# load the filtered SNP data
data <- read.table('Report-<project_code>_2row.csv', sep=',', header=F, row.names=NULL)

#assign colnames for ease of use
colnames(data) <- as.character(unlist(data[7,]))
# remove the top 7 rows and save for later
header <- data[1:7,1:ncol(data)]
data <- data[8:nrow(data),1:ncol(data)]

#create new column containing just the clone name (not clone name and SNP position, which is what CloneID is)
data$clone1 <- sapply(data$CloneID,gsub,pattern="\\|.*",replacement="",perl=TRUE)

#extract clone duplicates (same clone different SNP)
n_occur <- data.frame(table(data$clone1))
duplicates <- data[data$clone1 %in% n_occur$Var1[n_occur$Freq > 2],]
rm(n_occur)
#sort by clone name and call rate
duplicates <- duplicates[order(duplicates$clone1,duplicates$CallRate, decreasing=T),]
#extract the reference clone from each duplicate stack with the highest call rate (because they are sorted we can use 'duplicated' which always selects the first of any ambiguous values)
#these are the clones we want to keep!
duplicates2 <- duplicates[!duplicated(duplicates$clone1),]

#create new matching variable in duplicates and data based on cloneID to be used as a reference for filtering
duplicates$truncCID <- sapply(duplicates$CloneID,gsub,pattern=":.*",replacement="",perl=TRUE)
duplicates2$truncCID <- sapply(duplicates2$CloneID,gsub,pattern=":.*",replacement="",perl=TRUE)
duplicates2$truncCID2 <- sapply(duplicates2$CloneID,gsub,pattern=":.*",replacement="",perl=TRUE)
duplicates2$truncCID2 <- sapply(duplicates2$truncCID2,gsub,pattern="--",replacement="-",perl=TRUE)

#remove duplicates to keep from duplicate list
duplicates <- duplicates[! duplicates$truncCID %in% duplicates2$truncCID,]
duplicates <- duplicates[! duplicates$truncCID %in% duplicates2$truncCID2,]

#remove loci in duplicate list from orignial data
data$truncCID <- sapply(data$CloneID,gsub,pattern=":.*",replacement="",perl=TRUE)
data <- data[! data$truncCID %in% duplicates$truncCID,]
rm(duplicates,duplicates2)

#remove (now) superfluous 'helper' columns
data <- subset(data, select=-clone1)
data <- subset(data, select=-truncCID)
#add in header
data <- rbind(header,data)

# export data in same format that DArT_filter.pl imports and exports.
write.table(data, "Report-<project_code>_2row_cloneFiltered.csv",row.names=F, col.names=F, quote=F, sep=",")


