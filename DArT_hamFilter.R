# by Steven Myers
# 05.vii.2016
# 
# this script is for use in the DArT GBS processing pipeline
# after running 'DArT_filter_v2.pl' script it makes a list of loci to remove from the filtered DArT Rport based on a user set threshold for hamming distance (default <4) and call rate data, it then removes those loci.
# 
# input: '.ham_dist_values' and '_filtered.csv' output from 'DArT_filter_v2.pl' or 'DArT_filter2.pl' scripts
# 
# 
# A POTENTIAL PROBLEM WITH THIS METHOD IS THAT FOR LOCI WITH SMALL HAMMING DISTANCES TO MORE THAN ONE LOCUS, LOCI BEYOND THE FIRST MAY GET MISSED.
# --> A WORKAROUND FOR THIS WOULD BE TO RUN THE DArT_filter.pl again
# I have implemented this workaround, but after 7 passes there are still some clone pairs with hamming distances of 0!! However, these are in very small number and are likely to be removed by downstream processes (e.g. linkage analysis)
# I am confident that this routine is valuable, and technically sound; it's just a question of how easily it could be improved, and -- The most obvious way it could be improved is reducing the number of times required to run the perl script!
# 
# set the working directory
setwd("C:/path/to/data")
# increase the print limit (just in case)
options(max.print=99999999)  #increases the print capacity for demanding print requirements; required for next step

# set the minimum hamming distance below which loci will be removed
minHam = 5

# load the hamming distance data
hamList <- read.table("Report-<project_code>_2row_cloneFiltered.ham_dist_values", sep=",", col.names=FALSE, fill=FALSE, strip.white=TRUE) 

# load the filtered SNP data
filteredData <- read.table('Report-<project_code>_2row_cloneFiltered_filtered.csv', sep=',', header=TRUE, row.names=1) # whoops, I probably should've used row.names = F, but the script is written to fit this now (it just would've saved a step or two near the end)

hamData <- as.data.frame(sapply(hamList,gsub,pattern="\\).*",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
hamData <- as.data.frame(sapply(hamData,gsub,pattern=".*\\(",replacement="",perl=TRUE)) # replaces everything before last reverse bracket ")" with nothing ""

# these yielded undesired outcomes but kept for posterity
#hamData <- as.data.frame(sapply(hamList,gsub,pattern=".*\\(",replacement="")) # replaces everything before and including last forward bracket "(" with nothing ""
#hamData <- as.data.frame(sapply(hamData,gsub,pattern="\\)",replacement="")) # replaces all instances of reverse bracket ")" with nothing ""

#locus1 <- as.data.frame(sapply(hamList,gsub,pattern="\\|.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""; turns out the numeric value after |F| is part of the identifier
locus1 <- as.data.frame(sapply(hamList,gsub,pattern=":.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""

# I want to keep everything after 3rd ":" (replace everything before 3rd) or 
#kept for posterity
#locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(:.*)",replacement="",perl=TRUE)) # replaces everything after and including first colon ":" with nothing ""
#locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(.*:)",replacement="",perl=TRUE)) # replaces everything before and including last colon ":" with nothing ""
#locus2 <- as.data.frame(sapply(hamList,gsub,pattern="([^:])",replacement="",perl=TRUE)) # replaces everything except colons ":" with nothing ""
locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern="[^:]*:[^:]*:[^:]*:",replacement="",perl=TRUE)) # replaces everything up to and including the third colon ":" with nothing ""
#locus2 <- as.data.frame(sapply(locus2,gsub,pattern="\\|.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=":.*",replacement="")) # replaces everything after and including first colon ":" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=" ",replacement=""))   # replaces space in first position with nothing


# combine the three lists we generated into a single data frame
data <- cbind(locus1,locus2,hamData)
rm(locus1,locus2,hamData)

# format
colnames(data) <- c("locus1","locus2","hamData")
data$hamData <- as.numeric(data$hamData)


# no need to sort but kept for posterity
#data <- data[order(data$hamData)]   # sort by hamming distance in ascending order (i.e. smallest first)

# extract locus pairs with hamming distances less than user defined minimum hamming distance
shortList <- subset(data, data$hamData < minHam)  # this doesn't work

#prepare filteredData for use by adding a locus column to the end that matches locus IDs used in newList
filteredData$loci <- gsub(":.*","",row.names(filteredData))

# for each pair of loci identify the locus with the highest call rate and remove it from that pair in newList

# 1. make new column for call rates for locus1
shortList$locus1CR <- shortList$locus1 %in% filteredData$loci
r = 1

# 2. find matching locus "locus1 %in% loci" pull out associated call rate "filteredData$CallRate" and add to new column
repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,1],]  
  shortList[r,4] <- CR[1,4]
  
    r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

#do the same for locus 2
shortList$locus2CR <- shortList$locus2 %in% filteredData$loci
r = 1

repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,2],]  
  shortList[r,5] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

rm(CR,r,minHam)
# duplicate list
shortList1 <- shortList
shortList2 <- shortList

# from list 1 remove samples where locus1CR is greater than or equal to locus2CR
shortList1 <- subset(shortList1, shortList1$locus1CR >= shortList1$locus2CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus2 <- as.data.frame(factor(shortList1$locus2))
colnames(locus2) <- "locus"

shortList2 <- subset(shortList2, shortList2$locus2CR > shortList2$locus1CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus1 <- as.data.frame(factor(shortList2$locus1))
colnames(locus1) <- "locus"

#join the two lists
locus <- rbind(locus1,locus2)
rm(locus1,locus2,shortList1,shortList2,shortList,data,hamList)

#remove duplicates
locus <- unique(locus)

# remove all loci identified as sharing sequence similarity with other loci but having lower call rate from filteredData
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]
#now you have to remove their counterparts, as this is 2row_SNP data!!
locus <- as.data.frame(sapply(locus,gsub,pattern="-",replacement="--"))
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]

rm(locus)

# format data for export
filteredData <- subset(filteredData, select=-loci )  # remove helper column (loci) from filteredData
filteredData$CloneID <- row.names(filteredData)   # add CloneID column (because otherwise it won't have a header)
filteredData <- filteredData[c(ncol(filteredData),1:(ncol(filteredData)-1))]   # move CloneID to first column
names(filteredData) <- sub("X","",names(filteredData))   # remove stupid Xs from the headers

#not sure if this next part is even necessary now; it may have just been that I need "quote=F" in write table.
#filteredData2 <- read.table("Report-<project_code>.csv", sep=",") 
#filteredData2 <- filteredData2[-c(8:nrow(filteredData2)),]
#colnames(filteredData2) <- colnames(filteredData)   # give it the same column names as data
#filteredData2[c(nrow(filteredData2)),]<-colnames(filteredData2)   # fill the last row with column headers
#filteredData <- rbind(filteredData2,filteredData)   #append it to the top of the data
#rm(filteredData2)
# yep -- I didn't need this bit!

# export data in same format that DArT_filter.pl imports and exports.
write.table(filteredData, "Report-<project_code>_2Row_hamFiltered.csv",row.names=F, col.names=F, quote=F, sep=",")


###############################################################



# I've noticed that even after 7 cycles of this and the perl filter there are still SNPs called from the same clone!!! I should definitely remove duplicate clones in a first pass (and then filter with perl again before applying this filter!!)




# After running this filter it is necessary to pass through the perl script again, as due to the output of the hamming distances it makes it difficult to capture all at once.
# so run DArT_filter.pl again making sure you use the latest filtered version "Report-<project_code>_hamFiltered.csv" as the input.
# then run through this filter again -- below:

# set the minimum hamming distance below which loci will be removed
minHam = 5

# load the hamming distance data
hamList <- read.table("Report-<project_code>_2Row_hamFiltered.ham_dist_values", sep=",", col.names=FALSE, fill=FALSE, strip.white=TRUE) 

# load the filtered SNP data
filteredData <- read.table('Report-<project_code>_2Row_hamFiltered_filtered.csv', sep=',', header=TRUE, row.names=1) # whoops, I probably should've used row.names = F, but the script is written to fit this now (it just would've saved a step or two near the end)

hamData <- as.data.frame(sapply(hamList,gsub,pattern="\\).*",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
hamData <- as.data.frame(sapply(hamData,gsub,pattern=".*\\(",replacement="",perl=TRUE)) # replaces everything before last reverse bracket ")" with nothing ""

locus1 <- as.data.frame(sapply(hamList,gsub,pattern=":.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""

locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern="[^:]*:[^:]*:[^:]*:",replacement="",perl=TRUE)) # does something weird, need to investigate... It seems it works in some cases (where there is only one pair) but not others (where there is more than one pair)
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=":.*",replacement="")) # replaces everything after and including first pipe "|" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=" ",replacement=""))   # replaces space in first position with nothing

# combine the three lists we generated into a single data frame
data <- cbind(locus1,locus2,hamData)
rm(locus1,locus2,hamData)

# format
colnames(data) <- c("locus1","locus2","hamData")
data$hamData <- as.numeric(data$hamData)

# extract locus pairs with hamming distances less than user defined minimum hamming distance
shortList <- subset(data, data$hamData < minHam)  # this doesn't work

#prepare filteredData for use by adding a locus column to the end that matches locus IDs used in newList
filteredData$loci <- gsub(":.*","",row.names(filteredData))

# for each pair of loci identify the locus with the highest call rate and remove it from that pair in newList
# 1. make new column for call rates for locus1
shortList$locus1CR <- shortList$locus1 %in% filteredData$loci
r = 1

# 2. find matching locus "locus1 %in% loci" pull out associated call rate "filteredData$CallRate" and add to new column
repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,1],]  
  shortList[r,4] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

#do the same for locus 2
shortList$locus2CR <- shortList$locus2 %in% filteredData$loci
r = 1

repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,2],]  
  shortList[r,5] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

rm(CR,r,minHam)
# duplicate list
shortList1 <- shortList
shortList2 <- shortList

# from list 1 remove samples where locus1CR is greater than or equal to locus2CR
shortList1 <- subset(shortList1, shortList1$locus1CR >= shortList1$locus2CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus2 <- as.data.frame(factor(shortList1$locus2))
colnames(locus2) <- "locus"

shortList2 <- subset(shortList2, shortList2$locus2CR > shortList2$locus1CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus1 <- as.data.frame(factor(shortList2$locus1))
colnames(locus1) <- "locus"

#join the two lists
locus <- rbind(locus1,locus2)
rm(locus1,locus2,shortList1,shortList2,shortList,data,hamList)

#remove duplicates
locus <- unique(locus)

# remove all loci identified as sharing sequence similarity with other loci but having lower call rate from filteredData
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]
#now you have to remove their counterparts, as this is 2row_SNP data!!
locus <- as.data.frame(sapply(locus,gsub,pattern="-",replacement="--"))
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]

rm(locus)

# format data for export
filteredData <- subset(filteredData, select=-loci )  # remove helper column (loci) from filteredData
filteredData$CloneID <- row.names(filteredData)   # add CloneID column (because otherwise it won't have a header)
filteredData <- filteredData[c(ncol(filteredData),1:(ncol(filteredData)-1))]   # move CloneID to first column
names(filteredData) <- sub("X","",names(filteredData))   # remove stupid Xs from the headers

# export data in same format that DArT_filter.pl imports and exports.
write.table(filteredData, "Report-<project_code>_2Row_hamFiltered2.csv",row.names=F, col.names=F, quote=F, sep=",")



######################################################


# After running this filter it is necessary to pass through the perl script again, as due to the output of the hamming distances it makes it difficult to capture all at once.
# so run DArT_filter.pl again making sure you use the latest filtered version "Report-<project_code>_hamFiltered.csv" as the input.
# then run through this filter again -- below:

# set the minimum hamming distance below which loci will be removed
minHam = 5

# load the hamming distance data
hamList <- read.table("Report-<project_code>_2Row_hamFiltered2.ham_dist_values", sep=",", col.names=FALSE, fill=FALSE, strip.white=TRUE) 

# load the filtered SNP data
filteredData <- read.table('Report-<project_code>_2Row_hamFiltered2_filtered.csv', sep=',', header=TRUE, row.names=1) # whoops, I probably should've used row.names = F, but the script is written to fit this now (it just would've saved a step or two near the end)

hamData <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
hamData <- as.data.frame(sapply(hamData,gsub,pattern="(.*\\()",replacement="",perl=TRUE)) # replaces everything before last reverse bracket ")" with nothing ""

locus1 <- as.data.frame(sapply(hamList,gsub,pattern=":.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""

locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern="[^:]*:[^:]*:[^:]*:",replacement="",perl=TRUE)) # does something weird, need to investigate... It seems it works in some cases (where there is only one pair) but not others (where there is more than one pair)
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=":.*",replacement="")) # replaces everything after and including first pipe "|" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=" ",replacement=""))   # replaces space in first position with nothing

# combine the three lists we generated into a single data frame
data <- cbind(locus1,locus2,hamData)
rm(locus1,locus2,hamData)

# format
colnames(data) <- c("locus1","locus2","hamData")
data$hamData <- as.numeric(data$hamData)

# extract locus pairs with hamming distances less than user defined minimum hamming distance
shortList <- subset(data, data$hamData < minHam)  # this doesn't work

#prepare filteredData for use by adding a locus column to the end that matches locus IDs used in newList
filteredData$loci <- gsub(":.*","",row.names(filteredData))

# for each pair of loci identify the locus with the highest call rate and remove it from that pair in newList
# 1. make new column for call rates for locus1
shortList$locus1CR <- shortList$locus1 %in% filteredData$loci
r = 1

# 2. find matching locus "locus1 %in% loci" pull out associated call rate "filteredData$CallRate" and add to new column
repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,1],]  
  shortList[r,4] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

#do the same for locus 2
shortList$locus2CR <- shortList$locus2 %in% filteredData$loci
r = 1

repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,2],]  
  shortList[r,5] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

rm(CR,r,minHam)
# duplicate list
shortList1 <- shortList
shortList2 <- shortList

# from list 1 remove samples where locus1CR is greater than or equal to locus2CR
shortList1 <- subset(shortList1, shortList1$locus1CR >= shortList1$locus2CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus2 <- as.data.frame(factor(shortList1$locus2))
colnames(locus2) <- "locus"

shortList2 <- subset(shortList2, shortList2$locus2CR > shortList2$locus1CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus1 <- as.data.frame(factor(shortList2$locus1))
colnames(locus1) <- "locus"

#join the two lists
locus <- rbind(locus1,locus2)
rm(locus1,locus2,shortList1,shortList2,shortList,data,hamList)

#remove duplicates
locus <- unique(locus)

# remove all loci identified as sharing sequence similarity with other loci but having lower call rate from filteredData
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]
#now you have to remove their counterparts, as this is 2row_SNP data!!
locus <- as.data.frame(sapply(locus,gsub,pattern="-",replacement="--"))
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]

rm(locus)

# format data for export
filteredData <- subset(filteredData, select=-loci )  # remove helper column (loci) from filteredData
filteredData$CloneID <- row.names(filteredData)   # add CloneID column (because otherwise it won't have a header)
filteredData <- filteredData[c(ncol(filteredData),1:(ncol(filteredData)-1))]   # move CloneID to first column
names(filteredData) <- sub("X","",names(filteredData))   # remove stupid Xs from the headers

# export data in same format that DArT_filter.pl imports and exports.
write.table(filteredData, "Report-<project_code>_2Row_hamFiltered3.csv",row.names=F, col.names=F, quote=F, sep=",")



##############################################################

# After running this filter it is necessary to pass through the perl script again, as due to the output of the hamming distances it makes it difficult to capture all at once.
# so run DArT_filter.pl again making sure you use the latest filtered version "Report-<project_code>_hamFiltered.csv" as the input.
# then run through this filter again -- below:

# set the minimum hamming distance below which loci will be removed
minHam = 5

# load the hamming distance data
hamList <- read.table("Report-<project_code>_2Row_hamFiltered3.ham_dist_values", sep=",", col.names=FALSE, fill=FALSE, strip.white=TRUE) 

# load the filtered SNP data
filteredData <- read.table('Report-<project_code>_2Row_hamFiltered3_filtered.csv', sep=',', header=TRUE, row.names=1) # whoops, I probably should've used row.names = F, but the script is written to fit this now (it just would've saved a step or two near the end)

hamData <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
hamData <- as.data.frame(sapply(hamData,gsub,pattern="(.*\\()",replacement="",perl=TRUE)) # replaces everything before last reverse bracket ")" with nothing ""

locus1 <- as.data.frame(sapply(hamList,gsub,pattern=":.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""

locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern="[^:]*:[^:]*:[^:]*:",replacement="",perl=TRUE)) # does something weird, need to investigate... It seems it works in some cases (where there is only one pair) but not others (where there is more than one pair)
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=":.*",replacement="")) # replaces everything after and including first pipe "|" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=" ",replacement=""))   # replaces space in first position with nothing

# combine the three lists we generated into a single data frame
data <- cbind(locus1,locus2,hamData)
rm(locus1,locus2,hamData)

# format
colnames(data) <- c("locus1","locus2","hamData")
data$hamData <- as.numeric(data$hamData)

# extract locus pairs with hamming distances less than user defined minimum hamming distance
shortList <- subset(data, data$hamData < minHam)  # this doesn't work

#prepare filteredData for use by adding a locus column to the end that matches locus IDs used in newList
filteredData$loci <- gsub(":.*","",row.names(filteredData))

# for each pair of loci identify the locus with the highest call rate and remove it from that pair in newList
# 1. make new column for call rates for locus1
shortList$locus1CR <- shortList$locus1 %in% filteredData$loci
r = 1

# 2. find matching locus "locus1 %in% loci" pull out associated call rate "filteredData$CallRate" and add to new column
repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,1],]  
  shortList[r,4] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

#do the same for locus 2
shortList$locus2CR <- shortList$locus2 %in% filteredData$loci
r = 1

repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,2],]  
  shortList[r,5] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

rm(CR,r,minHam)
# duplicate list
shortList1 <- shortList
shortList2 <- shortList

# from list 1 remove samples where locus1CR is greater than or equal to locus2CR
shortList1 <- subset(shortList1, shortList1$locus1CR >= shortList1$locus2CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus2 <- as.data.frame(factor(shortList1$locus2))
colnames(locus2) <- "locus"

shortList2 <- subset(shortList2, shortList2$locus2CR > shortList2$locus1CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus1 <- as.data.frame(factor(shortList2$locus1))
colnames(locus1) <- "locus"

#join the two lists
locus <- rbind(locus1,locus2)
rm(locus1,locus2,shortList1,shortList2,shortList,data,hamList)

#remove duplicates
locus <- unique(locus)

# remove all loci identified as sharing sequence similarity with other loci but having lower call rate from filteredData
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]
#now you have to remove their counterparts, as this is 2row_SNP data!!
locus <- as.data.frame(sapply(locus,gsub,pattern="-",replacement="--"))
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]

rm(locus)

# format data for export
filteredData <- subset(filteredData, select=-loci )  # remove helper column (loci) from filteredData
filteredData$CloneID <- row.names(filteredData)   # add CloneID column (because otherwise it won't have a header)
filteredData <- filteredData[c(ncol(filteredData),1:(ncol(filteredData)-1))]   # move CloneID to first column
names(filteredData) <- sub("X","",names(filteredData))   # remove stupid Xs from the headers

# export data in same format that DArT_filter.pl imports and exports.
write.table(filteredData, "Report-<project_code>_2Row_hamFiltered4.csv",row.names=F, col.names=F, quote=F, sep=",")



#################

# After running this filter it is necessary to pass through the perl script again, as due to the output of the hamming distances it makes it difficult to capture all at once.
# so run DArT_filter.pl again making sure you use the latest filtered version "Report-<project_code>_hamFiltered.csv" as the input.
# then run through this filter again -- below:

# set the minimum hamming distance below which loci will be removed
minHam = 5

# load the hamming distance data
hamList <- read.table("Report-<project_code>_2Row_hamFiltered4.ham_dist_values", sep=",", col.names=FALSE, fill=FALSE, strip.white=TRUE) 

# load the filtered SNP data
filteredData <- read.table('Report-<project_code>_2Row_hamFiltered4_filtered.csv', sep=',', header=TRUE, row.names=1) # whoops, I probably should've used row.names = F, but the script is written to fit this now (it just would've saved a step or two near the end)

hamData <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
hamData <- as.data.frame(sapply(hamData,gsub,pattern="(.*\\()",replacement="",perl=TRUE)) # replaces everything before last reverse bracket ")" with nothing ""

locus1 <- as.data.frame(sapply(hamList,gsub,pattern=":.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""

locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern="[^:]*:[^:]*:[^:]*:",replacement="",perl=TRUE)) # does something weird, need to investigate... It seems it works in some cases (where there is only one pair) but not others (where there is more than one pair)
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=":.*",replacement="")) # replaces everything after and including first pipe "|" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=" ",replacement=""))   # replaces space in first position with nothing

# combine the three lists we generated into a single data frame
data <- cbind(locus1,locus2,hamData)
rm(locus1,locus2,hamData)

# format
colnames(data) <- c("locus1","locus2","hamData")
data$hamData <- as.numeric(data$hamData)

# extract locus pairs with hamming distances less than user defined minimum hamming distance
shortList <- subset(data, data$hamData < minHam)  # this doesn't work

#prepare filteredData for use by adding a locus column to the end that matches locus IDs used in newList
filteredData$loci <- gsub(":.*","",row.names(filteredData))

# for each pair of loci identify the locus with the highest call rate and remove it from that pair in newList
# 1. make new column for call rates for locus1
shortList$locus1CR <- shortList$locus1 %in% filteredData$loci
r = 1

# 2. find matching locus "locus1 %in% loci" pull out associated call rate "filteredData$CallRate" and add to new column
repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,1],]  
  shortList[r,4] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

#do the same for locus 2
shortList$locus2CR <- shortList$locus2 %in% filteredData$loci
r = 1

repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,2],]  
  shortList[r,5] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

rm(CR,r,minHam)
# duplicate list
shortList1 <- shortList
shortList2 <- shortList

# from list 1 remove samples where locus1CR is greater than or equal to locus2CR
shortList1 <- subset(shortList1, shortList1$locus1CR >= shortList1$locus2CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus2 <- as.data.frame(factor(shortList1$locus2))
colnames(locus2) <- "locus"

shortList2 <- subset(shortList2, shortList2$locus2CR > shortList2$locus1CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus1 <- as.data.frame(factor(shortList2$locus1))
colnames(locus1) <- "locus"

#join the two lists
locus <- rbind(locus1,locus2)
rm(locus1,locus2,shortList1,shortList2,shortList,data,hamList)

#remove duplicates
locus <- unique(locus)

# remove all loci identified as sharing sequence similarity with other loci but having lower call rate from filteredData
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]
#now you have to remove their counterparts, as this is 2row_SNP data!!
locus <- as.data.frame(sapply(locus,gsub,pattern="-",replacement="--"))
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]

rm(locus)

# format data for export
filteredData <- subset(filteredData, select=-loci )  # remove helper column (loci) from filteredData
filteredData$CloneID <- row.names(filteredData)   # add CloneID column (because otherwise it won't have a header)
filteredData <- filteredData[c(ncol(filteredData),1:(ncol(filteredData)-1))]   # move CloneID to first column
names(filteredData) <- sub("X","",names(filteredData))   # remove stupid Xs from the headers

# export data in same format that DArT_filter.pl imports and exports.
write.table(filteredData, "Report-<project_code>_2Row_hamFiltered5.csv",row.names=F, col.names=F, quote=F, sep=",")



#############

# After running this filter it is necessary to pass through the perl script again, as due to the output of the hamming distances it makes it difficult to capture all at once.
# so run DArT_filter.pl again making sure you use the latest filtered version "Report-<project_code>_hamFiltered.csv" as the input.
# then run through this filter again -- below:

# set the minimum hamming distance below which loci will be removed
minHam = 5

# load the hamming distance data
hamList <- read.table("Report-<project_code>_2Row_hamFiltered5.ham_dist_values", sep=",", col.names=FALSE, fill=FALSE, strip.white=TRUE) 

# load the filtered SNP data
filteredData <- read.table('Report-<project_code>_2Row_hamFiltered5_filtered.csv', sep=',', header=TRUE, row.names=1) # whoops, I probably should've used row.names = F, but the script is written to fit this now (it just would've saved a step or two near the end)

hamData <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
hamData <- as.data.frame(sapply(hamData,gsub,pattern="(.*\\()",replacement="",perl=TRUE)) # replaces everything before last reverse bracket ")" with nothing ""

locus1 <- as.data.frame(sapply(hamList,gsub,pattern=":.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""

locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern="[^:]*:[^:]*:[^:]*:",replacement="",perl=TRUE)) # does something weird, need to investigate... It seems it works in some cases (where there is only one pair) but not others (where there is more than one pair)
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=":.*",replacement="")) # replaces everything after and including first pipe "|" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=" ",replacement=""))   # replaces space in first position with nothing

# combine the three lists we generated into a single data frame
data <- cbind(locus1,locus2,hamData)
rm(locus1,locus2,hamData)

# format
colnames(data) <- c("locus1","locus2","hamData")
data$hamData <- as.numeric(data$hamData)

# extract locus pairs with hamming distances less than user defined minimum hamming distance
shortList <- subset(data, data$hamData < minHam)  # this doesn't work

#prepare filteredData for use by adding a locus column to the end that matches locus IDs used in newList
filteredData$loci <- gsub(":.*","",row.names(filteredData))

# for each pair of loci identify the locus with the highest call rate and remove it from that pair in newList
# 1. make new column for call rates for locus1
shortList$locus1CR <- shortList$locus1 %in% filteredData$loci
r = 1

# 2. find matching locus "locus1 %in% loci" pull out associated call rate "filteredData$CallRate" and add to new column
repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,1],]  
  shortList[r,4] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

#do the same for locus 2
shortList$locus2CR <- shortList$locus2 %in% filteredData$loci
r = 1

repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,2],]  
  shortList[r,5] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

rm(CR,r,minHam)
# duplicate list
shortList1 <- shortList
shortList2 <- shortList

# from list 1 remove samples where locus1CR is greater than or equal to locus2CR
shortList1 <- subset(shortList1, shortList1$locus1CR >= shortList1$locus2CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus2 <- as.data.frame(factor(shortList1$locus2))
colnames(locus2) <- "locus"

shortList2 <- subset(shortList2, shortList2$locus2CR > shortList2$locus1CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus1 <- as.data.frame(factor(shortList2$locus1))
colnames(locus1) <- "locus"

#join the two lists
locus <- rbind(locus1,locus2)
rm(locus1,locus2,shortList1,shortList2,shortList,data,hamList)

#remove duplicates
locus <- unique(locus)

# remove all loci identified as sharing sequence similarity with other loci but having lower call rate from filteredData
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]
#now you have to remove their counterparts, as this is 2row_SNP data!!
locus <- as.data.frame(sapply(locus,gsub,pattern="-",replacement="--"))
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]

rm(locus)

# format data for export
filteredData <- subset(filteredData, select=-loci )  # remove helper column (loci) from filteredData
filteredData$CloneID <- row.names(filteredData)   # add CloneID column (because otherwise it won't have a header)
filteredData <- filteredData[c(ncol(filteredData),1:(ncol(filteredData)-1))]   # move CloneID to first column
names(filteredData) <- sub("X","",names(filteredData))   # remove stupid Xs from the headers

# export data in same format that DArT_filter.pl imports and exports.
write.table(filteredData, "Report-<project_code>_2Row_hamFiltered6.csv",row.names=F, col.names=F, quote=F, sep=",")




################

# After running this filter it is necessary to pass through the perl script again, as due to the output of the hamming distances it makes it difficult to capture all at once.
# so run DArT_filter.pl again making sure you use the latest filtered version "Report-<project_code>_hamFiltered.csv" as the input.
# then run through this filter again -- below:

# set the minimum hamming distance below which loci will be removed
minHam = 5

# load the hamming distance data
hamList <- read.table("Report-<project_code>_2Row_hamFiltered6.ham_dist_values", sep=",", col.names=FALSE, fill=FALSE, strip.white=TRUE) 

# load the filtered SNP data
filteredData <- read.table('Report-<project_code>_2Row_hamFiltered6_filtered.csv', sep=',', header=TRUE, row.names=1) # whoops, I probably should've used row.names = F, but the script is written to fit this now (it just would've saved a step or two near the end)

hamData <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
hamData <- as.data.frame(sapply(hamData,gsub,pattern="(.*\\()",replacement="",perl=TRUE)) # replaces everything before last reverse bracket ")" with nothing ""

locus1 <- as.data.frame(sapply(hamList,gsub,pattern=":.*",replacement=""))   # replaces everything after and including first pipe "|" with nothing ""

locus2 <- as.data.frame(sapply(hamList,gsub,pattern="(\\).*)",replacement="",perl=TRUE)) # replaces everything after and including first reverse bracket ")" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern="[^:]*:[^:]*:[^:]*:",replacement="",perl=TRUE)) # does something weird, need to investigate... It seems it works in some cases (where there is only one pair) but not others (where there is more than one pair)
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=":.*",replacement="")) # replaces everything after and including first pipe "|" with nothing ""
locus2 <- as.data.frame(sapply(locus2,gsub,pattern=" ",replacement=""))   # replaces space in first position with nothing

# combine the three lists we generated into a single data frame
data <- cbind(locus1,locus2,hamData)
rm(locus1,locus2,hamData)

# format
colnames(data) <- c("locus1","locus2","hamData")
data$hamData <- as.numeric(data$hamData)

# extract locus pairs with hamming distances less than user defined minimum hamming distance
shortList <- subset(data, data$hamData < minHam)  # this doesn't work

#prepare filteredData for use by adding a locus column to the end that matches locus IDs used in newList
filteredData$loci <- gsub(":.*","",row.names(filteredData))

# for each pair of loci identify the locus with the highest call rate and remove it from that pair in newList
# 1. make new column for call rates for locus1
shortList$locus1CR <- shortList$locus1 %in% filteredData$loci
r = 1

# 2. find matching locus "locus1 %in% loci" pull out associated call rate "filteredData$CallRate" and add to new column
repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,1],]  
  shortList[r,4] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

#do the same for locus 2
shortList$locus2CR <- shortList$locus2 %in% filteredData$loci
r = 1

repeat{
  CR <- filteredData[filteredData$loci %in% shortList[r,2],]  
  shortList[r,5] <- CR[1,4]
  
  r = r+1
  if (nrow(shortList)+1 == r){
    break
  }
}

rm(CR,r,minHam)
# duplicate list
shortList1 <- shortList
shortList2 <- shortList

# from list 1 remove samples where locus1CR is greater than or equal to locus2CR
shortList1 <- subset(shortList1, shortList1$locus1CR >= shortList1$locus2CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus2 <- as.data.frame(factor(shortList1$locus2))
colnames(locus2) <- "locus"

shortList2 <- subset(shortList2, shortList2$locus2CR > shortList2$locus1CR)
# in these cases we want to ditch locus 2, so we extract a list of these names
locus1 <- as.data.frame(factor(shortList2$locus1))
colnames(locus1) <- "locus"

#join the two lists
locus <- rbind(locus1,locus2)
rm(locus1,locus2,shortList1,shortList2,shortList,data,hamList)

#remove duplicates
locus <- unique(locus)

# remove all loci identified as sharing sequence similarity with other loci but having lower call rate from filteredData
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]
#now you have to remove their counterparts, as this is 2row_SNP data!!
locus <- as.data.frame(sapply(locus,gsub,pattern="-",replacement="--"))
filteredData <- filteredData[! filteredData$loci %in% locus$locus,]

rm(locus)

# format data for export
filteredData <- subset(filteredData, select=-loci )  # remove helper column (loci) from filteredData
filteredData$CloneID <- row.names(filteredData)   # add CloneID column (because otherwise it won't have a header)
filteredData <- filteredData[c(ncol(filteredData),1:(ncol(filteredData)-1))]   # move CloneID to first column
names(filteredData) <- sub("X","",names(filteredData))   # remove stupid Xs from the headers

# export data in same format that DArT_filter.pl imports and exports.
write.table(filteredData, "Report-<project_code>_2Row_hamFiltered7.csv",row.names=F, col.names=F, quote=F, sep=",")



#you could keep doing this ad nauseum, due to the clunky solution used for this problem, but we're getting to the stage that we're only removing a trivial number of "loci" in each new pass, and theoretically if two "loci" should've been in the same cluster then they should be removed in the linkage test if they were missed here
#but either way there would be so few of them at this stage that it's unlikely they would be influencing 
#we can visualise 
filteredData <- read.table('Report-<project_code>_2Row_cloneFiltered_filtered.csv', sep=',', header=TRUE, row.names=1)
filteredData1 <- read.table('Report-<project_code>_2Row_hamFiltered_filtered.csv', sep=',', header=TRUE, row.names=1)
filteredData2 <- read.table('Report-<project_code>_2Row_hamFiltered2_filtered.csv', sep=',', header=TRUE, row.names=1)
filteredData3 <- read.table('Report-<project_code>_2Row_hamFiltered3_filtered.csv', sep=',', header=TRUE, row.names=1)
filteredData4 <- read.table('Report-<project_code>_2Row_hamFiltered4_filtered.csv', sep=',', header=TRUE, row.names=1)
filteredData5 <- read.table('Report-<project_code>_2Row_hamFiltered5_filtered.csv', sep=',', header=TRUE, row.names=1)
filteredData6 <- read.table('Report-<project_code>_2Row_hamFiltered6_filtered.csv', sep=',', header=TRUE, row.names=1)
filteredData7 <- read.table('Report-<project_code>_2Row_hamFiltered7.csv', sep=',', header=TRUE, row.names=1)
removal <- as.data.frame(c((nrow(filteredData)-nrow(filteredData1)),(nrow(filteredData1)-nrow(filteredData2)),(nrow(filteredData2)-nrow(filteredData3)),(nrow(filteredData3)-nrow(filteredData4)),(nrow(filteredData4)-nrow(filteredData5)),(nrow(filteredData5)-nrow(filteredData6)),(nrow(filteredData6)-nrow(filteredData7))))
removal$series <- row.names(removal)
colnames(removal)[1] <- "removal"
library(ggplot2)
ggplot(data=removal, aes(x=series, y=removal, group=1)) +
  geom_line() +
  geom_point()

rm(filteredData,filteredData1,filteredData2,filteredData3,filteredData4,filteredData5,filteredData6,filteredData7)
# so we've gone from removing ~7500 in the first pass, to maybe ~100 or so in the final pass, with a pattern of exponential decay.
# so this matches expectations and suggests only a few potential duplicate clones remain (I only wish I could write a better script to do what is needed, that could remove all the problem loci in a single pass)
# I just had an idea, maybe instead of only using the first column I could repeat the hamList and locus2 steps multiple times (adjusting them to use filter the appropriate subsequent data), and then bind them together to make a matrix for hamList and locus2. Then using the hamList after transposing and making a column for initial order I could sort by hamScore, make another column with sorted order and sort by initial order, then traspose back. If I copy the sorted order (now a row) to locus2 I can then sort (after transposing) to match the locus with the correct hamScore. --even if I only make the matrices with length 5, that's achieving 5x the work of one pass through this script.

# you may as well delete all the "in-between" files now to free up some space on your workstation, including "Report-<project_code>_2Row_cloneFiltered" and "Report-<project_code>_2Row_hamFiltered" through to "Report-<project_code>_2Row_hamFiltered6".

### NOW TO APPLY A FILTER TO THE 1ROW_SNP FILE TO MIRROR THE CHANGES TO THE 2ROW_SNP FILE READY FOR ANALYSIS PIPELINE ###


