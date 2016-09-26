# RADPipelineToNJTree.R
#
# written by Steven Myers 
# 26.ix.2016
#
# description: perform NJ tree analysis using 'ape' and 1. loci filtered for missing data and 2. neutral and unlinked loci
#
# input: various outputs from 'RAD_analysis_pipeline.R'. 
# 
# R dependencies: ape, adegenet, ggplot2, parallel


#############
### Steps ###

#1. load packages and data if necessary. 

#2. recode from minor and major alleles to SNP and Ref alleles

#2. format for ape (convert to IUPAC codes and concatenate)

#3. calculate distance matrix and generate NJ tree


setwd("C:/path/to/data/output/from/RAD_analysis_pipeline.R")
options(max.print=99999999) 


###############
###  STEP 1 ###    <<< load libraries, set directory, and load data >>>

#install.packages("ape")
#install.packages("adegenet")
#install.packages("ggplot2")


library(ape)
library(adegenet)
library(ggplot2)
library(parallel)


# load the data from STEP 5 of RAD_analysis_pipeline.R
x <-read.table("data_Filtered1row_RC_L15MD.csv", sep=',', header=T, row.names=1, stringsAsFactors=FALSE)
colnames(x) <- sub("X", "", colnames(x))




###############
### Step 2. ### <<< recode from minor and major alleles to SNP and Ref alleles >>>

y <- read.table("data_Filtered1row_recodedLoci.csv", sep=",", na = "NA", row.names = NULL, col.names = F) 

z <- as.data.frame(t(x))
z <- z[rownames(z) %in% y[,1],]

# recode the data
for(n in 1:nrow(z)){     
  recode <- (gsub(0,9,z[n,]))  #could use suppress warnings() here because it will give warning 'NAs introduced by coercion' in every case a change was made.
  recode <- (gsub(2,0,recode))
  z[n,]  <- (as.integer(gsub(9, 2, recode)))
}

# merge the recoded locus data and data for loci that didn't need recoding
x <- as.data.frame(t(x))
x <- rbind((x[!rownames(x) %in% rownames(z),]),z)

x <- as.data.frame(t(x))

rm(y,z,recode,n)




###############
### Step 3. ###    <<< format for 'ape' (this requires the original DArT Report file) >>>

#duplicate the data
xsampID <- rownames(x)    #first preserve sampleIDs
xlocID <- colnames(x)   #and locus IDs

xa <- x
xb <- x

xa = data.frame(sapply(xa, as.integer))   
xb = data.frame(sapply(xb, as.integer))   

xa[xa==0]<- 1
xa[is.na(xa)]<- 0
xb[xb==1]<- 2
xb[xb==0]<- 1
xb[is.na(xb)]<- 0

# reassign sample IDs to rownames in both duplicates
xa$sampID <- xsampID
xa<-xa[c(ncol(xa),1:(ncol(xa)-1))]
rownames(xa) = xa[,1]
xa<-subset(xa,select=-sampID)

xb$sampID <- xsampID
xb<-xb[c(ncol(xb),1:(ncol(xb)-1))]
rownames(xb) = xb[,1]
xb<-subset(xb,select=-sampID)

#merge the recoded data
xx <- cbind(xa,xb)
xlocID <- sub("X", "", colnames(xx))   #save new locus IDs

yy <- data.frame(sapply(xx, as.integer)) #save xx for later (need for making header elements)
yy$sampID <- xsampID   #reapply sample IDs
yy<-yy[c(ncol(yy),1:(ncol(yy)-1))]
rownames(yy) = yy[,1]
yy<-subset(yy,select=-sampID)

#sort columns by locus ID (by making rows and sorting by row, lol)
yy <- as.data.frame(t(yy))
yy <- yy[order(rownames(yy)),]
yy <- as.data.frame(t(yy))


# load the Filtered 1 Row DArT SNP data
SNP <- read.table('Report-<project_name>_1row_hamFiltered.csv', sep=',', header=F, row.names=NULL)

# remove first 7 rows
SNP <- SNP[8:nrow(SNP),1:ncol(SNP)]
# extract only columns 1 and 3 -- the information we need to reassign 1 & 2 allele calls with their corresponding base call
SNP <- SNP[,c(1,9)]
colnames(SNP) <- c("locus","SNP")
SNP$locus <- gsub("\\|.*","",SNP$locus )

SNP1 <- gsub(".*:","",SNP$SNP )
SNP1 <- gsub(">.*","",SNP1 )
SNP2 <- gsub(".*>","",SNP$SNP )

SNP <- cbind("locus" = SNP$locus,SNP1,SNP2)
SNP <- as.data.frame(SNP)
SNP$locus <- as.character(SNP$locus)
SNP$SNP1 <- as.character(SNP$SNP1)
SNP$SNP2 <- as.character(SNP$SNP2)

# transpose and append a locus IDs column
yy <- as.data.frame(t(yy))
yy <- cbind("locus" = row.names(yy),yy)
# fix locus names
yy$locus <- sub("\\.1","",yy$locus)
yy$locus <- sub("X","",yy$locus)

#trim 2Row locus list
SNP <- SNP[SNP$locus %in% yy$locus,]

# use a loop to convert the data from 1,2 to ATGC codes, and then from separate columns per allele per locus to IUPAC codes in a single column

# prepare elements for the loop
r=1
yy3 <- matrix(NA, nrow=1,ncol=(ncol(yy)-1))
yy3 <- as.data.frame(yy3)
yy3 <- cbind("locus" = NA, yy3)

# do the loop
repeat{
  
  # search the data for instances of the locus in row r in the locus list and subset to a new data frame
  yy2 <- yy[yy$locus %in% SNP[r,1],]    
  # convert 1s to the appropriate ATGC code value representing the major allele for that locus
  yy2[yy2==1] <- SNP$SNP1[r]
  # convert 2s to the appropriate ATGC code value representing the minor allele for that locus
  yy2[yy2==2] <- SNP$SNP2[r]
  # store the name of the locus that was just processed
  locus <- yy2[1,1]
  # transpose data for next process
  yy2 <- as.data.frame(t(yy2[2:ncol(yy2)]))
  # convert pairs of ATGC codes to IUPAC codes in a new column
  yy2$IUPAC[yy2[,1] == "A" & yy2[,2] == "A"] <- "A"
  yy2$IUPAC[yy2[,1] == "T" & yy2[,2] == "T"] <- "T"
  yy2$IUPAC[yy2[,1] == "G" & yy2[,2] == "G"] <- "G"
  yy2$IUPAC[yy2[,1] == "C" & yy2[,2] == "C"] <- "C"
  yy2$IUPAC[(yy2[,1] == "A" & yy2[,2] == "G") | (yy2[,1] == "G" & yy2[,2] == "A")] <- "R"
  yy2$IUPAC[(yy2[,1] == "T" & yy2[,2] == "C") | (yy2[,1] == "C" & yy2[,2] == "T")] <- "Y"
  yy2$IUPAC[(yy2[,1] == "G" & yy2[,2] == "C") | (yy2[,1] == "C" & yy2[,2] == "G")] <- "S"
  yy2$IUPAC[(yy2[,1] == "T" & yy2[,2] == "A") | (yy2[,1] == "A" & yy2[,2] == "T")] <- "W"
  yy2$IUPAC[(yy2[,1] == "T" & yy2[,2] == "G") | (yy2[,1] == "G" & yy2[,2] == "T")] <- "K"
  yy2$IUPAC[(yy2[,1] == "A" & yy2[,2] == "C") | (yy2[,1] == "C" & yy2[,2] == "A")] <- "M"
  yy2$IUPAC[yy2[,1] == 0 | yy2[,2] == 0] <- "N"    # we use "N" instead of "-" because we assume missing data is only missing due to lack of coverage (i.e. the locus exists but wasn't captured in sequencing) due to stachastic variation in various processes such as PCR rather than actually absent due to a cut-site mutation, in which case we would use "-" to indicate a gap in the alignment (not that I know these are treated differently in any analyses)
  
  # extract the new column and transpose
  yy2 <- as.data.frame(t(yy2[,3]))
  yy2 <- cbind(locus, yy2)
  # add the current locus being processed to a data frame containing the previously processed loci
  yy3 <- rbind(yy3,yy2)
  
  r=r+1
  
  if (r == nrow(SNP)+1){
    break
  }
}

# now for some final formatting of the data frame
# remove the blank row
yy3 <- yy3[-c(1),]
# transpose
yy3 <- as.data.frame(t(yy3))
# remove the row containing locus names
yy3 <- yy3[-c(1),]
# add sample IDs
rownames(yy3) <- rownames(x)


####### THIS PART NOT NECESSARY; NEEDS TO BE MODIFIED TO FIT DATA IF USED ########
# load up the PC1scores data to access sample info (alternatively you could do it the harder way, as we did to get the info in the PC1scores file, by loading in the masterfile)
y <- read.table("data_info.csv", sep=",", na = "NA", row.names = 1, header = T)
# sort both data frames by sample ID (rownames)
y <- y[order(rownames(y)),]
yy3 <- yy3[order(rownames(yy3)),]
identical(rownames(yy3),rownames(y))
# concatenate all columns into a single column
yy3 <- as.data.frame(do.call(paste, as.data.frame(yy3, stringsAsFactors=FALSE, sep="")))
# change sample labels in data to show sample.number and locality
rownames(yy3) <- paste(y$ABTC,y$locality,sep="_")
# remove spaces in row names (sample labels) to avoid problems with reading the data
rownames(yy3) <- gsub(" ","_",rownames(yy3))
####### END #########


# calculate number of samples and length of sequence (i.e. number of loci)
nSamp <- nrow(x)
nLoc <- (ncol(x))
# make element out of nSamp and nLoc, separated by a space
head <- paste(nSamp, nLoc, sep=" ")
# add new blank row with nSamp and nLoc, separated by a space, as the name
yy4 <- rbind(NA, yy3)
rownames(yy4)[1] <- head
#replace spaces in sequences with nothing (i.e. remove the spaces)
yy4[,1] <- gsub(" ","",yy4[,1])

# write to file (phylip file type)
write.table(yy4, "data_Filtered1row_RC_L15MD_IUPAC.phy", sep=" ", na = "", row.names = T, col.names = F, quote=F)

# tidy up the global environment
rm(xa,xb,xx,SNP,y,yy,yy2,yy3,yy4,head,nSamp,nLoc,SNP1,SNP2,xlocID,xsampID,r,locus)

